package CollagenAnalysis;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.process.ImageProcessor;
import ij.util.ArrayUtil;
import org.apache.commons.math3.linear.*;

import java.awt.*;
import java.awt.event.AWTEventListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Arrays;

import static CollagenAnalysis.ImageProcessingUtils.*;

public class EllipseFitting {

    private ArrayList<double[]> centroids;
    private int[] cluster_idx;
    private ArrayList<double[]> fibrilPixels = new ArrayList<>();
    public static final int DEGREES = 361; // 0 to 360 degrees

    private ArrayList<double[]> xEllipses = new ArrayList<>();
    private ArrayList<double[]> yEllipses = new ArrayList<>();

    public ArrayList<Polygon> ellipses = new ArrayList<>();

    private ArrayList<double[]> radius_pix = new ArrayList<>();

    public ArrayList<Double> area_pix = new ArrayList<>();

    private ArrayList<Double> angle_of_major_axis = new ArrayList<Double>();

    ArrayList<double[]> nonboundryCentroids = new ArrayList<>();
    private double[][] post_fibril;

    private ImagePlus imp;

    private ImageProcessor original;
    private ImageProcessor ip;
    private ImageCanvas canvas;

    public double[] componentProportion;

    public ArrayList<Ellipse> fibrilEllipses = new ArrayList<>();

    public EllipseFitting(ImagePlus imp, ArrayList<double[]> fibrilPixels, ArrayList<double[]> centroids, int[] cluster_idx, double[][] post_fibril, double[] componentProportion, ArrayList<double[]> nonBoundryCentroids){
        this.fibrilPixels = fibrilPixels;
        this.centroids = centroids;
        this.cluster_idx = cluster_idx;
        this.post_fibril = post_fibril;
        this.componentProportion = componentProportion;
        this.nonboundryCentroids = nonBoundryCentroids;
        this.ip = imp.getProcessor().duplicate().convertToColorProcessor();
        this.original = imp.getProcessor().duplicate().convertToColorProcessor();
        ip.snapshot();
        this.imp = new ImagePlus("Ellipse Fitting", ip);
        this.imp.show();
        setCanvasEvents();
    }
    private static double[] createThetaArray() {
        double[] theta = new double[361];
        for (int i = 0; i <= 360; i++) {
            theta[i] = Math.toRadians(i);
        }
        return theta;
    }

    public  void fitEllipses(ArrayList<Boolean> isBoundryCentroid) {
        //recalculate these when user deletes ellipse
        xEllipses.clear();
        yEllipses.clear();
        radius_pix.clear();
        angle_of_major_axis.clear();

        FibrilUtils util = new FibrilUtils();

        ImageProcessor test = ip.duplicate();

        util.calculateMuAndCovarianceForEachFibril(fibrilPixels, centroids, cluster_idx);
        util.validateCovariance(test);
        double[] areaForEachFibril = computeFibrilArea(post_fibril, 2, centroids.size());

        ArrayList<double[][]> covarianceForEachFibril = util.covarianceForEachFibril;
        //ArrayList<RealMatrix> covariance = util.covariance;
        ArrayList<double[]> mus = util.muForEachFibril;
        //remove prob centroids
        for(double[] pCentroids : util.problematicCentroids){
            centroids.remove(pCentroids);
            nonboundryCentroids.remove(pCentroids);
        }

        // Fit ellipses
        for (int i = 0; i < centroids.size(); i++) {
            Ellipse e = new Ellipse();
            double[] centroid = centroids.get(i);
            boolean notBoundry = nonboundryCentroids.contains(centroid);
            if(notBoundry){
                e.x = centroid[0];
                e.y = centroid[1];
                e.componentProportion = componentProportion[i];
                double[] mu = mus.get(i);
                e.mu = mu;
                RealMatrix covMatrix = new Array2DRowRealMatrix(covarianceForEachFibril.get(i));
                e.cov = covarianceForEachFibril.get(i);
                //RealMatrix covMatrix = new Array2DRowRealMatrix(covariance.get(i));
                EigenDecomposition eig = new EigenDecomposition(covMatrix);

                double[] eigenValues = eig.getRealEigenvalues();
                RealMatrix V = eig.getV();

                // Calculate major and minor radii
                double majorRadius = 2 * Math.sqrt(Math.max(eigenValues[0], eigenValues[1])); // semi-major radius
                double minorRadius = 2 * Math.sqrt(Math.min(eigenValues[0], eigenValues[1])); // semi-minor radius

                // Calculate the area of the ellipse
                double area = Math.PI * majorRadius * minorRadius;
                e.area = area;
                e.postProbArea = areaForEachFibril[i];
                e.majorRadius = majorRadius;
                e.minorRadius = minorRadius;
                e.aspectRatio = majorRadius / minorRadius;
                area_pix.add(area);


                // Store the radius information
                radius_pix.add(new double[]{majorRadius, minorRadius}); // Storing full diameter for reference



                double[] xEllipse = new double[DEGREES];
                double[] yEllipse = new double[DEGREES];

                for (int theta = 0; theta < DEGREES; theta++) {
                    //System.out.println("Degrees: " + theta);
                    double radians = Math.toRadians(theta);
                    double cosTheta = majorRadius * Math.cos(radians);
                    double sinTheta = minorRadius * Math.sin(radians);

                    //multiply RealMatrix by vector
                    double[] ellipsePoint = V.operate(new double[]{cosTheta, sinTheta});
                    xEllipse[theta] = ellipsePoint[0] + mu[0];
                    yEllipse[theta] = ellipsePoint[1] + mu[1];
                }

                // Store or display the ellipse coordinates
                xEllipses.add(xEllipse);
                e.xEllipse = xEllipse;
                yEllipses.add(yEllipse);
                e.yEllipse = yEllipse;


                // Extract the angle of the major axis
                double[] eigenVector = V.getColumnVector(eigenValues[0] > eigenValues[1] ? 0 : 1).toArray();
                double angle = Math.atan2(eigenVector[1], eigenVector[0]);
                e.angle = Math.toDegrees(angle);
                angle_of_major_axis.add(Math.toDegrees(angle));
                // Print the radii and area
                e.generateEllipsePolygon();
                fibrilEllipses.add(e);
                System.out.println("Ellipse " + i + ": Major Radius = " + majorRadius + ", Minor Radius = " + minorRadius + ", Area = " + area + ", Angle = " + angle);
               // System.out.println("NonBoundryCentroidCount: " + nonboundryCentroids.size());
                //System.out.println("all centroids: " + centroids.size());



//            System.out.printf("Ellipse %d coordinates (x, y):\n", i);
//            for (int j = 0; j < DEGREES; j++) {
//                System.out.printf("[%f, %f]\n", xEllipse[j], yEllipse[j]);
//            }
            }

        }
        drawEllipses();
    }

    private void generateEllipsePolygons(){
        // Iterate over all ellipses
        for (int i = 0; i < xEllipses.size(); i++) {
            double[] xCoords = xEllipses.get(i);
            double[] yCoords = yEllipses.get(i);

            if (xCoords.length != yCoords.length || xCoords.length != DEGREES) {
                throw new IllegalArgumentException("Each ellipse coordinate array must have the same length and match the expected number of degrees.");
            }

            int[] xPoints = new int[DEGREES];
            int[] yPoints = new int[DEGREES];

            // Prepare the coordinates for drawing
            for (int j = 0; j < DEGREES; j++) {
                xPoints[j] = (int) Math.round(xCoords[j]);
                yPoints[j] = (int) Math.round(yCoords[j]);
            }

            // Create a polygon from the ellipse points
            Polygon p = new Polygon(xPoints, yPoints, DEGREES);
            ellipses.add(p);
        }
    }

    private void drawEllipses(){
        ImageProcessor ip = imp.getProcessor();
        //draw ellipse polygon
        for(Ellipse e : fibrilEllipses){
            ip.setColor(Color.blue);
            ip.drawPolygon(e.ellipsePolygon);
            //draw centroid
            // Size of the maxima point to be drawn (radius of the circle around the maxima point)
            int pointSize = 2;
            int color = Color.red.getRGB();

            // Draw each maximum point as a red circle on the RGB image

            double x = e.x;
            double y = e.y;

            // Draw a circle or a larger point at (x, y) in red
            for (int dx = -pointSize; dx <= pointSize; dx++) {
                for (int dy = -pointSize; dy <= pointSize; dy++) {
                    if (dx * dx + dy * dy <= pointSize * pointSize) {
                        double newX = x + dx;
                        double newY = y + dy;
                        if (newX >= 0 && newX < ip.getWidth() && newY >= 0 && newY < ip.getHeight()) {
                            ip.set((int)newX, (int)newY, color);
                        }
                    }
                }
            }
        }
        imp.updateAndDraw();


    }
    private void setCanvasEvents(){
        canvas = imp.getCanvas();
        canvas.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.getButton() == MouseEvent.BUTTON1) { // Left click delete ellipse
                    int clickX = canvas.offScreenX(e.getX());
                    int clickY = canvas.offScreenY(e.getY());

                    int removedIndex = findAndRemoveNearestEllipse(clickX, clickY);
                    //do not reset image if user did not click in a ellipse to remove it
                    if (removedIndex != -1) {
                        resetImage();
                        drawEllipses();
                    }

                    e.consume();
                }
            }
        });

    }

    private int findAndRemoveNearestEllipse(int clickX, int clickY) {
        for (int i = 0; i < fibrilEllipses.size(); i++) {
            Ellipse e = fibrilEllipses.get(i);
            Polygon p = e.ellipsePolygon;
            //if the click is in the polygon and centroid is in ellipse
            if (p.contains(clickX, clickY) && p.contains(e.x, e.y)) {
                fibrilEllipses.remove(i);
                System.out.println("# of ellipses: " + fibrilEllipses.size());
                return i;
            }
        }
        return -1; // No ellipse contains the point
    }
    public void resetImage() {

        // Reset to the original processor
        imp.setProcessor(original.duplicate());

        // Update the display to show the reset image
        imp.updateAndDraw();
    }

    private double[] findNearestEllipseCentroid(int clickX, int clickY){
        double[] nearest = null;
        double minDistance = Double.MAX_VALUE;

        for(Ellipse e : fibrilEllipses){
            double distance = Math.sqrt(Math.pow(e.x - clickX, 2) + Math.pow(e.y - clickY, 2));
            if (distance < minDistance) {
                minDistance = distance;
                nearest = new double[]{e.x, e.y};
            }
        }
        return nearest;
    }



    public ArrayList<double[]> getxEllipses() {
        return xEllipses;
    }

    public ArrayList<double[]> getyEllipses() {
        return yEllipses;
    }

    public ArrayList<Double> getAngle_of_major_axis() {
        return angle_of_major_axis;
    }

    public ImagePlus getImage(){
        return imp;
    }



    public ArrayList<double[]> getRadius_pix() {
        return radius_pix;
    }
}
