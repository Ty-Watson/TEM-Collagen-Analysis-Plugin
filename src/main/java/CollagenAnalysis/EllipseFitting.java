package CollagenAnalysis;

import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.process.ImageProcessor;
import org.apache.commons.math3.linear.*;

import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;


import static CollagenAnalysis.Constants.SCALE;
import static CollagenAnalysis.Constants.pointSize;
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

    //scaled up image
    private ImagePlus imp;

    //scaled up image
    private ImageProcessor original;
    //scaled up processor
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
        setImageAndProcessors(imp);
        ip.snapshot();
        this.imp.show();
        setCanvasEvents();
    }
    private void setImageAndProcessors(ImagePlus imp){
        //show user scaled up image (original image they submitted, but we downsize image for faster processing)
        ImageProcessor imageProcessor = imp.getProcessor().duplicate().convertToColorProcessor();
        int h = imageProcessor.getHeight() * SCALE;
        int w = imageProcessor.getWidth() * SCALE;
        this.ip = imageProcessor.resize(w, h);
        this.original = ip.duplicate();
        this.imp = new ImagePlus("Ellipse Fitting", ip);
    }


    public  void fitEllipses() {
        //recalculate these when user deletes ellipse
        xEllipses.clear();
        yEllipses.clear();
        radius_pix.clear();
        angle_of_major_axis.clear();

        FibrilUtils util = new FibrilUtils();

        ImageProcessor test = ip.duplicate();

        util.calculateMuAndCovarianceForEachFibril(fibrilPixels, centroids, cluster_idx);
        util.validateCovariance(test);

        double[] areaForEachFibril = computeFibrilArea(post_fibril, centroids.size());

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

                // Calculate the area of the ellipse. using postProbArea in results
                double area = Math.PI * majorRadius * minorRadius;
                e.area = area; //not used in results. Just to compare ellipse area with fibril area (should be similar)
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

                //generate the ellipse that will be shown on the image for the current ellipse
                e.generateEllipsePolygon();
                fibrilEllipses.add(e);
                System.out.println("Ellipse " + i + ": Major Radius = " + majorRadius + ", Minor Radius = " + minorRadius + ", Area = " + area + ", Angle = " + angle);
            }

        }
        drawEllipses();
    }

    private void drawEllipses(){
        ImageProcessor ip = imp.getProcessor();
        //draw ellipse polygon
        for(Ellipse e : fibrilEllipses){
            ip.setColor(Color.blue);
            ip.drawPolygon(e.ellipsePolygon);
            //draw centroid
            // Size of the maxima point to be drawn (radius of the circle around the maxima point)

            int color = Color.red.getRGB();

            // Draw each maximum point as a red circle on the RGB image

            //have to scale up ellipse centroid coordinate to show user on full image
            double x = e.x * SCALE;
            double y = e.y * SCALE;

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
                    //do not need to scale down click because the ellipse polygon is generated for the full size image
                    //so we need to check if click is in the polygon to delete it
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
            if (p.contains(clickX, clickY) && p.contains(e.x * SCALE, e.y * SCALE)) { // have to scale up ellipse centroid coordinates because the polygon is generated on the full image
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
    public ImagePlus getImage(){
        return imp;
    }
}
