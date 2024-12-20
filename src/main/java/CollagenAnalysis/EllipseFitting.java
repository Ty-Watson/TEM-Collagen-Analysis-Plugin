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
import static CollagenAnalysis.ImageProcessingUtils.*;

public class EllipseFitting {

    private ArrayList<double[]> centroids;
    private int[] cluster_idx;
    private ArrayList<double[]> fibrilPixels = new ArrayList<>();
    public static final int DEGREES = 361; // 0 to 360 degrees
    ArrayList<double[]> nonboundryCentroids = new ArrayList<>();
    private double[][] post_fibril;
    //scaled up image
    private ImagePlus imp;

    //scaled up image
    private ImageProcessor original;
    //scaled up processor
    private ImageProcessor ip;
    private ImageCanvas canvas;
    public ArrayList<Ellipse> fibrilEllipses = new ArrayList<>();

    public EllipseFitting(ImagePlus imp, ArrayList<double[]> fibrilPixels, ArrayList<double[]> centroids, int[] cluster_idx, double[][] post_fibril, ArrayList<double[]> nonBoundryCentroids){
        this.fibrilPixels = fibrilPixels;
        this.centroids = centroids;
        this.cluster_idx = cluster_idx;
        this.post_fibril = post_fibril;
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
        FibrilUtils util = new FibrilUtils();

        ImageProcessor test = ip.duplicate();

        util.calculateMuAndCovarianceForEachFibril(fibrilPixels, centroids, cluster_idx);
        util.validateCovariance(test);

        double[] areaForEachFibril = computeFibrilArea(post_fibril, centroids.size());

        ArrayList<double[][]> covarianceForEachFibril = util.covarianceForEachFibril;
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
                e.area = area;
                e.postProbArea = areaForEachFibril[i];
                e.majorRadius = majorRadius;
                e.minorRadius = minorRadius;
                e.aspectRatio = majorRadius / minorRadius;

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

                // Store the ellipse coordinates
                e.xEllipse = xEllipse;
                e.yEllipse = yEllipse;


                // Extract the angle of the major axis
                double[] eigenVector = V.getColumnVector(eigenValues[0] > eigenValues[1] ? 0 : 1).toArray();
                double angle = Math.atan2(eigenVector[1], eigenVector[0]);
                e.angle = Math.toDegrees(angle);

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
            double x = e.centerX;
            double y = e.centerY;
            PointDrawer.drawCross(ip,x,y);
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
            if (p.contains(clickX, clickY) && p.contains(e.centerX, e.centerY)) { // dont have to scale center of ellipse coordinates because they were calculated from scaled up polygon already
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
