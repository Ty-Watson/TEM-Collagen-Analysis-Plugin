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

    private ArrayList<Polygon> ellipses = new ArrayList<>();

    private ArrayList<double[]> radius_pix = new ArrayList<>();

    private ArrayList<Double> angle_of_major_axis = new ArrayList<Double>();

    ArrayList<double[]> nonboundryCentroids = new ArrayList<>();
    private double[][] post_fibril;

    private ImagePlus imp;
    private ImageProcessor ip;
    private ImageCanvas canvas;

    public EllipseFitting(ImagePlus imp, ArrayList<double[]> fibrilPixels, ArrayList<double[]> centroids, int[] cluster_idx, double[][] post_fibril, ArrayList<double[]> nonBoundryCentroids){
        this.fibrilPixels = fibrilPixels;
        this.centroids = centroids;
        this.cluster_idx = cluster_idx;
        this.post_fibril = post_fibril;
        this.nonboundryCentroids = nonBoundryCentroids;
        ip = imp.getProcessor().duplicate().convertToColorProcessor();
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
        ip.snapshot();
        //recalculate these when user deletes ellipse
        xEllipses.clear();
        yEllipses.clear();
        radius_pix.clear();
        angle_of_major_axis.clear();

        FibrilUtils util = new FibrilUtils();

        util.calculateMuAndCovarianceForEachFibril(fibrilPixels, centroids, cluster_idx);

        ArrayList<double[][]> covarianceForEachFibril = util.covarianceForEachFibril;
        //ArrayList<RealMatrix> covariance = util.covariance;
        ArrayList<double[]> mus = util.muForEachFibril;

        // Fit ellipses
        for (int i = 0; i < centroids.size(); i++) {
            double[] centroid = centroids.get(i);
            boolean notBoundry = nonboundryCentroids.contains(centroid);
            if(notBoundry){
                double[] mu = mus.get(i);
                RealMatrix covMatrix = new Array2DRowRealMatrix(covarianceForEachFibril.get(i));
                //RealMatrix covMatrix = new Array2DRowRealMatrix(covariance.get(i));
                EigenDecomposition eig = new EigenDecomposition(covMatrix);

                double[] eigenValues = eig.getRealEigenvalues();
                RealMatrix V = eig.getV();

                double[] radiusPix = {
                        2 * Math.sqrt(Math.max(eigenValues[0], eigenValues[1])), //major radius
                        2 * Math.sqrt(Math.min(eigenValues[0], eigenValues[1]))  //minor radius
                };
                radius_pix.add(radiusPix);

                double[] xEllipse = new double[DEGREES];
                double[] yEllipse = new double[DEGREES];

                for (int theta = 0; theta < DEGREES; theta++) {
                    //System.out.println("Degrees: " + theta);
                    double radians = Math.toRadians(theta);
                    double cosTheta = radiusPix[0] * Math.cos(radians);
                    double sinTheta = radiusPix[1] * Math.sin(radians);

                    //multiply RealMatrix by vector
                    double[] ellipsePoint = V.operate(new double[]{cosTheta, sinTheta});
                    xEllipse[theta] = ellipsePoint[0] + mu[0];
                    yEllipse[theta] = ellipsePoint[1] + mu[1];
                }

                // Store or display the ellipse coordinates
                xEllipses.add(xEllipse);
                yEllipses.add(yEllipse);


                // Extract the angle of the major axis
                double[] eigenVector = V.getColumnVector(eigenValues[0] > eigenValues[1] ? 0 : 1).toArray();
                double angle = Math.atan2(eigenVector[1], eigenVector[0]);
                angle_of_major_axis.add(Math.toDegrees(angle));
                System.out.println("Ellipse " + i + ": Major Radius = " + radiusPix[0] + ", Minor Radius = " + radiusPix[1] + ", Angle = " + angle);
                System.out.println("NonBoundryCentroidCount: " + nonboundryCentroids.size());
                System.out.println("all centroids: " + centroids.size());



//            System.out.printf("Ellipse %d coordinates (x, y):\n", i);
//            for (int j = 0; j < DEGREES; j++) {
//                System.out.printf("[%f, %f]\n", xEllipse[j], yEllipse[j]);
//            }
            }

        }
        generateEllipsePolygons();
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
        OverlayManager.drawEllipsesOnImage(ip, ellipses);
        imp.updateAndDraw();
        OverlayManager.overlayCentroids(ip, nonboundryCentroids.toArray(new double[nonboundryCentroids.size()][]), Color.red.getRGB());
        imp.updateAndDraw();
    }

    //Ratio of major to minor radius
    public ArrayList<Double> getAspectRatios(){
        ArrayList<Double> ratios = new ArrayList<>();
        for (double[] radiusPix : radius_pix) {
            double ratio = radiusPix[0] / radiusPix[1];
            ratios.add(ratio);
        }
        return ratios;
    }

    private void setCanvasEvents(){
        canvas = imp.getCanvas();
        canvas.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.getButton() == MouseEvent.BUTTON1) { // Left click delete ellipse
                    int clickX = canvas.offScreenX(e.getX());
                    int clickY = canvas.offScreenY(e.getY());

                    // Find and remove the nearest centroid to the click
                    double[] nearestCentroid = findNearestCentroid(nonboundryCentroids, clickX, clickY);
                    if (nearestCentroid != null) {
                        int removedIndex = findAndRemoveNearestEllipse(clickX, clickY);
                        if (removedIndex != -1) {
                            //if click is not in ellipse, then dont remove centroid
                            nonboundryCentroids.remove(nearestCentroid);
                            ip.reset();
                            drawEllipses();
                        }
                    }
                }
            }
        });

    }

    private int findAndRemoveNearestEllipse(int clickX, int clickY) {
        for (int i = 0; i < ellipses.size(); i++) {
            Polygon p = ellipses.get(i);
            if (p.contains(clickX, clickY)) {
                xEllipses.remove(i);
                yEllipses.remove(i);
                radius_pix.remove(i);
                angle_of_major_axis.remove(i);
                ellipses.remove(i);
                return i;
            }
        }
        return -1; // No ellipse contains the point
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
