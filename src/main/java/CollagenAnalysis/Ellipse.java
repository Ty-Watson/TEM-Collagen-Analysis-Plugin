package CollagenAnalysis;

import java.awt.*;
import java.util.ArrayList;

import static CollagenAnalysis.Constants.SCALE;

public class Ellipse {
    double[] xEllipse;
    double[] yEllipse;

    public double x;

    public double y;

    public double aspectRatio;
    public double aspectRatio_nm;
    public double area;
    public double postProbArea;
    public double area_nm;
    public double angle;
    public double majorRadius;
    public double minorRadius;

    public double majorRadius_nm;
    public double minorRadius_nm;
    double[] mu;
    double[][] cov;

    double componentProportion;
    Polygon ellipsePolygon;


    public void generateEllipsePolygon() {

        int DEGREES = 361;
        int[] xPoints = new int[DEGREES];
        int[] yPoints = new int[DEGREES];

        // Prepare the coordinates for drawing
        for (int j = 0; j < DEGREES; j++) {
            //have to scale up the ellipse to show user on the full image
            xPoints[j] = (int) Math.round(xEllipse[j] * SCALE);
            yPoints[j] = (int) Math.round(yEllipse[j] * SCALE);
        }

        // Create a polygon from the ellipse points
        ellipsePolygon = new Polygon(xPoints, yPoints, DEGREES);

    }

    public void scale(int scale){
        // Scale center point of the ellipse
        x *= scale;
        y *= scale;

        // Scale the mean vector (mu)
        mu[0] *= scale;  // Scale x component of the mean
        mu[1] *= scale;  // Scale y component of the mean

        // Scale radii
        majorRadius *= scale;
        minorRadius *= scale;
        aspectRatio *= scale;

        // Area scales by the square of the scale factor
        area *= Math.pow(scale, 2);
        postProbArea *= Math.pow(scale,2);


        // Scale the ellipse points
        for (int i = 0; i < xEllipse.length; i++) {
            xEllipse[i] *= scale;
            yEllipse[i] *= scale;
        }

        // Scale the covariance matrix (cov)
        for (int i = 0; i < cov.length; i++) {
            for (int j = 0; j < cov[i].length; j++) {
                cov[i][j] *= Math.pow(scale, 2);
            }
        }
        componentProportion *= scale;

        // Update the polygon after scaling
        //generateEllipsePolygon();
    }
    public void convertToNM(double nanometers_over_pixels){
        //this is using major and minor radius to calculate area
        //area_nm = area * Math.pow(nanometers_over_pixels, 2);
        area_nm = postProbArea * Math.pow(nanometers_over_pixels, 2);
        majorRadius_nm = majorRadius * nanometers_over_pixels;
        minorRadius_nm = minorRadius * nanometers_over_pixels;
        aspectRatio_nm = aspectRatio * nanometers_over_pixels;
    }

}
