package CollagenAnalysis;

import java.awt.*;
import java.util.ArrayList;

import static CollagenAnalysis.Constants.SCALE;

public class Ellipse {
    double[] xEllipse;
    double[] yEllipse;

    public double x;
    //center of ellipse coordinates
    public double centerX;
    public double centerY;
    public double y;
    public double aspectRatio;
    public double area;
    public double area_nm;
    public double postProbArea;
    public double postProbArea_nm;
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
            //have to scale up the ellipse coordinates to show user on the full image
            xPoints[j] = (int) Math.round(xEllipse[j] * SCALE);
            yPoints[j] = (int) Math.round(yEllipse[j] * SCALE);
        }

        // Create a polygon from the ellipse points
        ellipsePolygon = new Polygon(xPoints, yPoints, DEGREES);
        //get center coordinates from generated ellipse
        centerX = ellipsePolygon.getBounds2D().getCenterX();
        centerY = ellipsePolygon.getBounds2D().getCenterY();


    }

    public void scale(int scale){
        // Scale radii
        majorRadius *= scale;
        minorRadius *= scale;

        // Area scales by the square of the scale factor
        area *= Math.pow(scale, 2);
        postProbArea *= Math.pow(scale,2);


        // Scale the ellipse points
        for (int i = 0; i < xEllipse.length; i++) {
            xEllipse[i] *= scale;
            yEllipse[i] *= scale;
        }
    }
    public void convertToNM(double nanometers_over_pixels){
        postProbArea_nm = postProbArea * Math.pow(nanometers_over_pixels, 2);
        //this is using major and minor radius to calculate area
        area_nm = area * Math.pow(nanometers_over_pixels, 2);
        majorRadius_nm = majorRadius * nanometers_over_pixels;
        minorRadius_nm = minorRadius * nanometers_over_pixels;
    }

}