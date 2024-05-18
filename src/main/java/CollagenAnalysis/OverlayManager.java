package CollagenAnalysis;


import de.alsclo.voronoi.Voronoi;
import de.alsclo.voronoi.graph.Graph;
import de.alsclo.voronoi.graph.Point;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.ArrayList;

import static CollagenAnalysis.EllipseFitting.DEGREES;

public  class OverlayManager {

    public static void overlayContourLines(ImageProcessor contour , double[][][] post_all, int color){
        int numFibrils = post_all.length;
        int height = contour.getHeight();
        int width = contour.getWidth();

        for (int m = 1; m < height - 1; m++) {
            for (int n = 1; n < width - 1; n++) {
                for (int i = 0; i < numFibrils; i++) {
                    double[] neighbors = {
                            post_all[i][m - 1][n - 1], post_all[i][m - 1][n], post_all[i][m - 1][n + 1],
                            post_all[i][m][n - 1], post_all[i][m][n], post_all[i][m][n + 1],
                            post_all[i][m + 1][n - 1], post_all[i][m + 1][n], post_all[i][m + 1][n + 1]
                    };

                    double min = Double.MAX_VALUE;
                    double max = Double.MIN_VALUE;

                    for (double prob : neighbors) {
                        if (prob < min) min = prob;
                        if (prob > max) max = prob;
                    }

                    if (min < 0.5 && max > 0.5) {
                        // Setting pixel color to blue if it's a boundary
                        contour.putPixel(n, m, color);
                    }
                }

            }
        }
    }

    public static void overlayVoronoi(ImageProcessor ip, Voronoi voronoi, int color){
        ip.setColor(color);
        Graph g = voronoi.getGraph();

        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            de.alsclo.voronoi.graph.Point a = e.getA().getLocation();
            Point b = e.getB().getLocation();
            ip.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });
    }
    public static void overlayCentroids(ImageProcessor ip, double[][] centroids, int color){


        // Size of the maxima point to be drawn (radius of the circle around the maxima point)
        int pointSize = 2;

        // Draw each maximum point as a red circle on the RGB image
        for (double[] point : centroids) {
            double x = point[0];
            double y = point[1];

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
    }
    public static void overlayPixels(ImageProcessor ip, double[][] pixels, int color){

        for (double[] pixel : pixels) {
            int x = (int) pixel[0];
            int y = (int) pixel[1];

            if (x >= 0 && x < ip.getWidth() && y >= 0 && y < ip.getHeight()) {
                ip.set(x, y, color);
            }
        }
    }
    public static void drawEllipsesOnImage(ImageProcessor ip, ArrayList<double[]> xEllipses, ArrayList<double[]> yEllipses, int nCentroid) {
        if (xEllipses.size() != yEllipses.size()) {
            throw new IllegalArgumentException("The number of x and y ellipse coordinate arrays must match.");
        }

        ip.setColor(Color.BLUE);

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
            ip.drawPolygon(p);
        }

    }
    public static void overlayColoredClusterAssignments(ImageProcessor clusterColorsProcessor,ArrayList<double[]> fibrilPixels, int[] clusterAssignments){
        Color[] clusterColors = new Color[] {
                new Color(255, 0, 0),      // Red
                new Color(0, 255, 0),      // Green
                new Color(0, 0, 255),      // Blue
                new Color(255, 255, 0),    // Yellow
                new Color(255, 0, 255),    // Magenta
                new Color(0, 255, 255),    // Cyan
                new Color(255, 165, 0),    // Orange
                new Color(128, 0, 128),    // Purple
                new Color(128, 128, 0),    // Olive
                new Color(255, 192, 203)   // Pink
        };
        for (int j = 0; j < fibrilPixels.size(); j++) {
            int x = (int) fibrilPixels.get(j)[0];
            int y = (int) fibrilPixels.get(j)[1];
            int clusterIndex = clusterAssignments[j];


            clusterColorsProcessor.setColor(clusterColors[clusterIndex % 10]);
            clusterColorsProcessor.drawPixel(x, y);
        }
    }


}
