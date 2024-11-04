package CollagenAnalysis;


import de.alsclo.voronoi.Voronoi;
import de.alsclo.voronoi.graph.Graph;
import de.alsclo.voronoi.graph.Point;
import ij.ImagePlus;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.process.*;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;

import static CollagenAnalysis.Constants.SCALE;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Graphics2D;

public  class OverlayManager {

    public static void overlayBinaryAttempts(ImageProcessor ip, ImageProcessor initalAttempt, ImageProcessor finalAttempt){
        int width = initalAttempt.getWidth();
        int height = initalAttempt.getHeight();
        for(int y = 0; y < height; y++){
            for(int x = 0; x < width; x++){
                float initialPixel = initalAttempt.getPixelValue(x, y);
                if(initialPixel == finalAttempt.getPixelValue(x,y)){
                    if(initialPixel == 255){
                        ip.setf(x, y, Color.white.getRGB());
                    }
                    else{
                        ip.setf(x, y, Color.black.getRGB());
                    }

                }
                else{
                    if(initialPixel == 0){
                        ip.setf(x,y, Color.magenta.getRGB());
                    }
                    else{
                        ip.setf(x, y, Color.green.getRGB());
                    }
                }
            }
        }
    }

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


    public static void overlayContourLines(ImageProcessor contour, double[][][] post_all, Color lineColor) {
        int numFibrils = post_all.length;
        int height = contour.getHeight();
        int width = contour.getWidth();

        // Use Graphics2D for smooth lines and anti-aliasing
        Graphics2D g2d = contour.getBufferedImage().createGraphics();
        g2d.setColor(lineColor);
        g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON);

        for (int i = 0; i < numFibrils; i++) {
            // ArrayList to store boundary points
            ArrayList<int[]> contourPoints = new ArrayList<>();

            for (int m = 1; m < height - 1; m++) {
                for (int n = 1; n < width - 1; n++) {
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
                        // This is a boundary pixel; add to contour points
                        contourPoints.add(new int[]{n, m});
                    }
                }
            }

            // Draw lines between adjacent contour points
            for (int j = 0; j < contourPoints.size() - 1; j++) {
                int[] startPoint = contourPoints.get(j);
                int[] endPoint = contourPoints.get(j + 1);

                g2d.drawLine(startPoint[0], startPoint[1], endPoint[0], endPoint[1]);
            }
        }

        g2d.dispose(); // Clean up
        // Create an ImagePlus to show the image
        ImagePlus imgWithContour = new ImagePlus("Image with Contour Lines", contour);
        imgWithContour.show(); // Display the image with the drawn contours
    }

    public static void overlayContourLines2(ImageProcessor contour, double[][][] post_all, int color) {
        int numFibrils = post_all.length;
        int height = contour.getHeight();
        int width = contour.getWidth();

        // Create an overlay to store the contour lines
        Overlay overlay = new Overlay();

        for (int i = 0; i < numFibrils; i++) {
            ArrayList<int[]> contourPoints = new ArrayList<>();

            for (int m = 1; m < height - 1; m++) {
                for (int n = 1; n < width - 1; n++) {
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
                        // Add to contour points if it's a boundary pixel
                        contourPoints.add(new int[]{n, m});
                    }
                }
            }

            if (!contourPoints.isEmpty()) {
                // Create a path for smoother contour lines
                GeneralPath path = new GeneralPath();
                int[] firstPoint = contourPoints.get(0);
                path.moveTo(firstPoint[0], firstPoint[1]);

                // Iterate over points and create a smooth path
                for (int j = 1; j < contourPoints.size(); j++) {
                    int[] point = contourPoints.get(j);
                    path.lineTo(point[0], point[1]);
                }

                // Optionally, close the path if it's a closed contour
                // path.closePath();

                // Add the path to the overlay as a shape
                ij.gui.ShapeRoi shapeRoi = new ij.gui.ShapeRoi(path);
                shapeRoi.setStrokeColor(new Color(color));
                shapeRoi.setStrokeWidth(1); // Set the stroke width for thin lines
                overlay.add(shapeRoi);
            }
        }

        // Apply the overlay to the image processor's ImagePlus and show the result
        ImagePlus imageWithOverlay = new ImagePlus("Contour Lines", contour);
        imageWithOverlay.setOverlay(overlay);
        imageWithOverlay.show(); // Display the image
    }


    public static void overlayVoronoi(ImageProcessor ip, Voronoi voronoi, int color){
        ip.setColor(color);
        Graph g = voronoi.getGraph();

        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            de.alsclo.voronoi.graph.Point a = e.getA().getLocation();
            Point b = e.getB().getLocation() ;
            ip.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });
    }

    public static void overlayExcludedRegions(ImageProcessor ip, boolean[][] exclusionMask){
        int width = ip.getWidth();
        int height = ip.getHeight();

//        if (exclusionMask.length != height || exclusionMask[0].length != width) {
//            throw new IllegalArgumentException("Exclusion mask size does not match image dimensions.");
//        }

        // Iterate over all pixels in the image
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Check if the pixel is in the excluded region
                if (exclusionMask[x][y]) {
                    // Set the pixel to black (RGB value: 0)
                    ip.putPixel(x, y, 0);  // Black pixel (for grayscale or RGB images)
                }
            }
        }
    }

    public static void overlayRedBinarization(ColorProcessor ip, ImageProcessor bp){
        for(int y = 0; y < ip.getHeight(); y++){
            for(int x = 0; x < ip.getWidth(); x++){
                double bpBValue = bp.getPixelValue(x,y);
                if(bpBValue == 0){ //fibril pixel

                    int[] rgbArray2 = new int[3];
                    ip.getPixel(x, y, rgbArray2);
                    int red2 = rgbArray2[0];
                    int green2 = rgbArray2[1];
                    int blue2 = rgbArray2[2];

                    int rgb = ip.getPixel(x, y);
                        // Extract the Red, Green, and Blue components using bitwise shifts and masks
                    int red = (rgb >> 16) & 0xFF;
                    int green = (rgb >> 8) & 0xFF;
                    int blue = rgb & 0xFF;

                    red = (red + 255) / 2;
//                    green = 0;
//                    blue = 0;
                    int[] rgbArray = new int[] {red, green, blue};
                    ip.putPixel(x, y, rgbArray);

                }
            }
        }

    }

    public static void overlayCentroids(ImageProcessor ip, double[][] centroids, int color){

        // Draw each maximum point as a red circle on the RGB image
        for (double[] point : centroids) {
            //have to scale centroids up to overlay on the full image because the original image is scaled down for faster processing
            double x = point[0] * SCALE;
            double y = point[1] * SCALE;

            PointDrawer.drawCross(ip, x, y, color);
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
