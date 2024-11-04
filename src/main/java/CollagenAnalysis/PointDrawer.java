package CollagenAnalysis;

import ij.process.ImageProcessor;

import static CollagenAnalysis.Constants.*;

public class PointDrawer {
    // Base resolution where the original values look good
    private static final int BASE_WIDTH = 1020;
    private static final int BASE_HEIGHT = 750;

    // Base values for crossSize and thickness
    private static final int BASE_CROSS_SIZE = 6;
    private static final int BASE_THICKNESS = 4;

    public static int crossSize;
    public static int thickness;
    // Method to adjust sizes based on the current image resolution
    public static void adjustSizes(int imageWidth, int imageHeight) {
        // Calculate scaling factors based on the image resolution
        double widthScale = (double) imageWidth / BASE_WIDTH;
        double heightScale = (double) imageHeight / BASE_HEIGHT;

        // Use the average of the scaling factors for uniformity
        double scaleFactor = (widthScale + heightScale) / 2;

        // Adjust crossSize and thickness
        crossSize = (int) Math.round(BASE_CROSS_SIZE * scaleFactor);
        thickness = (int) Math.round(BASE_THICKNESS * scaleFactor);
    }
    public static void drawCross(ImageProcessor ip, double x, double y) {
        // Draw vertical line with thickness
        for (int dx = -thickness / 2; dx <= thickness / 2; dx++) {
            for (int dy = -crossSize; dy <= crossSize; dy++) {
                double newX = x + dx;
                double newY = y + dy;
                if (newX >= 0 && newX < ip.getWidth() && newY >= 0 && newY < ip.getHeight()) {
                    ip.set((int)newX, (int)newY, pointColor);
                }
            }
        }

        // Draw horizontal line with thickness
        for (int dy = -thickness / 2; dy <= thickness / 2; dy++) {
            for (int dx = -crossSize; dx <= crossSize; dx++) {
                double newX = x + dx;
                double newY = y + dy;
                if (newX >= 0 && newX < ip.getWidth() && newY >= 0 && newY < ip.getHeight()) {
                    ip.set((int)newX, (int)newY, pointColor);
                }
            }
        }
    }
    public static void drawCross(ImageProcessor ip, double x, double y, int color) {
        // Draw vertical line with thickness
        for (int dx = -thickness / 2; dx <= thickness / 2; dx++) {
            for (int dy = -crossSize; dy <= crossSize; dy++) {
                double newX = x + dx;
                double newY = y + dy;
                if (newX >= 0 && newX < ip.getWidth() && newY >= 0 && newY < ip.getHeight()) {
                    ip.set((int)newX, (int)newY, color);
                }
            }
        }

        // Draw horizontal line with thickness
        for (int dy = -thickness / 2; dy <= thickness / 2; dy++) {
            for (int dx = -crossSize; dx <= crossSize; dx++) {
                double newX = x + dx;
                double newY = y + dy;
                if (newX >= 0 && newX < ip.getWidth() && newY >= 0 && newY < ip.getHeight()) {
                    ip.set((int)newX, (int)newY, color);
                }
            }
        }
    }





}
