package CollagenAnalysis;


import Jama.LUDecomposition;
import Jama.Matrix;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.plugin.filter.EDM;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.frame.RoiManager;
import ij.process.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import smile.classification.KNN;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;

public class ImageProcessingUtils {

    public static int findNearestEllipse(int clickX, int clickY, ArrayList<double[]> centroids) {
        double minDistance = Double.MAX_VALUE;
        int nearestEllipseIndex = -1;

        for (int i = 0; i < centroids.size(); i++) {
            double[] centroid = centroids.get(i);
            double distance = Math.sqrt(Math.pow(clickX - centroid[0], 2) + Math.pow(clickY - centroid[1], 2));
            if (distance < minDistance) {
                minDistance = distance;
                nearestEllipseIndex = i;
            }
        }

        return nearestEllipseIndex;
    }
    public static double[] findNearestCentroid(ArrayList<double[]> points, int x, int y) {
        double[] nearest = null;
        double minDistance = Double.MAX_VALUE;
        for (double[] point : points) {
            double distance = Math.sqrt(Math.pow(point[0] - x, 2) + Math.pow(point[1] - y, 2));
            if (distance < minDistance) {
                minDistance = distance;
                nearest = point;
            }
        }
        return nearest;
    }
    public static ArrayList<double[]> convertToArrayList(double[][] doubleArray) {
        ArrayList<double[]> arrayList = new ArrayList<>();
        for (double[] arr : doubleArray) {
            arrayList.add(arr);
        }
        return arrayList;
    }
    public static double[] computeFibrilArea(double[][] postFibril, int scale, int fibrilCount) {
        double[] areaPix2 = new double[fibrilCount];

        for (int i = 0; i < fibrilCount; i++) {
            for(int j = 0; j < postFibril[i].length; j++) {
                areaPix2[i] += postFibril[i][j];
            }
            areaPix2[i] *= Math.pow(scale, 2);
        }
        return areaPix2;
    }
    public static int[] assignEachFibrilPixelToFibril(int fibrilPixelsLength, int fibrilCentroidsLength, double[][] post_fibril, ArrayList<Boolean> isBoundryCentroid){
        int[] clusterAssignments = new int[fibrilPixelsLength];

        // Iterate over each fibril pixel to find the centroid with the maximum posterior probability
        for (int j = 0; j < fibrilPixelsLength; j++) {
            double maxProbability = -1; // Initialize with a value that will be lower than any probability
            int maxIndex = -1; // Initialize with a value that cannot be an index

            // Iterate over each centroid to find the max probability for the current pixel
            for (int i = 0; i < fibrilCentroidsLength; i++) {
                // Check if the current probability is greater than the max found so far

                if (post_fibril[i][j] > maxProbability) {
                    maxProbability = post_fibril[i][j]; // Update max probability
                    maxIndex = i; // Update index of the centroid with the max probability
                }

            }

            // Assign the pixel to the centroid with the highest probability
            clusterAssignments[j] = maxIndex;
        }
        return clusterAssignments;
    }
    public static ByteProcessor applyEuclideanDistanceMap(ByteProcessor bp) {
        EDM edm = new EDM();
        edm.toEDM(bp);
        return bp;
    }
    public static void drawBoundries(ImageProcessor ip, double[][] p){
        ImageProcessor copy = ip.duplicate();  // Create a copy of the original image.
        ImagePlus imp = new ImagePlus("Boundries");
        int width = ip.getWidth();
        int height = ip.getHeight();

        int numberOfFibrils = p.length;  // Assuming the number of fibrils is the length of the first dimension of p

        // Iterate through each pixel in the image excluding the border pixels
        for (int m = 1; m < height - 1; m++) {
            for (int n = 1; n < width - 1; n++) {
                int index = m * width + n; // Calculate the linear index of the pixel

                for (int i = 0; i < numberOfFibrils; i++) {  // Iterate over each fibril
                    double minProb = Double.MAX_VALUE;
                    double maxProb = Double.MIN_VALUE;

                    // Check the 3x3 neighborhood around the pixel
                    for (int dy = -1; dy <= 1; dy++) {
                        for (int dx = -1; dx <= 1; dx++) {
                            int neighborIndex = (m + dy) * width + (n + dx); // Linear index of the neighbor
                            double prob = p[i][neighborIndex];
                            if (prob < minProb) minProb = prob;
                            if (prob > maxProb) maxProb = prob;
                        }
                    }

                    // Change pixel color if it's a boundary
                    if (minProb < 0.5 && maxProb > 0.5) {
                        copy.putPixel(n, m, Color.BLUE.getRGB());  // Set boundary pixel to blue
                        break;  // Stop checking other fibrils once boundary is found for this pixel
                    }
                }
            }
        }

        // Update the ImagePlus object with the modified copy and refresh the display
        imp.setProcessor(copy);
        imp.updateAndDraw();
        imp.show();
    }
//    public  double[] computeIntensities(List<double[][]> fibrilPixels, double[] distances, int[] indices, int numCentroids) {
//        Map<Integer, List<Double>> centroidDistances = new HashMap<>();
//        Map<Integer, List<Double>> centroidIntensities = new HashMap<>();
//
//        // Collect distances and intensities for each centroid
//        for (int i = 0; i < indices.length; i++) {
//            centroidDistances.computeIfAbsent(indices[i], k -> new ArrayList<>()).add(distances[i]);
//            centroidIntensities.computeIfAbsent(indices[i], k -> new ArrayList<>()).add(fibrilPixels.get(i));
//        }
//
//        double[] charCellIntensity = new double[numCentroids];
//        // Calculate the mean intensity for pixels within the 5th percentile distance for each centroid
//        centroidDistances.forEach((centroid, distList) -> {
//            double fifthPercentile = percentile(distList, 5);
//            List<Double> intensities = centroidIntensities.get(centroid);
//            double sum = 0;
//            int count = 0;
//            for (int i = 0; i < distList.size(); i++) {
//                if (distList.get(i) < fifthPercentile) {
//                    sum += intensities.get(i);
//                    count++;
//                }
//            }
//            charCellIntensity[centroid] = (count > 0) ? sum / count : Double.NaN; // Compute mean or NaN if no pixels qualify
//        });
//
//        return charCellIntensity;
//    }

    //    private double percentile(List<Double> values, double percentile) {
//        Collections.sort(values);
//        int index = (int) Math.ceil(percentile / 100.0 * values.size()) - 1;
//        return values.get(Math.max(index, 0));
//    }
//}
    public static void extractContours(ImageProcessor ip, ImagePlus imp) {
        RoiManager roiManager = RoiManager.getRoiManager();
        roiManager.reset();

        // Wand tool to trace outlines
        Wand wand = new Wand(ip);

        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                if (ip.getPixel(x, y) == 255) { // Assuming foreground is white
                    wand.autoOutline(x, y); // Use the wand tool at the starting point (x, y)

                    // Create a PolygonRoi based on the Wand tool's outline
                    Roi roi = new PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.FREEROI);

                    // Add the ROI to the RoiManager
                    roiManager.addRoi(roi);

                    // Optionally, draw the contour on a new image or overlay
                    DrawROI(ip, roi, imp);
                }
            }
        }
    }

    // Method to draw ROI on the ImageProcessor (optional)
    public static void DrawROI(ImageProcessor ip, Roi roi, ImagePlus imp) {
        ip.setColor(Color.BLUE); // Set the color for drawing
        ip.draw(roi); // Draw the ROI on the image
        imp.updateAndDraw(); // Update the ImagePlus object to show the drawn ROI
    }

    public static ArrayList<int[]> thinData(ArrayList<int[]> originalList, int thinInc) {
        ArrayList<int[]> thinnedList = new ArrayList<>();

        for (int i = 0; i < originalList.size(); i += thinInc) {
            thinnedList.add(originalList.get(i));
        }

        return thinnedList;
    }
    public static ArrayList<double[]> getFibrilPixels(FloatProcessor bp) {
        ArrayList<double[]> cords = new ArrayList<>();
        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                if(bp.getPixel(x, y) == 0){
                    double[] c = new double[2];
                    c[0] = x;
                    c[1] = y;
                    cords.add(c);
                }
            }
        }
        return cords;
    }
    public static ArrayList<double[]> getAllPixels(ImageProcessor ip) {
        ArrayList<double[]> cords = new ArrayList<>();
        //double[][][] allPixels = new double[ip.getHeight()][ip.getWidth()][2];
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                double[] c = new double[2];
                c[0] = x;
                c[1] = y;
                cords.add(c);
                //allPixels[x][y] = c;
            }
        }
        return cords;
    }

    public static  int[] findNearestCentroidForFibrilPixel(ArrayList<double[]> centroids, ArrayList<double[]> fibrilPixels){
        int[] centroidLabels = new int[centroids.size()]; // Labels for each centroid
        for (int i = 0; i < centroidLabels.length; i++) {
            centroidLabels[i] = i; // Assign a unique label to each centroid
        }
        double[][] c = new double[centroids.size()][2];
        for(int i = 0; i < centroids.size(); i++){
            c[i][0] = centroids.get(i)[0];
            c[i][1] = centroids.get(i)[1];
        }
        KNN<double[]> knn = KNN.fit(c, centroidLabels, 1);
        int[] nearestCentroids = new int[fibrilPixels.size()];

        for (int i = 0; i < fibrilPixels.size(); i++) {
            nearestCentroids[i] = knn.predict(fibrilPixels.get(i));
        }
        return nearestCentroids;
    }


    public static ImagePlus smoothIntensityValues(ImageProcessor ip){
        double[][] thresholds = performRegression(ip);
        double[] th = new double[6];
        for(int i = 0; i < thresholds.length; i++){
            th[i] =thresholds[i][0];
        }

        ImageProcessor test = ip.duplicate();

        double sum = 0;
        double count = 0;
        //adjusting the filtered image intensities based on the local threshold values predicted by the bivariate quadratic model
        Integer y1 = ip.getHeight();
        Integer x1 = ip.getWidth();
        int index = 0;
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {

                double originalIntensity = ip.getPixelValue(x, y);
                sum += originalIntensity;
                count++;
                double value = f(th, x,y);
                test.putPixelValue(x,y, value);
                // Apply the adjustment formula, ensuring localThreshold is not 0 to avoid division by zero
                if (value != 0) {
                    double adjustedIntensity = Math.pow(originalIntensity / 255, Math.log(0.5) / Math.log(value / 255));
                    if(!Double.isNaN(adjustedIntensity))
                        ip.putPixelValue(x, y, adjustedIntensity * 255);
                }
            }
        }

        double avg = sum / count;

//        ImagePlus testimg = new ImagePlus("Test value",test);
//        testimg.show();

        ImagePlus filteredImg = new ImagePlus("Filtered Img",ip);
//        filteredImg.show();
        return filteredImg;
    }

    public static ByteProcessor binarizeUsingOtsu(FloatProcessor fp) {
        // Convert FloatProcessor to ByteProcessor for histogram calculation
        ByteProcessor bp = fp.duplicate().convertToByteProcessor();
       // ByteProcessor bp2 = bp.duplicate().convertToByteProcessor();

        // Calculate histogram
        int[] histogram = bp.getHistogram();

        // Calculate Otsu's threshold
        AutoThresholder thresholder = new AutoThresholder();
        int otsuThreshold = thresholder.getThreshold(AutoThresholder.Method.Otsu, histogram);

        // Apply Otsu's threshold to the FloatProcessor
        for (int y = 0; y <bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                float value =bp.getPixelValue(x, y);
                if(value == 0){
                    bp.setf(x, y, 255.0f);
                    continue;
                }

                if (value <= otsuThreshold) {
                    bp.setf(x, y, 0.0f); // Set pixels below the threshold to 0
                }
                else {
                    bp.setf(x, y, 255.0f); // Set pixels above the threshold to 255
                }
            }
        }
        return bp;
    }
    public static ByteProcessor binarizeUsingOtsu(FloatProcessor fp, boolean[][] exclusionMask) {
        // Convert FloatProcessor to ByteProcessor for histogram calculation
        ByteProcessor bp = fp.duplicate().convertToByteProcessor();

        // Calculate histogram, ignoring excluded regions
        int[] histogram = new int[256];
        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                if (exclusionMask[x][y]) {  // Skip excluded regions (mask value 255 for excluded)
                    continue;
                }
                int pixelValue = bp.get(x, y);
                histogram[pixelValue]++;
            }
        }

        // Calculate Otsu's threshold
        AutoThresholder thresholder = new AutoThresholder();
        int otsuThreshold = thresholder.getThreshold(AutoThresholder.Method.Otsu, histogram);

        // Apply Otsu's threshold, ignoring excluded regions
        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                if (exclusionMask[x][y]) {
                    bp.setf(x, y, 255.0f);// Skip excluded regions
                    continue;
                }

                float value = bp.getPixelValue(x, y);
                if (value == 0) {
                    bp.setf(x, y, 255.0f);
                    continue;
                }

                if (value <= otsuThreshold) {
                    bp.setf(x, y, 0.0f);  // Set pixels below the threshold to 0
                } else {
                    bp.setf(x, y, 255.0f);  // Set pixels above the threshold to 255
                }
            }
        }

        return bp;
    }



    public static double f(double[] coeffs, double u, double v) {
        if (coeffs == null || coeffs.length != 6) {
            throw new IllegalArgumentException("Coefficients array must be non-null and have exactly 6 elements.");
        }

        return coeffs[0] // C00
                + coeffs[1] * u // C10 * u
                + coeffs[2] * v // C01 * v
                + coeffs[3] * u * v // C11 * u * v
                + coeffs[4] * Math.pow(u, 2) // C20 * u^2
                + coeffs[5] * Math.pow(v, 2); // C02 * v^2
    }




    private static double[][] performRegression(ImageProcessor ip){
        double[] intensities = extractIntensities(ip);
        double[][] coordinates = extractCoordinates(ip);

        double[][] X = prepareDesignMatrix(coordinates);
        double[] Y = intensities;

//        RealMatrix m1 = MatrixUtils.createRealMatrix(X);
//        double[][] mt = new double[Y.length][];
//        for(int i = 0; i < Y.length; i++){
//            mt[i][0] = Y[i];
//        }
//        RealMatrix m2 = MatrixUtils.createRealMatrix(mt);
//        RealMatrix m = computePseudoInverse(m1);
//        RealMatrix p = m.multiply(m2);
        Matrix A = new Matrix(X);
        Matrix b = new Matrix(Y, Y.length);
        Matrix at = A.transpose();
        Matrix ata = at.times(A);
        Matrix atb = at.times(b);


        LUDecomposition lu = new LUDecomposition(ata);

        Matrix x = lu.solve(atb);
        double[][] resultingVector = x.getArrayCopy();

        return resultingVector;

    }

    private static double[] extractIntensities(ImageProcessor ip){
        // Calculate the total number of pixels
        int totalPixels = ip.getWidth() * ip.getHeight();

        // Initialize an array to hold the intensity values
        double[] intensities = new double[totalPixels];
        double[] intensities2 = new double[totalPixels];

        // Index for the intensities array
        int index = 0;
        int index2 = 0;

        /*byte[] pixels = (byte[])ip.getPixels();
        int w = ip.getWidth(), h = ip.getHeight();
        for (int j = 0; j < h; j++)
            for (int i = 0; i < w; i++) {
                // Java has no unsigned 8-bit data type, so we need to perform Boolean arithmetics
                int value = pixels[i + w * j] & 0xff;
                intensities2[index2++] = value;

            }*/
        // Iterate over all pixels in the image
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                // Get the intensity value of the current pixel
                //double value = ip.getPixelValue(x, y);
                float value = ip.getf(x,y);


                // Store the intensity value in the array
                intensities[index++] = (double)value;
            }
        }

        // Return the array containing all the intensity values
        return intensities;
    }

    private static double[][] extractCoordinates(ImageProcessor ip){


        int totalPixels = ip.getWidth() * ip.getHeight();

        // Initialize the array to hold the coordinates; each row has 2 columns for x and y coordinates
        double[][] coordinates = new double[totalPixels][2];


        int index = 0;


        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                // Store the current pixel's coordinates in the array
                coordinates[index][0] = x; // x-coordinate
                coordinates[index][1] = y; // y-coordinate
                index++;
            }
        }

        // Return the array of coordinates
        return coordinates;
    }

    private static double[][] prepareDesignMatrix(double[][] coordinates){
        int numRows = coordinates.length;
        // 6 columns for the intercept, x, y, x^2, xy, and y^2 terms
        double[][] X = new double[numRows][6];

        for (int i = 0; i < numRows; i++) {
            double x = coordinates[i][0];
            double y = coordinates[i][1];
            X[i][0] = 1;       // Intercept term
            X[i][1] = x;       // x term
            X[i][2] = y;       // y term
            X[i][3] = x * y;   // xy term
            X[i][4] = x * x;   // x^2 term
            X[i][5] = y * y;   // y^2 term
        }
        return X;
    }
    //
    public static void overlayBinaryOnGrayscale(ImagePlus originalImage, ImagePlus binaryImage) {
        // Convert the binary image to RGB
        ImagePlus binaryRGB = binaryImage.duplicate();
        binaryRGB.setTitle("Binarized Red");
        IJ.run(binaryRGB, "RGB Color", "");

        // Access the RGB processor
        ColorProcessor cp = (ColorProcessor) binaryRGB.getProcessor();

        // Iterate over all pixels
        int[] rgb = new int[3];
        for (int y = 0; y < cp.getHeight(); y++) {
            for (int x = 0; x < cp.getWidth(); x++) {
                float pixelValue = cp.getf(x, y);
                //cp.getPixel(x, y, rgb);

                // If the pixel in the binary image is foreground, set it to red
                if (pixelValue == 255) { // Check the blue channel for non-zero values
                    rgb[0] = 255; // Red
                    rgb[1] = 0;   // Green
                    rgb[2] = 0;   // Blue
                    cp.putPixel(x, y, rgb);
                } else {
                    // Make background pixels transparent
                    cp.putPixel(x, y, new int[]{0, 0, 0, 0});
                }
            }
        }

        // Create an overlay and add the RGB image
        Overlay overlay = new Overlay();
        overlay.add(new ImageRoi(0, 0, cp));
        originalImage.setOverlay(overlay);

        // Display the original image with the overlay
        originalImage.show();
        originalImage.updateAndDraw();
    }
    public void overlayBinaryOnGrayscale2(ImagePlus originalImage, ImagePlus binaryImage) {
        // Convert the original image to RGB to allow for color overlay
        ImagePlus originalRGB = originalImage.duplicate();
        IJ.run(originalRGB, "RGB Color", "");

        // Access the processors
        ImageProcessor ip = binaryImage.getProcessor().duplicate();
        //FloatProcessor fp = (FloatProcessor) binaryProcessor;
        ColorProcessor originalCP = (ColorProcessor) originalRGB.getProcessor();
        //float[] pixels = (float[]) fp.getPixels();
        // Iterate over all pixels in the binary image
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                float pixelValue = ip.getf(x, y);

                // If the pixel in the binary image is part of the foreground, set the corresponding pixel in the RGB image to red
                if (pixelValue == 255) { // Check for foreground in binary image
                    originalCP.putPixel(x, y, new int[]{181, 132, 132}); // Set to red
                }
                // If the binaryPixel is not 255, the original pixel is retained
            }
        }

        // Create a new ImagePlus to display the result
        ImagePlus overlayedImage = new ImagePlus("Original image overlayed with binarized image", originalCP);
        overlayedImage.show();
    }

    public static ImagePlus applyGaussianFiltering(ImageProcessor ipp){
        ImageProcessor ip = ipp.duplicate();
        // Calculate the sigma for the Gaussian filter
        double smallestDimension = Math.min(ip.getWidth(), ip.getHeight());
        double gsigma = smallestDimension * 0.003; // 0.3% of the smallest dimension

        // Apply Gaussian filter
        ip.blurGaussian(gsigma);

        //show gaus image
       ImagePlus impGaus = new ImagePlus("Gaus Image", ip);
//        impGaus.show();
        return impGaus;
    }

    public static ByteProcessor applyGaussianBlur(ByteProcessor bp, double sigma) {
        ImageProcessor fpFiltered = bp.duplicate().convertToFloatProcessor();
        fpFiltered.blurGaussian(sigma);
        return (ByteProcessor) fpFiltered.convertToByteProcessor();
    }

    public static ByteProcessor findMaxima(ByteProcessor bp, double tolerance) {
        MaximumFinder mf = new MaximumFinder();
        return mf.findMaxima(bp, tolerance, ImageProcessor.NO_THRESHOLD, MaximumFinder.SINGLE_POINTS, false, false);
    }

    public static int countNonZeroPixels(ByteProcessor bp) {
        int count = 0;
        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                if (bp.getPixel(x, y) != 0) count++;
            }
        }
        return count;
    }
    public static double[] derivePolynomial(double[] coefficients) {
        // The derivative of a polynomial of degree n has degree n-1
        double[] derivative = new double[coefficients.length - 1];

        // Compute the derivative coefficients
        for (int i = 1; i < coefficients.length; i++) {
            derivative[i - 1] = i * coefficients[i];
        }

        return derivative;
    }

    public static double[] squarePolynomial(double[] coefficients) {
        // The resulting polynomial degree will be twice that of the input polynomial minus 2
        // (since the derivative reduces the degree by 1, and squaring doubles it).
        double[] squared = new double[2 * coefficients.length - 3];

        // Multiply each coefficient by every other coefficient
        for (int i = 0; i < coefficients.length - 1; i++) {
            for (int j = 0; j < coefficients.length - 1; j++) {
                squared[i + j] += coefficients[i] * coefficients[j];
            }
        }

        return squared;
    }
    private static int[] generateColors(int numClusters) {
        int[] colors = new int[numClusters];
        for (int i = 0; i < numClusters; i++) {
            float hue = (float) i / numClusters; // Spread the hues across the spectrum
            float saturation = 0.6f + 0.4f * (i % 2); // Alternate saturation to increase distinction
            float brightness = 0.5f + 0.5f * ((i % 3) / 2.0f); // Vary brightness
            colors[i] = Color.HSBtoRGB(hue, saturation, brightness);
        }
        return colors;
    }

    public static void displayClusters(ImageProcessor ip, int[] cluster_idx) {
        int[] colors = generateColors(234); // Example colors

        int thinningFactor = 10;
        int width = ip.getWidth();
        int height = ip.getHeight();

        // Assuming clusterIdx corresponds to every 'thinningFactor' pixel linearly in a row-major order
        int thinnedWidth = width / thinningFactor;
        int thinnedHeight = height / thinningFactor;

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Calculate the thinned index
                int thinnedIndexX = x / thinningFactor;
                int thinnedIndexY = y / thinningFactor;

                // Check if the calculated index is within the bounds of the thinned data
                if (thinnedIndexY < thinnedHeight && thinnedIndexX < thinnedWidth) {
                    int index = thinnedIndexY * thinnedWidth + thinnedIndexX;  // Convert 2D index to 1D
                    if (index < cluster_idx.length) {  // Additional safety check
                        int colorIndex = cluster_idx[index];
                        ip.putPixel(x, y, colors[colorIndex % colors.length]);
                    }
                }
            }
        }
        //return ip;
    }
    public static void writeArrayToFile(double[][] data, String filePath) {
        try (ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filePath))) {
            out.writeObject(data);
        } catch (IOException e) {
            System.out.println("An error occurred while writing the array to a file.");
            e.printStackTrace();
        }
    }
    public static void writeArrayToFile(Object data, String filePath) {
        try (ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(filePath))) {
            out.writeObject(data);
        } catch (IOException e) {
            System.out.println("An error occurred while writing the array to a file.");
            e.printStackTrace();
        }
    }
    public static double[][] readArrayFromFile(String filePath) {
        double[][] data = null;
        try (ObjectInputStream in = new ObjectInputStream(new FileInputStream(filePath))) {
            data = (double[][]) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            System.out.println("An error occurred while reading the array from the file.");
            e.printStackTrace();
        }
        return data;
    }
    public static double[][][] readArrayFromFile2(String filePath) {
        double[][][] data = null;
        try (ObjectInputStream in = new ObjectInputStream(new FileInputStream(filePath))) {
            data = (double[][][]) in.readObject();
        } catch (IOException | ClassNotFoundException e) {
            System.out.println("An error occurred while reading the array from the file.");
            e.printStackTrace();
        }
        return data;
    }



    // Calculate the Cauchy weight for a single residual
//    public double cauchyWeight(double residual, double leverage, double[] allResiduals) {
//        double s = calculateRobustS(allResiduals);
//        double standardizedResidual = standardizedResidual(residual, s, leverage);
//        return 1.0 / (1.0 + Math.pow(standardizedResidual, 2));
//    }
//
//    // Standardize a single residual
//    private double standardizedResidual(double residual, double s, double leverage) {
//        return residual / (s * Math.sqrt(1.0 - leverage));
//    }
//
//    // Calculate the robust estimator of standard deviation (s)
//    private double calculateRobustS(double[] residuals) {
//        double medianResidual = median(residuals);
//        double[] absoluteDifferences = Arrays.stream(residuals)
//                .map(r -> Math.abs(r - medianResidual))
//                .toArray();
//        double mad = median(absoluteDifferences); // Median Absolute Deviation
//        return mad / 0.6745; // Approximation for F^{-1}_N(3/4)
//    }
//
//    // Utility method to calculate the median of an array
//    private double median(double[] values) {
//        double[] sortedValues = Arrays.copyOf(values, values.length);
//        Arrays.sort(sortedValues);
//        int middle = sortedValues.length / 2;
//        if (sortedValues.length % 2 == 0) {
//            return (sortedValues[middle - 1] + sortedValues[middle]) / 2.0;
//        } else {
//            return sortedValues[middle];
//        }
//
//    }
    public static double calculateDistance(double x1, double x2, double y1, double y2) {
        return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
    }
//
//    public double[][] computeCovarianceMatrix(double[][] points, double[] centroid) {
//        double varX = 0.0;
//        double varY = 0.0;
//        double covXY = 0.0;
//
//        for (double[] point : points) {
//            double dx = point[0] - centroid[0]; // Difference from centroid x
//            double dy = point[1] - centroid[1]; // Difference from centroid y
//
//            varX += dx * dx;
//            varY += dy * dy;
//            covXY += dx * dy;
//        }
//
//        int n = points.length;
//        varX /= n;
//        varY /= n;
//        covXY /= n;
//
//        return new double[][]{
//                {varX, covXY},
//                {covXY, varY}
//        };
//    }

}

