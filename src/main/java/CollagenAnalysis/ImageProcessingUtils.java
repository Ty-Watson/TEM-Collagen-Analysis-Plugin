package CollagenAnalysis;


import Jama.LUDecomposition;
import Jama.Matrix;
import ij.ImagePlus;
import ij.process.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;

public class ImageProcessingUtils {

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
            //areaPix2[i] *= Math.pow(scale, 2);
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

        // Index for the intensities array
        int index = 0;

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
    public static double calculateDistance(double x1, double x2, double y1, double y2) {
        return Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));
    }

}

