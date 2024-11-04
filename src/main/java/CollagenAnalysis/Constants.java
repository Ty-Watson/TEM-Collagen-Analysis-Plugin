package CollagenAnalysis;

import java.awt.*;

public  class Constants {
    //plugin constants
    public static boolean DEBUG_MODE = false;

    //EM constants
    public static int TOTAL_ITERATIONS = 20;
    public static int blockSize = 100;
    public static int maxCentroidThreshold = 450;


    //Image constants
    public static  int SCALE;
    public static int borderSize = 50;
    public static int pointColor = Color.red.getRGB();


    //csv constants
    public static String delimiter = ",";
    public static String newLine = "\n";
    // Header for Results_in_nm.csv
    public static String[] nmCsvHeader = new String[]{
            "True Area (nanometers^2)",
            "Ellipse Area (nanometers^2)",
            "Major Radius (nanometers)",
            "Minor Radius (nanometers)",
            "Aspect Ratio"
    };

    // Header for Results_in_pixels.csv
    public static String[] pixelsCsvHeader = new String[]{
            "True Area (pixels^2)",
            "Ellipse Area (pixels^2)",
            "Major Radius (pixels)",
            "Minor Radius (pixels)",
            "Angle (degrees)",
            "Ellipse Centroid X (pixels)",
            "Ellipse Centroid Y (pixels)"
    };

    // Header for Results_GaussianMixture.csv
    public static String[] gaussianMixtureCsvHeader = new String[]{
            "Mean X (pixels)",
            "Mean Y (pixels)",
            "Covariance XX (pixels^2)",
            "Covariance XY (pixels^2)",
            "Covariance YY (pixels^2)",
            "Component Proportion"
    };



    //GMM constants
    public static int thin_inc = 10;

    //Dialog constants
    public static String excludeRegionsInstructions = "Step 1: Click in the Extended Image.\n" +
            "Step 2: Draw a polygon then press the 'a' key to exclude the region from processing.\n" +
            "Step 3: Click ok once you are finished excluding regions";
}
