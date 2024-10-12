package CollagenAnalysis;

public  class Constants {
    //plugin constants
    public static boolean DEBUG_MODE = false;

    //Image constants
    public static  int SCALE;
    public static int borderSize = 50;
    public static int pointSize = 2; // This can be adjusted based on how big you want the maxima points to be
    //csv constants
    public static String delimiter = ",";
    public static String newLine = "\n";
    public static String[] pixelCsvHeader = new String[] {"Area (pixels^2)", "Major Radius (pixels)", "Minor Radius (pixels)",
            "Angle (degrees)", "Centroid X (pixels)", "Centroid Y (pixels)",
            "Covariance XX (pixels^2)", "Covariance XY (pixels^2)", "Covariance YY (pixels^2)",
            "Component Proportion"};
    public static String[] nmCsvHeader = new String[]{"Area (nanometers^2)", "Major Radius (nanometers)", "Minor Radius (nanometers)",
            "Aspect Ratio"};

    //GMM constants
    public static int thin_inc = 10;

    //Dialog constants
    public static String excludeRegionsInstructions = "Step 1: Click in the Extended Image.\n" +
            "Step 2: Draw a polygon then press the 'a' key to exclude the region from processing.\n" +
            "Step 3: Click ok once you are finished excluding regions";
}
