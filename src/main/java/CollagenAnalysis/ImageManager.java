package CollagenAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.Measurements;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.text.TextWindow;
import org.scijava.tool.Tool;

import java.awt.*;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import static CollagenAnalysis.Constants.*;
import static CollagenAnalysis.ImageProcessingUtils.calculateDistance;

public class ImageManager {
    private int scale = 1;
    private ImageProcessor ip;
    private ImageProcessor scaledUpOriginalProcessor;
    private int height;
    private int width;
    private double  conversionFactor;
    public boolean[][] exclusionMask;

    public ImageManager(ImagePlus imp){
        this.scaledUpOriginalProcessor = imp.getProcessor();
        this.ip = imp.getProcessor();
        this.height = ip.getHeight();
        this.width = ip.getWidth();
    }
    public void scaleDown(){
        while(height >= 1000 || width >= 1000){
            width /= 2;
            height /= 2;
            ip = ip.resize(width, height);
            scale *= 2;
        }
        //assign scale constant for the rest of plugin
        SCALE = scale;
    }
    public void scaleUp(ImagePlus imp){
        ImageProcessor imageProcessor = imp.getProcessor().duplicate();
        int h = imageProcessor.getHeight() * scale;
        int w = imageProcessor.getWidth() * scale;
        imp.setProcessor(imageProcessor.resize(w, h));
    }
    public void showScaledUp(ImagePlus imp){
        if(DEBUG_MODE){
            ImageProcessor imageProcessor = imp.getProcessor().duplicate();
            int h = imageProcessor.getHeight() * scale;
            int w = imageProcessor.getWidth() * scale;
            imageProcessor = imageProcessor.resize(w, h);
            ImagePlus newImage = new ImagePlus(imp.getTitle(), imageProcessor);
            newImage.show();
        }
    }
    public void showScaledUp(ImagePlus imp, boolean override){

        ImageProcessor imageProcessor = imp.getProcessor().duplicate();
        int h = imageProcessor.getHeight() * scale;
        int w = imageProcessor.getWidth() * scale;
        imageProcessor = imageProcessor.resize(w, h);
        ImagePlus newImage = new ImagePlus(imp.getTitle(), imageProcessor);
        newImage.show();
    }
    public int getScale(){
        return scale;
    }
    public ImageProcessor getProcessor(){
        return ip;
    }

    public void DrawScaleBar(ImagePlus imp){
        IJ.setTool("line");
        while(true){
            WaitForUserDialog wfud = new WaitForUserDialog("Action Required", "Draw a line for the scale bar then click ok");
            wfud.show();

            Roi roi = imp.getRoi();
            if (roi == null || !(roi instanceof Line)) {
                IJ.showMessage("Error", "Please draw a line ROI and try again.");
            }
            else{
                break;
            }
        }

        Line line = (Line) imp.getRoi();
        System.out.println("Line Coordinates: X1 = " + line.x1d + ", Y1 = " + line.y1d + ", X2 = " + line.x2d + ", Y2 = " + line.y2d);
        double lengthInPixels = calculateDistance(line.x1d, line.x2d, line.y1d, line.y2d);

        // Ask the user to enter the real-world length in nanometers
        GenericDialog gd = new GenericDialog("Enter Scale Bar Length");
        gd.addNumericField("Length (nm):", 100, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        double realWorldLengthInNm = gd.getNextNumber();
        if (Double.isNaN(realWorldLengthInNm) || realWorldLengthInNm <= 0) {
            IJ.showMessage("Error", "Invalid number entered.");
            return;
        }

        conversionFactor = realWorldLengthInNm / lengthInPixels;

    }
    public ImageProcessor excludeRegions(){
        if (!(ip instanceof ByteProcessor)) {
            ip = ip.convertToByteProcessor();
        }
        //ImageProcessor t = ip.get
        //Todo: ix returning this value it is not needed
        ImageProcessor processorWithExcludedRegions = ip.duplicate();

        ImageProcessor ipExt = scaledUpOriginalProcessor.createProcessor(scaledUpOriginalProcessor.getWidth()  + 2 * borderSize, scaledUpOriginalProcessor.getHeight()  + 2 * borderSize);
        ipExt.setColor(255);
        ipExt.fill();
        ipExt.insert(scaledUpOriginalProcessor, borderSize, borderSize);

        // Display extended image and set up key events
        ImagePlus impExt = new ImagePlus("Extended Image", ipExt);
        RoiManager roiManager = new RoiManager();
        impExt.show();
        ImageCanvas canvas = impExt.getCanvas();
        canvas.addKeyListener(new KeyAdapter() {
            @Override
            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_A) {
                    addPolygonToRoiManager(impExt, roiManager);
                }
            }
        });
        IJ.setTool("polygon");
        IJ.log("Excluding regions instructions: " + excludeRegionsInstructions);
        WaitForUserDialog wfud = new WaitForUserDialog("Action Required", excludeRegionsInstructions);
        wfud.show();

        //generate a mask where coord = true if its in an excluded region
        //scaled down image exclusion mask
        exclusionMask = new boolean[ip.getWidth()][ip.getHeight()];
        // Get ROIs and create the exclusion mask
        Roi[] rois = roiManager.getRoisAsArray();
        for (Roi roi : rois) {
            // Shift the ROI coordinates back to the original image coordinate system
            roi = adjustAndScaleDownRoi(roi, borderSize, scale);
            for (int y = 0; y < ip.getHeight(); y++) {
                for (int x = 0; x < ip.getWidth(); x++) {
                    if (roi.contains(x,y)) {
                        exclusionMask[x][y] = true;
                    }
                }
            }
        }

        //Finalize editing
        impExt.changes = false;
        impExt.close();
        IJ.setTool("hand");

        return processorWithExcludedRegions;
    }
    // Helper method to adjust for the border and scale down the ROI since the exclusion mask is for the scaled down image and
    //the user is excluding regions on the full image
    private Roi adjustAndScaleDownRoi(Roi roi, int borderSize, double scale) {
        Polygon polygon = roi.getPolygon();
        int[] xPoints = polygon.xpoints;
        int[] yPoints = polygon.ypoints;

        int[] adjustedXPoints = new int[xPoints.length];
        int[] adjustedYPoints = new int[yPoints.length];

        // Adjust each point in the ROI to account for the border and scale down
        for (int i = 0; i < xPoints.length; i++) {
            // Subtract the border size first, then scale down the coordinates
            adjustedXPoints[i] = (int) ((xPoints[i] - borderSize) / scale);
            adjustedYPoints[i] = (int) ((yPoints[i] - borderSize) / scale);
        }

        // Return a new ROI with the adjusted and scaled-down points
        return new PolygonRoi(adjustedXPoints, adjustedYPoints, xPoints.length, Roi.POLYGON);
    }

    // Add the polygon to the RoiManager when 'A' is pressed
    private void addPolygonToRoiManager(ImagePlus imp, RoiManager roiManager) {
        Roi roi = imp.getRoi();
        ImageProcessor ip = imp.getProcessor();
        if (roi != null && roi.getType() == Roi.POLYGON) {
            roiManager.addRoi(roi);
            IJ.log("Region excluded");
            ip.setColor(255);
            ip.fill(roi);
            imp.updateAndDraw();
            imp.killRoi(); // Remove the drawn polygon after adding it
        } else {
            IJ.log("No valid polygon to add.");
        }
    }
    public double getConversionFactor() {
        return conversionFactor;
    }
}
