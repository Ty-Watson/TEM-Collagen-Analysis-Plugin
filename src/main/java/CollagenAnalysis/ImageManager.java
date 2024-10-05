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

import static CollagenAnalysis.ImageProcessingUtils.calculateDistance;

public class ImageManager {
    private int scale = 1;
    private ImageProcessor ip;
    private int height;
    private int width;
    private double  conversionFactor;

    public boolean[][] exclusionMask;
    private boolean exitLoop = false;

    public ImageManager(ImagePlus imp){
        this.ip = imp.getProcessor();
        this.height = ip.getHeight();
        this.width = ip.getWidth();
        checkProcessor(imp);
    }
    public void scaleDown(){
        while(height >= 1000 || width >= 1000){
            width /= 2;
            height /= 2;
            ip = ip.resize(width, height);
            scale *= 2;
        }
    }
    public void scaleUp(ImagePlus imp){
        ImageProcessor imageProcessor = imp.getProcessor().duplicate();
        //drawExcludedRegions(imageProcessor);
        int h = imageProcessor.getHeight() * scale;
        int w = imageProcessor.getWidth() * scale;
        imp.setProcessor(imageProcessor.resize(w, h));
    }
    public void showScaledUp(ImagePlus imp){
        ImageProcessor imageProcessor = imp.getProcessor().duplicate();
        //drawExcludedRegions(imageProcessor);
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
        //TODO make a loop if fail first time
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
        //double lengthInPixels = line.getLength();
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
        //GlobalKeyListener.addGlobalKeyListener();
        if (!(ip instanceof ByteProcessor)) {
            ip = ip.convertToByteProcessor();
        }
        //ImageProcessor t = ip.get
        ImageProcessor processorWithExcludedRegions = ip.duplicate();
        int borderSize = 50;
        ImageProcessor ipExt = ip.createProcessor(ip.getWidth() + 2 * borderSize, ip.getHeight() + 2 * borderSize);
        ipExt.setColor(255);
        ipExt.fill();
        ipExt.insert(ip, borderSize, borderSize);

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
                } else if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                    exitLoop = true;
                    IJ.log("Exiting exclusion loop...");
                }
            }
        });

        String s = "Draw a polygon then press the 'a' key to exclude the region from processing. Press 'Escape' key to finish and exit";
        //IJ.log("Draw a polygon then press the 'a' key to exclude the region from processing. Press 'Escape' key to finish and exit");
        TextWindow textWindow = new TextWindow("Action Required", s, 800, 100);
        // User interaction for drawing polygons and excluding regions
        IJ.setTool("polygon");
        // Loop until Escape is pressed ot exit out of excluding regions
        while (!exitLoop) {
            IJ.setTool("polygon");
            IJ.wait(1000); // Small delay to prevent the loop from consuming too much CPU

//            Roi roi = impExt.getRoi();
//            if (roi != null) {
//                // User has drawn a polygon, wait for input (handled by key events)
//                IJ.log("Polygon drawn, press 'A' to add or 'Escape' to finish.");
//            }
        }

        //generate a mask where coord = true if its in an excluded region
        exclusionMask = new boolean[ip.getWidth()][ip.getHeight()];
        // Get ROIs and create the exclusion mask
        Roi[] rois = roiManager.getRoisAsArray();
        for (Roi roi : rois) {
            // Shift the ROI coordinates back to the original image coordinate system
            roi.setLocation(roi.getBounds().x - borderSize, roi.getBounds().y - borderSize);
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
    public void drawExcludedRegions(ImageProcessor ip){
        for (int y = 0; y < ip.getHeight(); y++) {
            for (int x = 0; x < ip.getWidth(); x++) {
                if (exclusionMask[x][y]) {
                    ip.putPixelValue(x, y, 0);
                }
            }
        }
    }
    private void checkProcessor(ImagePlus imp){
//        if(!(ip instanceof ByteProcessor)){
//            ip = ip.convertToByte(true);
//            imp.setProcessor(ip);
//        }
    }

    public double getConversionFactor() {
        return conversionFactor;
    }
}
