package CollagenAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Line;
import ij.gui.Roi;
import ij.gui.WaitForUserDialog;
import ij.measure.Measurements;
import ij.plugin.frame.RoiManager;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.*;

import static CollagenAnalysis.ImageProcessingUtils.calculateDistance;

public class ImageManager {
    private int scale = 1;
    private ImageProcessor ip;
    private int height;
    private int width;

    private double  conversionFactor;

    public ImageManager(ImageProcessor ip){
        this.ip = ip;
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
    }
    public int getScale(){
        return scale;
    }
    public ImageProcessor getProcessor(){
        return ip;
    }

    public void DrawScaleBar(ImagePlus imp){
        //TODO make a loop if fail first time
        WaitForUserDialog wfud = new WaitForUserDialog("Action Required", "Draw a line and click OK.");
        wfud.show();

        Roi roi = imp.getRoi();
        if (roi == null || !(roi instanceof Line)) {
            IJ.showMessage("Error", "Please draw a line ROI and try again.");
            return;
        }
        Line line = (Line) roi;
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
        if (!(ip instanceof ByteProcessor)) {
            ip = ip.convertToByteProcessor();
        }
        //ImageProcessor t = ip.get
        ImageProcessor processorWithExcludedRegions = ip.duplicate();
        int borderSize = 50;
        ImageProcessor ipExt = ip.createProcessor(ip.getWidth() + 2 * borderSize, ip.getHeight() + 2 * borderSize);
        ipExt.setColor(255); // White for 8-bit, adjust for other types
        ipExt.fill();
        ipExt.insert(ip, borderSize, borderSize);

        // Display extended image
        ImagePlus impExt = new ImagePlus("Extended Image", ipExt);
        impExt.show();

        // Instructions for the user
        GenericDialog gd = new GenericDialog("Interactive Image Editing");
        gd.addMessage("Draw polygons to exclude regions.\nDouble-click to finish.");
        gd.showDialog();
//        if (gd.wasCanceled()) {
//            return;
//        }

        // User interaction for drawing polygons and excluding regions
        boolean isEmptyPolygon = false;
        RoiManager roiManager = new RoiManager();


//        while (!isEmptyPolygon) {
//
////            WaitForUserDialog wfud = new WaitForUserDialog("Action Required", "Draw a polygon and click OK.");
////            wfud.show();
////
////            if (wfud.escPressed()) break; // Exit on escape
//
//            Roi roi = impExt.getRoi();
//            roi.getPolygon();
//            if (roi != null && roi.getType() == Roi.POLYGON) {
//                roiManager.addRoi(roi);
//                ImageStatistics stats = impExt.getStatistics(Measurements.AREA);
//                if (stats.area <= 0) {
//                    isEmptyPolygon = true;
//                } else {
//                    ipExt.setColor(255); // Set excluded regions to NaN, adjust for image type
//                    ipExt.fill(roi);
//                    impExt.updateAndDraw();
//                }
//                impExt.killRoi(); // Remove the drawn polygon
//            }
//        }
        ImageProcessor ipCropped = ipExt.duplicate(); // Duplicate to preserve the extended image
        ipCropped.setRoi(borderSize, borderSize, ip.getWidth(), ip.getHeight());
        ImageProcessor croppedImg = ipCropped.crop();

        for (int x = 0; x < croppedImg.getWidth(); x++) {
            for (int y = 0; y < croppedImg.getHeight(); y++) {
                if (croppedImg.getPixel(x, y) == 255) { // Check if the pixel is part of an excluded region
                    processorWithExcludedRegions.putPixel(x, y, 255); //set to white in original image
                }
            }
        }
        //Finalize editing
        impExt.changes = false;
        impExt.close();
        return processorWithExcludedRegions;
    }

    public double getConversionFactor() {
        return conversionFactor;
    }
}