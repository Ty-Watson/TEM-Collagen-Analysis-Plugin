/*
 * To the extent possible under law, the ImageJ developers have waived
 * all copyright and related or neighboring rights to this tutorial code.
 *
 * See the CC0 1.0 Universal license for details:
 *     http://creativecommons.org/publicdomain/zero/1.0/
 */

package CollagenAnalysis;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.gui.WaitForUserDialog;
import ij.io.FileSaver;
import ij.plugin.filter.PlugInFilter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.text.TextWindow;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVPrinter;
import smile.neighbor.KDTree;
import smile.neighbor.Neighbor;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import static CollagenAnalysis.ImageProcessingUtils.*;
import static java.awt.Color.cyan;

/**

 *
 * @author Ty Watson
 */
public class CollagenAnalysisPlugin implements PlugInFilter {
    protected ImagePlus imp;
    // plugin parameters
    public double value;
    public String name;
    private double[] nearestCentroidDistances;

    @Override
    public int setup(String arg, ImagePlus image) {
        if (arg.equals("about")) {
            //showAbout();
            return DONE;
        }

        imp = image;
        return DOES_8G | DOES_16 | DOES_32 | DOES_RGB;
    }

    @Override
    public void run(ImageProcessor preProcessProcessor) {


        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.showMessage("No image open");
            return;
        }

        //setting up post processing items
        String imageName = IJ.getImage().getTitle();
        ArrayList<ImagePlus> imagesToSave = new ArrayList<>();



        //scale down image to managable size
        ImageManager imageManager = new ImageManager(imp);
        imageManager.DrawScaleBar(imp);
        imageManager.scaleDown();
        ImageProcessor ip = imageManager.getProcessor();
        ImagePlus originalImgScaled = new ImagePlus("OrignalImgScaled", ip);
        //originalImgScaled.show();

        // allow users to exclude regions from processing on image
        ImageProcessor processorWithExcludedRegions = imageManager.excludeRegions();

        //orignal image with excluded regions
        ImageProcessor originalImgWithExcludedRegions = originalImgScaled.getProcessor().duplicate();
        OverlayManager.overlayExcludedRegions(originalImgWithExcludedRegions, imageManager.exclusionMask);
        ImagePlus originalImgWithExcludedRegionsImg = new ImagePlus("Original image with excluded regions", originalImgWithExcludedRegions);
        imageManager.showScaledUp(originalImgWithExcludedRegionsImg);
        imagesToSave.add(originalImgWithExcludedRegionsImg);

        //[*]Filtering, thresholding, and binarization[*]

        //apply gaussian filtering to original img
        ImagePlus imgAfterGausFiltering = applyGaussianFiltering(processorWithExcludedRegions);
        ImageProcessor ip2 = imgAfterGausFiltering.getProcessor().duplicate().convertToFloatProcessor();

        //smooth intensity values for binarization step
        ImagePlus imageAfterSmoothing = smoothIntensityValues(ip2);
        //binarize the smoothed img using Otsu optimal threshold
        FloatProcessor ip3 = imageAfterSmoothing.getProcessor().duplicate().convertToFloatProcessor();
        ByteProcessor binary1 = binarizeUsingOtsu(ip3, imageManager.exclusionMask); //set exluded regions to white
        ImagePlus binaryImg = new ImagePlus("Initial binary attempt",binary1);
        binaryImg.show();
        //imageManager.showScaledUp(binaryImg);

        //after binarization step, we can find centroids
        CentroidManager cm = new CentroidManager(binaryImg.getProcessor().duplicate());
        ArrayList<double[]> centroids = cm.findCentroids();
        ArrayList<double[]> initialCentroidCopy = (ArrayList<double[]>) centroids.clone();
        //set up image for a voronoi diagram
        VoronoiManager vm = new VoronoiManager(originalImgScaled, imageManager.exclusionMask);
        vm.setCentroids(centroids);
        //draw voronoi diagram based on initial guess of the centroids for the fibrils
        vm.initDrawVoronoi();
        vm.showImg();
        //user can add and remove centroids and the voronoi diagram will adjust
        WaitForUserDialog wfud = new WaitForUserDialog("Action Required", "Right click to remove centroids. Left click to add centroids. Click ok when finished.");
        wfud.show();
        if (wfud.escPressed()) return;

        //get non boundry centroids for later processing
        ArrayList<double[]> allPixels = getAllPixels(ip2);
        ArrayList<double[]> nonBoundryCentroids = vm.identifyBoundaryFibrils(imageManager.exclusionMask, allPixels);
        ArrayList<Boolean> isBoundryCentroid = vm.identifyBoundaryFibrilsBoolean(imageManager.exclusionMask, allPixels);

        //completed voronoi diagram with user approved centroids
        double[][] convertedCentroids = vm.getCentroids().toArray(new double[0][]);

        // Create a KD-tree for centroids
        //Todo: need to have a way to delete all of the centroids and mu and cov assosciated if it is aboundry centroid. No need to run Gmm on data we are not going to use
        KDTree<double[]> kdTree = KDTree.of(convertedCentroids);
        double[] nearestCentroidDistancesAll = new double[allPixels.size()];
        int[] nearestCentroidIndicesAll = new int[allPixels.size()];
        // Perform the nearest neighbor search for every pixel to each fibril
        for (int i = 0; i < allPixels.size(); i++) {
            Neighbor<double[], double[]> nearest = kdTree.nearest(allPixels.get(i));
            nearestCentroidIndicesAll[i] = nearest.index;
            nearestCentroidDistancesAll[i] = nearest.distance;
        }

        //set up interpolation to create a smooother image for second binary attempt
        vm.SetNearestCentroidDistances(nearestCentroidDistancesAll);
        vm.SetNearestCentroidIndices(nearestCentroidIndicesAll);
        ImagePlus smooth2 = vm.Interpolation(imageAfterSmoothing.getProcessor(), convertToArrayList(convertedCentroids));
        //smooth2.show();

        //binarize the smoothed img using Otsu optimal threshold
        FloatProcessor ip4 = smooth2.getProcessor().duplicate().convertToFloatProcessor();
        ImageProcessor binary2 = binarizeUsingOtsu(ip4, imageManager.exclusionMask);
        ImagePlus binaryImg2 = new ImagePlus("Final binary attempt", binary2);
        binaryImg2.show();

        //red binary image over original image
        ImageProcessor binary2dup = binary2.duplicate();
        ColorProcessor redBinaryProcessor = originalImgScaled.getProcessor().duplicate().convertToColorProcessor();
        OverlayManager.overlayRedBinarization(redBinaryProcessor, binary2dup);
        ImagePlus redBinaryImg = new ImagePlus("Original image overlayed with final binarized image", redBinaryProcessor);
        imageManager.showScaledUp(redBinaryImg);
        imagesToSave.add(redBinaryImg);


        //binarization overlay that shows attempt one vs attempt with colors
        ImageProcessor binaryOverlay = originalImgScaled.getProcessor().duplicate().convertToColorProcessor();
        OverlayManager.overlayBinaryAttempts(binaryOverlay, binary1, binary2);
        ImagePlus binaryOverlayImg = new ImagePlus("Binarization Overlay", binaryOverlay);
        String binaryOverlayDescription = "Binarization Overlay Description: Green: non-fibrillar to fibrillar; magenta: fibrillar to non-fibrillar; black/ white: unchanged.";
        // Show the description in a TextWindow
        //TextWindow textWindow1 = new TextWindow("Binarization Overlay Description", binaryOverlayDescription, 800, 100);
        IJ.log(binaryOverlayDescription);
        imageManager.showScaledUp(binaryOverlayImg);
        imagesToSave.add(binaryOverlayImg);

        //TODO:  abstract this code into its own class. if centroids over certain number increase thin_inc
        //TODO: need a way to somehow increase decrease this based on centroids number and image size to maximaze results relative to time
        int thin_inc = 10;

        //get fibril pixels from 2 binarization attempt
        ArrayList<double[]> fibrilPixels = getFibrilPixels(binaryImg2.getProcessor().convertToFloatProcessor());
        KDTree<double[]> kdTree2 = KDTree.of(convertedCentroids);
        double[] nearestCentroidDistances = new double[fibrilPixels.size()];
        //associate each fibril pixel with a centroid.
        //The index represents each fibril pixel and the value represents the index in centroids array.
        int[] nearestCentroidIndices = new int[fibrilPixels.size()];
        // Perform the nearest neighbor search
        for (int i = 0; i < fibrilPixels.size(); i++) {
            Neighbor<double[], double[]> nearest = kdTree2.nearest(fibrilPixels.get(i));
            nearestCentroidIndices[i] = nearest.index;
            nearestCentroidDistances[i] = nearest.distance;
        }
        double[][] fibrilPixelsThinned = new double[(fibrilPixels.size() / thin_inc) + 1][2];
        double[] centroidsIdxThinned = new double[(nearestCentroidIndices.length / thin_inc) + 1];
        for (int i = 0, j = 0; i < fibrilPixels.size(); i += thin_inc, j++) {
            fibrilPixelsThinned[j] = fibrilPixels.get(i);
            centroidsIdxThinned[j] = nearestCentroidIndices[i];
        }

        //test to see what pixels are being processed by em
        ImageProcessor pixelsForGmm = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
        OverlayManager.overlayPixels(pixelsForGmm,fibrilPixelsThinned, Color.red.getRGB());
        ImagePlus t = new ImagePlus("Pixels for gmm", pixelsForGmm);
        //t.show();

        //EM Stem
        GaussianMixtureModel gmm = new GaussianMixtureModel(convertedCentroids.length, fibrilPixelsThinned, centroidsIdxThinned);
        gmm.EM(false, 50);

        ImageProcessor boundryIp = imgAfterGausFiltering.getProcessor();

        //this is only used to calculate the contour lines but do we need to do this really?
        //double[][][] post_all = gmm.calculatePosterProbabilityMatrixAll(allPixels, boundryIp.getHeight(), boundryIp.getWidth());

        //calculate posterior probabilities that each pixel belongs to a fibril
        double[][] post_fibril = gmm.calculatePosterProbabilityMatrixFibril(fibrilPixels, boundryIp.getHeight(), boundryIp.getWidth());

        IJ.showStatus("Assigning Clusters");
        IJ.log("Assigning Clusters");

        //assign pixels to a fibril based on probabilites that each pixel belongs to a centroid
        int[] clusterAssignments = assignEachFibrilPixelToFibril(fibrilPixels.size(), convertedCentroids.length, post_fibril, isBoundryCentroid);

        //processor setup
        ImageProcessor contour = originalImgScaled.getProcessor().duplicate().convertToColorProcessor();
        ImageProcessor contour_lines = originalImgScaled.getProcessor().duplicate().convertToColorProcessor();
        ImageProcessor intiVoronoi = originalImgScaled.getProcessor().duplicate().convertToColorProcessor();

        //test to show cluster assignments
//        ImageProcessor clusterColorsProcessor = imgAfterGausFiltering.getProcessor().duplicate();
//        OverlayManager.overlayColoredClusterAssignments(clusterColorsProcessor, fibrilPixels, clusterAssignments);
//        new ImagePlus("Clusters check", clusterColorsProcessor).show();


        //{*}IMAGE PREP
        IJ.showStatus("Preparing Overlayed Images");
        IJ.log("Preparing Images");

        OverlayManager.overlayCentroids(contour, gmm.getMus(), Color.red.getRGB()); // gmm centroids red
        OverlayManager.overlayCentroids(contour, convertedCentroids, Color.magenta.getRGB()); // voronoi centroids magenta
        OverlayManager.overlayVoronoi(contour, vm.voronoi, cyan.getRGB());
        //OverlayManager.overlayContourLines(contour, post_all, Color.blue.getRGB());

        //OverlayManager.overlayContourLines(contour_lines, post_all, Color.blue.getRGB());
        OverlayManager.overlayCentroids(contour_lines, gmm.getMus(), Color.red.getRGB()); // gmm centroids red


        OverlayManager.overlayVoronoi(intiVoronoi, vm.init_voronoi, cyan.getRGB()); //overlay corrected voronoi over initial guess of voronoi
        OverlayManager.overlayCentroids(intiVoronoi, initialCentroidCopy.toArray(new double[0][]), cyan.getRGB());
        OverlayManager.overlayVoronoi(intiVoronoi, vm.voronoi, Color.blue.getRGB()); //overlay corrected voronoi over initial guess of voronoi
        OverlayManager.overlayCentroids(intiVoronoi, convertedCentroids, Color.red.getRGB()); // corrected voronoi centroids

        //shows contour overlays not needed
        ImagePlus overlayImage = new ImagePlus("Overlays", contour);

        ImagePlus contourLinesImage = new ImagePlus("Contour Lines", contour_lines);
         //contourLinesImage.show();
        ImagePlus correctedVoronoiImage = new ImagePlus("Corrected Voronoi", intiVoronoi);
        //imageManager.showScaledUp(contourLinesImage);
        //imageManager.showScaledUp(overlayImage);
        imageManager.showScaledUp(correctedVoronoiImage);
        //imagesToSave.add(overlayImage);
        //imagesToSave.add(contourLinesImage);
        imagesToSave.add(correctedVoronoiImage);


        IJ.showStatus("Fitting Ellipses");
        IJ.log("Fitting Ellipses");

        //gmm centroids or user corrected centroids?
        //TODO: figure this out, test both. right now using vm centroids
        ArrayList<double[]> centroidsFromGmm = convertToArrayList(gmm.getMus());

        EllipseFitting ef = new EllipseFitting(originalImgScaled, fibrilPixels, vm.getCentroids() , clusterAssignments, post_fibril, gmm.getComponentProportions() ,  nonBoundryCentroids);
        ef.fitEllipses(isBoundryCentroid);

        //user can delete ellipses that did not perform well
        WaitForUserDialog wfud2 = new WaitForUserDialog("Action Required", "Left click to remove ellipses. Once finished, click ok.");
        wfud2.show();
        if (wfud2.escPressed()) return;
        imagesToSave.add(ef.getImage());

        //{*}POST PROCESSING
        IJ.showStatus("Post Processing...");

        //conversion from pixels to nm
        double nanometers_over_pixels = imageManager.getConversionFactor();
        //Number times original image got scaled down by half
        int scale = imageManager.getScale();

        //Ellipse data
        ArrayList<Ellipse> ellipses = ef.fibrilEllipses;

        //ellipse scaling
        for(Ellipse e : ellipses){
            e.scale(scale);
        }

        //gmm scaling
        //double[] area_pix2 = computeFibrilArea(post_fibril, scale, nonBoundryCentroids.size());

        ArrayList<Double> areaFromEllipse_pix2 = ef.area_pix;
        for(int i = 0; i < areaFromEllipse_pix2.size(); i++){
            areaFromEllipse_pix2.set(i, areaFromEllipse_pix2.get(i) * Math.pow(scale, 2));
        }

        // converting from pixels to nm
        for(Ellipse e : ellipses){
            e.convertToNM(nanometers_over_pixels);
        }

        //get image name without file extension
        int dotIndex = imageName.lastIndexOf('.');
        if (dotIndex > 0 && dotIndex <= imageName.length() - 2) {
            imageName = imageName.substring(0, dotIndex);
        }
        // Get the path of the image currently open in ImageJ
        //String path = IJ.getDirectory("Choose a directory to create");
        JFileChooser chooser = new JFileChooser();
        chooser.setDialogTitle("Choose a directory to create");
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);

        int returnValue = chooser.showOpenDialog(null);

        String path;
        if (returnValue == JFileChooser.APPROVE_OPTION) {
            File selectedDirectory = chooser.getSelectedFile();
            System.out.println("Selected directory: " + selectedDirectory.getAbsolutePath());
            path = selectedDirectory.getAbsolutePath() + File.separator;
        } else {
            System.out.println("No directory selected.");
            path = null;
        }

        Path newFolderPath;

        if (path != null) {
            createDirectory(path + imageName);
            newFolderPath = Paths.get(path + imageName);
        } else {
            IJ.showMessage("No directory chosen to save images");
            newFolderPath = null;
        }

        // Generate csv files
        try{
            File Results_In_Pixels = new File(newFolderPath.toString() + File.separator + "Results_In_Pixels.csv");
            if(Results_In_Pixels.delete()){
                Results_In_Pixels.createNewFile();
            }
            FileWriter out = new FileWriter(Results_In_Pixels);
            CSVPrinter printer = new CSVPrinter(out, CSVFormat.DEFAULT.withHeader(
                    "Area (pixels^2)", "Major Radius (pixels)", "Minor Radius (pixels)",
                    "Angle (degrees)", "Centroid X (pixels)", "Centroid Y (pixels)",
                    "Covariance XX (pixels^2)", "Covariance XY (pixels^2)", "Covariance YY (pixels^2)",
                    "Component Proportion"));

            //TODO: if remove problematic centroids, need to use that size for this loop. Need to use ellipse size
            for(Ellipse e : ellipses){
                printer.printRecord(formatField(e.postProbArea),
                        formatField(e.majorRadius),
                        formatField(e.minorRadius),
                        formatField(e.angle),
                        formatField(e.x),
                        formatField(e.y),
                        formatField(e.cov[0][0]),
                        formatField(e.cov[0][1]),
                        formatField(e.cov[1][1]),
                        formatField(e.componentProportion));
            }
            out.close();

        }catch (Exception e){
            e.printStackTrace();
        }

        try{
            File Results_In_Nm = new File(newFolderPath.toString() + File.separator + "Results_In_Nm.csv");
            if(Results_In_Nm.delete()){
                Results_In_Nm.createNewFile();
            }
            FileWriter out2 = new FileWriter(Results_In_Nm);
            CSVPrinter printer = new CSVPrinter(out2, CSVFormat.DEFAULT.withHeader(
                    "Area (nanometers^2)", "Major Radius (nanometers)", "Minor Radius (nanometers)",
                    "Aspect Ratio"));

            //TODO: need to use final ellipse size
            for(Ellipse e : ellipses){
                printer.printRecord(formatField(e.area_nm),
                        formatField(e.majorRadius_nm),
                        formatField(e.minorRadius_nm),
                        formatField(e.aspectRatio_nm));
            }
            out2.close();

        }catch (Exception e){
            e.printStackTrace();
        }


        //save images to path specified from user
        for(ImagePlus image: imagesToSave){
            imageManager.scaleUp(image);
            String fileName = image.getTitle();
            String savePath = newFolderPath.resolve(fileName + ".tiff").toString();
            //String savePath = newFolderPath + fileName + ".tiff";
            FileSaver fileSaver = new FileSaver(image);
            if(fileSaver.saveAsTiff(savePath)){
                System.out.println("Saved " + fileName);
            }
            else{
                System.out.println("Save failed");
            }
        }


        //DONE
        IJ.showStatus("Processing done");
        IJ.log("Processing done. Check directory specified for results");
        System.out.println("Processing done");
    }

    private static String formatField(double value) {
        return String.format("%-28.8f", value);
    }

    public void createDirectory(String path) {
        File directory = new File(path);

        if (!directory.exists()) {
            try {
                boolean isCreated = directory.mkdirs();

                if (isCreated) {
                    IJ.showMessage("Directory created successfully: " + path);
                } else {
                    IJ.showMessage("Failed to create directory: " + path);
                }
            } catch (SecurityException e) {
                IJ.showMessage("Permission denied: Unable to create directory. " + e.getMessage());
            }
        } else {
            IJ.showMessage("Directory already exists: " + path);
        }
    }


    /**
     * Main method for debugging.
     *
     * For debugging, it is convenient to have a method that starts ImageJ, loads
     * an image and calls the plugin, e.g. after setting breakpoints.
     *
     * @param args unused
     */
    public static void main(String[] args) throws Exception {
        // set the plugins.dir property to make the plugin appear in the Plugins menu
        // see: https://stackoverflow.com/a/7060464/1207769
        Class<?> clazz = CollagenAnalysisPlugin.class;
        java.net.URL url = clazz.getProtectionDomain().getCodeSource().getLocation();
        java.io.File file = new java.io.File(url.toURI());
        System.setProperty("plugins.dir", file.getAbsolutePath());

        // start ImageJ
        new ImageJ();

        // open the Clown sample
        //TEM_4.jpg
        //L3915_3c_BA190R-1.tif
        ImagePlus image = IJ.openImage("/Users/tywatson/development/repos/Rego-Imagej-pluggin/testimages/Adventitia_RIRII.tif");
        //ImagePlus image = IJ.openImage("C:\\Users\\Ty Watson\\Downloads\\testimages\\Adventitia_RIRII.tif");
        image.show();
        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}
