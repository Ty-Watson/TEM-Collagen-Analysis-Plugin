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
import ij.gui.GenericDialog;
import ij.gui.WaitForUserDialog;
import ij.plugin.filter.PlugInFilter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import smile.neighbor.KDTree;
import smile.neighbor.Neighbor;

import java.awt.*;
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

    // image property members
    private int width;
    private int height;

    // plugin parameters
    public double value;
    public String name;

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


        ImageManager imageManager = new ImageManager(preProcessProcessor);
        imageManager.DrawScaleBar(imp);
        imageManager.scaleDown();
        ImageProcessor ip = imageManager.getProcessor();
        ImagePlus originalImgScaled = new ImagePlus("OrignalImgScaled", ip);
        originalImgScaled.show();
//        Img img = ij.convert().convert(originalImgScaled, Img.class);
//        ij.ui().show(img);

        // Create extended image with white border
        ImageProcessor processorWithExcludedRegions = imageManager.excludeRegions();

        ImagePlus imageAfterUserExclude = new ImagePlus("Image after excluding regions ", processorWithExcludedRegions);
        imageAfterUserExclude.show();


        //IJ.showMessage("Interactive Image Editing", "Editing complete.");
        //toolbar.setTool(Toolbar.HAND);

        //[*]Filtering, thresholding, and binarization[*]

        //apply gaussian filtering to original img
        ImagePlus imgAfterGausFiltering = applyGaussianFiltering(processorWithExcludedRegions);

        ImageProcessor ip2 = processorWithExcludedRegions.duplicate().convertToFloatProcessor();

        //smooth intensity values for binarization step
        ImagePlus imageAfterSmoothing = smoothIntensityValues(ip2);

        FloatProcessor ip3 = ip2.duplicate().convertToFloatProcessor();

        //binarize the smoothed img using Otsu optimal threshold
        ImagePlus binaryImg = binarizeUsingOtsu(ip3);

//        for (int x = 0; x < ip3.getWidth(); x++) {
//            for (int y = 0; y < ip3.getHeight(); y++) {
//                if (ip3.getPixel(x, y) == 0) { // Check if the pixel is part of an excluded region
//                    ip3.putPixel(x, y, 255); // Set excluded pixels to 0 in the original image which turns them black
//                }
//            }
//        }

        //overlay the binarized img(in red) with the original
        //overlayBinaryOnGrayscale2(originalImgScaled, binaryImg);

        //FloatProcessor fp = (FloatProcessor)ip3.duplicate();


        CentroidManager cm = new CentroidManager(binaryImg.getProcessor().duplicate());
        ArrayList<double[]> centroids = cm.findCentroids();

        ArrayList<double[]> initialCentroidCopy = (ArrayList<double[]>) centroids.clone();





        VoronoiManager vm = new VoronoiManager(imgAfterGausFiltering);
        vm.setCentroids(centroids);

        vm.initDrawVoronoi();
        vm.showImg();

        WaitForUserDialog wfud = new WaitForUserDialog("Action Required", "edit centroids");
        wfud.show();
        if (wfud.escPressed()) return;

        ArrayList<double[]> allPixels = getAllPixels(ip2);

        double[][] convertedCentroids = vm.getCentroids().toArray(new double[0][]);
        writeArrayToFile(convertedCentroids, "centroids.Text");
        // Create a KD-tree for centroids
        KDTree<double[]> kdTree = KDTree.of(convertedCentroids);
        double[] nearestCentroidDistancesAll = new double[allPixels.size()];
        int[] nearestCentroidIndicesAll = new int[allPixels.size()];
        // Perform the nearest neighbor search
        for (int i = 0; i < allPixels.size(); i++) {
            Neighbor<double[], double[]> nearest = kdTree.nearest(allPixels.get(i));
            nearestCentroidIndicesAll[i] = nearest.index;
            nearestCentroidDistancesAll[i] = nearest.distance;
        }

        vm.SetNearestCentroidDistances(nearestCentroidDistancesAll);
        vm.SetNearestCentroidIndices(nearestCentroidIndicesAll);


        ImagePlus smooth2 = vm.Interpolation(imageAfterSmoothing.getProcessor());
        smooth2.show();

        FloatProcessor ip4 = smooth2.getProcessor().duplicate().convertToFloatProcessor();

        //binarize the smoothed img using Otsu optimal threshold
        ImagePlus binaryImg2 = binarizeUsingOtsu(ip4);

        //TODO:  abstract this code into its own class. if centroids over certain number increase thin_inc
        int thin_inc = 20;
        double regularizationValue = 0.1;
        //get all fibril pixels


        //need to get from 2 binarization attempt
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

        ImageProcessor pixelsForGmm = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
        OverlayManager.overlayPixels(pixelsForGmm,fibrilPixelsThinned, Color.red.getRGB());
        ImagePlus t = new ImagePlus("Pixels for gmm", pixelsForGmm);
        t.show();

        writeArrayToFile(centroidsIdxThinned, "c_idx_thinned.Text");

        GaussianMixtureModel gmm = new GaussianMixtureModel(convertedCentroids.length, fibrilPixelsThinned, centroidsIdxThinned);
        gmm.EM(false, 50);
        double[][] p = gmm.getPosteriorProbabilityMatrix();
        //writeArrayToFile(p, "p.Text");
        //double[][] p = readArrayFromFile("p.Text");


        ImageProcessor boundryIp = imgAfterGausFiltering.getProcessor();

        double[][][] post_all = gmm.calculatePosterProbabilityMatrixAll(allPixels, boundryIp.getHeight(), boundryIp.getWidth());
        double[][] post_fibril = gmm.calculatePosterProbabilityMatrixFibril(fibrilPixels, boundryIp.getHeight(), boundryIp.getWidth());
        //writeArrayToFile(post_all, "p_all.Text");
        //double[][][] post_all = readArrayFromFile2("p_all.Text");

        //drawBoundries(boundryIp, post_all);

        IJ.showStatus("Assigning Clusters");
        int[] clusterAssignments = assignEachFibrilPixelToFibril(fibrilPixels.size(), convertedCentroids.length, post_fibril);


        ImageProcessor contour = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
        ImageProcessor contour_lines = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
        ImageProcessor intiVoronoi = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
//        ImageProcessor clusterColorsProcessor = imgAfterGausFiltering.getProcessor().duplicate();
//
//        OverlayManager.overlayColoredClusterAssignments(clusterColorsProcessor, fibrilPixels, clusterAssignments);
//
//        new ImagePlus("Clusters check", clusterColorsProcessor).show();


        IJ.showStatus("Preparing Overlayed Images");

        OverlayManager.overlayCentroids(contour, gmm.getMus(), Color.red.getRGB()); // gmm centroids red
        OverlayManager.overlayCentroids(contour, convertedCentroids, Color.magenta.getRGB()); // voronoi centroids magenta
        OverlayManager.overlayVoronoi(contour, vm.voronoi, cyan.getRGB());
        OverlayManager.overlayContourLines(contour, post_all, Color.blue.getRGB());

        OverlayManager.overlayContourLines(contour_lines, post_all, Color.blue.getRGB());
        OverlayManager.overlayCentroids(contour_lines, gmm.getMus(), Color.red.getRGB()); // gmm centroids red


        OverlayManager.overlayVoronoi(intiVoronoi, vm.init_voronoi, cyan.getRGB()); //overlay corrected voronoi over initial guess of voronoi
        OverlayManager.overlayCentroids(intiVoronoi, initialCentroidCopy.toArray(new double[0][]), cyan.getRGB());
        OverlayManager.overlayVoronoi(intiVoronoi, vm.voronoi, Color.blue.getRGB()); //overlay corrected voronoi over initial guess of voronoi
        OverlayManager.overlayCentroids(intiVoronoi, convertedCentroids, Color.red.getRGB()); // corrected voronoi centroids



        ImagePlus resultImage = new ImagePlus("Overlays", contour);
        ImagePlus resultImage2 = new ImagePlus("Contour Lines", contour_lines);
        ImagePlus resultImage3 = new ImagePlus("Corrected Voronoi", intiVoronoi);
        resultImage.show();
        resultImage2.show();
        resultImage3.show();

        IJ.showStatus("Fitting Ellipses");
        EllipseFitting ef = new EllipseFitting();
        ef.fitEllipses(fibrilPixels, convertedCentroids.length, clusterAssignments);
        ImageProcessor ellipseP = imgAfterGausFiltering.getProcessor().duplicate().convertToColorProcessor();
        OverlayManager.drawEllipsesOnImage(ellipseP, ef.getxEllipses(), ef.getyEllipses(), convertedCentroids.length);
        OverlayManager.overlayCentroids(ellipseP, convertedCentroids, Color.red.getRGB());
        ImagePlus ei = new ImagePlus("Ellipse", ellipseP);
        ei.show();


        //POST PROCESSING
        IJ.showStatus("Post Processing...");

        double nanometers_over_pixels = imageManager.getConversionFactor();
        int scale = imageManager.getScale();

        ArrayList<double[]> radiusPix = ef.getRadius_pix();
        ArrayList<double[]> x_ellipse = ef.getxEllipses();
        ArrayList<double[]> y_ellipse = ef.getyEllipses();
        double[][] mus = gmm.getMus();


        radiusPix.forEach(arr -> Arrays.setAll(arr, i -> arr[i] * scale));
        x_ellipse.forEach(arr -> Arrays.setAll(arr, i -> arr[i] * scale));
        y_ellipse.forEach(arr -> Arrays.setAll(arr, i -> arr[i] * scale));
        Arrays.stream(mus).forEach(arr -> Arrays.setAll(arr, i -> arr[i] * scale));


        double[] area_pix2 = computeFibrilArea(post_fibril, scale, mus.length);
        Arrays.setAll(area_pix2, i -> area_pix2[i] * Math.pow(nanometers_over_pixels, 2));



    }

    private boolean showDialog() {
        GenericDialog gd = new GenericDialog("Process pixels");

        // default value is 0.00, 2 digits right of the decimal point
        gd.addNumericField("value", 0.00, 2);
        gd.addStringField("name", "John");

        gd.showDialog();
        if (gd.wasCanceled())
            return false;

        // get entered values
        value = gd.getNextNumber();
        name = gd.getNextString();

        return true;
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
        ImagePlus image = IJ.openImage("/Users/tywatson/development/repos/Rego-Imagej-Pluggin/testimages/Adventitia_RIRII.tif");
        image.show();
        // run the plugin
        IJ.runPlugIn(clazz.getName(), "");
    }
}
