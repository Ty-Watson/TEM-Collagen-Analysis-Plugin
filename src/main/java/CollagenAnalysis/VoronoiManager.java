package CollagenAnalysis;

import de.alsclo.voronoi.Voronoi;
import de.alsclo.voronoi.graph.Edge;
import de.alsclo.voronoi.graph.Graph;
import de.alsclo.voronoi.graph.Point;
import de.alsclo.voronoi.graph.Vertex;
import ij.ImagePlus;
import ij.gui.ImageCanvas;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import java.awt.*;
import java.awt.event.AWTEventListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.util.stream.Stream;

import static CollagenAnalysis.Constants.SCALE;
import static CollagenAnalysis.Constants.pointSize;
import static CollagenAnalysis.ImageProcessingUtils.findNearestCentroid;


public class VoronoiManager {
    public VoronoiManager(ImagePlus imp, boolean[][] exclusionMask){

        ip = imp.getProcessor().duplicate();
        setScaleUpImage();
        this.exclusionMask = exclusionMask;
        applyExclusionMask();
        this.scaledUpImage.show();

        setCanvasEvents();
    }
    public ImagePlus scaledUpImage;
    //have to have processor here for scaled up image because ip.reset only works with instatiated processor rather than imageplus.getprocessor()
    //so in order to have the image reset and draw the updated changes, have to call imageprocessor.reset
    public ImageProcessor scaledUpProcessor;
    private ImageProcessor ip;
    //user corrected voronoi diagram
    public Voronoi voronoi;
    public boolean[][] exclusionMask;
    //inital voronoi diagram generated
    public Voronoi init_voronoi;

    private ArrayList<double[]> centroids;

    private List<Point> voronoi_centroids;
    private  ImageCanvas canvas;
    private HashMap<Point, Set<Vertex>> cells = new HashMap<>();
    double[] nearestCentroidDistances;
    int[] nearestCentroidIndices;
    List<int[]> _allpixels = new ArrayList<>();

    public ImagePlus Interpolation(ImageProcessor ip, ArrayList<double[]> centroids){
        ImageProcessor p = ip.duplicate();
        ImagePlus n = new ImagePlus("Filtered image before Interpolation", ip);
        //n.show();
        double[] charIntensities = generateCharCellIntensity(ip, centroids);
        Interpolator interpolator = new Interpolator(centroids, charIntensities);
        ImageProcessor testp = p.duplicate();
        ImagePlus test = new ImagePlus("interpolatedIntesity", testp);

        for (int y = 0; y < p.getHeight(); y++) {
            for (int x = 0; x < p.getWidth(); x++) {
                double originalIntensity = p.getPixelValue(x, y);
                double interpolatedIntensity = interpolator.interpolate(x, y);
                testp.putPixelValue(x,y,interpolatedIntensity);
                double transformedIntensity = getTransformedIntensity(interpolatedIntensity, originalIntensity);

                p.putPixelValue(x, y, transformedIntensity);
            }
        }
        //test.show();
        ImagePlus i = new ImagePlus("after interpolation", p);
        return i;


    }

    private static double getTransformedIntensity(double interpolatedIntensity, double originalIntensity) {
        double transformedIntensity;


        if (interpolatedIntensity == 0 || Double.isNaN(interpolatedIntensity)) {
            //System.out.printf("[*] Interpolation failed for pixel (%d, %d) \n", x, y);
            transformedIntensity = originalIntensity;  // Avoid division by zero or log(0)
        } else {
            if(originalIntensity == interpolatedIntensity){
                transformedIntensity = 0.5 * 255;
            }
            else{
                transformedIntensity = Math.pow(originalIntensity / 255, Math.log(0.5) / Math.log(interpolatedIntensity / 255 )) * 255;
            }

        }
        return transformedIntensity;
    }

    public double[] generateCharCellIntensity(ImageProcessor ip, ArrayList<double[]> centroids){

        int width = ip.getWidth();
        int height = ip.getHeight();


        List<Double> pixelIntensities = new ArrayList<>();
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double intensity = ip.getPixelValue(x, y);
                pixelIntensities.add(intensity);
                _allpixels.add(new int[]{x, y});
            }
        }

        //double[] charCellIntensity = computeIntensities(pixelIntensities, nearestCentroidDistances, nearestCentroidIndices, centroids.size());
        double[] charCellIntensity2 = computeIntensities2(pixelIntensities, nearestCentroidDistances, nearestCentroidIndices, ip, centroids);


        //charCellIntensity now contains the characteristic intensity for each cell
        return charCellIntensity2;
    }
    public void generateCells(){
        Stream<Edge> edges = voronoi.getGraph().edgeStream();
        for (Iterator<Edge> it = edges.iterator(); it.hasNext(); ) {
            Edge edge = it.next();
            if (edge.getA() != null && edge.getB() != null) {
                // Scale down the points because voronoi is being generated on full image
                Point scaledA = new Point(edge.getA().getLocation().x / SCALE, edge.getA().getLocation().y / SCALE);
                Point scaledB = new Point(edge.getB().getLocation().x / SCALE, edge.getB().getLocation().y / SCALE);

                //scale down these values because the voronoi edge sites are coming from full image
                Point site1 = new Point(edge.getSite1().x / 2, edge.getSite1().y / 2);
                Point site2 = new Point(edge.getSite2().x / 2, edge.getSite2().y / 2);

                cells.computeIfAbsent(site1, k -> new HashSet<>()).add(new Vertex(scaledA));
                cells.computeIfAbsent(site1, k -> new HashSet<>()).add(new Vertex(scaledB));
                cells.computeIfAbsent(site2, k -> new HashSet<>()).add(new Vertex(scaledA));
                cells.computeIfAbsent(site2, k -> new HashSet<>()).add(new Vertex(scaledB));
            }
        }

        List<Point> centroids_ = new ArrayList<>();

        cells.forEach((point, vertices) -> {
            // System.out.println("Cell centered at: " + point);
            centroids_.add(point);
            //System.out.println("Vertices:");
            //vertices.forEach(vertex -> System.out.println(vertex));
        });
        voronoi_centroids = centroids_;
    }

    public void drawVoronoi(){
        scaledUpProcessor.snapshot();
        Collection<Point> points = new ArrayList<>();
        for( double[] centroid : centroids){
            de.alsclo.voronoi.graph.Point point = new de.alsclo.voronoi.graph.Point(centroid[0] * SCALE, centroid[1] * SCALE);
            points.add(point);
        }
        voronoi = new Voronoi(points);
        Graph g = voronoi.getGraph();
        generateCells();


        //ImageProcessor ip = imp.getProcessor().duplicate();

        scaledUpProcessor.setColor(Color.BLUE);


        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            Point a = e.getA().getLocation();
            Point b = e.getB().getLocation();
            scaledUpProcessor.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });
        scaledUpImage.updateAndDraw();
        overlayMaxima(Color.red.getRGB());
        applyExclusionMask();
    }
    public void initDrawVoronoi(){
        scaledUpProcessor.snapshot();
        Collection<Point> points = new ArrayList<>();
        for( double[] centroid : centroids){
            Point point = new Point(centroid[0] * SCALE, centroid[1] * SCALE);
            points.add(point);
        }
        voronoi = new Voronoi(points);
        init_voronoi = voronoi;
        Graph g = voronoi.getGraph();
        generateCells();


        //ImageProcessor ip = imp.getProcessor().duplicate();

        scaledUpProcessor.setColor(Color.BLUE);


        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            Point a = e.getA().getLocation();
            Point b = e.getB().getLocation();
           scaledUpProcessor.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });

        //imageManager.drawExcludedRegions(ip);

        scaledUpImage.updateAndDraw();
        overlayMaxima(Color.red.getRGB());
        applyExclusionMask();
    }
    private void overlayMaxima(int color) {
        // Convert the binarized image to RGB to overlay red points
        ImageProcessor d = scaledUpProcessor.duplicate();
        ColorProcessor colorProcessor = d.convertToColorProcessor();

        // Define the red color
        //int red = (255 << 16); // RGB value for red

        // Draw each maximum point as a red circle on the RGB image
        for (double[] point : centroids) {
            double x = point[0] * SCALE;
            double y = point[1] * SCALE;

            // Draw a circle or a larger point at (x, y) in red
            for (int dx = -pointSize; dx <= pointSize; dx++) {
                for (int dy = -pointSize; dy <= pointSize; dy++) {
                    if (dx * dx + dy * dy <= pointSize * pointSize) {
                        double newX = x + dx;
                        double newY = y + dy;
                        if (newX >= 0 && newX < colorProcessor.getWidth() && newY >= 0 && newY < colorProcessor.getHeight()) {
                            colorProcessor.set((int)newX, (int)newY, color);
                        }
                    }
                }
            }
        }

        // Display the result
        scaledUpImage.setProcessor(colorProcessor);
    }
    public void setCentroids(ArrayList<double[]> centroids){
        this.centroids = centroids;
    }
    public void showImg(){
        scaledUpImage.show();
    }
    public void closeImg(){scaledUpImage.close();}
    public void setScaleUpImage(){
        ImageProcessor imageProcessor = ip.duplicate();
        int h = imageProcessor.getHeight() * SCALE;
        int w = imageProcessor.getWidth() * SCALE;
        scaledUpProcessor = imageProcessor.resize(w, h);
        scaledUpImage = new ImagePlus("Voronoi Image", scaledUpProcessor);
    }

    private void setCanvasEvents(){
        canvas = scaledUpImage.getCanvas();

        // Add an AWTEventListener to globally intercept mouse events
        Toolkit.getDefaultToolkit().addAWTEventListener(new AWTEventListener() {
            @Override
            public void eventDispatched(AWTEvent event) {
                if (event instanceof MouseEvent) {
                    MouseEvent mouseEvent = (MouseEvent) event;
                    if (mouseEvent.getID() == MouseEvent.MOUSE_PRESSED && mouseEvent.getButton() == MouseEvent.BUTTON3) {
                        mouseEvent.consume(); // Consume the event to prevent context menu

                        // Check if the event occurred on our canvas
                        if (mouseEvent.getSource() == canvas) {
                            int clickX = canvas.offScreenX(mouseEvent.getX()) / SCALE;
                            int clickY = canvas.offScreenY(mouseEvent.getY()) / SCALE;


                            // Find and remove the nearest centroid to the click
                            double[] nearestCentroid = findNearestCentroid(centroids, clickX, clickY);
                            if (nearestCentroid != null) {
                                centroids.remove(nearestCentroid);
                                scaledUpProcessor.reset();
                                drawVoronoi();
                            }
                        }
                    }
                }
            }
        }, AWTEvent.MOUSE_EVENT_MASK);

        canvas.addMouseListener(new MouseAdapter() {
            @Override
            public void mousePressed(MouseEvent e) {
                if (e.getButton() == MouseEvent.BUTTON1) { // Left click add centroid
                    int clickX = canvas.offScreenX(e.getX()) / SCALE;
                    int clickY = canvas.offScreenY(e.getY()) / SCALE;



                    //cant add centroid in exclusion mask
                    if(exclusionMask[clickX][clickY]){
                        return;
                    }


                    double[] newCentroid = new double[]{ clickX, clickY};
                    centroids.add(newCentroid);
                    scaledUpProcessor.reset();
                    drawVoronoi();
                }
            }
        });

    }
    public double[] computeIntensities2(List<Double> pixelIntensities, double[] nearest_centroid_dist_all, int[] nearest_centroid_idx_all, ImageProcessor ip, ArrayList<double[]> centroids) {
// The result array
        double[] char_cell_intensity = new double[centroids.size()];
        Arrays.fill(char_cell_intensity, Double.NaN);

        ImageProcessor test = ip.duplicate();
        ImagePlus testimg = new ImagePlus("White pixels represent char cell intensity for each centroid", test);

        for (int i = 0; i < char_cell_intensity.length; i++) {
            // Create a list for collecting intensities
            List<Double> intensities = new ArrayList<>();

            // Find the 5th percentile distance for the current centroid
            List<Double> distances = new ArrayList<>();
            for (int j = 0; j < nearest_centroid_idx_all.length; j++) {
                if (nearest_centroid_idx_all[j] == i) {
                    distances.add(nearest_centroid_dist_all[j]);
                }
            }
            double fifthPercentileDistance = quantile(distances.stream().mapToDouble(d -> d).toArray(), 0.05);

            // Collect intensities for pixels closer than the 5th percentile distance
            for (int j = 0; j < nearest_centroid_idx_all.length; j++) {
                if (nearest_centroid_idx_all[j] == i && nearest_centroid_dist_all[j] < fifthPercentileDistance) {
                    // Assuming you can get the intensity from I_filt_tform based on the pixel's index j
                    //test for showing white pixels representing char pixels in original image
                    int[] pixel = _allpixels.get(j);

                    test.putPixelValue(pixel[0], pixel[1], 255);
                    intensities.add(pixelIntensities.get(j));
                }
            }

            // Calculate the mean intensity for the current centroid
            char_cell_intensity[i] = mean(intensities);
        }
        //testimg.show();
        return char_cell_intensity;
    }
    // Assuming you have a method to calculate the quantile
    public static double quantile(double[] data, double quantile) {
        Arrays.sort(data);
        int index = (int) Math.ceil(quantile * data.length);
        return data[index];
    }

    // Assuming you have a method to calculate the mean of a list
    public static double mean(List<Double> list) {
        return list.stream()
                .mapToDouble(Double::doubleValue)
                .average()
                .orElse(Double.NaN); // handle empty list case
    }
    public void SetNearestCentroidIndices(int[] nearestCentroidIndices){
        this.nearestCentroidIndices = nearestCentroidIndices;
    }
    public void SetNearestCentroidDistances(double[] nearestCentroidDistances){
        this.nearestCentroidDistances = nearestCentroidDistances;
    }
    public ArrayList<double[]> getCentroids(){
        return centroids;
    }

    public ArrayList<double[]> identifyBoundaryFibrils(boolean[][] exclusionMask) {
        ArrayList<double[]> nonBoundryCentroids = new ArrayList<>();
        for(double[] centroid: centroids){
            Point centroidPont = new Point(centroid[0], centroid[1]);
            Set<Vertex> cellVertices = cells.get(centroidPont);
            if(cellVertices != null){
                boolean isBoundary = false;

                for(Vertex v : cellVertices){
                    Point locationPoint = v.getLocation();
                    double[] location = new double[]{locationPoint.x, locationPoint.y};
                    if (location[0] < 0 || location[1] < 0 || location[0] >= ip.getWidth() || location[1] >= ip.getHeight()) {
                        isBoundary = true;
                        break;
                    }
                    else if(exclusionMask[(int)location[0]][(int)location[1]]){
                        isBoundary = true;
                        break;
                    }
                }
                if(!isBoundary){
                    nonBoundryCentroids.add(centroid);
                }
            }
            else{
                System.out.printf("Cant get vertices for centroids: %f %f \n", centroid[0], centroid[1]);
            }

        }
        return nonBoundryCentroids;

    }
    public ArrayList<Boolean> identifyBoundaryFibrilsBoolean(boolean[][] exclusionMask) {

        ArrayList<Boolean> isBoundaryCentroid = new ArrayList<Boolean>();
        for(double[] centroid: centroids){
            Point centroidPont = new Point(centroid[0], centroid[1]);
            Set<Vertex> cellVertices = cells.get(centroidPont);
            if(cellVertices != null){
                boolean isBoundary = false;

                for(Vertex v : cellVertices){
                    Point locationPoint = v.getLocation();
                    double[] location = new double[]{locationPoint.x, locationPoint.y};
                    if (location[0] < 0 || location[1] < 0 || location[0] >= ip.getWidth() || location[1] >= ip.getHeight()) {
                        isBoundary = true;
                        break;
                    }
                    else if(exclusionMask[(int)location[0]][(int)location[1]]){
                        isBoundary = true;
                        break;
                    }
                }
                isBoundaryCentroid.add(isBoundary);
            }
            else{
                System.out.printf("Cant get vertices for centroids: %f %f \n", centroid[0], centroid[1]);
            }


        }
        return isBoundaryCentroid;
    }
    public  void applyExclusionMask() {
        int width = scaledUpProcessor.getWidth();
        int height = scaledUpProcessor.getHeight();

        // Iterate over all pixels in the image
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Check if the pixel is in the excluded region
                // Map the full-size pixel coordinates to the scaled-down coordinates
                int scaledX = x / SCALE;
                int scaledY = y / SCALE;
                if (exclusionMask[scaledX][scaledY]) {
                    // Set the pixel to black (RGB value: 0)
                    scaledUpProcessor.putPixel(x, y, 0);  // Black pixel (for grayscale or RGB images)
                }
            }
        }
    }
}
