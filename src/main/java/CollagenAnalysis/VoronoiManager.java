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

import smile.classification.KNN;
import smile.neighbor.KDTree;
import smile.neighbor.Neighbor;

import java.awt.*;
import java.awt.event.AWTEventListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Path2D;
import java.util.*;
import java.util.List;
import java.util.stream.Stream;

import static CollagenAnalysis.ImageProcessingUtils.findNearestCentroid;


public class VoronoiManager {
    public VoronoiManager(ImagePlus imp, boolean[][] exclusionMask){

        ip = imp.getProcessor().duplicate();
        this.exclusionMask = exclusionMask;
        applyExclusionMask();
        this.imp = new ImagePlus("Voronoi", ip);
        this.imp.show();

        setCanvasEvents();
    }
    private ImagePlus imp;
    private ImageProcessor ip;
    public Voronoi voronoi;
    public boolean[][] exclusionMask;
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
                cells.computeIfAbsent(edge.getSite1(), k -> new HashSet<>()).add(edge.getA());
                cells.computeIfAbsent(edge.getSite1(), k -> new HashSet<>()).add(edge.getB());
                cells.computeIfAbsent(edge.getSite2(), k -> new HashSet<>()).add(edge.getA());
                cells.computeIfAbsent(edge.getSite2(), k -> new HashSet<>()).add(edge.getB());
            }
        }

        List<Point> centroids = new ArrayList<>();

        cells.forEach((point, vertices) -> {
            // System.out.println("Cell centered at: " + point);
            centroids.add(point);
            //System.out.println("Vertices:");
            //vertices.forEach(vertex -> System.out.println(vertex));
        });
        voronoi_centroids = centroids;
    }

    public void drawVoronoi(){
        ip.snapshot();
        Collection<Point> points = new ArrayList<>();
        for( double[] centroid : centroids){
            de.alsclo.voronoi.graph.Point point = new de.alsclo.voronoi.graph.Point(centroid[0], centroid[1]);
            points.add(point);
        }
        voronoi = new Voronoi(points);
        Graph g = voronoi.getGraph();
        generateCells();


        //ImageProcessor ip = imp.getProcessor().duplicate();

        ip.setColor(Color.BLUE);


        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            Point a = e.getA().getLocation();
            Point b = e.getB().getLocation();
            ip.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });
        imp.updateAndDraw();
        overlayMaxima(centroids, Color.red.getRGB());
        applyExclusionMask();
    }
    public void initDrawVoronoi(){
        ip.snapshot();
        Collection<Point> points = new ArrayList<>();
        for( double[] centroid : centroids){
            Point point = new Point(centroid[0], centroid[1]);
            points.add(point);
        }
        voronoi = new Voronoi(points);
        init_voronoi = voronoi;
        Graph g = voronoi.getGraph();
        generateCells();


        //ImageProcessor ip = imp.getProcessor().duplicate();

        ip.setColor(Color.BLUE);


        // Use forEach to iterate through the stream and draw each edge
        g.edgeStream().filter(e -> e.getA() != null && e.getB() != null).forEach(e -> {
            Point a = e.getA().getLocation();
            Point b = e.getB().getLocation();
            ip.drawLine((int) a.x,  (int)a.y, (int) b.x, (int) b.y);
        });

        //imageManager.drawExcludedRegions(ip);

        imp.updateAndDraw();
        overlayMaxima(centroids, Color.red.getRGB());
        applyExclusionMask();
    }
    private void overlayMaxima(List<double[]> centroids, int color) {
        // Convert the binarized image to RGB to overlay red points
        ImageProcessor d = ip.duplicate();
        ColorProcessor colorProcessor = d.convertToColorProcessor();

        // Define the red color
        //int red = (255 << 16); // RGB value for red

        // Size of the maxima point to be drawn (radius of the circle around the maxima point)
        int pointSize = 2; // This can be adjusted based on how big you want the maxima points to be

        // Draw each maximum point as a red circle on the RGB image
        for (double[] point : centroids) {
            double x = point[0];
            double y = point[1];

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
        imp.setProcessor(colorProcessor);
    }
    public void setCentroids(ArrayList<double[]> centroids){
        this.centroids = centroids;
    }
    public void showImg(){
        imp.show();
    }

    private void setCanvasEvents(){
        canvas = imp.getCanvas();

//        canvas.addMouseListener(new MouseAdapter() {
//            @Override
//            public void mouseClicked(MouseEvent e) {
//
//                if(e.getButton() == MouseEvent.BUTTON1){   //left click add centroid
//                    int clickX = canvas.offScreenX(e.getX());
//                    int clickY = canvas.offScreenY(e.getY());
//
//                    double[] newCentroid = new double[]{ clickX, clickY};
//
//                    centroids.add(newCentroid);
//                    ip.reset();
//                    drawVoronoi();
//                }
//                else if(e.getButton() == MouseEvent.BUTTON3){ //right click remove centroid
//                    int clickX = canvas.offScreenX(e.getX());
//                    int clickY = canvas.offScreenY(e.getY());
//
//                    // Find and remove the nearest centroid to the click
//                    double[] nearestCentroid = findNearestCentroid(centroids, clickX, clickY);
//                    centroids.remove(nearestCentroid);
//                    ip.reset();
//                    drawVoronoi();
//                }
//
//
//            }
//
//        });
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
                            int clickX = canvas.offScreenX(mouseEvent.getX());
                            int clickY = canvas.offScreenY(mouseEvent.getY());


                            // Find and remove the nearest centroid to the click
                            double[] nearestCentroid = findNearestCentroid(centroids, clickX, clickY);
                            if (nearestCentroid != null) {
                                centroids.remove(nearestCentroid);
                                ip.reset();
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
                    int clickX = canvas.offScreenX(e.getX());
                    int clickY = canvas.offScreenY(e.getY());

                    //cant add centroid in exclusion mask
                    if(exclusionMask[clickX][clickY]){
                        return;
                    }


                    double[] newCentroid = new double[]{ clickX, clickY};
                    centroids.add(newCentroid);
                    ip.reset();
                    drawVoronoi();
                }
            }
        });

    }


    private boolean belongsToCell(Point pixel, Point centroid){
        Set<Vertex> verticesForVoronoiCell = cells.get(centroid);
        List<Point> vertices = new ArrayList<>();
        Path2D.Double path = new Path2D.Double();

        for(Vertex vertex : verticesForVoronoiCell){
            Point point = vertex.getLocation();
            vertices.add(point);
        }
        path.moveTo(vertices.get(0).x, vertices.get(0).y);
        vertices.stream().skip(1).forEach(v -> path.lineTo(v.x, v.y));
        path.closePath(); //
        return path.contains(pixel.x, pixel.y);
    }
    public double[] computeIntensities(List<Double> pixelIntensities, double[] distances, int[] indices, int numCentroids) {
        Map<Integer, List<Double>> centroidDistances = new HashMap<>();
        Map<Integer, List<Double>> centroidIntensities = new HashMap<>();

        // Collect distances and intensities for each centroid
        for (int i = 0; i < indices.length; i++) {
            centroidDistances.computeIfAbsent(indices[i], k -> new ArrayList<>()).add(distances[i]);
            centroidIntensities.computeIfAbsent(indices[i], k -> new ArrayList<>()).add(pixelIntensities.get(i));
        }

        double[] charCellIntensity = new double[numCentroids];
        // Calculate the mean intensity for pixels within the 5th percentile distance for each centroid
        centroidDistances.forEach((centroid, distList) -> {
            double fifthPercentile = percentile(distList, .05);
            List<Double> intensities = centroidIntensities.get(centroid);
            double sum = 0;
            int count = 0;
            for (int i = 0; i < distList.size(); i++) {
                if (distList.get(i) < fifthPercentile) {
                    sum += intensities.get(i);
                    count++;
                }
            }
            charCellIntensity[centroid] = (count > 0) ? sum / count : Double.NaN; // Compute mean
        });

        return charCellIntensity;
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

    private double percentile(List<Double> values, double percentile) {
        Collections.sort(values);
        int index = (int) Math.ceil(percentile
                * values.size());
        return values.get(index - 1);
    }

    private double distance(Point p1, Point p2) {
        return Math.sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y));
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

    public ArrayList<double[]> identifyBoundaryFibrils(boolean[][] exclusionMask, ArrayList<double[]> pixels) {
//        List<Vertex> voronoiVertices = new ArrayList<>();
//        List<List<Integer>> voronoiCells = new ArrayList<>();
//
//        // Extract Voronoi vertices and cells from the Voronoi diagram
//        for (Point centroid : voronoi_centroids) {
//            Set<Vertex> cellVertices = cells.get(centroid);
//            List<Integer> cellIndices = new ArrayList<>();
//            for (Vertex vertex : cellVertices) {
//                voronoiVertices.add(vertex);
//                cellIndices.add(voronoiVertices.size() - 1);
//            }
//            voronoiCells.add(cellIndices);
//        }
//
//        // Perform nearest neighbor search for Voronoi vertices
//        int[] nearestPixelIndices = nearestNeighborSearch(pixels, voronoiVertices);


        boolean[] isBoundaryCentroid = new boolean[voronoi_centroids.size()];
        ArrayList<double[]> nonBoundryCentroids = new ArrayList<>();
        for(double[] centroid: centroids){
            Point centroidPont = new Point(centroid[0], centroid[1]);
            Set<Vertex> cellVertices = cells.get(centroidPont);
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
        return nonBoundryCentroids;
//        for (int i = 0; i < voronoiCells.size(); i++) {
//            List<Integer> cell = voronoiCells.get(i);
//            boolean isBoundary = false;
//
//            for (int idx : cell) {
//                Vertex vertex = voronoiVertices.get(idx);
//                Point locationPoint = vertex.getLocation();
//                double[] location = new double[]{locationPoint.x, locationPoint.y};
//                if (location[0] < 0 || location[1] < 0 || location[0] >= ip.getWidth() || location[1] >= ip.getHeight()) {
//                    isBoundary = true;
//                    break;
//                }
//
////                // Check exclusion mask
////                int nearestPixelIndex = nearestPixelIndices[idx];
////                double[] nearestPixel = pixels.get(nearestPixelIndex);
////                if (exclusionMask[(int) nearestPixel[1]][(int) nearestPixel[0]]) {
////                    isBoundary = true;
////                    break;
////                }
//            }
//
//            isBoundaryCentroid[i] = isBoundary;
//        }
//
//        // Output results
//        for (int i = 0; i < isBoundaryCentroid.length; i++) {
//            System.out.println("Centroid " + i + " is boundary: " + isBoundaryCentroid[i]);
//        }
        //return isBoundaryCentroid;
    }
    public ArrayList<Boolean> identifyBoundaryFibrilsBoolean(boolean[][] exclusionMask, ArrayList<double[]> pixels) {

        ArrayList<Boolean> isBoundaryCentroid = new ArrayList<Boolean>();
        for(Point centroid: voronoi_centroids){
            Set<Vertex> cellVertices = cells.get(centroid);
            boolean isBoundary = false;

            for(Vertex v : cellVertices){
                Point locationPoint = v.getLocation();
                double[] location = new double[]{locationPoint.x, locationPoint.y};
                if (location[0] < 0 || location[1] < 0 || location[0] >= ip.getWidth() || location[1] >= ip.getHeight()) {
                    isBoundary = true;
                    break;
                }
            }
           isBoundaryCentroid.add(isBoundary);
        }
        return isBoundaryCentroid;
//
    }
    private int[] nearestNeighborSearch(ArrayList<double[]> pixels, List<Vertex> voronoiVertices) {


        // Convert the list of pixels to an array for KDTree
        double[][] pixelArray = pixels.toArray(new double[0][]);

        // Create a KD-tree for the pixels
        KDTree<double[]> pixelTree = KDTree.of(pixelArray);

        // Perform the nearest neighbor search for Voronoi vertices
        int[] nearestPixelIndices = new int[voronoiVertices.size()];
        for (int i = 0; i < voronoiVertices.size(); i++) {
            Point vertexLocationPoint = voronoiVertices.get(i).getLocation();
            double[] vertexLocation = new double[]{vertexLocationPoint.x, vertexLocationPoint.y};
            Neighbor<double[], double[]> nearest = pixelTree.nearest(vertexLocation);
            nearestPixelIndices[i] = nearest.index;
        }

        // Output or use the nearest pixel indices
        for (int i = 0; i < nearestPixelIndices.length; i++) {
            System.out.println("Voronoi Vertex " + i + " nearest pixel index: " + nearestPixelIndices[i]);
        }
        return nearestPixelIndices;
    }

    public  void applyExclusionMask() {
        int width = ip.getWidth();
        int height = ip.getHeight();

//        if (exclusionMask.length != height || exclusionMask[0].length != width) {
//            throw new IllegalArgumentException("Exclusion mask size does not match image dimensions.");
//        }

        // Iterate over all pixels in the image
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                // Check if the pixel is in the excluded region
                if (exclusionMask[x][y]) {
                    // Set the pixel to black (RGB value: 0)
                    ip.putPixel(x, y, 0);  // Black pixel (for grayscale or RGB images)
                }
            }
        }
    }

}
