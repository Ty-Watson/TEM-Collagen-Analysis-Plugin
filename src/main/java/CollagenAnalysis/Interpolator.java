package CollagenAnalysis;


import org.tinfour.common.Vertex;
import org.tinfour.interpolation.NaturalNeighborInterpolator;
import org.tinfour.standard.IncrementalTin;
import smile.neighbor.KDTree;
import smile.neighbor.Neighbor;

import java.util.ArrayList;
import java.util.List;

public class Interpolator {
    private List<Vertex> vertices = new ArrayList<Vertex>();
    private IncrementalTin tin;

    NaturalNeighborInterpolator interpolator;

    double[][] _centroids;

    KDTree<double[]> kdTree2;

    double[] _charIntensity;

    public Interpolator(List<double[]> centroids, double[] charIntensity){
        generateCentroidsForInterpolator(centroids, charIntensity);
        tin = new IncrementalTin();
        tin.add(vertices, null);
        interpolator = new NaturalNeighborInterpolator(tin);
        kdTree2 = KDTree.of(_centroids);
    }
    public double interpolate(double queryX, double queryY){
        double interpolatedValue = interpolator.interpolate(queryX, queryY, null);
        if(Double.isNaN(interpolatedValue)){
            Neighbor<double[], double[]> nearest = kdTree2.nearest(new double[]{queryX, queryY});
            return _charIntensity[nearest.index];

        }
        else{
            return interpolatedValue;
        }
    }

    private void generateCentroidsForInterpolator(List<double[]> centroids, double[] charIntensity){
        _centroids = centroids.toArray(new double[0][]);
        _charIntensity = charIntensity;
        for(int i = 0; i < centroids.size(); i++){
            Vertex v = new Vertex(centroids.get(i)[0], centroids.get(i)[1], charIntensity[i]);
            vertices.add(v);
        }
    }

}
