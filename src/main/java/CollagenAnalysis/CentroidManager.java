package CollagenAnalysis;

import ij.ImagePlus;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.EDM;
import ij.plugin.filter.MaximumFinder;
import ij.plugin.filter.ParticleAnalyzer;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.solvers.LaguerreSolver;
import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import java.util.ArrayList;
import java.util.List;
import static CollagenAnalysis.ImageProcessingUtils.countNonZeroPixels;
import static CollagenAnalysis.ImageProcessingUtils.derivePolynomial;
import static CollagenAnalysis.ImageProcessingUtils.squarePolynomial;

public class CentroidManager {

    final double[] sigma = new double[20];

    public CentroidManager(ImageProcessor ip){
        this.ip = ip;
        generateSigmaArray();
    }
    ImageProcessor ip;


    public ArrayList<double[]> findCentroids(){
        double optimalSigma = findOptimalSigma();
        // double optimalSigma = 5;
        ByteProcessor bp = ip.duplicate().convertToByteProcessor();

        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                double value = bp.getPixelValue(x,y) / 255;
                bp.putPixelValue(x, y, 1 - value);
            }
        }

        EDM edm = new EDM();
        edm.toEDM(bp);

        //new ImagePlus("Euclidean distance map",bp);

        ImageProcessor fpFiltered  = bp.duplicate().convertToFloatProcessor();
        fpFiltered.blurGaussian(optimalSigma);

        MaximumFinder mf = new MaximumFinder();
        ByteProcessor maxima = mf.findMaxima(fpFiltered, 0.1, ImageProcessor.NO_THRESHOLD, MaximumFinder.SINGLE_POINTS, true, false);

        ResultsTable rt = new ResultsTable();
        int measurements = Measurements.CENTROID;

        // Explicitly set the threshold to 0-0 to treat only pixels with value 0 as foreground
        maxima.setThreshold(0, 0, ImageProcessor.NO_LUT_UPDATE);  // Setting threshold to 0-0

        ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE, measurements, rt, 0, Double.POSITIVE_INFINITY);
        int numberOfMaxima = countNonZeroPixels(maxima);

        ImagePlus t;
        pa.analyze(t = new ImagePlus(numberOfMaxima + " centroids at sigma " + optimalSigma, maxima));

        ArrayList<double[]> centroids = getCentroidsCords(maxima);

        return centroids;

    }
    private double findOptimalSigma(){
        double[] NumberOfCentroidsForEachSigmaIndex = findCentroidsForEachSigmaValue();
        // Fit a 5th-degree polynomial
        PolynomialCurveFitter fitter = PolynomialCurveFitter.create(5);
        WeightedObservedPoints obs = new WeightedObservedPoints();
        for (int i = 0; i < sigma.length; i++) {
            obs.add(sigma[i], NumberOfCentroidsForEachSigmaIndex[i]);
        }
        double[] originalCoefficients = fitter.fit(obs.toList());



        double[] derivativeCoefficients = derivePolynomial(originalCoefficients);
        double[] squaredDerivativeCoefficients = squarePolynomial(derivativeCoefficients);

        PolynomialFunction derivativeSquared = new PolynomialFunction(squaredDerivativeCoefficients);

        // Use a solver to find the roots
        LaguerreSolver solver = new LaguerreSolver();
        Complex[] complexRoots = solver.solveAllComplex(derivativeSquared.getCoefficients(), 0);

        List<Double> realRootsList = new ArrayList<>();
        for (Complex root : complexRoots) {
            //if (Math.abs(root.getImaginary()) < 1e-6) { // 1e-6 is a small tolerance value
            realRootsList.add(root.getReal());
            //}
        }

// Convert the list of real roots to a double array
        double[] realRoots = new double[realRootsList.size()];
        for (int i = 0; i < realRootsList.size(); i++) {
            realRoots[i] = realRootsList.get(i);
        }


        // Filter for positive real roots
        double sigmaOpt = Double.MAX_VALUE;
        for (double root : realRoots) {
            if (root > 0 && !Double.isNaN(root) && root < sigmaOpt) {
                sigmaOpt = root;
            }
        }
        System.out.println("Optimal Sigma:  " + sigmaOpt);
        return sigmaOpt;
    }
    private double[] findCentroidsForEachSigmaValue(){
        ByteProcessor bp2 = ip.duplicate().convertToByteProcessor();

        for (int y = 0; y < bp2.getHeight(); y++) {
            for (int x = 0; x < bp2.getWidth(); x++) {
                double value = bp2.getPixelValue(x,y) / 255;
                bp2.putPixelValue(x, y, 1 - value);
            }
        }

        EDM edm = new EDM();
        edm.toEDM(bp2);
        double[] NumberOfCentroidsForEachSigmaIndex = new double[sigma.length];
        for (int i = 0; i < sigma.length; i++) {
            ImageProcessor fpFiltered  = bp2.duplicate().convertToFloatProcessor();
            fpFiltered.blurGaussian(sigma[i]);

            MaximumFinder mf = new MaximumFinder();
            ByteProcessor maxima = mf.findMaxima(fpFiltered, 0.1, ImageProcessor.NO_THRESHOLD, MaximumFinder.SINGLE_POINTS, true, false);
            ResultsTable rt = new ResultsTable();
            int measurements = Measurements.CENTROID;

            // Explicitly set the threshold to 0-0 to treat only pixels with value 0 as foreground
            maxima.setThreshold(0, 0, ImageProcessor.NO_LUT_UPDATE);  // Setting threshold to 0-0
            ParticleAnalyzer pa = new ParticleAnalyzer(ParticleAnalyzer.SHOW_NONE, measurements, rt, 0, Double.POSITIVE_INFINITY);
            int numberOfMaxima = countNonZeroPixels(maxima);
            NumberOfCentroidsForEachSigmaIndex[i] = numberOfMaxima;
            pa.analyze(new ImagePlus(numberOfMaxima + " maxima at sigma " + sigma[i], maxima));

        }
        return NumberOfCentroidsForEachSigmaIndex;
    }
    private ArrayList<double[]> getCentroidsCords(ByteProcessor bp) {
        ArrayList<double[]> cords = new ArrayList<>();
        for (int y = 0; y < bp.getHeight(); y++) {
            for (int x = 0; x < bp.getWidth(); x++) {
                if(bp.getPixel(x, y) != 0){
                    double[] c = new double[2];
                    c[0] = x;
                    c[1] = y;
                    cords.add(c);
                }
            }
        }
        return cords;
    }

    private void generateSigmaArray(){
        for (int i = 0; i < sigma.length; i++) {
            sigma[i] = 0.5 * (i + 1);
        }
    }

}
