package CollagenAnalysis;

import ij.process.ImageProcessor;
import java.util.ArrayList;

public class FibrilUtils {
    private final double regularizationValue = 0.1;

    public ArrayList<double[]> muForEachFibril = new ArrayList<>();
    public ArrayList<double[][]> covarianceForEachFibril = new ArrayList<>();

    public ArrayList<double[]> problematicCentroids = new ArrayList<>();
    public ArrayList<double[]> _centroids = new ArrayList<>();
    //public ArrayList<RealMatrix> covariance = new ArrayList<>();
    public FibrilUtils(){

    }

    public void validateCovariance(ImageProcessor ip){
        for(int i = 0; i < covarianceForEachFibril.size(); i++){
            double[][] matrix = covarianceForEachFibril.get(i);
            boolean hasProblem = false;

            for(int row = 0; row < matrix.length; row++){
                for(int col = 0; col < matrix[row].length; col++){
                    double value = matrix[row][col];
                    if (Double.isNaN(value) || Double.isInfinite(value)) {
                        hasProblem = true;
                        break;
                    }
                }
                if(hasProblem)break;
            }
            if(hasProblem){
                problematicCentroids.add(_centroids.get(i));
                //remove problematic centroid;
                muForEachFibril.remove(i);
                covarianceForEachFibril.remove(i);
            }
        }
        if (problematicCentroids.isEmpty()) {
            System.out.println("No problematic centroids found.");
        } else {
            //generateImageWithProblematicCentroids(ip);
            System.out.println("Problematic centroids at coordinates:");
            for (double[] centroid : problematicCentroids) {
                System.out.println("Centroid: (" + centroid[0] + ", " + centroid[1] + ")");
            }
        }
    }

    public void calculateMuAndCovarianceForEachFibril(ArrayList<double[]> fibrilPixels, ArrayList<double[]> centroids, int[] cluster_idx){
        _centroids = centroids;
        ArrayList<double[]> fibrilPixelsForThisCentroid = new ArrayList<>();
        int[] numberOfPixelsForEachCentroid = new int[centroids.size()];

        ArrayList<double[][]> covariance = new ArrayList<>();
        for(int i = 0; i < centroids.size(); i++) {
            // Reset the temporary storage for the next centroid
            fibrilPixelsForThisCentroid.clear();
            for (int j = 0; j < fibrilPixels.size(); j++) {
                //Loop through all fibril pixels and find the ones that belong to this centroid.
                if (cluster_idx[j] == i) {
                    //count number of pixels for this centroid
                    numberOfPixelsForEachCentroid[i] += 1;
                    //add the fibril pixel coordinates to this centroid
                    fibrilPixelsForThisCentroid.add(fibrilPixels.get(j));
                }


            }
            //  sums for the current centroid
            double xsum = 0;
            double ysum = 0;

            // Sum up x and y coordinates for all pixels associated with the current centroid
            for (double[] pixel : fibrilPixelsForThisCentroid) {
                xsum += pixel[0];
                ysum += pixel[1];
            }

            // Calculate and store the mean (mu) for the current centroid
            double[] mu = new double[2];
            mu[0] = xsum / numberOfPixelsForEachCentroid[i];
            mu[1] = ysum / numberOfPixelsForEachCentroid[i];
            muForEachFibril.add(mu);

            // Variables for calculating the variances and covariance
            double varx = 0;
            double vary = 0;
            double covxy = 0;

            // Calculate variances and covariance for the current centroid
            for (double[] pixel : fibrilPixelsForThisCentroid) {
                double dx = pixel[0] - mu[0];
                double dy = pixel[1] - mu[1];
                varx += dx * dx;
                vary += dy * dy;
                covxy += dx * dy;
            }
            // Normalize the variances and covariance by the number of pixels
            varx = varx / numberOfPixelsForEachCentroid[i];
            vary = vary / numberOfPixelsForEachCentroid[i];
            covxy = covxy / numberOfPixelsForEachCentroid[i];

            // Store the covariance matrix for the current centroid
            double[][] cov = new double[2][2];
            cov[0][0] = varx;
            cov[0][1] = covxy;
            cov[1][0] = covxy;
            cov[1][1] = vary;

            // Add regularization value to the diagonal elements
            cov[0][0] += regularizationValue;
            cov[1][1] += regularizationValue;

            covarianceForEachFibril.add(cov);
        }
    }
}
