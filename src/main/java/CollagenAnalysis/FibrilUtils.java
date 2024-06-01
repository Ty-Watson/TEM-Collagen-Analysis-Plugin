package CollagenAnalysis;

import java.util.ArrayList;

public class FibrilUtils {
    public ArrayList<double[]> muForEachFibril = new ArrayList<>();
    public ArrayList<double[][]> covarianceForEachFibril = new ArrayList<>();
    public FibrilUtils(){

    }

    public void calculateMuAndCovarianceForEachFibril(ArrayList<double[]> fibrilPixels, ArrayList<double[]> centroids, int[] cluster_idx){
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
//            for(int t = 0; t < numberOfPixelsForEachCentroid[i]; t++){
//
//                varx  += Math.pow(fibrilPixelsForThisCentroid.get(t)[0] - mus.get(i)[0], 2);
//                vary  += Math.pow(fibrilPixelsForThisCentroid.get(t)[1] - mus.get(i)[1], 2);
//                covxy  += (fibrilPixelsForThisCentroid.get(t)[0] - mus.get(i)[0]) *  (fibrilPixelsForThisCentroid.get(t)[1] - mus.get(i)[1]);
//            }
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

            covarianceForEachFibril.add(cov);
        }
    }


}
