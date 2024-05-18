package CollagenAnalysis;

import org.apache.commons.math3.linear.*;

import java.util.ArrayList;

public class EllipseFitting {
    private static final double[] theta = createThetaArray();
    public static final int DEGREES = 361; // 0 to 360 degrees

    private ArrayList<double[]> xEllipses = new ArrayList<>();
    private ArrayList<double[]> yEllipses = new ArrayList<>();

    private ArrayList<double[]> radius_pix = new ArrayList<>();

    private ArrayList<double[]> mus = new ArrayList<>();


    private static double[] createThetaArray() {
        double[] theta = new double[361];
        for (int i = 0; i <= 360; i++) {
            theta[i] = Math.toRadians(i);
        }
        return theta;
    }

    public  void fitEllipses(ArrayList<double[]> fibrilPixels, int centroidLength, int[] cluster_idx/*boolean[] isBoundaryCentroid*/) {
        ArrayList<double[]> fibrilPixelsForThisCentroid = new ArrayList<>();
        int[] numberOfPixelsForEachCentroid = new int[centroidLength];

        ArrayList<double[][]> covariance = new ArrayList<>();
        for(int i = 0; i < centroidLength; i++) {
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
            mus.add(mu);

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

            covariance.add(cov);
        }



        // Fit ellipses
        for (int i = 0; i < centroidLength; i++) {
            double[] mu = mus.get(i);
            RealMatrix covMatrix = new Array2DRowRealMatrix(covariance.get(i));
            EigenDecomposition eig = new EigenDecomposition(covMatrix);

            double[] eigenValues = eig.getRealEigenvalues();
            RealMatrix V = eig.getV();

            double[] radiusPix = {
                    2 * Math.sqrt(Math.max(eigenValues[0], eigenValues[1])),
                    2 * Math.sqrt(Math.min(eigenValues[0], eigenValues[1]))
            };
            radius_pix.add(radiusPix);

            double[] xEllipse = new double[DEGREES];
            double[] yEllipse = new double[DEGREES];

            for (int theta = 0; theta < DEGREES; theta++) {
                //System.out.println("Degrees: " + theta);
                double radians = Math.toRadians(theta);
                double cosTheta = radiusPix[0] * Math.cos(radians);
                double sinTheta = radiusPix[1] * Math.sin(radians);

                //multiply RealMatrix by vector
                double[] ellipsePoint = V.operate(new double[]{cosTheta, sinTheta});
                xEllipse[theta] = ellipsePoint[0] + mu[0];
                yEllipse[theta] = ellipsePoint[1] + mu[1];
            }

            // Store or display the ellipse coordinates
            xEllipses.add(xEllipse);
            yEllipses.add(yEllipse);


//            System.out.printf("Ellipse %d coordinates (x, y):\n", i);
//            for (int j = 0; j < DEGREES; j++) {
//                System.out.printf("[%f, %f]\n", xEllipse[j], yEllipse[j]);
//            }
        }
    }

    public ArrayList<double[]> getxEllipses() {
        return xEllipses;
    }

    public ArrayList<double[]> getyEllipses() {
        return yEllipses;
    }
}
