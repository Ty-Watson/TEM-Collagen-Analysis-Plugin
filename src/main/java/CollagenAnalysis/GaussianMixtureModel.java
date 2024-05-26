package CollagenAnalysis;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ProgressBar;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.scijava.ui.StatusBar;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

public class GaussianMixtureModel {
    private final double regularizationValue = 0.1;

    private final int centroidLength;
    private double[][] centroids;
    private final double[][] fibrilPixelsThinned;
    private double[] centroidsIdxThinned;
    private double[][] mus;
    private double[][][] covariance;

    // Initialize a matrix to store the probability density function (PDF) values
    public double[][] pdf_values;
    // Initialize a matrix to store the posterior probabilities (responsibilities) that each pixel belongs to every centroid
    private double[][] p;

    public double[] componentProportions;

    private int[] numberOfPixelsForEachCentroid;
    private ArrayList<double[]> fibrilPixelsForThisCentroid = new ArrayList<double[]>();

    public GaussianMixtureModel(int centroidLength, double[][] fibrilPixelsThinned, double[] centroidIdx){
        this.centroidLength = centroidLength;
        this.fibrilPixelsThinned = fibrilPixelsThinned;
        this.centroidsIdxThinned = centroidIdx;
        pdf_values = new double[this.centroidLength][fibrilPixelsThinned.length];
        p = new double[this.centroidLength][fibrilPixelsThinned.length];
        numberOfPixelsForEachCentroid = new int[this.centroidLength];
        componentProportions = new double[this.centroidLength];
        mus = new double[this.centroidLength][2];
        covariance = new double[this.centroidLength][2][2];
    }
    public void EM(boolean combined, int totalIterations){
        TimeKeeper totalTimer = new TimeKeeper();
        totalTimer.start();

        initMuAndCov();
        initPdfValuesForEachFibrilPixel();
        initPosteriorProbabilityMatrix();
        // First phase of EM algorithm: Adjust mu and ComponentProportion, keeping Sigma fixed
        //Expectation Step (E-step): Assigns "responsibilities" to each data point based on the current parameters of the Gaussians, reflecting the probability that a data point belongs to each cluster.
        //Maximization Step (M-step): Updates the parameters of the Gaussians (mean, covariance, and mixture coefficients) to maximize the likelihood of the data given these responsibilities.
        TimeKeeper firstPhase = new TimeKeeper();
        firstPhase.start();

        if(combined){
            for(int k = 0; k < totalIterations; k++){
                System.out.println("iteration step " + k);
                //progressBar.setStatus("Fitting Gaussian mixture model");

                // E-step: Recalculate responsibilities
                //use the initial guess of prob matrix first because component proportions have not been calculated yet
                if(k != 0)
                    calculatePosteriorProbabilityMatrix();
                // M step: Accumulate sums for updating mu
                mStepCombined(k, totalIterations);
                // Update the progress bar after each iteration
                //progressBar.setProgress(k, totalIterations);
            }
        }
        else{
            for(int k = 0; k < totalIterations / 2; k++){
                System.out.println("First iteration step " + k);
                //progressBar.setStatus("Fitting Gaussian mixture model");
                IJ.showStatus("Fitting Gaussian mixture model");

                // E-step: Recalculate responsibilities
                //use the initial guess of prob matrix first because component proportions have not been calculated yet
                if(k != 0)
                    calculatePosteriorProbabilityMatrix();
                // M step: Accumulate sums for updating mu
                mStep(k, totalIterations, true);
                // Update the progress bar after each iteration
                //progressBar.setProgress(k, totalIterations);
                IJ.showProgress(k, totalIterations);

            }
            firstPhase.stop();
            System.out.println("First phase duration: " + firstPhase.getElapsedTimeFormatted());

            // Second phase of EM algorithm: Adjust Sigma and ComponentProportion, keeping mu fixed
            for(int k = totalIterations / 2; k < totalIterations + 1; k++) {
                System.out.println("Second iteration step " + k);
                // Reset sums for each iteration
                double[] sumP = new double[centroidLength];

                // E-step: Recalculate responsibilities
                calculatePosteriorProbabilityMatrix();

                // M-step: Update covariance matrices
                mStep(k, totalIterations, false);

                // Update the progress bar after each iteration
                //progressBar.setProgress(k, totalIterations);
                IJ.showProgress(k, totalIterations);
            }
        }

        // Complete the progress bar
        //progressBar.setStatus("Fitting Complete");
        IJ.showStatus("Fitting Gaussian Mixture Model Complete");
        totalTimer.stop();
        System.out.println("Total EM time: " + totalTimer.getElapsedTimeFormatted());
        //calculatePosteriorProbabilityMatrix();
    }

    public void mStep(int k, int totalIterations,  boolean adjustMu){
        double[] sumP = new double[centroidLength];
        double[] sumPxX = new double[centroidLength];
        double[] sumPxY = new double[centroidLength];


        // M step: Accumulate sums for updating mu
        IntStream.range(0, centroidLength).parallel().forEach(i -> {
            double[][] cov = new double[2][2];
//            double[] muForCentroid = mus.get(i);
//            double[][] covForCentroid = covariance.get(i);

            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                sumP[i] += p[i][j];
                if(adjustMu){
                    sumPxX[i] += p[i][j] * fibrilPixelsThinned[j][0];
                    sumPxY[i] += p[i][j] * fibrilPixelsThinned[j][1];
                }
                else{
                    double dx = fibrilPixelsThinned[j][0] - mus[i][0];
                    double dy = fibrilPixelsThinned[j][1] - mus[i][1];

                    cov[0][0] += p[i][j] * dx * dx; // Variance in X
                    cov[1][1] += p[i][j] * dy * dy; // Variance in Y
                    cov[0][1] += p[i][j] * dx * dy; // Covariance XY
                    cov[1][0] = cov[0][1]; // Symmetric
                }
            }
            componentProportions[i] = sumP[i] / fibrilPixelsThinned.length;
            // Update mu for each centroid
            if(adjustMu){
                double[] updatedMu = new double[]{sumPxX[i] / sumP[i], sumPxY[i] / sumP[i]};
                mus[i] = updatedMu;

                if(k == totalIterations / 2 -1)
                    System.out.printf("Updated Mu for centroid %d: [%f, %f] \n", i, updatedMu[0], updatedMu[1]);
            }
            else{
                // Add regularization value to the diagonal elements
                cov[0][0] += regularizationValue;
                cov[1][1] += regularizationValue;

                // Normalize and update the covariance matrix for this centroid
                cov[0][0] /= sumP[i];
                cov[1][1] /= sumP[i];
                cov[0][1] /= sumP[i];
                cov[1][0] = cov[0][1];
                covariance[i] =  cov;
                if(k == totalIterations)
                    System.out.printf("updated Covariance for centroid %d: [%f, %f][%f, %f] \n", i, cov[0][0], cov[0][1], cov[1][0], cov[1][1]);

            }

            //update pdf values
            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                pdf_values[i][j] = mvnpdf(fibrilPixelsThinned[j], mus[i], covariance[i]);
            }
        });




    }

    public void mStepCombined(int k, int totalIterations){
        double[] sumP = new double[centroidLength];
        double[] sumPxX = new double[centroidLength];
        double[] sumPxY = new double[centroidLength];


        // M step: Accumulate sums for updating mu
        for(int i = 0; i < centroidLength; i++){

            double[][] cov = new double[2][2];
//            double[] muForCentroid = mus.get(i);
//            double[][] covForCentroid = covariance.get(i);

            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                sumP[i] += p[i][j];

                sumPxX[i] += p[i][j] * fibrilPixelsThinned[j][0];
                sumPxY[i] += p[i][j] * fibrilPixelsThinned[j][1];


                double dx = fibrilPixelsThinned[j][0] - mus[i][0];
                double dy = fibrilPixelsThinned[j][1] - mus[i][1];

                cov[0][0] += p[i][j] * dx * dx; // Variance in X
                cov[1][1] += p[i][j] * dy * dy; // Variance in Y
                cov[0][1] += p[i][j] * dx * dy; // Covariance XY
                cov[1][0] = cov[0][1]; // Symmetric

            }
            componentProportions[i] = sumP[i] / fibrilPixelsThinned.length;
            // Update mu for each centroid

            double[] updatedMu = new double[]{sumPxX[i] / sumP[i], sumPxY[i] / sumP[i]};
            mus[i] = updatedMu;

            if(k == totalIterations -1)
                System.out.printf("Updated Mu for centroid %d: [%f, %f] \n", i, updatedMu[0], updatedMu[1]);


            // Add regularization value to the diagonal elements
            cov[0][0] += regularizationValue;
            cov[1][1] += regularizationValue;

            // Normalize and update the covariance matrix for this centroid
            cov[0][0] /= sumP[i];
            cov[1][1] /= sumP[i];
            cov[0][1] /= sumP[i];
            cov[1][0] = cov[0][1];
            covariance[i] = cov;
            if(k == totalIterations - 1)
                System.out.printf("updated Covariance for centroid %d: [%f, %f][%f, %f] \n", i, cov[0][0], cov[0][1], cov[1][0], cov[1][1]);



            //update pdf values
            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                pdf_values[i][j] = mvnpdf(fibrilPixelsThinned[j], mus[i], covariance[i]);
            }
        }
    }


    private void initMuAndCov(){

        for(int i = 0; i < centroidLength; i++){
            // Reset the temporary storage for the next centroid
            fibrilPixelsForThisCentroid.clear();
            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                //Loop through all fibril pixels and find the ones that belong to this centroid.
                if(centroidsIdxThinned[j] == i){
                    //count number of pixels for this centroid
                    numberOfPixelsForEachCentroid[i] += 1;
                    //add the fibril pixel coordinates to this centroid
                    fibrilPixelsForThisCentroid.add(fibrilPixelsThinned[j]);
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
            mus[i] = mu;

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
            for(double[] pixel : fibrilPixelsForThisCentroid){
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

            System.out.printf("Initial covxy, varx, vary guess for centroid %d: %f, %f, %f \n", i, covxy, varx, vary);
            covariance[i] = cov;
        }
    }

    //parallel proccessing
    private void initPdfValuesForEachFibrilPixel(){
        // Calculate the PDF values for each thinned fibril pixel across all centroids
        IntStream.range(0 , centroidLength).parallel().forEach( i -> {
            double[] muForCentroid = mus[i];
            double[][] covForCentroid = covariance[i];
            for(int j = 0; j < fibrilPixelsThinned.length; j++){
                pdf_values[i][j] = mvnpdf(fibrilPixelsThinned[j], muForCentroid,covForCentroid);
            }
        });
//        for(int i = 0; i < centroidLength; i++){
//            double[] muForCentroid = mus.get(i);
//            double[][] covForCentroid = covariance.get(i);
//            for(int j = 0; j < fibrilPixelsThinned.length; j++){
//                pdf_values[i][j] = mvnpdf(fibrilPixelsThinned[j], muForCentroid,covForCentroid);
//            }
//        }
    }
    private void initPosteriorProbabilityMatrix(){
        for(int j = 0; j < fibrilPixelsThinned.length; j++){
            double total = 0;
            // Calculate the denominator for the responsibility formula (sum over all centroids)
            for(int i = 0; i < centroidLength; i++){
                total += numberOfPixelsForEachCentroid[i] *  pdf_values[i][j];
            }
            // Calculate and store the responsibilities
            for(int i = 0 ; i < centroidLength; i++){
                if(total != 0)
                    p[i][j] = numberOfPixelsForEachCentroid[i] * pdf_values[i][j] / total;
            }
        }
    }
    public double[][][] calculatePosterProbabilityMatrixAll(ArrayList<double[]> allPixels, int height, int width){
        long startTime = System.currentTimeMillis();  // Start timing
        //progressBar.setStatus("Calculating posterior probability for all image pixels");
        IJ.showStatus("Calculating posterior probability for all image pixels");
        int numCentroids = mus.length;
        double[][][] posteriors = new double[numCentroids][height][width];
        double[] totalProbabilities = new double[allPixels.size()];
        double[][] pdf_values = new double[centroidLength][allPixels.size()];

        // Parallel computation of PDF values
        IntStream.range(0, centroidLength).parallel().forEach(i -> {
            double[] mu = mus[i];
            double[][] sigma = covariance[i];
            for (int j = 0; j < allPixels.size(); j++) {
                pdf_values[i][j] = mvnpdf(allPixels.get(j), mu, sigma);
            }
        });

        // Calculate the total probability for each pixel over all centroids
        //Arrays.parallelSetAll(totalProbabilities, i -> 0.0); // Initialize with zeros using parallel operation

        IntStream.range(0, numCentroids).parallel().forEach(i -> {
            double[] mu = mus[i];
            double[][] sigma = covariance[i];
            for (int j = 0; j < allPixels.size(); j++) {
                //double pdf = mvnpdf(allPixels.get(j), mu, sigma);
                int y = j / width;
                int x = j % width;

                //not sure what synchronized does remove this
                synchronized (posteriors) {
                    posteriors[i][y][x] = componentProportions[i] * pdf_values[i][j];
                    totalProbabilities[j] += componentProportions[i] * pdf_values[i][j];
                }
            }
        });

        // Normalize the probabilities to get posteriors
        for (int i = 0; i < numCentroids; i++) {
            for (int j = 0; j < allPixels.size(); j++) {
                int y = j / width;
                int x = j % width;
                if (totalProbabilities[j] != 0) {
                    posteriors[i][y][x] = componentProportions[i] * pdf_values[i][j] / totalProbabilities[j];
                }
            }
        }
        //progressBar.setStatus("Calculating posterior probability complete");
        long endTime = System.currentTimeMillis();  // End timing
        long duration = endTime - startTime;

        // Print the elapsed time in minutes and seconds
        long minutes = (duration / 1000) / 60;
        long seconds = (duration / 1000) % 60;
        if (minutes > 0) {
            System.out.printf("Posterior Prob Time: %d minutes, %d seconds\n", minutes, seconds);
        } else {
            System.out.printf("Posterior Prob Time: %d seconds\n", seconds);
        }

        return posteriors;

    }
    public double[][] calculatePosterProbabilityMatrixFibril(ArrayList<double[]> allPixels, int height, int width){
        long startTime = System.currentTimeMillis();  // Start timing
        //progressBar.setStatus("Calculating posterior probability for all fibril pixels");
        IJ.showStatus("Calculating posterior probability for all fibril pixels");
        int numCentroids = mus.length;
        double[][] posteriors = new double[numCentroids][allPixels.size()];
        double[] totalProbabilities = new double[allPixels.size()];
        double[][] pdf_values = new double[centroidLength][allPixels.size()];

        // Parallel computation of PDF values
        IntStream.range(0, centroidLength).parallel().forEach(i -> {
            double[] mu = mus[i];
            double[][] sigma = covariance[i];
            for (int j = 0; j < allPixels.size(); j++) {
                pdf_values[i][j] = mvnpdf(allPixels.get(j), mu, sigma);
            }
        });

        // Calculate the total probability for each pixel over all centroids
        for (int i = 0; i < numCentroids; i++) {
            for (int j = 0; j < allPixels.size(); j++) {
                posteriors[i][j] = componentProportions[i] * pdf_values[i][j];
                totalProbabilities[j] += posteriors[i][j];
            }
        }

        // Normalize the probabilities to get posteriors
        for (int i = 0; i < numCentroids; i++) {
            for (int j = 0; j < allPixels.size(); j++) {
                if (totalProbabilities[j] != 0) {
                    posteriors[i][j] = componentProportions[i] * pdf_values[i][j] / totalProbabilities[j];
                }
            }
        }

        //progressBar.setStatus("Calculating posterior probability for all fibril pixels complete");
        IJ.showStatus("Calculating posterior probability for all fibril pixels complete");
        long endTime = System.currentTimeMillis();  // End timing
        long duration = endTime - startTime;

        // Print the elapsed time in minutes and seconds
        long minutes = (duration / 1000) / 60;
        long seconds = (duration / 1000) % 60;
        if (minutes > 0) {
            System.out.printf("Posterior Prob Time for all fibril pixels: %d minutes, %d seconds\n", minutes, seconds);
        } else {
            System.out.printf("Posterior Prob Time for all fibril pixels: %d seconds\n", seconds);
        }

        return posteriors;

    }
    private void calculatePosteriorProbabilityMatrix(){
        for(int j = 0; j < fibrilPixelsThinned.length; j++) {
            double total = 0;
            // Calculate the denominator for the responsibility formula (sum over all centroids)
            for(int i = 0; i < centroidLength; i++) {
                total += componentProportions[i] * pdf_values[i][j]; // Total probability for normalization
            }
            // Update responsibilities
            for(int i = 0 ; i < centroidLength; i++) {
                if(total > 0){
                    p[i][j] = componentProportions[i] * pdf_values[i][j] / total; // Normalized responsibility
                }
                else{
                    p[i][j] = 0;
                }

            }
        }
    }
    private static double mvnpdf(double[] x, double[] mu, double[][] sigma) {
        int k = mu.length;
        RealMatrix xMatrix = new Array2DRowRealMatrix(x);
        RealMatrix muMatrix = new Array2DRowRealMatrix(mu);
        RealMatrix sigmaMatrix = new Array2DRowRealMatrix(sigma);

        // Calculate determinant and inverse of sigma
        double sigmaDet = new org.apache.commons.math3.linear.LUDecomposition(sigmaMatrix).getDeterminant();
        RealMatrix sigmaInverse = new org.apache.commons.math3.linear.LUDecomposition(sigmaMatrix).getSolver().getInverse();

        // Calculate the log of the normalization factor
        double logNormalizationFactor = -0.5 * (k * Math.log(2 * Math.PI) + Math.log(sigmaDet));

        // Calculate the exponent part
        RealMatrix diff = xMatrix.subtract(muMatrix);
        double exponent = -0.5 * diff.transpose().multiply(sigmaInverse).multiply(diff).getEntry(0, 0);

        // Calculate the log of the PDF value
        double logPdfValue = logNormalizationFactor + exponent;

        // Return the PDF in the log domain or exponentiate if necessary, checking for underflow
        double result = Math.exp(logPdfValue);
        return result;

//        try{
//            MultivariateNormalDistribution md = new MultivariateNormalDistribution(mu, sigma);
//            double pdf_value = md.density(x);
//            return pdf_value;
//        } catch(Exception e){
//            System.out.printf("[*] pdf calculation failed for mu: [%f, %f] and covariance: [[%f, %f][%f, %f]]", mu[0], mu[1], sigma[0][0], sigma[0][1], sigma[1][0], sigma[1][1]);
//            return  0;
//        }

    }

    public double[][] getPosteriorProbabilityMatrix() {
        return p;
    }
    public double[][] getMus(){return mus;}

    public double[][][] getCovariance() {
        return covariance;
    }

    public double[] getComponentProportions() {
        return componentProportions;
    }
}
