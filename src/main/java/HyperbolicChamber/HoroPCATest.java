package HyperbolicChamber;

import static HyperbolicChamber.HoroPCAValidator.printTotalTime;
import java.util.List;

/**
 *
 * @author sean phillips
 */
public class HoroPCATest {
    public static void main(String[] args) {
        int inputDim = 100;
        int targetDim = 3;
        int maxIterations = 300;
        int numClusters = 4;
        int pointsPerCluster = 50;
        int seed = 42;
        double learningRate = 0.01;
        double tolerance = 1e-8;
        System.out.println("Generating synthetic clustered data.");
        long startTime = System.nanoTime();
        // Step 1: Generate synthetic clustered data on the Poincaré ball
        PoincareBallFactory factory = new PoincareBallFactory(inputDim, seed);
//        List<VectorN> isotropicDataset = factory.generateIsotropicClusteredData(
//                numClusters, pointsPerCluster, 0.05);
//        // Step 2: Train HoroPCA model
//        HoroPCA horo = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
//            HoroPCA.InitializationStrategy.KMEANS_PLUS_PLUS, 
//                HyperbolicUtils.hyperbolicSquaredDist, seed);
//        horo.verbose = true;
//        // Step 3: Transform and inverse transform for reconstruction
//        List<VectorN> transformed = horo.fitTransform(isotropicDataset);
//        List<VectorN> reconstructed = transformed.stream().map(v -> {
//            VectorN r = horo.inverseTransform(v, maxIterations, learningRate);
//            return r;
//        }).toList();
//        // Step 4: Print explained variance and reconstruction error
//        double reconError = HoroPCAValidator.reconstructionError(isotropicDataset, reconstructed);
//        System.out.printf("Mean Reconstruction Error: %.6f\n", reconError);


        // Repeat for anisotropic set
        double[] anisotropicSpread = new double[inputDim];
        for (int i = 0; i < inputDim; i++) {
            anisotropicSpread[i] = Math.pow(0.5, i);  // strong in lower dims
        }        
        List<VectorN> anisotropicDataset = factory.generateAnisotropicClusteredData(
                numClusters, pointsPerCluster, anisotropicSpread);  
        printTotalTime(startTime);

        System.out.println("Fitting HoroPCA model and transforming data...");
        startTime = System.nanoTime();
        HoroPCA horo2 = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
            HoroPCA.InitializationStrategy.KMEANS_PLUS_PLUS, 
                HyperbolicUtils.hyperbolicSquaredDist, seed);
        horo2.verbose = true;
        List<VectorN> transformed2  = horo2.fitTransform(anisotropicDataset);
        printTotalTime(startTime);
        
        System.out.println("Reconstructing transformed data using inverseTransform()...");
        startTime = System.nanoTime();
        List<VectorN> reconstructed2 = transformed2.stream()
            .peek(t -> { System.out.println("Reconstructing:" + t); })
            .map(v -> { return horo2.inverseTransform(v, maxIterations, learningRate, true); })
            .toList();
        double reconError2 = HoroPCAValidator.reconstructionError(anisotropicDataset, reconstructed2);
        printTotalTime(startTime);
        System.out.printf("Mean Reconstruction Error: %.6f\n", reconError2);        
    }
}