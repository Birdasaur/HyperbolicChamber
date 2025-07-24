package HyperbolicChamber;

import java.util.List;

/**
 *
 * @author sean phillips
 */
public class HoroPCATest {
    public static void main(String[] args) {
        int inputDim = 10;
        int targetDim = 3;
        int maxIterations = 500;
        int numClusters = 4;
        int pointsPerCluster = 50;
        int seed = 42;
        double learningRate = 0.005;
        double tolerance = 1e-6;

        // Step 1: Generate synthetic clustered data on the Poincaré ball
        PoincareBallFactory factory = new PoincareBallFactory(inputDim, seed);
        List<VectorN> dataset = factory.generateClusteredData(numClusters, pointsPerCluster, 0.05);

        // Step 2: Train HoroPCA model
        HoroPCA horo = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
            HoroPCA.InitializationStrategy.RANDOM_UNIT_SPHERE, 
                HyperbolicUtils.hyperbolicSquaredDist, seed);
        horo.verbose = true;
        
        // Step 3: Transform and inverse transform for reconstruction
        List<VectorN> transformed = horo.fitTransform(dataset);
        List<VectorN> reconstructed = transformed.stream().map(v -> {
            VectorN r = horo.inverseTransform(v, maxIterations, learningRate);
            return r;
        }).toList();

        // Step 4: Print explained variance and reconstruction error
        HoroPCAValidator.printExplainedVariance(horo, transformed);
        double reconError = HoroPCAValidator.reconstructionError(dataset, reconstructed);
        System.out.printf("Mean Reconstruction Error: %.6f\n", reconError);
    }
    
}
