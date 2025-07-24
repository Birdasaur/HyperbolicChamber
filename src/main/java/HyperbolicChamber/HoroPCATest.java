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
        double learningRate = 0.001;
        double tolerance = 1e-8;

        // Step 1: Generate synthetic clustered data on the Poincaré ball
        PoincareBallFactory factory = new PoincareBallFactory(inputDim, seed);
        List<VectorN> isotropicDataset = factory.generateIsotropicClusteredData(
                numClusters, pointsPerCluster, 0.05);
        // Step 2: Train HoroPCA model
        HoroPCA horo = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
            HoroPCA.InitializationStrategy.RANDOM_UNIT_SPHERE, 
                HyperbolicUtils.hyperbolicSquaredDist, seed);
        horo.verbose = true;
        // Step 3: Transform and inverse transform for reconstruction
        List<VectorN> transformed = horo.fitTransform(isotropicDataset);
        List<VectorN> reconstructed = transformed.stream().map(v -> {
            VectorN r = horo.inverseTransform(v, maxIterations, learningRate);
            return r;
        }).toList();
        // Step 4: Print explained variance and reconstruction error
        double reconError = HoroPCAValidator.reconstructionError(isotropicDataset, reconstructed);
        System.out.printf("Mean Reconstruction Error: %.6f\n", reconError);

        // Repeat for anisotropic set
//        double[] scales = new double[] { 1.0, 0.6, 0.2, 
//            0.08, 0.06, 0.04, 0.02, 0.006, 0.004, 0.002 };   
        double[] anisotropicSpread = new double[inputDim];
        for (int i = 0; i < inputDim; i++) {
            anisotropicSpread[i] = Math.pow(0.5, i);  // strong in lower dims
        }        
        List<VectorN> anisotropicDataset = factory.generateAnisotropicClusteredData(
                numClusters, pointsPerCluster, anisotropicSpread);  
        HoroPCA horo2 = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
            HoroPCA.InitializationStrategy.RANDOM_UNIT_SPHERE, 
                HyperbolicUtils.hyperbolicSquaredDist, seed);
        horo2.verbose = true;
        List<VectorN> transformed2  = horo2.fitTransform(anisotropicDataset);
        List<VectorN> reconstructed2 = transformed2.stream().map(v -> {
            VectorN r = horo2.inverseTransform(v, maxIterations, learningRate);
            return r;
        }).toList();
        double reconError2 = HoroPCAValidator.reconstructionError(anisotropicDataset, reconstructed2);
        System.out.printf("Mean Reconstruction Error: %.6f\n", reconError2);        
        
    }
    
}
