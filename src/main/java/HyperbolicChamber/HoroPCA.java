package HyperbolicChamber;

/**
 *
 * @author Sean Phillips
 */
import static HyperbolicChamber.HoroPCAValidator.printExplainedVariance;
import static HyperbolicChamber.HoroPCAValidator.printVariancePerDimension;
import static HyperbolicChamber.HyperbolicUtils.hyperbolicSquaredDist;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class HoroPCA {

    public enum InitializationStrategy {
        RANDOM_UNIT_SPHERE,
        KMEANS_PLUS_PLUS
    }

    private final int numComponents;
    private final int maxIterations;
    private final double learningRate;
    private final double tolerance;
    private final InitializationStrategy initStrategy;
    private final DistanceFunction distanceFunction; // used for KMeans++ if needed
    private final int kmeansSeed;
    private List<VectorN> projectionDirections; // learned ξ vectors
    public boolean verbose = false;
    public boolean verboseOptimization = false;
    
    public HoroPCA() {
        this(3, 1000, 0.05, 1e-8, InitializationStrategy.KMEANS_PLUS_PLUS, hyperbolicSquaredDist, 42);
    }

    public HoroPCA(int numComponents,
            int maxIterations,
            double learningRate,
            double tolerance,
            InitializationStrategy initStrategy,
            DistanceFunction distanceFunction,
            int kmeansSeed) {
        this.numComponents = numComponents;
        this.maxIterations = maxIterations;
        this.learningRate = learningRate;
        this.tolerance = tolerance;
        this.initStrategy = initStrategy;
        this.distanceFunction = distanceFunction;
        this.kmeansSeed = kmeansSeed;
        this.projectionDirections = new ArrayList<>();
    }

    /**
     * Fit the HoroPCA model to the hyperbolic data by finding k principal
     * horospherical directions.
     *
     * @param data List of points in hyperbolic space (assumed to be in Poincaré
     * ball model)
     * @param k Number of principal directions to extract
     */
    public void fit(List<VectorN> data, int k) {
        projectionDirections.clear();
        if (verbose) {
            System.out.println("[DEBUG] Starting HoroPCA fit...");
        }

        // STEP 1: Variance of raw data (for baseline)
        if (verbose) {
            System.out.println("[DEBUG] Raw data variance:");
            printVariancePerDimension(data);
        }

        // STEP 2: Initialize directions
        List<VectorN> initialDirections = switch (initStrategy) {
            case RANDOM_UNIT_SPHERE ->
                OptimizationUtils.generateRandomUnitVectors(k, data.get(0).dimension());
            case KMEANS_PLUS_PLUS ->
                OptimizationUtils.selectInitialDirectionsWithKMeans(
                data, kmeansSeed, hyperbolicSquaredDist, k);
        };

        // Optional: Check variance in initial projections
        if (verbose) {
            System.out.println("[DEBUG] Initial horospherical projections:");
            printExplainedVariance(data, initialDirections);
        }
        // STEP 3: Optimize each direction
        for (int i = 0; i < k; i++) {
            VectorN direction = OptimizationUtils.optimizeDirection(
                    initialDirections.get(i), data, learningRate,
                    maxIterations, tolerance, verboseOptimization);
            projectionDirections.add(direction);
        }

        // STEP 4: Final projections
        if (verbose) {
            System.out.println("[DEBUG] Final horospherical projections:");
            printExplainedVariance(data, projectionDirections);
        }
    }

    /**
     * Transform the data into k-dimensional Euclidean coordinates using
     * Busemann projections.
     *
     * @param data Hyperbolic data points
     * @return List of projected vectors in R^k
     */
    public List<VectorN> transform(List<VectorN> data) {
        List<VectorN> result = new ArrayList<>();
        for (VectorN point : data) {
            double[] projection = new double[projectionDirections.size()];
            for (int i = 0; i < projectionDirections.size(); i++) {
                projection[i] = HyperbolicUtils.busemannProjection(point, projectionDirections.get(i));
            }
            result.add(new VectorN(projection));
        }
        return result;
    }

    public List<VectorN> fitTransform(List<VectorN> data) {
        return fitTransform(data, numComponents);
    }

    /**
     * Convenience method to fit and transform the data.
     */
    public List<VectorN> fitTransform(List<VectorN> data, int k) {
        fit(data, k);
        return transform(data);
    }

    /**
     * Attempt to reconstruct a point in hyperbolic space (Poincaré ball) from
     * its projected coordinates.
     *
     * @param projection A k-dimensional vector (output of transform)
     * @param maxIters Maximum number of optimization iterations
     * @param stepSize Gradient descent learning rate
     * @return Approximate point in original space (Poincaré ball)
     */
    public VectorN inverseTransform(VectorN projection, int maxIters, double stepSize) {
        int dim = projectionDirections.get(0).dimension();
        VectorN x = OptimizationUtils.randomPointInPoincareBall(dim);

        for (int iter = 0; iter < maxIters; iter++) {
            double[] grad = new double[dim];

            for (int i = 0; i < projectionDirections.size(); i++) {
                VectorN e = projectionDirections.get(i); //e == ξ
                double z_i = projection.get(i);
                double β_i = HyperbolicUtils.busemannProjection(x, e);
                double error = β_i - z_i;

                VectorN grad_i = HyperbolicUtils.busemannGradient(x, e).scale(2 * error);
                for (int j = 0; j < dim; j++) {
                    grad[j] += grad_i.get(j);
                }
            }

            VectorN gradient = new VectorN(grad).normalize().scale(stepSize);
            VectorN updated = x.subtract(gradient);

            // Project back into ball if needed
            if (updated.norm() >= 1.0) {
                updated = updated.normalize().scale(0.999);
            }

            x = updated;
        }

        return x;
    }

    /**
     * Computes the explained variance along each horospherical direction ξ.
     *
     * @param projectedData The output of the transform() method (n × k).
     * @return An array of variances, one per component.
     */
    public double[] explainedVariance(List<VectorN> transformedData) {
        if (transformedData.isEmpty()) {
            return new double[0];
        }

        int k = transformedData.get(0).dimension();  // number of components
        int n = transformedData.size();              // number of data points
        double[] variances = new double[k];

        for (int j = 0; j < k; j++) {
            double sum = 0.0;
            double sumSq = 0.0;

            for (VectorN v : transformedData) {
                double x = v.get(j);
                sum += x;
                sumSq += x * x;
            }

            double mean = sum / n;
            variances[j] = (sumSq / n) - (mean * mean);  // population variance
        }

        return variances;
    }

    /**
     * Computes the percentage of total variance explained by each horospherical
     * component.
     *
     * @param projectedData The output of the transform() method (n × k).
     * @return An array of percentages summing to 1.0 (or 100.0 if scaled).
     */

    public double[] explainedVarianceRatio(List<VectorN> transformedData) {
        double[] variances = explainedVariance(transformedData);
        double total = Arrays.stream(variances).sum();
        if (total == 0.0) {
            return new double[variances.length];
        }

        double[] ratio = new double[variances.length];
        for (int i = 0; i < variances.length; i++) {
            ratio[i] = variances[i] / total;
        }
        return ratio;
    }

    /**
     * @return The list of learned principal horospherical directions.
     */
    public List<VectorN> getPrincipalDirections() {
        return new ArrayList<>(projectionDirections);
    }


}
