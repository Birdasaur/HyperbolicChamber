package HyperbolicChamber;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.Collectors;

public class OptimizationUtils {

    /**
     * Perform gradient ascent to find the direction ξ that maximizes Busemann variance.
     */
public static VectorN optimizeDirection(VectorN initialDirection,
                                        List<VectorN> data,
                                        double initialLearningRate,
                                        int maxIterations,
                                        double tolerance,
                                        boolean verbose) {

    VectorN direction = initialDirection.normalize();
    double prevVariance = computeBusemannVariance(direction, data);

    double learningRate = initialLearningRate;
    double decayRate = 0.98;           // Multiplied each 50 iterations
    int decayInterval = 50;
    double maxGradientNorm = 10.0;

    if (verbose) {
        System.out.printf("Initial Busemann Variance: %.6f\n", prevVariance);
    }

    for (int iter = 1; iter <= maxIterations; iter++) {
        VectorN gradient = computeBusemannGradient(direction, data);

        // Gradient clipping (safe step control)
        double gradNorm = gradient.norm();
        if (gradNorm > maxGradientNorm) {
            gradient = gradient.scale(maxGradientNorm / gradNorm);
        }

        // Learning rate decay
        if (iter % decayInterval == 0) {
            learningRate *= decayRate;
        }

        // Gradient ascent update with projection back to unit sphere
        direction = direction.add(gradient.scale(learningRate)).normalize();

        double variance = computeBusemannVariance(direction, data);
        double delta = Math.abs(variance - prevVariance);

        if (verbose) {
            System.out.printf("Iter %3d | Variance: %.6f | Δ: %.8f | GradNorm: %.4f | LR: %.6f\n",
                    iter, variance, delta, gradNorm, learningRate);
        }

        if (delta < tolerance) {
            if (verbose) {
                System.out.println("Converged: Δ variance below tolerance.");
            }
            break;
        }

        prevVariance = variance;
    }

    return direction;
}

public static VectorN randomPointInPoincareBall(int dim) {
    double[] coords = new double[dim];

    // Step 1: Sample direction from standard normal distribution
    double normSq = 0.0;
    for (int i = 0; i < dim; i++) {
        double val = ThreadLocalRandom.current().nextGaussian();
        coords[i] = val;
        normSq += val * val;
    }

    // Step 2: Normalize to unit vector
    double norm = Math.sqrt(normSq);
    for (int i = 0; i < dim; i++) {
        coords[i] /= norm;
    }

    // Step 3: Sample radius using inverse CDF for uniformity inside hypersphere
    double u = ThreadLocalRandom.current().nextDouble();  // in [0,1)
    double radius = Math.pow(u, 1.0 / dim) * 0.999;  // safely inside unit ball

    for (int i = 0; i < dim; i++) {
        coords[i] *= radius;
    }

    return new VectorN(coords);
}

    /**
     * Compute the variance of Busemann projections along a direction.
     */
    private static double computeBusemannVariance(VectorN direction, List<VectorN> data) {
        List<Double> projections = data.stream()
                .map(p -> HyperbolicUtils.busemannProjection(p, direction))
                .collect(Collectors.toList());

        double mean = projections.stream().mapToDouble(Double::doubleValue).average().orElse(0.0);
        return projections.stream()
                .mapToDouble(v -> (v - mean) * (v - mean))
                .average()
                .orElse(0.0);
    }

    /**
     * Compute gradient of Busemann variance w.r.t. direction.
     */
    private static VectorN computeBusemannGradient(VectorN direction, List<VectorN> data) {
        int dim = direction.dimension();
        double[] grad = new double[dim];

        // Mean of projections
        double[] projections = new double[data.size()];
        for (int i = 0; i < data.size(); i++) {
            projections[i] = HyperbolicUtils.busemannProjection(data.get(i), direction);
        }
        double mean = Arrays.stream(projections).average().orElse(0.0);

        // Gradient: 2 * sum_i ((proj_i - mean) * ∇proj_i)
        for (int i = 0; i < data.size(); i++) {
            VectorN point = data.get(i);
            double diff = projections[i] - mean;
            VectorN grad_i = HyperbolicUtils.busemannProjectionGradientDirection(point, direction);
            for (int j = 0; j < dim; j++) {
                grad[j] += 2 * diff * grad_i.get(j);
            }
        }

        return new VectorN(grad);
    }

    /**
     * Generate a list of random unit vectors on the sphere.
     */
    public static List<VectorN> generateRandomUnitVectors(int count, int dimension) {
        Random random = new Random();
        List<VectorN> vectors = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            double[] components = new double[dimension];
            double norm = 0;
            for (int j = 0; j < dimension; j++) {
                double val = random.nextGaussian();
                components[j] = val;
                norm += val * val;
            }
            norm = Math.sqrt(norm);
            for (int j = 0; j < dimension; j++) {
                components[j] /= norm;
            }
            vectors.add(new VectorN(components));
        }
        return vectors;
    }

    /**
     * Select k smart directions using KMeans++ centers normalized to the unit sphere.
     */
    public static List<VectorN> selectInitialDirectionsWithKMeans(List<VectorN> data,
                                                                   int k,
                                                                   DistanceFunction distanceFunction,
                                                                   int seed) {
        // Convert List<VectorN> to double[][] for KMeans++
        int n = data.size();
        int dim = data.get(0).dimension();
        double[][] dataArray = new double[n][dim];
        for (int i = 0; i < n; i++) {
            dataArray[i] = data.get(i).toArray();
        }

        // Use your custom KMeans++ initializer
        KMeansPlusPlus kmeans = new KMeansPlusPlus(dataArray, k, seed, distanceFunction);
        double[][] centroidArray = kmeans.getCentroids();

        // Convert centroids to normalized VectorN
        List<VectorN> centroids = new ArrayList<>();
        for (double[] c : centroidArray) {
            VectorN vec = new VectorN(c).normalize();
            centroids.add(vec);
        }

        return centroids;
    }

    /**
     * Basic Euclidean distance squared fallback.
     */
    public static double euclideanDistanceSquared(VectorN a, VectorN b) {
        double sum = 0;
        for (int i = 0; i < a.dimension(); i++) {
            double diff = a.get(i) - b.get(i);
            sum += diff * diff;
        }
        return sum;
    }
}