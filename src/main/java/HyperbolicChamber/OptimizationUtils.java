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
                                            double learningRate,
                                            int maxIterations,
                                            double tolerance) {
        VectorN direction = initialDirection;
        double prevVariance = computeBusemannVariance(direction, data);

        for (int iter = 0; iter < maxIterations; iter++) {
            VectorN gradient = computeBusemannGradient(direction, data);
            direction = direction.add(gradient.scale(learningRate)).normalize();

            double variance = computeBusemannVariance(direction, data);
            if (Math.abs(variance - prevVariance) < tolerance) break;

            prevVariance = variance;
        }

        return direction;
    }
    public static VectorN randomPointInPoincareBall(int dim) {
        VectorN v;
        do {
            double[] coords = new double[dim];
            for (int i = 0; i < dim; i++) {
                coords[i] = ThreadLocalRandom.current().nextDouble(-1.0, 1.0);
            }
            v = new VectorN(coords);
        } while (v.norm() >= 1.0);
        return v;
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
            VectorN grad_i = HyperbolicUtils.busemannProjectionGradient(point, direction);
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