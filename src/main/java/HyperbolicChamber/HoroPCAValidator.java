package HyperbolicChamber;

import static HyperbolicChamber.HyperbolicUtils.horosphericalProjection;
import java.util.List;

/**
 *
 * @author sean phillips
 */
public class HoroPCAValidator {

    public static void printExplainedVariance(List<VectorN> data, List<VectorN> directions) {
        int k = directions.size();
        double[] variances = new double[k];
        double total = 0;

        for (VectorN x : data) {
            for (int i = 0; i < k; i++) {
                double p = horosphericalProjection(x, directions.get(i));
                variances[i] += p * p;
            }
        }

        for (int i = 0; i < k; i++) {
            variances[i] /= data.size();
            total += variances[i];
        }

        System.out.println("Explained Variance per Direction:");
        for (int i = 0; i < k; i++) {
            double percent = 100.0 * variances[i] / total;
            System.out.printf("  Dir %d: %.6f (%.2f%%)%n", i, variances[i], percent);
        }
    }

    public static void printVariancePerDimension(List<VectorN> data) {
        if (data.isEmpty()) {
            return;
        }

        int dim = data.get(0).dimension();
        double[] means = new double[dim];
        double[] variances = new double[dim];
        int n = data.size();

        // Compute means
        for (VectorN vec : data) {
            for (int i = 0; i < dim; i++) {
                means[i] += vec.get(i);
            }
        }
        for (int i = 0; i < dim; i++) {
            means[i] /= n;
        }

        // Compute variances
        for (VectorN vec : data) {
            for (int i = 0; i < dim; i++) {
                double diff = vec.get(i) - means[i];
                variances[i] += diff * diff;
            }
        }
        for (int i = 0; i < dim; i++) {
            variances[i] /= n;
        }

        // Print result
        System.out.println("Raw Variance per Dimension:");
        for (int i = 0; i < dim; i++) {
            System.out.printf("  Dim %d: %.6f%n", i, variances[i]);
        }
    }

    public static double reconstructionError(List<VectorN> original, List<VectorN> reconstructed) {
        double totalError = 0.0;
        int count = Math.min(original.size(), reconstructed.size());

        for (int i = 0; i < count; i++) {
            VectorN diff = original.get(i).subtract(reconstructed.get(i));
            totalError += diff.normSq();
        }

        return totalError / count;
    }
    public static void printTotalTime(long startTime) {
        System.out.println(totalTimeString(startTime));
    }

    public static String totalTimeString(long startTime) {
        long estimatedTime = System.nanoTime() - startTime;
        long totalNanos = estimatedTime;
        long s = totalNanos / 1000000000;
        totalNanos -= s * 1000000000;
        long ms = totalNanos / 1000000;
        totalNanos -= ms * 1000000;

        long us = totalNanos / 1000;
        totalNanos -= us * 1000;
        return "Total elapsed time: " + s + ":s:" + ms + ":ms:" + us + ":us:" + totalNanos + ":ns";
    }    
}
