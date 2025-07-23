package HyperbolicChamber;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class KMeansPlusPlus {

    private final double[][] centroids;          // Selected centroids
    private final List<Integer> centroidIndices; // Indices of selected centroids in original dataset

    /**
     * Constructs the KMeans++ initializer with custom distance.
     *
     * @param data Input dataset (n x d)
     * @param k Number of centroids to select
     * @param seed Random seed for reproducibility
     * @param distFunc Distance function to use for computing distances
     */
    public KMeansPlusPlus(double[][] data, int k, int seed, DistanceFunction distFunc) {
        Random rand = new Random(seed);
        int n = data.length;
        int d = data[0].length;

        centroids = new double[k][d];
        centroidIndices = new ArrayList<>();

        // Step 1: Choose first centroid randomly
        int firstIndex = rand.nextInt(n);
        centroids[0] = Arrays.copyOf(data[firstIndex], d);
        centroidIndices.add(firstIndex);

        // Step 2: Choose remaining centroids with weighted probability
        for (int i = 1; i < k; i++) {
            double[] distances = new double[n];

            // Compute minimum squared distance to any existing centroid
            for (int j = 0; j < n; j++) {
                double minDist = Double.MAX_VALUE;
                for (int c = 0; c < i; c++) {
                    double dist = distFunc.distance(data[j], centroids[c]);
                    if (dist < minDist) {
                        minDist = dist;
                    }
                }
                distances[j] = minDist;
            }

            // Weighted random selection
            double sum = Arrays.stream(distances).sum();
            double r = rand.nextDouble() * sum;
            double cumulative = 0;

            for (int j = 0; j < n; j++) {
                cumulative += distances[j];
                if (cumulative >= r) {
                    centroids[i] = Arrays.copyOf(data[j], d);
                    centroidIndices.add(j);
                    break;
                }
            }
        }
    }

    /**
     * Returns the selected centroids.
     *
     * @return Array of k centroids (k x d)
     */
    public double[][] getCentroids() {
        return centroids;
    }

    /**
     * Returns the indices of the selected centroids in the original dataset.
     *
     * @return List of selected indices (size k)
     */
    public List<Integer> getCentroidIndices() {
        return centroidIndices;
    }
}

