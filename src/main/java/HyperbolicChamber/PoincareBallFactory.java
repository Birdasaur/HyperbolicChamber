package HyperbolicChamber;

/**
 *
 * @author Sean Phillips
 */
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class PoincareBallFactory {

    private final int dim;
    private final Random random;

    public PoincareBallFactory(int dim, long seed) {
        this.dim = dim;
        this.random = new Random(seed);
    }
    //Generates clustered data with per-dimension spread
    public List<VectorN> generateAnisotropicClusteredData(int numClusters, int pointsPerCluster, double[] clusterSpreadPerDim) {
        List<VectorN> data = new ArrayList<>();
        for (int c = 0; c < numClusters; c++) {
            VectorN center = randomUnitVector().scale(0.3); // small center radius
            for (int i = 0; i < pointsPerCluster; i++) {
                VectorN noise = randomAnisotropicVector(clusterSpreadPerDim);
                VectorN point = center.add(noise);

                // Avoid projection: ensure the point stays within ball
                double norm = Math.sqrt(point.normSq());
                if (norm >= 1.0) {
                    point = point.scale(0.99 / norm); // manual safety
                }
                data.add(point);
            }
        }
        return data;
    }

    public List<VectorN> generateIsotropicClusteredData(int numClusters, int pointsPerCluster, double clusterSpread) {
        List<VectorN> data = new ArrayList<>();
        for (int c = 0; c < numClusters; c++) {
            VectorN center = randomVector(0.3); // Cluster centers
            for (int i = 0; i < pointsPerCluster; i++) {
                VectorN offset = randomVector(clusterSpread);
                VectorN point = center.add(offset);
                data.add(projectToPoincareBall(point));
            }
        }
        return data;
    }

    public List<VectorN> generateUniformSpherePoints(int count) {
        List<VectorN> points = new ArrayList<>();
        for (int i = 0; i < count; i++) {
            double[] vec = new double[dim];
            double norm = 0;
            for (int d = 0; d < dim; d++) {
                double val = random.nextGaussian();
                vec[d] = val;
                norm += val * val;
            }
            norm = Math.sqrt(norm);
            for (int d = 0; d < dim; d++) {
                vec[d] /= norm; // Normalize to unit sphere
            }
            points.add(new VectorN(vec));
        }
        return points;
    }

    public VectorN randomUnitVector() {
        double[] vec = new double[dim];
        double norm = 0;
        for (int i = 0; i < dim; i++) {
            double v = random.nextGaussian();
            vec[i] = v;
            norm += v * v;
        }
        norm = Math.sqrt(norm);
        for (int i = 0; i < dim; i++) {
            vec[i] /= norm;
        }
        return new VectorN(vec);
    }

    public VectorN randomVector(double scale) {
        double[] values = new double[dim];
        for (int i = 0; i < dim; i++) {
            values[i] = random.nextGaussian() * scale;
        }
        return new VectorN(values);
    }

    public VectorN randomAnisotropicVector(double[] scales) {
        double[] values = new double[dim];
        for (int i = 0; i < dim; i++) {
            values[i] = random.nextGaussian() * scales[i];
        }
        return new VectorN(values);
    }

    public VectorN projectToPoincareBall(VectorN v) {
        double norm = Math.sqrt(v.normSq());
        double scale = norm >= 1.0 ? 0.99 / norm : 1.0;
        return v.scale(scale);
    }
}