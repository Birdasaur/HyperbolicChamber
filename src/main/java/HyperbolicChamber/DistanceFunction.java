package HyperbolicChamber;

public interface DistanceFunction {
    /**
     * Compute distance between two points.
     * @param a first point as double array
     * @param b second point as double array
     * @return distance (non-negative)
     */
    double distance(double[] a, double[] b);
}