package HyperbolicChamber;

/**
 *
 * @author Sean Phillips
 */
import java.util.Arrays;

public class VectorN {
    private final double[] values;

    public VectorN(int dimension) {
        this.values = new double[dimension];
    }

    public VectorN(double[] values) {
        this.values = Arrays.copyOf(values, values.length);
    }

    public int dimension() {
        return values.length;
    }

    public double get(int i) {
        return values[i];
    }

    public void set(int i, double value) {
        values[i] = value;
    }

    public double[] toArray() {
        return Arrays.copyOf(values, values.length);
    }

    public VectorN add(VectorN other) {
        checkSameDimension(other);
        double[] result = new double[dimension()];
        for (int i = 0; i < dimension(); i++) {
            result[i] = this.values[i] + other.values[i];
        }
        return new VectorN(result);
    }

    public VectorN subtract(VectorN other) {
        checkSameDimension(other);
        double[] result = new double[dimension()];
        for (int i = 0; i < dimension(); i++) {
            result[i] = this.values[i] - other.values[i];
        }
        return new VectorN(result);
    }

    public VectorN scale(double scalar) {
        double[] result = new double[dimension()];
        for (int i = 0; i < dimension(); i++) {
            result[i] = this.values[i] * scalar;
        }
        return new VectorN(result);
    }

    public double dot(VectorN other) {
        checkSameDimension(other);
        double sum = 0.0;
        for (int i = 0; i < dimension(); i++) {
            sum += this.values[i] * other.values[i];
        }
        return sum;
    }

    public double norm() {
        return Math.sqrt(dot(this));
    }
    /**
     * Returns the squared Euclidean norm of this vector.
     * Equivalent to dot(this, this).
     */
    public double normSq() {
        double sum = 0.0;
        for (double v : values) {
            sum += v * v;
        }
        return sum;
    }

    public VectorN normalize() {
        double norm = norm();
        if (norm == 0) throw new ArithmeticException("Cannot normalize a zero vector.");
        return scale(1.0 / norm);
    }

    public VectorN copy() {
        return new VectorN(this.values);
    }

    private void checkSameDimension(VectorN other) {
        if (this.dimension() != other.dimension()) {
            throw new IllegalArgumentException("Vector dimensions do not match.");
        }
    }

    @Override
    public String toString() {
        return Arrays.toString(values);
    }
}
