package HyperbolicChamber;

/**
 *
 * @author Sean Phillips
 */
public class HyperbolicUtils {

    // Define hyperbolic squared distance
    public static DistanceFunction hyperbolicSquaredDist = (a, b) -> {
        VectorN va  = new VectorN(a);
        VectorN vb = new VectorN(b);
        double d = poincareDistance(va, vb);
        return d * d;
    };

    /**
     * Compute the Busemann function B_xi(x) in the Poincaré ball model. xi must
     * be on the boundary (||xi|| = 1).
     */
    public static double busemannProjection(VectorN x, VectorN xi) {
        double normX = x.norm();
        double dot = x.dot(xi);

        double numerator = 1 + normX * normX - 2 * dot;
        double denominator = 1 - normX * normX;

        if (denominator <= 0) {
            throw new ArithmeticException("Point is at or beyond the Poincaré boundary.");
        }

        return Math.log(numerator / denominator);
    }

    public static VectorN busemannProjectionGradientDirection(VectorN x, VectorN xi) {
        double numerator = 1 + x.normSq() - 2 * x.dot(xi);
        if (numerator == 0) {
            numerator = 1e-10; // avoid div by zero
        }
        VectorN grad = x.scale(-2.0 / numerator);
        return grad;
    }

    /**
     * Gradient of the Busemann function B_xi(x) with respect to x in the
     * Poincaré ball model. xi must be on the boundary (||xi|| = 1).
     */
    public static VectorN busemannProjectionGradient(VectorN x, VectorN xi) {
        double xNormSquared = x.dot(x);
        double dotXxi = x.dot(xi);

        double denom = 1 - xNormSquared;
        if (denom <= 0.0) {
            throw new ArithmeticException("Point is at or beyond the Poincaré boundary.");
        }

        // Compute gradient of the log(numerator / denominator)
        VectorN gradNumerator = x.scale(2.0).subtract(xi.scale(2.0));
        VectorN gradDenominator = x.scale(-2.0);

        double numerator = 1 + xNormSquared - 2 * dotXxi;
        double denominator = denom;

        VectorN grad = gradNumerator.scale(1.0 / numerator)
                .subtract(gradDenominator.scale(1.0 / denominator));

        return grad;
    }

    public static VectorN busemannGradient(VectorN x, VectorN direction) {
        VectorN u = direction.normalize();
        double normSq = x.normSq();
        double denom = 1.0 - normSq;
        if (denom <= 0) {
            denom = 1e-6;
        }

        // ∇_x β(x, ξ) ≈ (2 / (1 - ||x||²)) * (x - ⟨x, u⟩ * u)
        double dot = x.dot(u);
        VectorN proj = u.scale(dot);
        VectorN diff = x.subtract(proj);

        return diff.scale(2.0 / denom);
    }

    public static double horosphericalProjection(VectorN x, VectorN v) {
        double dot = x.dot(v);  // dot product ⟨x, v⟩
        double value = 1.0 - dot;

        // Guard against numerical instability
        if (value <= 1e-12) {
            value = 1e-12;
        }

        return Math.log(2.0) - Math.log(value);
    }

    //Inverse versions of some math functions
    public static double acosh(double x) {
        if (x < 1.0) {
            throw new IllegalArgumentException("acosh is undefined for x < 1");
        }
        return Math.log(x + Math.sqrt(x * x - 1));
    }

    public static double atanh(double x) {
        if (x <= -1.0 || x >= 1.0) {
            throw new IllegalArgumentException("atanh is only defined for -1 < x < 1");
        }
        return 0.5 * Math.log((1 + x) / (1 - x));
    }

    //vector operations delegate to VectorN
    public static double norm(VectorN v) {
        return v.norm();
    }

    public static double dot(VectorN a, VectorN b) {
        return a.dot(b);
    }

    public static boolean isInsidePoincareBall(VectorN v) {
        return v.norm() < 1.0;
    }

    public static VectorN normalize(VectorN v) {
        return v.normalize();
    }

    /**
     * Hyperbolic distance in the Poincaré ball model.
     */
    public static double poincareDistance(VectorN x, VectorN y) {
        VectorN diff = x.subtract(y);
        double normDiffSquared = diff.dot(diff);
        double normX = x.norm();
        double normY = y.norm();

        double denom = (1 - normX * normX) * (1 - normY * normY);
        double argument = 1 + 2 * normDiffSquared / denom;

        if (argument < 1.0) {
            argument = 1.0;
        }

        return acosh(argument);
    }

    // Möbius operations in the Poincaré ball model 
    /**
     * Möbius addition: x ⊕ y
     */
    public static VectorN mobiusAdd(VectorN x, VectorN y) {
        double xNormSq = x.dot(x);
        double yNormSq = y.dot(y);
        double xyDot = x.dot(y);

        double denominator = 1 + 2 * xyDot + xNormSq * yNormSq;
        if (denominator == 0.0) {
            throw new ArithmeticException("Möbius addition division by zero");
        }

        VectorN term1 = x.scale(1 + 2 * xyDot + yNormSq);
        VectorN term2 = y.scale(1 - xNormSq);

        return term1.add(term2).scale(1.0 / denominator);
    }

    /**
     * Möbius scalar multiplication: r ⊗ x
     */
    public static VectorN mobiusScalarMultiply(double r, VectorN x) {
        double normX = x.norm();
        if (normX == 0.0) {
            return x.copy();
        }

        double scale = Math.tanh(r * atanh(normX)) / normX;
        return x.scale(scale);
    }

    /**
     * Möbius exponential map: exp_x(v)
     */
    public static VectorN mobiusExpMap(VectorN x, VectorN v) {
        double normV = v.norm();
        if (normV == 0.0) {
            return x.copy();
        }

        double lambda = 2.0 / (1.0 - x.dot(x));
        double factor = Math.tanh(lambda * normV / 2.0) / normV;
        VectorN mapped = v.scale(factor);
        return mobiusAdd(x, mapped);
    }

    /**
     * Möbius logarithmic map: log_x(y)
     */
    public static VectorN mobiusLogMap(VectorN x, VectorN y) {
        VectorN diff = mobiusAdd(y, x.scale(-1));
        double normDiff = diff.norm();
        if (normDiff == 0.0) {
            return new VectorN(x.dimension());
        }

        double lambda = 2.0 / (1.0 - x.dot(x));
        double factor = atanh(normDiff) * 2.0 / (lambda * normDiff);
        return diff.scale(factor);
    }
}
