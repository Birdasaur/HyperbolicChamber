package HyperbolicChamber;

/**
 *
 * @author sean phillips
 */
import java.util.List;

public class HoroPCAValidator {

public static void printExplainedVariance(HoroPCA model, List<VectorN> transformedData) {
    double[] ev = model.explainedVariance(transformedData);
    double[] ratio = model.explainedVarianceRatio(transformedData);

    System.out.println("Explained Variance per Component:");
    for (int i = 0; i < ev.length; i++) {
        System.out.printf("  [%d]: %.9f (%.2f%%)\n", i, ev[i], ratio[i] * 100.0);
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
}

