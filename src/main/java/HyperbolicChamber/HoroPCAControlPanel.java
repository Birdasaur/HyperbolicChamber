package HyperbolicChamber;

import HyperbolicChamber.HoroPCA.InitializationStrategy;
import static HyperbolicChamber.HoroPCAValidator.printTotalTime;
import java.util.List;
import javafx.scene.layout.VBox;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Group;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ComboBox;
import javafx.scene.control.Control;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Sphere;
import javafx.util.StringConverter;

/**
 *
 * @author Sean Phillips
 */
public class HoroPCAControlPanel extends VBox {

    // Controls
    private final Spinner<Integer> inputDimSpinner = createIntegerSpinner(1, 10000, 100);
    private final Spinner<Integer> targetDimSpinner = createIntegerSpinner(1, 100, 3);
    private final Spinner<Integer> maxIterationsSpinner = createIntegerSpinner(1, 10000, 300);
    private final Spinner<Integer> numClustersSpinner = createIntegerSpinner(1, 1000, 4);
    private final Spinner<Integer> pointsPerClusterSpinner = createIntegerSpinner(1, 10000, 50);
    private final Spinner<Double> scaleSpinner = createDoubleSpinner(1, 1000.0, 50.0, 1);
    private final Spinner<Double> learningRateSpinner = createDoubleSpinner(0.0001, 1.0, 0.01, 0.0001);
    private final Spinner<Double> toleranceSpinner = createDoubleSpinner(1e-10, 1e-2, 1e-8, 1e-10);

    private final ComboBox<InitializationStrategy> initStrategyCombo = new ComboBox<>();
    private final CheckBox verboseCheckBox = new CheckBox("Verbose Output");
    private final CheckBox verboseOptCheckBox = new CheckBox("Verbose Optimization");

    Group dataGroup;
    
    public HoroPCAControlPanel(Group dataGroup) {
        this.dataGroup = dataGroup;
        setSpacing(10);
        setPadding(new Insets(10));
        setMinWidth(200);

        getChildren().addAll(
            labeled("Input Dimension", inputDimSpinner),
            labeled("Target Dimension", targetDimSpinner),
            labeled("Max Iterations", maxIterationsSpinner),
            labeled("Clusters", numClustersSpinner),
            labeled("Points per Cluster", pointsPerClusterSpinner),
            labeled("Learning Rate", learningRateSpinner),
            labeled("Tolerance", toleranceSpinner),
            labeled("Output Scale", scaleSpinner),            
            labeled("Initialization Strategy", initStrategyCombo),
            verboseCheckBox,
            verboseOptCheckBox
        );

        // Set up combo box with enum values
        initStrategyCombo.getItems().setAll(InitializationStrategy.values());
        initStrategyCombo.getSelectionModel().selectFirst();

        // Run button
        Button runButton = new Button("Run HoroPCA");
        runButton.setOnAction(e -> runHoroPCA());
        getChildren().add(runButton);
        Button clearButton = new Button("Clear Output");
        clearButton.setOnAction(e -> dataGroup.getChildren().clear());
        getChildren().add(clearButton);
        
    }

    private VBox labeled(String labelText, Control control) {
        Label label = new Label(labelText);
//        label.setMinWidth(180);
        VBox box = new VBox(5, label, control);
        box.setAlignment(Pos.TOP_LEFT);
        return box;
    }

    private Spinner<Integer> createIntegerSpinner(int min, int max, int initial) {
        Spinner<Integer> spinner = new Spinner<>(min, max, initial);
        spinner.setEditable(true);
        return spinner;
    }

    private Spinner<Double> createDoubleSpinner(double min, double max, double initial, double step) {
        Spinner<Double> spinner = new Spinner<>(min, max, initial, step);
        spinner.setEditable(true);
        spinner.getValueFactory().setConverter(new StringConverter<>() {
            @Override public String toString(Double value) {
                return String.format("%.10f", value);
            }
            @Override public Double fromString(String s) {
                try {
                    return Double.parseDouble(s);
                } catch (NumberFormatException e) {
                    return spinner.getValue(); // fallback
                }
            }
        });
        return spinner;
    }

    /**
     * Gathers all current parameter values and runs HoroPCA.
     * Replace this stub with your actual invocation logic.
     */
    private void runHoroPCA() {
    double scale = 50;
    double radius = 1;
    int divisions = 8;
    
        int inputDim = inputDimSpinner.getValue();
        int targetDim = targetDimSpinner.getValue();
        int maxIterations = maxIterationsSpinner.getValue();
        int numClusters = numClustersSpinner.getValue();
        int pointsPerCluster = pointsPerClusterSpinner.getValue();
        double learningRate = learningRateSpinner.getValue();
        double tolerance = toleranceSpinner.getValue();
        InitializationStrategy strategy = initStrategyCombo.getValue();
        boolean verbose = verboseCheckBox.isSelected();
        boolean verboseOptimization = verboseOptCheckBox.isSelected();
        int seed = 42;
        
        System.out.printf("""
            Running HoroPCA with parameters:
              inputDim=%d, targetDim=%d, maxIterations=%d
              clusters=%d, pointsPerCluster=%d
              learningRate=%.10f, tolerance=%.1e
              initStrategy=%s
              verbose=%b, verboseOptimization=%b
            """,
            inputDim, targetDim, maxIterations,
            numClusters, pointsPerCluster,
            learningRate, tolerance,
            strategy, verbose, verboseOptimization
        );

        System.out.println("Generating synthetic clustered data.");
        long startTime = System.nanoTime();
        // Step 1: Generate synthetic clustered data on the Poincaré ball
        PoincareBallFactory factory = new PoincareBallFactory(inputDim, seed);
        // generate anisotropic set
        double[] anisotropicSpread = new double[inputDim];
        for (int i = 0; i < inputDim; i++) {
            anisotropicSpread[i] = Math.pow(0.5, i);  // strong in lower dims
        }        
        List<VectorN> anisotropicDataset = factory.generateAnisotropicClusteredData(
                numClusters, pointsPerCluster, anisotropicSpread);  
        printTotalTime(startTime);

        System.out.println("Fitting HoroPCA model and transforming data...");
        startTime = System.nanoTime();
        HoroPCA horo2 = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
            HoroPCA.InitializationStrategy.KMEANS_PLUS_PLUS, 
                HyperbolicUtils.hyperbolicSquaredDist, seed);
        List<VectorN> transformed2  = horo2.fitTransform(anisotropicDataset);
        printTotalTime(startTime);

        //make some spheres
        dataGroup.getChildren().addAll(
            transformed2.stream()
            .map(v -> {
                Sphere sphere = new Sphere(radius, divisions);
                sphere.setTranslateX(v.get(0)*scale);
                sphere.setTranslateY(-v.get(1)*scale);
                sphere.setTranslateZ(v.get(2)*scale);
                PhongMaterial phong = new PhongMaterial(Color.CYAN);
                sphere.setMaterial(phong);
                return sphere;
            })
            .toList()
        );    
    }
}

