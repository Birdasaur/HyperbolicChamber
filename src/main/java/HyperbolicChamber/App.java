package HyperbolicChamber;

import static HyperbolicChamber.HoroPCAValidator.printTotalTime;
import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Point3D;
import javafx.scene.AmbientLight;
import javafx.scene.Group;
import javafx.scene.PerspectiveCamera;
import javafx.scene.Scene;
import javafx.scene.SceneAntialiasing;
import javafx.scene.SubScene;
import javafx.scene.input.KeyCode;
import javafx.scene.input.MouseEvent;
import javafx.scene.input.ScrollEvent;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.CornerRadii;
import javafx.scene.layout.StackPane;
import javafx.scene.paint.Color;
import javafx.scene.paint.PhongMaterial;
import javafx.scene.shape.Sphere;
import javafx.stage.Stage;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class App extends Application {

    PerspectiveCamera camera = new PerspectiveCamera(true);
    public Group sceneRoot = new Group();
    public SubScene subScene;
    public CameraTransformer cameraTransform = new CameraTransformer();
    private double cameraDistance = -500;
    private final double sceneWidth = 4000;
    private final double sceneHeight = 4000;

    private double mousePosX;
    private double mousePosY;
    private double mouseOldX;
    private double mouseOldY;
    private double mouseDeltaX;
    private double mouseDeltaY;

    ArrayList<Point3D> positions;
    Random rando = new Random();
    double scale = 50;
    int totalPoints = 10;
    double radius = 1;
    int divisions = 8;
    Group dataGroup;
    
    @Override
    public void start(Stage primaryStage) throws Exception {
        dataGroup = new Group();
        subScene = new SubScene(sceneRoot, sceneWidth, sceneHeight, true, SceneAntialiasing.BALANCED);
        //Start Tracking mouse movements only when a button is pressed
        subScene.setOnMousePressed((MouseEvent me) -> {
            mousePosX = me.getSceneX();
            mousePosY = me.getSceneY();
            mouseOldX = me.getSceneX();
            mouseOldY = me.getSceneY();
        });
        subScene.setOnMouseDragged((MouseEvent me) -> mouseDragCamera(me));
        subScene.setOnScroll((ScrollEvent event) -> {
            double modifier = 2.0;
            double modifierFactor = 0.1;

            if (event.isControlDown()) {
                modifier = 1;
            }
            if (event.isShiftDown()) {
                modifier = 50.0;
            }
            double z = camera.getTranslateZ();
            double newZ = z + event.getDeltaY() * modifierFactor * modifier;
            camera.setTranslateZ(newZ);
        });
        StackPane stackPane = new StackPane(subScene);
        subScene.widthProperty().bind(stackPane.widthProperty());
        subScene.heightProperty().bind(stackPane.heightProperty());
        subScene.setFill(Color.BLACK);

        camera = new PerspectiveCamera(true);
        //setup camera transform for rotational support
        cameraTransform.setTranslate(0, 0, 0);
        cameraTransform.getChildren().add(camera);
        camera.setNearClip(0.1);
        camera.setFarClip(100000.0);
        camera.setTranslateZ(cameraDistance);
        cameraTransform.ry.setAngle(-45.0);
        cameraTransform.rx.setAngle(-10.0);
        subScene.setCamera(camera);

        AmbientLight ambientLight = new AmbientLight(Color.WHITE);

        Sphere sphereOrigin = new Sphere(5);
        sphereOrigin.setMaterial(new PhongMaterial(Color.GRAY.deriveColor(1, 1, 1, 0.1)));
        
        
        Sphere sphereX = new Sphere(5);
        sphereX.setTranslateX(scale);
        sphereX.setMaterial(new PhongMaterial(Color.RED));

        Sphere sphereY = new Sphere(5);
        sphereY.setTranslateY(-scale);
        sphereY.setMaterial(new PhongMaterial(Color.GREEN));

        Sphere sphereZ = new Sphere(5);
        sphereZ.setTranslateZ(scale);
        sphereZ.setMaterial(new PhongMaterial(Color.BLUE));

        sceneRoot.getChildren().addAll(cameraTransform, ambientLight, 
            sphereOrigin, sphereX, sphereY, sphereZ, dataGroup);

        subScene.setOnKeyPressed(event -> {
            //What key did the user press?
            KeyCode keycode = event.getCode();

            double change = 10.0;
            //Add shift modifier to simulate "Running Speed"
            if (event.isShiftDown()) {
                change = 100.0;
            }

            //Zoom controls
            if (keycode == KeyCode.W) {
                camera.setTranslateZ(camera.getTranslateZ() + change);
            }
            if (keycode == KeyCode.S) {
                camera.setTranslateZ(camera.getTranslateZ() - change);
            }
            //Strafe controls
            if (keycode == KeyCode.A) {
                camera.setTranslateX(camera.getTranslateX() - change);
            }
            if (keycode == KeyCode.D) {
                camera.setTranslateX(camera.getTranslateX() + change);
            }

            if (keycode == KeyCode.SPACE) {
                camera.setTranslateY(camera.getTranslateY() - change);
            }
            if (keycode == KeyCode.C) {
                camera.setTranslateY(camera.getTranslateY() + change);
            }
        });

        BorderPane bpOilSpill = new BorderPane(stackPane);
        HoroPCAControlPanel controls = new HoroPCAControlPanel(dataGroup);
        bpOilSpill.setLeft(controls);
//        stackPane.getChildren().clear();
//        stackPane.getChildren().addAll(bpOilSpill);
//        stackPane.setPadding(new Insets(10));
//        stackPane.setBackground(new Background(new BackgroundFill(Color.rgb(255, 255, 255), CornerRadii.EMPTY, Insets.EMPTY)));
        Scene scene = new Scene(bpOilSpill, 1000, 1000);
        scene.setOnMouseEntered(event -> subScene.requestFocus());

        primaryStage.setTitle("Hyperbolic Chamber");
        primaryStage.setScene(scene);
        primaryStage.show();

        genEuclid2Hyper();
//        runHoroPCA();        
    }
//    private void runHoroPCA() {
//        int inputDim = 100;
//        int targetDim = 3;
//        int maxIterations = 300;
//        int numClusters = 4;
//        int pointsPerCluster = 50;
//        int seed = 42;
//        double learningRate = 0.01;
//        double tolerance = 1e-8;
//        System.out.println("Generating synthetic clustered data.");
//        long startTime = System.nanoTime();
//        // Step 1: Generate synthetic clustered data on the Poincaré ball
//        PoincareBallFactory factory = new PoincareBallFactory(inputDim, seed);
//        // generate anisotropic set
//        double[] anisotropicSpread = new double[inputDim];
//        for (int i = 0; i < inputDim; i++) {
//            anisotropicSpread[i] = Math.pow(0.5, i);  // strong in lower dims
//        }        
//        List<VectorN> anisotropicDataset = factory.generateAnisotropicClusteredData(
//                numClusters, pointsPerCluster, anisotropicSpread);  
//        printTotalTime(startTime);
//
//        System.out.println("Fitting HoroPCA model and transforming data...");
//        startTime = System.nanoTime();
//        HoroPCA horo2 = new HoroPCA(targetDim, maxIterations, learningRate, tolerance, 
//            HoroPCA.InitializationStrategy.KMEANS_PLUS_PLUS, 
//                HyperbolicUtils.hyperbolicSquaredDist, seed);
//        List<VectorN> transformed2  = horo2.fitTransform(anisotropicDataset);
//        printTotalTime(startTime);
//
//        //make some spheres
//        dataGroup.getChildren().addAll(
//            transformed2.stream()
//            .map(v -> {
//                Sphere sphere = new Sphere(radius, divisions);
//                sphere.setTranslateX(v.get(0)*scale);
//                sphere.setTranslateY(-v.get(1)*scale);
//                sphere.setTranslateZ(v.get(2)*scale);
//                PhongMaterial phong = new PhongMaterial(Color.CYAN);
//                sphere.setMaterial(phong);
//                return sphere;
//            })
//            .toList()
//        );        
//    }
    private void genEuclid2Hyper() {
        int pointCount = 1000;
        positions = new ArrayList<>(pointCount);
        
        //generate some random positions
        double randomScale = 1; //scale/2.0;
        double alpha = 0.333;
        for (int i = 0; i < pointCount; i++) {
            Point3D p3D = new Point3D(
//                    rando.nextDouble() * scale,
//                    rando.nextDouble() * -scale,
//                    rando.nextDouble() * scale);
//                rando.nextGaussian() * randomScale,
//                rando.nextGaussian() * randomScale,
//                rando.nextGaussian() * randomScale);
                rando.nextExponential() * randomScale,
                rando.nextExponential() * randomScale,
                rando.nextExponential() * randomScale);

            double [] hyper = 
//                toPoincareBall3D(
                toHyperboloid(
                new double[] { p3D.getX(), p3D.getY(), p3D.getZ() }
//                ,alpha
                )
            ;            
            
            positions.add(p3D);
            
            Sphere sphere = new Sphere(radius, divisions);
            sphere.setTranslateX(hyper[0]*scale);
            sphere.setTranslateY(-hyper[1]*scale);
            sphere.setTranslateZ(hyper[2]*scale);
            sceneRoot.getChildren().add(sphere);
        }           
    }
    
    /**
     * Converts an N-dimensional Euclidean vector to a 3D point on the Poincare hyperboloid model.
     *
     * @param euclideanVector N-dimensional Euclidean input vector
     * @return double[3] representing the (x0, x1, x2) on the hyperboloid
     */
    public static double[] toHyperboloid(double[] euclideanVector) {
        if (euclideanVector == null || euclideanVector.length < 2) {
            throw new IllegalArgumentException("Input must be at least 2-dimensional.");
        }

        // Project to 2D plane (take first two dimensions)
        double x = euclideanVector[0];
        double y = euclideanVector[1];

        // Compute radius
        double r = Math.sqrt(x * x + y * y);

        // Handle zero vector
        if (r == 0) {
            return new double[] {1.0, 0.0, 0.0}; // Origin of hyperboloid
        }

        // Compute angle
        double theta = Math.atan2(y, x);

        // Compute hyperboloid coordinates
        double x0 = Math.sqrt(1 + r * r);
        double x1 = r * Math.cos(theta);
        double x2 = r * Math.sin(theta);

        return new double[] {x0, x1, x2};
    }    
    /**
     * Converts an N-dimensional Euclidean vector to a 4D hyperbolic sphere (3D spatial output + time).
     * Returns a point on the upper sheet of a hyperboloid: -x0^2 + x1^2 + x2^2 + x3^2 = -R^2
     *
     * @param euclideanVector N-dimensional input
     * @param R Hyperbolic radius (default is 1.0)
     * @return double[4] representing (x0, x1, x2, x3)
     */
    public static double[] toHyperbolicSphere(double[] euclideanVector, double R) {
        if (euclideanVector == null || euclideanVector.length < 3) {
            throw new IllegalArgumentException("Input must be at least 3-dimensional.");
        }

        // Use first 3 components (or you can apply PCA here)
        double x1 = euclideanVector[0];
        double x2 = euclideanVector[1];
        double x3 = euclideanVector[2];

        // Compute r^2 = x1^2 + x2^2 + x3^2
        double r2 = x1 * x1 + x2 * x2 + x3 * x3;

        // Compute x0 to satisfy hyperboloid equation: -x0^2 + r^2 = -R^2 => x0 = sqrt(R^2 + r^2)
        double x0 = Math.sqrt(R * R + r2);

        return new double[] {x0, x1, x2, x3};
    }

    // Overloaded method with R = 1
    public static double[] toHyperbolicSphere(double[] euclideanVector) {
        return toHyperbolicSphere(euclideanVector, 1.0);
    }
    /**
     * Converts an N-dimensional Euclidean vector to a point in the 3D Poincare ball model.
     * Reduces to 3D and projects inside the unit ball.
     *
     * @param euclideanVector N-dimensional input
     * @param alpha Scaling factor to control "squash" toward the unit ball (typical: 0.5–2)
     * @return double[3] representing (x, y, z) inside the Poincare ball
     */
    public static double[] toPoincareBall3D(double[] euclideanVector, double alpha) {
        if (euclideanVector == null || euclideanVector.length < 3) {
            throw new IllegalArgumentException("Input must be at least 3-dimensional.");
        }

        // Use first 3 components for spatial projection
        double x = euclideanVector[0];
        double y = euclideanVector[1];
        double z = euclideanVector[2];

        double norm = Math.sqrt(x * x + y * y + z * z);

        // Handle zero vector case
        if (norm == 0) {
            return new double[] {0.0, 0.0, 0.0};
        }

        // Apply tanh squashing to ensure inside unit ball
        double scale = Math.tanh(alpha * norm) / norm;

        return new double[] {
            x * scale,
            y * scale,
            z * scale
        };
    }

    // Overloaded with default alpha = 1.0
    public static double[] toPoincareBall3D(double[] euclideanVector) {
        return toPoincareBall3D(euclideanVector, 1.0);
    }
    
    /**
     * Converts an N-dimensional Euclidean vector to a 3D point in the Klein model of hyperbolic space.
     *
     * @param euclideanVector N-dimensional input vector
     * @return double[3] representing (x, y, z) inside the Klein ball
     */
    public static double[] toKleinModel3D(double[] euclideanVector) {
        if (euclideanVector == null || euclideanVector.length < 3) {
            throw new IllegalArgumentException("Input must be at least 3-dimensional.");
        }

        // Use first 3 components
        double x = euclideanVector[0];
        double y = euclideanVector[1];
        double z = euclideanVector[2];

        double norm = Math.sqrt(x * x + y * y + z * z);

        // Handle the zero vector (maps to center of ball)
        if (norm == 0.0) {
            return new double[] {0.0, 0.0, 0.0};
        }

        // Scale to ensure the point is inside the unit ball using the Klein model formula
        double scale = 1.0 / (1.0 + norm);

        return new double[] {
            x * scale,
            y * scale,
            z * scale
        };
    }    

    private void mouseDragCamera(MouseEvent me) {
        mouseOldX = mousePosX;
        mouseOldY = mousePosY;
        mousePosX = me.getSceneX();
        mousePosY = me.getSceneY();
        mouseDeltaX = (mousePosX - mouseOldX);
        mouseDeltaY = (mousePosY - mouseOldY);
        double modifier = 1.0;
        double modifierFactor = 0.1;

        if (me.isControlDown()) {
            modifier = 0.1;

        }
        if (me.isShiftDown()) {
            modifier = 10.0;
        }
        if (me.isPrimaryButtonDown()) {
            if (me.isAltDown()) { //roll
                cameraTransform.rz.setAngle(
                        ((cameraTransform.rz.getAngle() + mouseDeltaX * modifierFactor * modifier * 2.0) % 360 + 540) % 360 - 180); // +
            } else {
                cameraTransform.ry.setAngle(
                        ((cameraTransform.ry.getAngle() + mouseDeltaX * modifierFactor * modifier * 2.0) % 360 + 540) % 360 - 180); // +
                cameraTransform.rx.setAngle(
                        ((cameraTransform.rx.getAngle() - mouseDeltaY * modifierFactor * modifier * 2.0) % 360 + 540) % 360 - 180);
            }
        } else if (me.isMiddleButtonDown()) {
            cameraTransform.t.setX(cameraTransform.t.getX() + mouseDeltaX * modifierFactor * modifier * 0.3); // -
            cameraTransform.t.setY(cameraTransform.t.getY() + mouseDeltaY * modifierFactor * modifier * 0.3); // -
        }
    }

    public static void main(String[] args) {
        Application.launch(App.class, args);
    }
}
