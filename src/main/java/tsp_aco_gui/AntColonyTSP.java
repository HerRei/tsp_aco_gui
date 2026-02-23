package tsp_aco_gui;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A small demonstration program that applies Ant Colony Optimisation (ACO)
 * to solve the Travelling Salesman Problem (TSP).  The program offers a
 * graphical user interface implemented using Swing.  It allows the user to
 * configure the core parameters of the ACO algorithm (number of ants,
 * number of iterations, pheromone influence alpha, heuristic influence beta,
 * pheromone evaporation rate rho and pheromone deposit constant Q) and
 * visualises the resulting route.  In addition to randomly generated
 * cities the program provides a predefined "NYC" example featuring
 * several well-known locations in New York City.  One of the direct
 * connections between two NYC locations can be flagged as closed to
 * simulate a roadblock: this shows how a single blocked edge forces the
 * ants to adapt their tours.
 */
public class AntColonyTSP extends JFrame {
    // Constants defining the geographic bounding box used to project the
    // NYC example onto the background map.  These values were chosen
    // manually to align well with the silhouette of Manhattan on the
    // provided map image.  They roughly correspond to the southern and
    // northern extremes of the island and the western and eastern
    // longitude limits.  If you use a different map image these
    // constants may need to be adjusted.
    // Bounding box for the New York City location map.  These values
    // correspond to the definition used by the Wikipedia location map
    // module (USA_New_York_City_location_map.svg), which is an
    // equirectangular projection covering the five boroughs.  Using
    // this bounding box ensures that latitude/longitude coordinates
    // project correctly onto the provided map image.  See
    // https://en.wikipedia.org/wiki/Module:Location_map/data/USA_New_York_City
    // for details: top=40.92, bottom=40.49, left=-74.27, right=-73.68.
    private static final double NYC_LAT_MIN = 40.49;
    private static final double NYC_LAT_MAX = 40.92;
    private static final double NYC_LON_MIN = -74.27;
    private static final double NYC_LON_MAX = -73.68;

    /**
     * Precomputed relative coordinates (x,y) for the NYC example.  Each entry
     * corresponds to one of the eight locations defined in loadNYCCities().
     * The values specify the fractional position within the background map
     * image: x in [0,1] from left to right, y in [0,1] from top to bottom.
     * These coordinates were derived manually to align with the provided
     * map of New York City.  Using fixed positions avoids issues with
     * inconsistent latitude/longitude scaling on different displays or
     * images.  If you replace the map image, you may need to update
     * these values accordingly.
     */
    /**
     * Updated relative positions for the NYC example.  These values were
     * computed by projecting the latitude/longitude of each landmark
     * into a bounding box covering the Manhattan area (roughly
     * 40.68–40.82°N, −74.02–−73.92°E).  The x coordinate increases to
     * the east and the y coordinate increases downwards after
     * normalising and inverting the latitude.  Using this smaller
     * bounding box ensures that all eight landmarks fall within the
     * red silhouette of Manhattan on the supplied map.  If you wish to
     * adjust the overlay further simply tweak these numbers.
     */
    private static final double[][] NYC_RELATIVE_COORDS = {
            {0.345, 0.443}, // Times Square
            {0.535, 0.277}, // Central Park
            {0.343, 0.511}, // Empire State
            {0.231, 0.814}, // Brooklyn Bridge
            {0.000, 0.934}, // Statue of Liberty (clamped to left edge)
            {0.114, 0.814}, // Wall Street
            {0.428, 0.481}, // Grand Central
            {0.066, 0.766}  // One World Trade
    };
    // GUI fields
    private final DrawingPanel drawingPanel;
    private final JTextField antsField;
    private final JTextField iterationsField;
    private final JTextField alphaField;
    private final JTextField betaField;
    private final JTextField evaporationField;
    private final JTextField qField;
    private final JComboBox<String> datasetCombo;
    private final JCheckBox blockCheck;
    private final JLabel statusLabel;
    // Text field to specify number of cities for random instances.  This
    // value is ignored when the NYC dataset is selected.
    private final JTextField cityCountField;
    private City[] cities;
    private double[][] distances;
    private boolean[][] blocked;
    private int[] bestTour;
    private List<int[]> currentTours;
    // Array of colours used to draw each ant's tour.  This array is
    // initialised when a simulation starts and contains one entry per
    // ant.  If there are more ants than colours, the colours will be
    // reused in a cyclic manner.
    private Color[] antColors;
    // flag indicating whether a simulation is currently running
    private volatile boolean simulationRunning = false;
    // Slider controlling the delay between iterations in milliseconds.  The
    // user can move this during a run to speed up or slow down the
    // animation.  See runSimulation() for how it is used.  The
    // corresponding animationDelay field is updated via a change
    // listener attached to this slider.
    private final JSlider speedSlider;
    // Current delay between iterations.  Declared volatile so that
    // updates from the UI thread are visible to the background
    // simulation thread.
    private volatile int animationDelay = 200;
    /**
     * Background map image for the NYC example.  This image is loaded from
     * the resources folder when the frame is constructed.  The file is
     * derived from the Wikimedia Commons file "New York City location
     * Manhattan.svg" and is licensed under the Creative Commons
     * Attribution 3.0 Unported licence.  The image highlights Manhattan in
     * red and provides geographical context for the NYC example.  If the
     * image cannot be found on the classpath it will be attempted to load
     * from the relative path "resources/nyc_map.png".
     */
    private BufferedImage nycMapImage;
    public AntColonyTSP() {
        super("TSP with Ant Colony Optimisation");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout());
        drawingPanel = new DrawingPanel();
        drawingPanel.setPreferredSize(new Dimension(800, 600));
        add(drawingPanel, BorderLayout.CENTER);

        // Attempt to load the NYC map image from the file system first.  By
        // prioritising the local file before classpath resources we reduce
        // the likelihood of a missing map when running from within an IDE
        // such as VS Code.  The image is expected to reside in a
        // "resources" folder alongside the compiled classes.  If it
        // cannot be found there we try alternative locations and finally
        // fall back to loading from the classpath.
        nycMapImage = null;
        String[] fileNames = {
                "resources/nyc_map.png",
                "tsp_aco_gui/resources/nyc_map.png",
                "nyc_map.png"
        };
        for (String name : fileNames) {
            if (nycMapImage != null) break;
            java.io.File f = new java.io.File(name);
            if (f.exists()) {
                try {
                    nycMapImage = ImageIO.read(f);
                } catch (IOException ignore) {
                    // ignore and try next
                }
            }
        }
        if (nycMapImage == null) {
            // load from classpath if not found on disk
            String[] resourceNames = {
                    "/nyc_map.png",
                    "/resources/nyc_map.png",
                    "/tsp_aco_gui/resources/nyc_map.png",
                    "nyc_map.png",
                    "resources/nyc_map.png",
                    "tsp_aco_gui/resources/nyc_map.png"
            };
            for (String res : resourceNames) {
                if (nycMapImage != null) break;
                try {
                    InputStream in = AntColonyTSP.class.getResourceAsStream(res);
                    if (in != null) {
                        nycMapImage = ImageIO.read(in);
                    }
                } catch (IOException ignore) {
                    // ignore and try next
                }
            }
        }
        if (nycMapImage == null) {
            System.err.println("Warning: could not load NYC map image from any known location. The NYC example will be drawn without a map.");
        }
        // Controls
        JPanel controlPanel = new JPanel();
        controlPanel.setLayout(new GridLayout(0, 2));
        antsField = new JTextField("50", 5);
        iterationsField = new JTextField("100", 5);
        alphaField = new JTextField("1.0", 5);
        betaField = new JTextField("5.0", 5);
        evaporationField = new JTextField("0.5", 5);
        qField = new JTextField("100.0", 5);
        datasetCombo = new JComboBox<>(new String[]{"Random", "NYC Example"});
        blockCheck = new JCheckBox("Block road (NYC)");
        // field to allow the user to set the number of cities for random
        // instances.  Default to 10 cities.  It will be enabled only
        // when the random dataset is selected.
        cityCountField = new JTextField("10", 5);
        // Initially the random dataset is selected, so enable the
        // city count field and disable the road block checkbox because
        // there is no fixed road network in random instances.
        cityCountField.setEnabled(true);
        blockCheck.setEnabled(false);

        // Create a slider to control the delay between iterations.  The
        // minimum delay of 0 ms gives essentially instantaneous updates,
        // while the maximum delay of 1000 ms produces a very slow
        // animation.  The initial value of 200 ms matches the original
        // fixed delay.
        speedSlider = new JSlider(JSlider.HORIZONTAL, 0, 1000, animationDelay);
        speedSlider.setMajorTickSpacing(250);
        speedSlider.setMinorTickSpacing(50);
        speedSlider.setPaintTicks(true);
        speedSlider.setPaintLabels(true);
        // Update the animation delay whenever the slider is adjusted.  This
        // method is called on the event dispatch thread, so using a
        // volatile field ensures visibility to the background thread.
        speedSlider.addChangeListener(e -> {
            animationDelay = speedSlider.getValue();
        });
        controlPanel.add(new JLabel("Number of ants:"));
        controlPanel.add(antsField);
        controlPanel.add(new JLabel("Iterations:"));
        controlPanel.add(iterationsField);
        controlPanel.add(new JLabel("Alpha (pheromone importance):"));
        controlPanel.add(alphaField);
        controlPanel.add(new JLabel("Beta (heuristic importance):"));
        controlPanel.add(betaField);
        controlPanel.add(new JLabel("Evaporation rate (rho):"));
        controlPanel.add(evaporationField);
        controlPanel.add(new JLabel("Q (pheromone deposit):"));
        controlPanel.add(qField);
        controlPanel.add(new JLabel("Dataset:"));
        controlPanel.add(datasetCombo);
        // number of cities field for random dataset
        controlPanel.add(new JLabel("Cities:"));
        controlPanel.add(cityCountField);
        controlPanel.add(blockCheck);
        // Slider to adjust animation speed (delay between iterations)
        controlPanel.add(new JLabel("Delay (ms):"));
        controlPanel.add(speedSlider);
        JButton runButton = new JButton("Run simulation");
        controlPanel.add(runButton);
        JButton newRandomButton = new JButton("Generate Random");
        controlPanel.add(newRandomButton);
        // Status label spans two columns
        statusLabel = new JLabel("Ready");
        controlPanel.add(statusLabel);
        controlPanel.add(new JLabel(""));
        add(controlPanel, BorderLayout.SOUTH);
        // Action listeners
        newRandomButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                // parse the number of cities from the user-input field.
                int n;
                try {
                    n = Integer.parseInt(cityCountField.getText().trim());
                    if (n < 2) {
                        throw new NumberFormatException("Number of cities must be >= 2");
                    }
                } catch (NumberFormatException ex) {
                    JOptionPane.showMessageDialog(AntColonyTSP.this,
                            "Invalid number of cities: " + ex.getMessage(),
                            "Error", JOptionPane.ERROR_MESSAGE);
                    n = 10;
                    cityCountField.setText(Integer.toString(n));
                }
                loadRandomCities(n);
                bestTour = null;
                drawingPanel.repaint();
            }
        });
        datasetCombo.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                String sel = (String) datasetCombo.getSelectedItem();
                if ("NYC Example".equals(sel)) {
                    loadNYCCities();
                    // Disable city count field since the number of cities is fixed
                    cityCountField.setEnabled(false);
                    // Enable road block checkbox for the NYC example
                    blockCheck.setEnabled(true);
                } else {
                    // Random dataset selected. Enable city count field and disable road block
                    cityCountField.setEnabled(true);
                    blockCheck.setEnabled(false);
                    // Determine how many cities to load
                    int n;
                    try {
                        n = Integer.parseInt(cityCountField.getText().trim());
                        if (n < 2) {
                            throw new NumberFormatException("Number of cities must be >= 2");
                        }
                    } catch (NumberFormatException ex2) {
                        JOptionPane.showMessageDialog(AntColonyTSP.this,
                                "Invalid number of cities: " + ex2.getMessage(),
                                "Error", JOptionPane.ERROR_MESSAGE);
                        n = 10;
                        cityCountField.setText(Integer.toString(n));
                    }
                    loadRandomCities(n);
                }
                bestTour = null;
                drawingPanel.repaint();
            }
        });
        runButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if (!simulationRunning) {
                    runSimulation();
                }
            }
        });
        // initial dataset: load the random instance using the value from the city count field.
        int initialCities;
        try {
            initialCities = Integer.parseInt(cityCountField.getText().trim());
            if (initialCities < 2) {
                initialCities = 10;
            }
        } catch (NumberFormatException ex) {
            initialCities = 10;
        }
        loadRandomCities(initialCities);
        pack();
        setLocationRelativeTo(null);
        setVisible(true);
    }

    static void main(String[] args) {
        SwingUtilities.invokeLater(new Runnable() {
            @Override
            public void run() {
                new AntColonyTSP();
            }
        });
    }

    /**
     * Generates a random TSP instance with the given number of cities.
     */
    private void loadRandomCities(int n) {
        Random r = new Random();
        cities = new City[n];
        for (int i = 0; i < n; i++) {
            double x = 50 + r.nextDouble() * (drawingPanel.getWidth() - 100);
            double y = 50 + r.nextDouble() * (drawingPanel.getHeight() - 100);
            cities[i] = new City("C" + i, x, y);
        }
        computeDistances();
        blocked = new boolean[n][n];
    }

    /**
     * Loads the NYC example with fixed coordinates and names.  Distances
     * are computed using the haversine formula to approximate the road
     * network.  A specific edge (Times Square to Central Park) can be
     * blocked using the checkbox.
     */
    private void loadNYCCities() {
        cities = new City[]{
                new City("Times Square", 40.7580, -73.9855),
                new City("Central Park", 40.7812, -73.9665),
                new City("Empire State", 40.7484, -73.9857),
                new City("Brooklyn Bridge", 40.7061, -73.9969),
                new City("Statue of Liberty", 40.6892, -74.0445),
                new City("Wall Street", 40.7060, -74.0086),
                new City("Grand Central", 40.7527, -73.9772),
                new City("One World Trade", 40.7127, -74.0134)
        };
        computeDistancesHaversine();
        blocked = new boolean[cities.length][cities.length];
        // optionally block the edge between Times Square (0) and Central Park (1)
        // the checkbox state will be evaluated at run time
    }

    /**
     * Computes Euclidean distances between cities for the random instance.
     */
    private void computeDistances() {
        int n = cities.length;
        distances = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    distances[i][j] = 0;
                } else {
                    double dx = cities[i].x - cities[j].x;
                    double dy = cities[i].y - cities[j].y;
                    distances[i][j] = Math.hypot(dx, dy);
                }
            }
        }
    }

    /**
     * Computes geodesic distances using the haversine formula for NYC.
     */
    private void computeDistancesHaversine() {
        int n = cities.length;
        distances = new double[n][n];
        double R = 6371.0; // Earth radius in km
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    distances[i][j] = 0;
                } else {
                    double lat1 = Math.toRadians(cities[i].x);
                    double lon1 = Math.toRadians(cities[i].y);
                    double lat2 = Math.toRadians(cities[j].x);
                    double lon2 = Math.toRadians(cities[j].y);
                    double dlat = lat2 - lat1;
                    double dlon = lon2 - lon1;
                    double a = Math.sin(dlat / 2) * Math.sin(dlat / 2) +
                            Math.cos(lat1) * Math.cos(lat2) *
                                    Math.sin(dlon / 2) * Math.sin(dlon / 2);
                    double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
                    distances[i][j] = R * c; // kilometres
                }
            }
        }
    }

    /**
     * Runs the ACO algorithm in simulation mode.  A solver instance is
     * created from the current parameter settings and the selected
     * dataset (random or NYC).  Each iteration constructs tours for
     * every ant, updates pheromones and the best tour, and repaints
     * the drawing panel to visualise progress.  A small delay between
     * iterations allows the user to observe the process.  All work is
     * performed on a background thread to keep the UI responsive.
     */
    private void runSimulation() {
        // parse parameters
        final int ants;
        final int iterations;
        final double alpha;
        final double beta;
        final double evap;
        final double Q;
        try {
            ants = Integer.parseInt(antsField.getText().trim());
            iterations = Integer.parseInt(iterationsField.getText().trim());
            alpha = Double.parseDouble(alphaField.getText().trim());
            beta = Double.parseDouble(betaField.getText().trim());
            evap = Double.parseDouble(evaporationField.getText().trim());
            Q = Double.parseDouble(qField.getText().trim());
        } catch (NumberFormatException ex) {
            JOptionPane.showMessageDialog(this, "Invalid parameter: " + ex.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
            return;
        }
        // prepare blocked edges for NYC
        if ("NYC Example".equals(datasetCombo.getSelectedItem())) {
            for (int i = 0; i < blocked.length; i++) {
                for (int j = 0; j < blocked.length; j++) {
                    blocked[i][j] = false;
                }
            }
            if (blockCheck.isSelected()) {
                blocked[0][1] = true;
                blocked[1][0] = true;
            }
        }
        simulationRunning = true;
        currentTours = null;
        statusLabel.setText("Running...");
        // disable input fields during simulation
        antsField.setEnabled(false);
        iterationsField.setEnabled(false);
        alphaField.setEnabled(false);
        betaField.setEnabled(false);
        evaporationField.setEnabled(false);
        qField.setEnabled(false);
        datasetCombo.setEnabled(false);
        blockCheck.setEnabled(false);
        cityCountField.setEnabled(false);
        // create solver instance
        final ACO solver = new ACO(ants, iterations, alpha, beta, evap, Q, distances, blocked);

        // Generate a distinct colour palette for each ant.  We use
        // uniformly spaced hues in the HSB colour space to ensure
        // clearly distinguishable colours.  Saturation and brightness
        // values are set to 0.7 for vivid pastel tones.  The array
        // length equals the number of ants specified by the user.
        antColors = new Color[ants];
        for (int i = 0; i < ants; i++) {
            float hue = i / (float) ants;
            antColors[i] = Color.getHSBColor(hue, 0.7f, 0.8f);
        }
        bestTour = null;
        // run on new thread
        new Thread(() -> {
            for (int iter = 0; iter < iterations && simulationRunning; iter++) {
                // perform one iteration and capture tours
                currentTours = solver.iterateOnce();
                bestTour = solver.bestTour;
                // update UI on EDT
                int iterDisplay = iter + 1;
                double lengthDisplay = solver.bestLength;
                SwingUtilities.invokeLater(() -> {
                    drawingPanel.repaint();
                    statusLabel.setText(String.format("Iteration %d / %d, best length: %.2f", iterDisplay, iterations, lengthDisplay));
                });
                // delay to visualise progress.  The current animationDelay
                // value is read on each iteration; it can be modified via
                // the speed slider while the simulation is running.
                try {
                    Thread.sleep(animationDelay);
                } catch (InterruptedException e) {
                    // ignore
                }
            }
            // simulation finished
            SwingUtilities.invokeLater(() -> {
                simulationRunning = false;
                statusLabel.setText(String.format("Finished. Best length: %.2f", solver.bestLength));
                // re-enable controls
                antsField.setEnabled(true);
                iterationsField.setEnabled(true);
                alphaField.setEnabled(true);
                betaField.setEnabled(true);
                evaporationField.setEnabled(true);
                qField.setEnabled(true);
                datasetCombo.setEnabled(true);
                blockCheck.setEnabled(true);
                // re-enable city count field only when the random dataset is selected
                if ("Random".equals(datasetCombo.getSelectedItem())) {
                    cityCountField.setEnabled(true);
                    // disable road block checkbox for random dataset
                    blockCheck.setEnabled(false);
                }
                JOptionPane.showMessageDialog(AntColonyTSP.this,
                        String.format("Best tour length: %.2f", solver.bestLength),
                        "ACO result", JOptionPane.INFORMATION_MESSAGE);
            });
        }).start();
    }

    /**
         * Represents a city with a name and geographic coordinates.  The
         * x/y values are interpreted as screen coordinates when drawing
         * random instances.  For the NYC example they correspond to
         * latitude and longitude values (degrees) which are rescaled for
         * display.
         */
        record City(String name, double x, double y) {
    }

    /**
     * Implements the core Ant Colony Optimisation algorithm for the
     * travelling salesman problem.  The algorithm maintains a matrix of
     * pheromone values on every edge, a population of ants that
     * construct tours, and updates pheromones based on tour quality.
     */
    static class ACO {
        final int numAnts;
        final int maxIterations;
        final double alpha;
        final double beta;
        final double evaporation;
        final double Q;
        final double[][] distances;
        final boolean[][] blocked; // edges flagged as blocked
        final int numCities;
        // random number generator used during construction of tours.  Note that
        // this instance is not used within parallel regions; thread-local
        // random is employed there for safety.
        final Random rng = new Random();
        double[][] pheromones;
        int[] bestTour;
        double bestLength = Double.MAX_VALUE;

        ACO(int numAnts, int maxIterations, double alpha, double beta,
            double evaporation, double Q, double[][] distances,
            boolean[][] blocked) {
            this.numAnts = numAnts;
            this.maxIterations = maxIterations;
            this.alpha = alpha;
            this.beta = beta;
            this.evaporation = evaporation;
            this.Q = Q;
            this.distances = distances;
            this.blocked = blocked;
            this.numCities = distances.length;
            this.pheromones = new double[numCities][numCities];
            // initialise pheromones with a small constant
            for (int i = 0; i < numCities; i++) {
                for (int j = 0; j < numCities; j++) {
                    pheromones[i][j] = 1.0;
                }
            }
        }

        /**
         * Executes the ACO optimisation until completion and returns the best tour found.
         * This method is retained for compatibility but no longer used in the simulation.
         */
        public int[] solve() {
            for (int iter = 0; iter < maxIterations; iter++) {
                iterateOnce();
            }
            return bestTour;
        }

        /**
         * Performs a single iteration of the ACO algorithm.  It constructs a tour
         * for every ant, updates the global pheromone matrix, updates the
         * best solution if a shorter tour is found, and returns the list of
         * tours constructed in this iteration.  This method allows step-by-step
         * simulation of the algorithm.
         *
         * @return a list of tours, where each tour is represented as an array of
         * city indices.
         */
        public List<int[]> iterateOnce() {
            // Arrays to hold tours and lengths for each ant.  Using fixed
            // arrays avoids synchronization overhead on collections.  Each
            // ant builds its tour independently in a parallel stream.
            int[][] tours = new int[numAnts][numCities];
            double[] lengths = new double[numAnts];
            // Construct tours in parallel.  We avoid using the shared rng
            // within the parallel region; instead we rely on ThreadLocalRandom
            // for thread-safe random number generation.  The constructTour()
            // method does not modify shared state and therefore can be
            // safely invoked concurrently.
            java.util.stream.IntStream.range(0, numAnts).parallel().forEach(i -> {
                int start = java.util.concurrent.ThreadLocalRandom.current().nextInt(numCities);
                int[] tour = constructTour(start);
                double len = computeTourLength(tour);
                tours[i] = tour;
                lengths[i] = len;
            });
            // Update the global best tour and length in a single-threaded
            // section to avoid race conditions.  We identify the best tour
            // from this iteration and compare it to the historical best.
            int bestIndex = -1;
            double bestLenThisIter = Double.MAX_VALUE;
            for (int i = 0; i < numAnts; i++) {
                if (lengths[i] < bestLenThisIter) {
                    bestLenThisIter = lengths[i];
                    bestIndex = i;
                }
            }
            if (bestIndex >= 0 && bestLenThisIter < bestLength) {
                bestLength = bestLenThisIter;
                bestTour = tours[bestIndex].clone();
            }
            // Convert arrays to lists for the pheromone update and for the
            // caller to draw the ant paths.  Lists are used by the
            // existing updatePheromones() method.
            List<int[]> antTours = new ArrayList<>(numAnts);
            List<Double> antLengths = new ArrayList<>(numAnts);
            for (int i = 0; i < numAnts; i++) {
                antTours.add(tours[i]);
                antLengths.add(lengths[i]);
            }
            // Update pheromones based on the tours constructed in this
            // iteration.  This method remains sequential to avoid
            // concurrent writes to the pheromone matrix.
            updatePheromones(antTours, antLengths);
            return antTours;
        }

        /**
         * Constructs a tour for a single ant starting at the given city.
         */
        private int[] constructTour(int startCity) {
            int[] tour = new int[numCities];
            boolean[] visited = new boolean[numCities];
            tour[0] = startCity;
            visited[startCity] = true;
            for (int step = 1; step < numCities; step++) {
                int current = tour[step - 1];
                int next = selectNextCity(current, visited);
                // if no valid next city is found (all are blocked), return
                // incomplete path to avoid infinite loops
                if (next == -1) {
                    // fill remaining with an arbitrary order of unvisited
                    int idx = step;
                    for (int c = 0; c < numCities; c++) {
                        if (!visited[c]) {
                            tour[idx++] = c;
                            visited[c] = true;
                        }
                    }
                    break;
                }
                tour[step] = next;
                visited[next] = true;
            }
            return tour;
        }

        /**
         * Chooses the next city using a roulette wheel based on pheromone and
         * heuristic information.  Returns -1 if no unvisited city is reachable.
         */
        private int selectNextCity(int current, boolean[] visited) {
            double[] probs = new double[numCities];
            double sum = 0.0;
            for (int j = 0; j < numCities; j++) {
                if (!visited[j] && !blocked[current][j]) {
                    double tau = Math.pow(pheromones[current][j], alpha);
                    double eta = Math.pow(1.0 / (distances[current][j] + 1e-6), beta);
                    probs[j] = tau * eta;
                    sum += probs[j];
                } else {
                    probs[j] = 0.0;
                }
            }
            if (sum == 0) {
                return -1;
            }
            // Use ThreadLocalRandom for thread-safe random number generation
            double r = java.util.concurrent.ThreadLocalRandom.current().nextDouble() * sum;
            double cumulative = 0.0;
            for (int j = 0; j < numCities; j++) {
                cumulative += probs[j];
                if (r <= cumulative) {
                    return j;
                }
            }
            // fallback
            for (int j = 0; j < numCities; j++) {
                if (probs[j] > 0) return j;
            }
            return -1;
        }

        /**
         * Computes the length of a full tour (including returning to the start).
         */
        private double computeTourLength(int[] tour) {
            double length = 0.0;
            // accumulate distance between consecutive cities
            for (int i = 0; i < numCities - 1; i++) {
                int a = tour[i];
                int b = tour[i + 1];
                // If this edge is blocked, treat the tour as invalid by
                // returning a very large cost.  This should never happen
                // during construction because selectNextCity avoids blocked
                // edges, but we check for completeness.
                if (blocked[a][b]) {
                    return Double.POSITIVE_INFINITY;
                }
                length += distances[a][b];
            }
            // handle the return to the start city; if the closing edge is
            // blocked, assign an infinite cost so that ants cannot use
            // blocked roads to complete the cycle
            int last = tour[numCities - 1];
            int first = tour[0];
            if (blocked[last][first]) {
                return Double.POSITIVE_INFINITY;
            }
            length += distances[last][first];
            return length;
        }

        /**
         * Applies pheromone evaporation and deposits pheromones based on each
         * ant's tour length.
         */
        private void updatePheromones(List<int[]> tours, List<Double> lengths) {
            // Evaporation
            for (int i = 0; i < numCities; i++) {
                for (int j = 0; j < numCities; j++) {
                    pheromones[i][j] = (1.0 - evaporation) * pheromones[i][j];
                }
            }
            // Deposit
            for (int k = 0; k < tours.size(); k++) {
                int[] tour = tours.get(k);
                double len = lengths.get(k);
                double deposit = Q / len;
                for (int idx = 0; idx < numCities - 1; idx++) {
                    int i = tour[idx];
                    int j = tour[idx + 1];
                    pheromones[i][j] += deposit;
                    pheromones[j][i] += deposit;
                }
                // include returning to start
                int i = tour[numCities - 1];
                int j = tour[0];
                pheromones[i][j] += deposit;
                pheromones[j][i] += deposit;
            }
        }
    }

    /**
     * Panel responsible for drawing the cities and the best tour.
     */
    class DrawingPanel extends JPanel {
        private final int CITY_RADIUS = 8;

        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);
            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            if (cities == null) return;

            // Draw the NYC map in the background for the NYC example.  The
            // map is scaled to maintain aspect ratio and cover the entire
            // drawing area.  This is done before computing city positions
            // so that the cities and tours are drawn on top of the map.
            if ("NYC Example".equals(datasetCombo.getSelectedItem()) && nycMapImage != null) {
                int panelWidth = getWidth();
                int panelHeight = getHeight();
                int imgWidth = nycMapImage.getWidth();
                int imgHeight = nycMapImage.getHeight();
                double imgRatio = (double) imgWidth / imgHeight;
                double panelRatio = (double) panelWidth / panelHeight;
                int drawWidth;
                int drawHeight;
                if (panelRatio > imgRatio) {
                    // Panel is relatively wider than the image: scale by height
                    drawHeight = panelHeight;
                    drawWidth = (int) (panelHeight * imgRatio);
                } else {
                    // Panel is relatively taller: scale by width
                    drawWidth = panelWidth;
                    drawHeight = (int) (panelWidth / imgRatio);
                }
                // Center the image
                int xOffset = (panelWidth - drawWidth) / 2;
                int yOffset = (panelHeight - drawHeight) / 2;
                g2.drawImage(nycMapImage, xOffset, yOffset, drawWidth, drawHeight, null);
            }
            /*
             * Compute scaling for city coordinates. For randomly generated
             * instances we simply normalise the x/y values to the panel
             * dimensions with a margin. For the NYC example we instead
             * normalise the latitude/longitude values relative to the
             * geographic bounding box of the cities and map them into
             * the area where the NYC map image is drawn. This ensures
             * the points and labels align correctly with the map. The
             * coordinate system of the image has its origin at the top
             * left, x increasing to the right and y increasing downwards.
             */
            Point[] positions = new Point[cities.length];
            if ("NYC Example".equals(datasetCombo.getSelectedItem()) && nycMapImage != null) {
                // If we have precomputed relative positions for the NYC
                // example and the number of cities matches, use these
                // coordinates.  This provides a more accurate overlay on
                // the map image than computing from latitude/longitude
                // values at runtime.  If the array length does not match,
                // fall back to the bounding-box projection below.
                if (NYC_RELATIVE_COORDS.length == cities.length) {
                    int panelWidth = getWidth();
                    int panelHeight = getHeight();
                    int imgW = nycMapImage.getWidth();
                    int imgH = nycMapImage.getHeight();
                    double imgRatio = (double) imgW / imgH;
                    double panelRatio = (double) panelWidth / panelHeight;
                    int drawW;
                    int drawH;
                    if (panelRatio > imgRatio) {
                        drawH = panelHeight;
                        drawW = (int) (panelHeight * imgRatio);
                    } else {
                        drawW = panelWidth;
                        drawH = (int) (panelWidth / imgRatio);
                    }
                    int xOff = (panelWidth - drawW) / 2;
                    int yOff = (panelHeight - drawH) / 2;
                    for (int i = 0; i < cities.length; i++) {
                        double relX = NYC_RELATIVE_COORDS[i][0];
                        double relY = NYC_RELATIVE_COORDS[i][1];
                        int px = (int) (xOff + relX * drawW);
                        int py = (int) (yOff + relY * drawH);
                        positions[i] = new Point(px, py);
                    }
                } else {
                    // Use a fixed geographic bounding box to project the
                    // NYC coordinates onto the map image.  The constants
                    // are defined at the top of the class.  Clamp values
                    // to [0,1] so that outliers stay within the drawable
                    // area.
                    double latRange = NYC_LAT_MAX - NYC_LAT_MIN;
                    double lonRange = NYC_LON_MAX - NYC_LON_MIN;
                    int panelWidth = getWidth();
                    int panelHeight = getHeight();
                    int imgW = nycMapImage.getWidth();
                    int imgH = nycMapImage.getHeight();
                    double imgRatio = (double) imgW / imgH;
                    double panelRatio = (double) panelWidth / panelHeight;
                    int drawW;
                    int drawH;
                    if (panelRatio > imgRatio) {
                        drawH = panelHeight;
                        drawW = (int) (panelHeight * imgRatio);
                    } else {
                        drawW = panelWidth;
                        drawH = (int) (panelWidth / imgRatio);
                    }
                    int xOff = (panelWidth - drawW) / 2;
                    int yOff = (panelHeight - drawH) / 2;
                    for (int i = 0; i < cities.length; i++) {
                        double lat = cities[i].x;
                        double lon = cities[i].y;
                        double lonNorm = (lon - NYC_LON_MIN) / lonRange;
                        if (lonNorm < 0.0) lonNorm = 0.0;
                        if (lonNorm > 1.0) lonNorm = 1.0;
                        double latNorm = (lat - NYC_LAT_MIN) / latRange;
                        if (latNorm < 0.0) latNorm = 0.0;
                        if (latNorm > 1.0) latNorm = 1.0;
                        double invLatNorm = 1.0 - latNorm;
                        int px = (int) (xOff + lonNorm * drawW);
                        int py = (int) (yOff + invLatNorm * drawH);
                        positions[i] = new Point(px, py);
                    }
                }
            } else {
                // random instance: simple normalisation into panel with margins
                double minX = Double.MAX_VALUE, maxX = Double.MIN_VALUE;
                double minY = Double.MAX_VALUE, maxY = Double.MIN_VALUE;
                for (City c : cities) {
                    if (c.x < minX) minX = c.x;
                    if (c.x > maxX) maxX = c.x;
                    if (c.y < minY) minY = c.y;
                    if (c.y > maxY) maxY = c.y;
                }
                double panelW = getWidth() - 100;
                double panelH = getHeight() - 100;
                for (int i = 0; i < cities.length; i++) {
                    double normX = (cities[i].x - minX) / (maxX - minX);
                    double normY = (cities[i].y - minY) / (maxY - minY);
                    int px = (int) (50 + normX * panelW);
                    int py = (int) (50 + normY * panelH);
                    positions[i] = new Point(px, py);
                }
            }
            // draw edges for reference
            g2.setColor(new Color(220, 220, 220));
            for (int i = 0; i < cities.length; i++) {
                for (int j = i + 1; j < cities.length; j++) {
                    // highlight blocked road segments with a thicker red line
                    if (blocked != null && blocked[i][j]) {
                        // Save current stroke to restore later
                        Stroke oldStroke = g2.getStroke();
                        // Draw a thick red line to emphasise the road block
                        g2.setStroke(new BasicStroke(4.0f));
                        g2.setColor(Color.RED);
                        g2.drawLine(positions[i].x, positions[i].y, positions[j].x, positions[j].y);
                        // Restore previous stroke and colour for subsequent edges
                        g2.setStroke(oldStroke);
                        g2.setColor(new Color(200, 200, 200));
                    } else {
                        g2.drawLine(positions[i].x, positions[i].y, positions[j].x, positions[j].y);
                    }
                }
            }
            // draw best tour
            // draw ant tours from current iteration.  Each ant is
            // rendered with its own colour to make the individual tours
            // distinguishable.  If the antColours array has fewer
            // entries than the number of ants, colours will be reused.
            if (simulationRunning && currentTours != null) {
                g2.setStroke(new BasicStroke(1));
                // iterate through tours with index so we can assign
                // colours per ant
                for (int antIndex = 0; antIndex < currentTours.size(); antIndex++) {
                    int[] tour = currentTours.get(antIndex);
                    // choose colour for this ant.  If antColours is
                    // unavailable fall back to a default light gray.
                    Color c;
                    if (antColors != null && antColors.length > 0) {
                        c = antColors[antIndex % antColors.length];
                        // apply some transparency so the ant tours
                        // remain visible without overpowering other
                        // drawings.  Alpha value 150 gives ~60%
                        // opacity.
                        c = new Color(c.getRed(), c.getGreen(), c.getBlue(), 150);
                    } else {
                        c = new Color(190, 190, 190, 150);
                    }
                    g2.setColor(c);
                    for (int i = 0; i < tour.length - 1; i++) {
                        int a = tour[i];
                        int b = tour[i + 1];
                        g2.drawLine(positions[a].x, positions[a].y, positions[b].x, positions[b].y);
                    }
                    // close tour
                    int last = tour[tour.length - 1];
                    int first = tour[0];
                    g2.drawLine(positions[last].x, positions[last].y, positions[first].x, positions[first].y);
                }
            }
            // draw best tour in blue
            if (bestTour != null) {
                g2.setStroke(new BasicStroke(3));
                g2.setColor(Color.BLUE);
                for (int i = 0; i < bestTour.length - 1; i++) {
                    int a = bestTour[i];
                    int b = bestTour[i + 1];
                    g2.drawLine(positions[a].x, positions[a].y, positions[b].x, positions[b].y);
                }
                // return to start
                int last = bestTour[bestTour.length - 1];
                int first = bestTour[0];
                g2.drawLine(positions[last].x, positions[last].y, positions[first].x, positions[first].y);
            }
            // draw cities
            for (int i = 0; i < cities.length; i++) {
                Point p = positions[i];
                Shape circle = new Ellipse2D.Double(p.x - CITY_RADIUS / 2.0, p.y - CITY_RADIUS / 2.0, CITY_RADIUS, CITY_RADIUS);
                g2.setColor(Color.BLACK);
                g2.fill(circle);
                g2.drawString(cities[i].name, p.x + CITY_RADIUS, p.y - CITY_RADIUS);
            }
        }
    }
}