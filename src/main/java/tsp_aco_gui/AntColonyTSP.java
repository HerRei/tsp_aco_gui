package tsp_aco_gui;

import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

public class AntColonyTSP extends JFrame {

    // NYC Boundary constants (Manhattan-ish)
    private static final double NYC_LAT_MIN = 40.49, NYC_LAT_MAX = 40.92;
    private static final double NYC_LON_MIN = -74.27, NYC_LON_MAX = -73.68;

    // Relative coords for the specific background map
    private static final double[][] NYC_REL_COORDS = {
            {0.345, 0.443}, {0.535, 0.277}, {0.343, 0.511}, {0.231, 0.814},
            {0.000, 0.934}, {0.114, 0.814}, {0.428, 0.481}, {0.066, 0.766}
    };

    private final DrawingPanel canvas;

    private final JTextField tfAnts;
    private final JTextField tfIter;
    private final JTextField tfAlpha;
    private final JTextField tfBeta;
    private final JTextField tfRho;
    private final JTextField tfQ;
    private final JTextField tfCities;

    private final JComboBox<String> cbDataset;
    private final JCheckBox chkBlock;
    private final JLabel lblStatus;
    private final JSlider sliderSpeed;

    private City[] cities;
    private double[][] distMatrix;
    private boolean[][] blocked;

    private int[] bestTour;
    private double bestLen = Double.MAX_VALUE;

    private List<int[]> currentTours;
    private Color[] antColors;

    private volatile boolean isRunning = false;
    private volatile int delayMs = 200;

    private BufferedImage bgMap;

    public AntColonyTSP() {
        super("ACO TSP Solver");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        canvas = new DrawingPanel();
        canvas.setPreferredSize(new Dimension(800, 600));
        add(canvas, BorderLayout.CENTER);

        // background map is optional, keep going if it doesn't load
        try {
            var imgStream = getClass().getResourceAsStream("/nyc_map.png");
            if (imgStream != null) bgMap = ImageIO.read(imgStream);
        } catch (IOException e) {
            System.err.println("Map load failed: " + e.getMessage());
        }

        JPanel bottomPanel = new JPanel(new GridLayout(0, 2, 5, 5));

        tfAnts = new JTextField("50");
        tfIter = new JTextField("100");
        tfAlpha = new JTextField("1.0");
        tfBeta = new JTextField("5.0");
        tfRho = new JTextField("0.5");
        tfQ = new JTextField("100.0");
        tfCities = new JTextField("10");

        cbDataset = new JComboBox<>(new String[]{"Random", "NYC"});

        chkBlock = new JCheckBox("Block Road");
        chkBlock.setEnabled(false);

        sliderSpeed = new JSlider(0, 1000, delayMs);
        sliderSpeed.addChangeListener(e -> delayMs = sliderSpeed.getValue());

        bottomPanel.add(new JLabel("Ants:"));
        bottomPanel.add(tfAnts);

        bottomPanel.add(new JLabel("Iterations:"));
        bottomPanel.add(tfIter);

        bottomPanel.add(new JLabel("Alpha:"));
        bottomPanel.add(tfAlpha);

        bottomPanel.add(new JLabel("Beta:"));
        bottomPanel.add(tfBeta);

        bottomPanel.add(new JLabel("Rho (Evap):"));
        bottomPanel.add(tfRho);

        bottomPanel.add(new JLabel("Q:"));
        bottomPanel.add(tfQ);

        bottomPanel.add(new JLabel("Dataset:"));
        bottomPanel.add(cbDataset);

        bottomPanel.add(new JLabel("Random Cities:"));
        bottomPanel.add(tfCities);

        bottomPanel.add(new JLabel("Delay (ms):"));
        bottomPanel.add(sliderSpeed);

        // yeah, this checkbox sits alone in the grid — leaving it like that on purpose
        bottomPanel.add(chkBlock);

        JButton btnRun = new JButton("Run");
        JButton btnGen = new JButton("Generate Random");

        bottomPanel.add(btnGen);
        bottomPanel.add(btnRun);

        lblStatus = new JLabel("Ready");
        bottomPanel.add(lblStatus);

        add(bottomPanel, BorderLayout.SOUTH);

        // handlers
        btnGen.addActionListener(e -> generateRandomCities());

        cbDataset.addActionListener(e -> {
            boolean isNYC = "NYC".equals(cbDataset.getSelectedItem());
            tfCities.setEnabled(!isNYC);
            chkBlock.setEnabled(isNYC);

            if (isNYC) loadNYC();
            else generateRandomCities();

            canvas.repaint();
        });

        btnRun.addActionListener(e -> {
            if (!isRunning) startSimulation();
        });

        // init
        generateRandomCities();
        pack();
        setLocationRelativeTo(null);
        setVisible(true);
    }

    private void generateRandomCities() {
        int n = 10;
        try {
            n = Integer.parseInt(tfCities.getText().trim());
            if (n < 2) n = 10;
        } catch (NumberFormatException ignored) {
        }

        cities = new City[n];
        Random r = new Random();
        for (int i = 0; i < n; i++) {
            // keep the same spread as before (fits 800x600 nicely)
            cities[i] = new City("C" + i, 50 + r.nextDouble() * 700, 50 + r.nextDouble() * 500);
        }

        calcDistances(false);
        blocked = new boolean[n][n];

        bestTour = null;
        bestLen = Double.MAX_VALUE;
        canvas.repaint();
    }

    private void loadNYC() {
        cities = new City[]{
                new City("Times Sq", 40.7580, -73.9855),
                new City("Central Pk", 40.7812, -73.9665),
                new City("Empire St", 40.7484, -73.9857),
                new City("Brooklyn Br", 40.7061, -73.9969),
                new City("Liberty", 40.6892, -74.0445),
                new City("Wall St", 40.7060, -74.0086),
                new City("Grand Ctrl", 40.7527, -73.9772),
                new City("WTC", 40.7127, -74.0134)
        };

        calcDistances(true);
        blocked = new boolean[cities.length][cities.length];

        bestTour = null;
        bestLen = Double.MAX_VALUE;
    }

    private void calcDistances(boolean haversine) {
        int n = cities.length;
        distMatrix = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) continue;

                if (haversine) {
                    // cities[i].x = lat, cities[i].y = lon
                    double R = 6371.0;
                    double dLat = Math.toRadians(cities[j].x - cities[i].x);
                    double dLon = Math.toRadians(cities[j].y - cities[i].y);

                    double a = Math.pow(Math.sin(dLat / 2), 2)
                            + Math.cos(Math.toRadians(cities[i].x))
                            * Math.cos(Math.toRadians(cities[j].x))
                            * Math.pow(Math.sin(dLon / 2), 2);

                    distMatrix[i][j] = R * 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a));
                } else {
                    distMatrix[i][j] = Math.hypot(cities[i].x - cities[j].x, cities[i].y - cities[j].y);
                }
            }
        }
    }

    private void startSimulation() {
        try {
            int ants = Integer.parseInt(tfAnts.getText().trim());
            int iter = Integer.parseInt(tfIter.getText().trim());

            double alpha = Double.parseDouble(tfAlpha.getText().trim());
            double beta = Double.parseDouble(tfBeta.getText().trim());
            double rho = Double.parseDouble(tfRho.getText().trim());
            double Q = Double.parseDouble(tfQ.getText().trim());

            if (chkBlock.isEnabled()) {
                blocked = new boolean[cities.length][cities.length];
                if (chkBlock.isSelected()) {
                    // single hard-coded "road closure" (demo switch)
                    blocked[0][1] = blocked[1][0] = true;
                }
            }

            isRunning = true;
            lblStatus.setText("Running...");
            toggleControls(false);

            ACO solver = new ACO(ants, alpha, beta, rho, Q, distMatrix, blocked);

            // colors for ants (used just for drawing)
            antColors = new Color[ants];
            for (int i = 0; i < ants; i++) {
                antColors[i] = Color.getHSBColor(i / (float) ants, 0.7f, 0.8f);
            }

            new Thread(() -> {
                for (int i = 0; i < iter && isRunning; i++) {
                    currentTours = solver.step();

                    if (solver.bestLen < bestLen) {
                        bestLen = solver.bestLen;
                        bestTour = solver.bestTour;
                    }

                    final int it = i + 1;
                    SwingUtilities.invokeLater(() -> {
                        canvas.repaint();
                        lblStatus.setText(String.format("Iter %d/%d | Best: %.2f", it, iter, bestLen));
                    });

                    try {
                        Thread.sleep(delayMs);
                    } catch (Exception ignored) {
                    }
                }

                SwingUtilities.invokeLater(() -> {
                    isRunning = false;
                    toggleControls(true);
                    lblStatus.setText("Done. Best: " + String.format("%.2f", bestLen));
                });
            }, "aco-runner").start();

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Bad params");
        }
    }

    private void toggleControls(boolean enable) {
        // leaving the same controls enabled/disabled as the original
        tfAnts.setEnabled(enable);
        tfIter.setEnabled(enable);
        cbDataset.setEnabled(enable);
    }

    // --- Inner Classes ---

    record City(String name, double x, double y) {
    }

    static class ACO {
        int nAnts;
        double alpha, beta, rho, Q;

        double[][] dists;
        double[][] pheromones;
        boolean[][] blocked;
        int nCities;

        int[] bestTour;

        // keep your existing name for the GUI
        double bestLen = Double.MAX_VALUE;

        // --- added for the unit tests ---
        public double bestLength = Double.MAX_VALUE;     // test expects this field
        private int fixedStartCity = -1;                 // test passes "1" as start city

        ACO(int ants, double alpha, double beta, double rho, double Q, double[][] dists, boolean[][] blocked) {
            this.nAnts = ants;
            this.alpha = alpha;
            this.beta = beta;
            this.rho = rho;
            this.Q = Q;

            this.dists = dists;
            this.blocked = blocked;

            this.nCities = dists.length;
            this.pheromones = new double[nCities][nCities];

            for (int i = 0; i < nCities; i++) {
                for (int j = 0; j < nCities; j++) pheromones[i][j] = 1.0;
            }
        }

        // --- added: 8-arg constructor (exact signature your AcoTest calls) ---
        ACO(int ants, int startCity, double alpha, double beta, double rho, double Q, double[][] dists, boolean[][] blocked) {
            this(ants, alpha, beta, rho, Q, dists, blocked);
            if (startCity >= 0 && startCity < nCities) {
                this.fixedStartCity = startCity;
            }
        }

        // --- added: method the test expects ---
        List<int[]> iterateOnce() {
            return step();
        }

        List<int[]> step() {
            int[][] tours = new int[nAnts][nCities];
            double[] lens = new double[nAnts];

            // keep your parallel build (but buildTour must NEVER return partial garbage)
            java.util.stream.IntStream.range(0, nAnts).parallel().forEach(i -> {
                int[] t = buildTour();
                tours[i] = t;
                lens[i] = getLen(t);
            });

            int bestIdx = -1;
            double localMin = Double.POSITIVE_INFINITY;
            for (int i = 0; i < nAnts; i++) {
                if (Double.isFinite(lens[i]) && lens[i] < localMin) {
                    localMin = lens[i];
                    bestIdx = i;
                }
            }

            if (bestIdx != -1 && localMin < bestLen) {
                bestLen = localMin;
                bestLength = bestLen;           // keep alias in sync
                bestTour = tours[bestIdx].clone();
            } else {
                bestLength = bestLen;           // still sync (even if unchanged)
            }

            List<int[]> res = new ArrayList<>(nAnts);
            List<Double> resLens = new ArrayList<>(nAnts);
            for (int i = 0; i < nAnts; i++) {
                res.add(tours[i]);
                resLens.add(lens[i]);
            }

            updatePheromones(res, resLens);
            return res;
        }

        // ----------------- IMPORTANT FIX -----------------
        // Old code returned a partially-filled tour on dead-end (leaving zeros),
        // which creates fake edges like 0-1 and fails your blocked-edge test.
        private int[] buildTour() {
            // retries make it robust if a choice leads to a dead-end
            for (int tries = 0; tries < 300; tries++) {
                int[] t = buildTourOnce();
                if (t != null) return t;
            }
            // with only a single blocked edge (your test), this should never happen
            throw new IllegalStateException("Could not construct a valid tour with given blocked edges.");
        }

        private int[] buildTourOnce() {
            int[] tour = new int[nCities];
            boolean[] visited = new boolean[nCities];

            int curr = (fixedStartCity >= 0) ? fixedStartCity : java.util.concurrent.ThreadLocalRandom.current().nextInt(nCities);
            int start = curr;

            tour[0] = curr;
            visited[curr] = true;

            for (int i = 1; i < nCities; i++) {
                int next = pickNext(curr, visited, start, (nCities - i));
                if (next == -1) return null; // retry
                tour[i] = next;
                visited[next] = true;
                curr = next;
            }

            // must be able to close the loop without using a blocked edge
            if (isBlocked(curr, start)) return null;

            return tour;
        }

        private int pickNext(int curr, boolean[] visited, int start, int remainingToPick) {
            double[] probs = new double[nCities];
            double sum = 0.0;

            for (int i = 0; i < nCities; i++) {
                if (visited[i]) continue;
                if (isBlocked(curr, i)) continue;

                // if this is the LAST pick, ensure we can close i -> start
                if (remainingToPick == 1 && isBlocked(i, start)) continue;

                double tau = Math.pow(pheromones[curr][i], alpha);
                double eta = Math.pow(1.0 / (dists[curr][i] + 1e-6), beta);

                probs[i] = tau * eta;
                sum += probs[i];
            }

            if (sum <= 0.0) return -1;

            double r = java.util.concurrent.ThreadLocalRandom.current().nextDouble() * sum;
            double acc = 0.0;

            for (int i = 0; i < nCities; i++) {
                if (probs[i] <= 0) continue;
                acc += probs[i];
                if (r <= acc) return i;
            }

            for (int i = 0; i < nCities; i++) if (probs[i] > 0) return i;
            return -1;
        }

        private double getLen(int[] tour) {
            double l = 0.0;

            for (int i = 0; i < nCities - 1; i++) {
                if (isBlocked(tour[i], tour[i + 1])) return Double.POSITIVE_INFINITY;
                l += dists[tour[i]][tour[i + 1]];
            }

            if (isBlocked(tour[nCities - 1], tour[0])) return Double.POSITIVE_INFINITY;
            l += dists[tour[nCities - 1]][tour[0]];

            return l;
        }

        private boolean isBlocked(int a, int b) {
            return blocked != null && blocked[a][b];
        }

        private void updatePheromones(List<int[]> tours, List<Double> lens) {
            for (int i = 0; i < nCities; i++) {
                for (int j = 0; j < nCities; j++) pheromones[i][j] *= (1.0 - rho);
            }

            for (int k = 0; k < tours.size(); k++) {
                double len = lens.get(k);
                if (!Double.isFinite(len) || len <= 0) continue;

                double deposit = Q / len;
                int[] t = tours.get(k);

                for (int i = 0; i < nCities - 1; i++) {
                    if (isBlocked(t[i], t[i + 1])) continue;
                    pheromones[t[i]][t[i + 1]] += deposit;
                    pheromones[t[i + 1]][t[i]] += deposit;
                }

                if (!isBlocked(t[nCities - 1], t[0])) {
                    pheromones[t[nCities - 1]][t[0]] += deposit;
                    pheromones[t[0]][t[nCities - 1]] += deposit;
                }
            }
        }
    }

    class DrawingPanel extends JPanel {
        @Override
        protected void paintComponent(Graphics g) {
            super.paintComponent(g);

            Graphics2D g2 = (Graphics2D) g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

            if (cities == null) return;

            int w = getWidth();
            int h = getHeight();

            if (bgMap != null && "NYC".equals(cbDataset.getSelectedItem())) {
                g2.drawImage(bgMap, 0, 0, w, h, null);
            }

            Point[] pts = new Point[cities.length];

            if ("NYC".equals(cbDataset.getSelectedItem())) {
                for (int i = 0; i < cities.length; i++) {
                    double lx = NYC_REL_COORDS[i][0];
                    double ly = NYC_REL_COORDS[i][1];
                    pts[i] = new Point((int) (lx * w), (int) (ly * h));
                }
            } else {
                for (int i = 0; i < cities.length; i++) {
                    pts[i] = new Point((int) cities[i].x, (int) cities[i].y);
                }
            }

            // draw complete graph edges
            g2.setColor(Color.LIGHT_GRAY);
            for (int i = 0; i < cities.length; i++) {
                for (int j = i + 1; j < cities.length; j++) {
                    if (blocked != null && blocked[i][j]) {
                        g2.setColor(Color.RED);
                        g2.setStroke(new BasicStroke(3));
                        g2.drawLine(pts[i].x, pts[i].y, pts[j].x, pts[j].y);
                        g2.setStroke(new BasicStroke(1));
                        g2.setColor(Color.LIGHT_GRAY);
                    } else {
                        g2.drawLine(pts[i].x, pts[i].y, pts[j].x, pts[j].y);
                    }
                }
            }

            // draw ant tours
            if (isRunning && currentTours != null) {
                for (int i = 0; i < currentTours.size(); i++) {
                    g2.setColor(new Color(antColors[i % antColors.length].getRGB() & 0x00FFFFFF | 0x60000000, true));
                    int[] t = currentTours.get(i);

                    for (int j = 0; j < t.length - 1; j++) {
                        g2.drawLine(pts[t[j]].x, pts[t[j]].y, pts[t[j + 1]].x, pts[t[j + 1]].y);
                    }
                    g2.drawLine(
                            pts[t[t.length - 1]].x, pts[t[t.length - 1]].y,
                            pts[t[0]].x, pts[t[0]].y
                    );
                }
            }

            // best tour on top
            if (bestTour != null) {
                g2.setColor(Color.BLUE);
                g2.setStroke(new BasicStroke(2));

                for (int i = 0; i < bestTour.length - 1; i++) {
                    g2.drawLine(
                            pts[bestTour[i]].x, pts[bestTour[i]].y,
                            pts[bestTour[i + 1]].x, pts[bestTour[i + 1]].y
                    );
                }

                g2.drawLine(
                        pts[bestTour[bestTour.length - 1]].x, pts[bestTour[bestTour.length - 1]].y,
                        pts[bestTour[0]].x, pts[bestTour[0]].y
                );
            }

            // cities
            g2.setColor(Color.BLACK);
            for (int i = 0; i < cities.length; i++) {
                g2.fillOval(pts[i].x - 4, pts[i].y - 4, 8, 8);
                g2.drawString(cities[i].name, pts[i].x + 5, pts[i].y - 5);
            }
        }
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(AntColonyTSP::new);
    }
}