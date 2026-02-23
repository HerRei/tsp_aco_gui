package tsp_aco_gui;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class AcoTest {

    private static final int NUM_CITIES = 6;
    private static final int NUM_ANTS = 20;

    private double[][] distances;
    private boolean[][] blocked;
    private AntColonyTSP.ACO solver;

    @BeforeEach
    void setUp() {
        distances = new double[NUM_CITIES][NUM_CITIES];
        for (int i = 0; i < NUM_CITIES; i++) {
            for (int j = 0; j < NUM_CITIES; j++) {
                distances[i][j] = (i == j) ? 0.0 : (Math.abs(i - j) + 1) * 10.0;
            }
        }

        blocked = new boolean[NUM_CITIES][NUM_CITIES];
        solver = new AntColonyTSP.ACO(NUM_ANTS, 1, 1.0, 2.0, 0.5, 100.0, distances, blocked);
    }

    @Test
    void iterateOnceReturnsOneTourPerAnt() {
        List<?> tours = solver.iterateOnce();
        assertEquals(NUM_ANTS, tours.size());
    }

    @Test
    void toursVisitEachCityExactlyOnce() {
        List<?> tours = solver.iterateOnce();
        int[] tour = (int[]) tours.get(0);

        assertEquals(NUM_CITIES, tour.length);

        Set<Integer> seen = new HashSet<>();
        for (int city : tour) {
            assertTrue(city >= 0 && city < NUM_CITIES, "City index out of bounds: " + city);
            seen.add(city);
        }
        assertEquals(NUM_CITIES, seen.size(), "Tour should contain each city exactly once");
    }

    @Test
    void blockedEdgeIsNotUsedInTours() {
        blocked[0][1] = true;
        blocked[1][0] = true;
        solver = new AntColonyTSP.ACO(NUM_ANTS, 1, 1.0, 2.0, 0.5, 100.0, distances, blocked);

        List<?> tours = solver.iterateOnce();
        for (Object o : tours) {
            int[] tour = (int[]) o;
            assertFalse(usesEdge(tour, 0, 1), "Tour used blocked edge (0-1)");
        }
    }

    @Test
    void bestLengthGetsUpdated() {
        solver.iterateOnce();
        assertNotNull(solver.bestTour);
        assertTrue(Double.isFinite(solver.bestLength));
        assertTrue(solver.bestLength < Double.MAX_VALUE);
    }

    private static boolean usesEdge(int[] tour, int a, int b) {
        for (int i = 0; i < tour.length - 1; i++) {
            if (isEdge(tour[i], tour[i + 1], a, b)) return true;
        }
        return isEdge(tour[tour.length - 1], tour[0], a, b);
    }

    private static boolean isEdge(int u, int v, int a, int b) {
        return (u == a && v == b) || (u == b && v == a);
    }
}
