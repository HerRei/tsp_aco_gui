A small Java Swing app that visualizes an Ant Colony Optimization (ACO) approach to the Traveling Salesman Problem (TSP).  
It can generate random city sets or use a small NYC landmark dataset (optionally drawn over a background map).

Features
- Random dataset**: generate `N` cities and solve/visualize the TSP tour.
- NYC dataset**: fixed landmarks + optional `nyc_map.png` background.
- Live visualization**: shows current ant tours (semi-transparent) and the best tour found so far.
- Adjustable parameters**: ants, iterations, alpha, beta, rho (evaporation), Q, and animation delay.
- Road block demo** (NYC only): toggle a single blocked edge to see how it affects solutions.

Project layout
This is a Maven project:

- src/main/java/tsp_aco_gui/AntColonyTSP.java – GUI + solver classes
- src/main/resources/nyc_map.png – optional NYC background image
- src/test/java/... – unit tests (JUnit 5)

## Requirements
- **Java 21+**
- **Maven 3.9+**

Check:
bash
java -version
mvn -v
