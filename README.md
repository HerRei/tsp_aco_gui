This was a small project for a final assignment, at Unibasel. It was about Meta-Heuristic Algorithms. I built this to have a tool to experiment and show my findings. Most of the online tools i found for ACO were very performance constraint so this was the main point of concern. 

Features
- Random dataset**: generate *N* cities and solve/visualize the TSP tour.
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
