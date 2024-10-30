import java.util.ArrayList;
import java.util.List;

public class GrafoCoexpresion {
    private final List<List<Integer>> grafo;

    public GrafoCoexpresion(double[][] matrizCoexpresion, double minCorrelacion) {
        int n = matrizCoexpresion.length;
        grafo = new ArrayList<>();

        for (int i = 0; i < n; i++) {
            grafo.add(new ArrayList<>());
            for (int j = 0; j < n; j++) {
                if (i != j && matrizCoexpresion[i][j] >= minCorrelacion) {
                    grafo.get(i).add(j);
                }
            }
        }
    }

    public List<List<Integer>> getGrafo() {
        return grafo;
    }
}
