import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class App {
    public static void main(String[] args) throws IOException {
        // Crear una instancia de MatrizCoexpresion y cargar datos
        MatrizCoexpresion matriz = new MatrizCoexpresion("data/LiverFemale3600.csv");

        // Imprimir la matriz de coexpresión
        matriz.imprimirMatriz(2, 4);

        // Generar el grafo de coexpresión
        GrafoCoexpresion grafoCoexpresion = new GrafoCoexpresion(matriz.getMatrizCoexpresion(), 0.8);
        List<List<Integer>> grafo = grafoCoexpresion.getGrafo();

        // Calcular métricas del grafo
        int[] grados = calcularGrado(grafo);
        int[] distribucionGrados = distribucionDeGrados(grados);
        double densidad = calcularDensidad(grafo);
        double[] coeficientesAgrupamiento = coeficienteAgrupamiento(grafo);
        double espectroAgrupamiento = calcularEspectroAgrupamiento(coeficientesAgrupamiento);

        // Generar archivo para Cytoscape
        generarArchivoCytoscape(grafo, matriz.getIdentificadores(), matriz.getMatrizCoexpresion(), "coexpresion_network.txt");

        // Imprimir resultados
        System.out.println("Grados de cada vértice: " + java.util.Arrays.toString(grados));
        System.out.println("Distribución de grados: " + java.util.Arrays.toString(distribucionGrados));
        System.out.println("Densidad total: " + densidad);
        System.out.println("Coeficientes de agrupamiento: " + java.util.Arrays.toString(coeficientesAgrupamiento));
        System.out.println("Espectro de agrupamiento: " + espectroAgrupamiento);
    }

    public static int[] calcularGrado(List<List<Integer>> grafo) {
        int n = grafo.size();
        int[] grados = new int[n];

        for (int i = 0; i < n; i++) {
            grados[i] = grafo.get(i).size();
        }
        return grados;
    }

    public static int[] distribucionDeGrados(int[] grados) {
        int maxGrado = 0;

        for (int grado : grados) {
            if (grado > maxGrado) {
                maxGrado = grado;
            }
        }

        int[] distribucion = new int[maxGrado + 1];

        for (int grado : grados) {
            distribucion[grado]++;
        }
        return distribucion;
    }

    public static double calcularDensidad(List<List<Integer>> grafo) {
        int n = grafo.size();
        int aristas = 0;

        for (List<Integer> adyacentes : grafo) {
            aristas += adyacentes.size();
        }

        return (double) aristas / (n * (n - 1));
    }

    public static double[] coeficienteAgrupamiento(List<List<Integer>> grafo) {
        int n = grafo.size();
        double[] coeficientes = new double[n];

        for (int i = 0; i < n; i++) {
            List<Integer> vecinos = grafo.get(i);
            int k = vecinos.size();
            if (k < 2) {
                coeficientes[i] = 0.0;
                continue;
            }

            int enlacesVecinos = 0;
            for (int j = 0; j < k; j++) {
                for (int l = j + 1; l < k; l++) {
                    if (grafo.get(vecinos.get(j)).contains(vecinos.get(l))) {
                        enlacesVecinos++;
                    }
                }
            }

            coeficientes[i] = (2.0 * enlacesVecinos) / (k * (k - 1));
        }

        return coeficientes;
    }

    public static double calcularEspectroAgrupamiento(double[] coeficientes) {
        double sumaCoeficientes = 0.0;
        for (double coef : coeficientes) {
            sumaCoeficientes += coef;
        }
        return sumaCoeficientes / coeficientes.length;
    }

    public static void generarArchivoCytoscape(List<List<Integer>> grafo, List<String> identificadores, double[][] matrizCoexpresion, String archivo) throws IOException {
        FileWriter writer = new FileWriter(archivo);

        writer.write("Id gen 1\tId gen 2\tCoeficiente de correlación\n");

        for (int i = 0; i < grafo.size(); i++) {
            for (int j : grafo.get(i)) {
                if (i < j) { // Evitar duplicar conexiones en grafo no dirigido
                    writer.write(identificadores.get(i) + "\t" + identificadores.get(j) + "\t" + matrizCoexpresion[i][j] + "\n");
                }
            }
        }

        writer.close();
    }
}
