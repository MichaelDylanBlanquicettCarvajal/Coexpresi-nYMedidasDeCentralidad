import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MatrizCoexpresion {
    
    private List<String> identificadores;
    private double[][] matrizCoexpresion;

    public MatrizCoexpresion(String archivo) throws IOException {
        this.identificadores = new ArrayList<>();
        List<double[]> datosExpresion = new ArrayList<>();

        BufferedReader br = new BufferedReader(new FileReader(archivo));
        String linea;

        // Leer la primera línea para los identificadores
        if ((linea = br.readLine()) != null) {
            String[] cabecera = linea.split(",");
            for (int i = 8; i < cabecera.length; i++) {
                identificadores.add(cabecera[i]);
            }
        }

        // Leer los datos de expresión
        List<String[]> valoresFiltrados = new ArrayList<>();
        while ((linea = br.readLine()) != null) {
            String[] valores = linea.split(",");
            valoresFiltrados.add(valores);
        }
        br.close();

        // Inicializar promedios
        int numColumnas = valoresFiltrados.get(0).length - 8;
        double[] sumas = new double[numColumnas];
        int[] conteoValidos = new int[numColumnas];

        // Calcular sumas y conteo de valores válidos para calcular promedios
        for (String[] valores : valoresFiltrados) {
            double[] expresionGen = new double[numColumnas];
            for (int i = 8; i < valores.length; i++) {
                if (!valores[i].equals("NA")) {
                    double valor = Double.parseDouble(valores[i]);
                    expresionGen[i - 8] = valor;
                    sumas[i - 8] += valor;
                    conteoValidos[i - 8]++;
                }
            }
            datosExpresion.add(expresionGen);
        }

        // Calcular promedios
        double[] promedios = new double[numColumnas];
        for (int i = 0; i < numColumnas; i++) {
            promedios[i] = sumas[i] / conteoValidos[i];
        }

        // Reemplazar "NA" con promedios
        for (String[] valores : valoresFiltrados) {
            double[] expresionGen = new double[numColumnas];
            for (int i = 8; i < valores.length; i++) {
                if (valores[i].equals("NA")) {
                    expresionGen[i - 8] = promedios[i - 8];
                } else {
                    expresionGen[i - 8] = Double.parseDouble(valores[i]);
                }
            }
            datosExpresion.add(expresionGen);
        }

        // Convertir la lista en una matriz
        int numFilas = datosExpresion.size();
        double[][] matrizDatos = new double[numFilas][numColumnas];

        for (int i = 0; i < numFilas; i++) {
            matrizDatos[i] = datosExpresion.get(i);
        }

        // Calcular la matriz de coexpresión con Pearson
        matrizCoexpresion = new double[numColumnas][numColumnas];
        for (int i = 0; i < numColumnas; i++) {
            for (int j = i; j < numColumnas; j++) {
                double correlacion = calcularCorrelacionPearson(obtenerColumna(matrizDatos, i), obtenerColumna(matrizDatos, j));
                matrizCoexpresion[i][j] = correlacion;
                matrizCoexpresion[j][i] = correlacion; // Matriz simétrica
            }
        }
    }

    private double[] obtenerColumna(double[][] matriz, int columna) {
        double[] col = new double[matriz.length];
        for (int i = 0; i < matriz.length; i++) {
            col[i] = matriz[i][columna];
        }
        return col;
    }

    private double calcularCorrelacionPearson(double[] x, double[] y) {
        int n = x.length;
        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
        for (int i = 0; i < n; i++) {
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i] * y[i];
            sumX2 += x[i] * x[i];
            sumY2 += y[i] * y[i];
        }
        double numerador = n * sumXY - sumX * sumY;
        double denominador = Math.sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));
        return numerador / denominador;
    }


    public void imprimirMatriz(int x, int y) {
        int numIdentificadores = identificadores.size();
    
        // Encabezado
        System.out.printf("====| %-30s ====%n", String.join(" | ", identificadores.subList(0, x)));
        System.out.printf("----| %-30s ----%n", String.join(" | ", identificadores.subList(numIdentificadores - y, numIdentificadores)));
    
        // Esquina superior izquierda
        imprimirEsquina(0, x, 0, x, identificadores, matrizCoexpresion);
    
        System.out.println("----|");
    
        // Esquina inferior derecha
        imprimirEsquina(numIdentificadores - y, numIdentificadores, numIdentificadores - y, numIdentificadores, identificadores, matrizCoexpresion);
    }
    
    private void imprimirEsquina(int startRow, int endRow, int startCol, int endCol, List<String> identificadores, double[][] matrizCoexpresion) {
        for (int i = startRow; i < endRow; i++) {
            System.out.printf("%-5s |", identificadores.get(i));
            for (int j = startCol; j < endCol; j++) {
                System.out.printf("\t%.2f", matrizCoexpresion[i][j]);
            }
            System.out.println();
        }
    }
    
    
    

    public static void main(String[] args) throws IOException {
            MatrizCoexpresion matriz = new MatrizCoexpresion("data/LiverFemale3600.csv");
            matriz.imprimirMatriz(2,4);

    }

    public List<String> getIdentificadores() {
        return identificadores;
    }

    public double[][] getMatrizCoexpresion() {
        return matrizCoexpresion;
    }

}
