package tfm;

import us.hebi.matlab.mat.ejml.Mat5Ejml;
import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.types.MatFile;
import java.io.FileWriter;
import java.io.IOException;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.SingularOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import java.util.Arrays;
import java.util.Random;
import org.apache.commons.lang3.ArrayUtils;

public class App {

    /* @param double [] R - a vector
     * @param double [] E with nulls values
     * @param Integer Count - count of how many non null values has E vector
     * @return double [] vector with R + E (without nulls values)
     */
    public static double[] getVectorWithoutNulls(double[] R, double[] E, Integer Count) {
        double[] vectorWithoutNulls = new double[R.length + Count];
        if (Count == 0) {
            vectorWithoutNulls = R;
        } else {
            double[] vectorAux = new double[Count];
            for (int i = R.length - 1; i < R.length + Count; i++) {
                for (int j = 0; j < Count; j++) {
                    vectorAux[j] = E[Count - 1];
                }
            }
            vectorWithoutNulls = ArrayUtils.addAll(R, vectorAux);
        }
        return vectorWithoutNulls;
    }
/*
 * @param DMatrixRMaj matrixMixed 
 * @param double numFila
 * @return DmatrixRMaj v_p2 = mixed-mixed(R(1),:) , return the matrix v_p2 with the indicated row in numFila with 0.0 values
 */
    public static DMatrixRMaj getv_p2Matrix(DMatrixRMaj matrixMixed, double numFila) {
        double[][] fila = new double[matrixMixed.getNumRows()][matrixMixed.getNumCols()];
        for (int i = 0; i < matrixMixed.getNumRows(); i++) {
            for (int j = 0; j < matrixMixed.getNumCols(); j++) {
                fila[i][j] = matrixMixed.get(i, j) - matrixMixed.get((int) numFila, j);
            }
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(fila);
        return nuevaMatriz;
    }



    /*
 * @param DMatrixRMaj matrixNullSpace 
 * @param Integer numColumna
 * @return DmatrixRMaj with only the column of matrixNullSpace indicated by numColumna
     */
    public static DMatrixRMaj getColumna(DMatrixRMaj matrixNullspace, Integer numColumna) {
        int numRows = matrixNullspace.numRows; // Obtiene el número de filas de la matriz
        double[][] primeraColumna = new double[numRows][1]; // Crea un arreglo para almacenar la primera columna
        for (int i = 0; i < numRows; i++) {
            primeraColumna[i][0] = matrixNullspace.get(i, 0); // Obtiene el elemento de la primera columna y la i-ésima
                                                              // fila
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(primeraColumna);
        return nuevaMatriz;
    }

        /*
 * @param DMatrixRMaj matrixNullSpace 
 * @param double [] E
 * @return DmatrixRMaj with only the columns of matrixNullSpace indicated in the E vector
     */
    public static DMatrixRMaj getColumnas(DMatrixRMaj matrixNullspace, double[] E) {
        double[][] fila = new double[matrixNullspace.getNumRows()][E.length];
        for (int j = 0; j < E.length; j++) {
            for (int i = 0; i < matrixNullspace.getNumRows(); i++) {
                fila[i][j] = matrixNullspace.get(i, (int) E[j]);
            }
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(fila);
        return nuevaMatriz;
    }

    /*
     * @param Integer Np
     * @return double[] vector of Np size with random values between 0 and Np
     */
    public static double[] getRandomPermutation(Integer Np) {
        double[] arr = new double[Np];
        for (int i = 0; i < arr.length; i++) {
            arr[i] = i + 1;
        }

        Random rnd = new Random();
        for (int i = arr.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            double temp = arr[index];
            arr[index] = arr[i];
            arr[i] = temp;
        }
        return arr;
    }
    /*
     * @param Double [] R vector
     * @return double[] vector with the three first values of R
     */
    public static double[] getThreeFirsts(double[] R) {
        double[] arr = new double[3];
        for (int i = 0; i < 3; i++) {
            arr[i] = R[i];
        }
        return arr;
    }

    /**
     * @param double R
     * @param DMatrixRMaj matrixMixed
     * @return DMatrixRMaj with the R values rows of matrixMixed
     */
    public static DMatrixRMaj getSubMatrix(double[] R, DMatrixRMaj matrixMixed) {
        double[][] filasSeleccionadas = new double[R.length][matrixMixed.numCols];
        for (int i = 0; i < R.length; i++) {
            for (int j = 0; j < matrixMixed.getNumCols(); j++) {
                filasSeleccionadas[i][j] = matrixMixed.get((((int) R[i])), j);
            }
        }
        // Creamos una nueva matriz con las filas seleccionadas
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(filasSeleccionadas);
        return nuevaMatriz;
    }
    /**
     * @param DMatrixRMaj a
     * 
     * @return double with the Index when the max value of the DmatrixRMaj is. 
     */
    public static double elementMax_Index(DMatrixRMaj a) {
        final int size = a.getNumElements();
        double index = 0;
        double max = a.get(0);
        for (int i = 1; i < size; i++) {
            double val = a.get(i);
            if (val >= max) {
                max = val;
                index = i;
            }
        }

        return index;
    }


    public static void main(String[] args) throws IOException {
        // Lectura del fichero//
        MatFile mat = Mat5.readFromFile("/home/panda/TFM/TFM/first-version/src/main/java/tfm/DATA5.mat");
        System.out.println(mat);
        DMatrixRMaj matrixMixed = Mat5Ejml.convert(mat.getArray("mixed"), new DMatrixRMaj(0, 0));
        double valueC = mat.getMatrix("C").getDouble(0);

        // Set parameteres//
        int intC = (int) valueC;
        Integer Np = (Integer) matrixMixed.numRows;
        Double Nb = (double) matrixMixed.numCols;
        System.out.println("NP = " + Np);
        System.out.println("Nb = " + Nb);

        // Creación de la permutación
        // R = R(1:3);
        double[] R = getRandomPermutation(Np);
        double[] R_3 = getThreeFirsts(R);
        System.out.println(Arrays.toString(R_3));

        // E = [];
        double[] E = new double[5];
        Integer ContadorE = 0;
        // for i = 1:c
        for (int i = 0; i < intC; i++) {
            DMatrixRMaj v_p2 = null;
            DMatrixRMaj newMatrix = null;
            // R = [R,E];
            R_3 = getVectorWithoutNulls(R_3, E, ContadorE);
            // if length(E)<3 %||length(E)==3
            if (ContadorE < 3) {
                System.out.println(Arrays.toString(R_3));
                // R = R(end-2:end);
                R_3 = Arrays.copyOfRange(R_3, R_3.length - 3, R_3.length);
                System.out.println("RaUx = " + Arrays.toString(R_3));

                // newM = mixed(R,:);
                newMatrix = getSubMatrix(R_3, matrixMixed);

                // v_p2 = mixed-mixed(R(1),:);
                v_p2 = getv_p2Matrix(matrixMixed, R_3[0]);
            } else {

                // newM = mixed(E,:);
                // v_p2 = mixed-mixed(E(1),:);
                newMatrix = getSubMatrix(E, matrixMixed);
                v_p2 = getv_p2Matrix(matrixMixed, E[0]);
            }

            // r =rank(newM);
            int rank = MatrixFeatures_DDRM.rank(newMatrix);
            System.out.println("Rank =" + rank);

            // Gen_sol = null(newM,r);
            SingularValueDecomposition_F64<DMatrixRMaj> svd = DecompositionFactory_DDRM.svd(newMatrix.numRows,
                    newMatrix.numCols, true, true, false);

            svd.decompose(newMatrix);
            // Obtener la matriz S de valores singulares
            System.out.println(newMatrix.numCols + "," + newMatrix.numRows);
            DMatrixRMaj nullspace = new DMatrixRMaj(newMatrix.getNumCols(), newMatrix.getNumCols() - (rank));
            System.out.println("nullspaceIndices: " + nullspace.numCols + nullspace.numRows);
            nullspace = SingularOps_DDRM.nullSpace(svd, nullspace, rank);
            DMatrixRMaj columna = getColumna(nullspace, 0);
            // *******************************************************************************
            // [V_p_123,E_lab] = max(abs(v_p2*n)/norm(n)); donde n es columna y v_p2 pues
            // v_p2
            // *******************************
            DMatrixRMaj C = new DMatrixRMaj(v_p2.numRows, columna.numCols);
            // v_p2*n
            CommonOps_DDRM.mult(v_p2, columna, C);
            // System.out.println(C.numCols + "," + C.numRows);
            // System.out.println(C);
            // absoluto de C
            CommonOps_DDRM.abs(C);
            // Double test = NormOps_DDRM.normF(columna);
            // normalización de la matriz columna
            Double norm = NormOps_DDRM.normP2(columna);
            System.out.println(norm);
            // abs(v_p2*n)/norm(n)
            DMatrixRMaj Caux = new DMatrixRMaj(v_p2.numRows, columna.numCols);
            CommonOps_DDRM.divide(C, norm, Caux);
            // System.out.println(Caux);
            double V_p_123 = CommonOps_DDRM.elementMax(Caux);
            // max (abs(v_p2*n)/norm(n))
            double maxIndex = elementMax_Index(Caux);
            System.out.println("Columnas: " + Caux.getNumCols() + " /Filas: " + Caux.getNumRows());
            System.out.println("El mayor valor obtenido es : " + V_p_123 + " en indice: " + maxIndex);
            if (V_p_123 > 0.00001) {
                // E(:,i) = E_lab;
                E[i] = maxIndex;
                ContadorE = ContadorE + 1;
                System.out.println("Contador : " + ContadorE);
            }
        }

        // mixed = mixed';

        DMatrixRMaj mixedTransColumn = new DMatrixRMaj(matrixMixed.numCols, matrixMixed.numRows);
        CommonOps_DDRM.transpose(matrixMixed, mixedTransColumn);
        //System.out.println(mixedTransColumn);
        mixedTransColumn = getColumnas(mixedTransColumn, E);
        // Nx = 58;Ny = 58;Nb = 188;
        FileWriter escribir = new FileWriter("endmembers.txt", false);

        for (int j = 0; j < mixedTransColumn.getNumCols(); j++) {
            for (int i = 0; i < mixedTransColumn.getNumRows(); i++) {
                if (i < 187) {
                    escribir.write(mixedTransColumn.get(i, j) + ",");

                } else {
                    escribir.write(mixedTransColumn.get(i, j) + "");

                }
            }
            escribir.write("\n");
        }

        escribir.close();
        System.out.println("");
    }
}
