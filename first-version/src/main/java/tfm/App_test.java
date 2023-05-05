package tfm;

import us.hebi.matlab.mat.ejml.Mat5Ejml;
import tfm.performance_metrics;
/**
 * Hello world!
 *
 */

import us.hebi.matlab.mat.format.Mat5;
import us.hebi.matlab.mat.format.Mat5Reader;
import us.hebi.matlab.mat.format.Mat5Serializable;
import us.hebi.matlab.mat.format.Mat5Type;
import us.hebi.matlab.mat.types.AbstractArray;
import us.hebi.matlab.mat.types.Char;
import us.hebi.matlab.mat.types.MatFile;
import us.hebi.matlab.mat.types.MatlabType;
import us.hebi.matlab.mat.types.Matrix;
import us.hebi.matlab.mat.types.Sink;
import us.hebi.matlab.mat.types.Sinks;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.fixed.CommonOps_DDF2;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.MatrixFeatures_DDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.SingularOps_DDRM;
import org.ejml.dense.row.SingularOps_FDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_FDRM;
import org.ejml.dense.row.NormOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.SingularValueDecomposition;
import org.ejml.interfaces.decomposition.SingularValueDecomposition_F64;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;

import static us.hebi.matlab.mat.format.Mat5.*;
import static us.hebi.matlab.mat.format.Mat5WriteUtil.*;
import static us.hebi.matlab.mat.util.Bytes.*;
import static us.hebi.matlab.mat.util.Preconditions.*;

import java.util.Arrays;
import java.util.Random;
import org.apache.commons.lang3.ArrayUtils;

public class App_test {

    public static double[] getVectorWithoutNulls(double[] R, double[] E, Integer Count) {
        double[] vectorWithoutNulls = new double[R.length + Count];
        System.out.println("R.LEGHT ES: " + R.length);
        if (Count == 0) {
            vectorWithoutNulls = R;
        } else {
            double[] vectorAux = new double[Count];
            for (int i = R.length - 1; i < (R.length - 1) + Count; i++) {
                for (int j = 0; j < Count; j++) {
                    vectorAux[j] = E[Count - 1];
                }
            }
            vectorWithoutNulls = ArrayUtils.addAll(R, vectorAux);
        }
        return vectorWithoutNulls;
    }

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

    public static DMatrixRMaj getFila(DMatrixRMaj matrixMixed, Integer numFila) {
        double[][] fila = new double[1][matrixMixed.getNumCols()];
        for (int i = 0; i < matrixMixed.getNumCols(); i++) {
            fila[0][i] = matrixMixed.get(numFila, i);
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(fila);
        return nuevaMatriz;
    }

    // public static DMatrixRMaj getColumna(DMatrixRMaj matrixNullspace, Integer
    // numColumna) {
    // double[][] fila = new double[matrixNullspace.getNumRows()][1];
    // for (int i = 0; i < matrixNullspace.getNumRows() -1; i++) {
    // fila[i][0] = (double)matrixNullspace.get(numColumna, i);
    // }
    // DMatrixRMaj nuevaMatriz = new DMatrixRMaj(fila);
    // return nuevaMatriz;
    // }
    public static DMatrixRMaj getColumna(DMatrixRMaj matrixNullspace, Integer numColumna) {
        int numRows = matrixNullspace.numRows; // Obtiene el número de filas de la matriz
        double[][] primeraColumna = new double[numRows][1]; // Crea un arreglo para almacenar la primera columna
        for (int i = 0; i < numRows; i++) {
            primeraColumna[i][0] = matrixNullspace.get(i, 0); // Obtiene el elemento de la primera columna y la i-ésima fila
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(primeraColumna);
        return nuevaMatriz;
    }

    public static DMatrixRMaj getColumnas(DMatrixRMaj matrixNullspace, double[] E) {
        double[][] fila = new double[matrixNullspace.getNumRows()][E.length];
        for (int j = 0; j < E.length - 1; j++) {
            for (int i = 0; i < matrixNullspace.getNumRows() - 1; i++) {
                fila[i][0] = matrixNullspace.get(i, (int) E[j]);
            }
        }
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(fila);
        return nuevaMatriz;
    }

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

    public static double[] getThreeFirsts(double[] R) {
        double[] arr = new double[3];
        for (int i = 0; i < 3; i++) {
            arr[i] = R[i];
        }
        return arr;
    }

    /**
     * @param R
     * @param matrixMixed
     * @return
     */
    public static DMatrixRMaj getSubMatrix(double[] R, DMatrixRMaj matrixMixed) {
        double[][] filasSeleccionadas = new double[R.length][matrixMixed.numCols];
        for (int i = 0; i < R.length; i++) {
            for (int j = 0; j < matrixMixed.getNumCols(); j++) {
                double x = matrixMixed.get((((int) R[i])), j);
                filasSeleccionadas[i][j] = matrixMixed.get((((int) R[i])), j);
            }
        }
        // Creamos una nueva matriz con las filas seleccionadas
        DMatrixRMaj nuevaMatriz = new DMatrixRMaj(filasSeleccionadas);
        return nuevaMatriz;
    }

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

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        // Lectura del fichero//

        MatFile mat = Mat5.readFromFile("/home/panda/TFM/TFM/first-version/src/main/java/tfm/DATA5.mat");
        System.out.println(mat);
        System.out.println("Hola, probemos a ver porque MIEDO xd");
        DMatrixRMaj matrixA = Mat5Ejml.convert(mat.getArray("A"), new DMatrixRMaj(0, 0));
        DMatrixRMaj matrixABF = Mat5Ejml.convert(mat.getArray("abf"), new DMatrixRMaj(0, 0));
        ;
        DMatrixRMaj matrixMixed = Mat5Ejml.convert(mat.getArray("mixed"), new DMatrixRMaj(0, 0));
        ;
        double valueC = mat.getMatrix("C").getDouble(0);

        // Set parameteres//
        int intC = (int) valueC;
        Integer Np = (Integer) matrixMixed.numRows;
        Double Nb = (double) matrixMixed.numCols;
        System.out.println("NP = " + Np);
        System.out.println("Nb = " + Nb);
        double[] R = getRandomPermutation(Np);
        double[] R_3 = getThreeFirsts(R);
        R_3[0] = 0;
        R_3[1] = 0;
        R_3[2] = 0;

        System.out.println(Arrays.toString(R_3));
        double[] E = new double[5];
        double[] Eaux = new double[1];
        double[] aux = null;
        Integer ContadorE = 0;
        for (int i = 0; i < intC; i++) {
            DMatrixRMaj v_p2 = null;
            DMatrixRMaj newMatrix = null;
            // aux = ArrayUtils.addAll(R_3, E);
            R_3 = getVectorWithoutNulls(R_3, E, ContadorE);
            if (ContadorE < 3) {
                System.out.println(Arrays.toString(R_3));
                // System.out.println(E.length);
                R_3 = Arrays.copyOfRange(R_3, R_3.length - 3, R_3.length);
                System.out.println("RaUx = " + Arrays.toString(R_3));
                newMatrix = getSubMatrix(R_3, matrixMixed);
                // System.out.println("NewMatrix is:" + newMatrix);
                v_p2 = getv_p2Matrix(matrixMixed, R_3[0]);
                if (ContadorE == 0) {
                    System.out.println(v_p2);
                }
            } else {
                newMatrix = getSubMatrix(E, matrixMixed);
                v_p2 = getv_p2Matrix(matrixMixed, E[0]);
            }
            // System.out.println(v_p2);
            // System.out.println("HERE" + newMatrix.numCols + "" + newMatrix.numRows);
            int rank = MatrixFeatures_DDRM.rank(newMatrix);
            System.out.println("Rank =" + rank);
            SingularValueDecomposition_F64<DMatrixRMaj> svd = DecompositionFactory_DDRM.svd(newMatrix.numRows,
                    newMatrix.numCols, true, true, false);
            SingularValueDecomposition<DMatrixRMaj> svdTest = DecompositionFactory_DDRM.svd(newMatrix.numRows,
                    newMatrix.numCols, true, true, false);
            svd.decompose(newMatrix);
            // Obtener la matriz S de valores singulares
            svdTest.decompose(newMatrix);
            System.out.println(newMatrix.numCols + "," + newMatrix.numRows);
            // calcular epsilon de la máquina (valor mas pequeño de un double o algo asi...)
            // double eps2 = Math.ulp(1.0);
            // System.out.println(eps2);
            DMatrixRMaj nullspace = new DMatrixRMaj(newMatrix.getNumCols(), newMatrix.getNumCols() - (rank));
            System.out.println("nullspaceIndices: " + nullspace.numCols + nullspace.numRows);
            nullspace = SingularOps_DDRM.nullSpace(svd, nullspace, rank);
            // nullspace = SingularOps_DDRM.nullspaceQRP(nullspace, rank);
            // nullspace =
            // System.out.println(nullspace);
            // N = la primera columna de nullspace, es decir, n = Gen_sol(:,1);
            DMatrixRMaj columna = getColumna(nullspace, 0);
            // System.out.println(columna);
            DMatrixRMaj C = new DMatrixRMaj(v_p2.numRows, columna.numCols);
            // *******************************************************************************
            // [V_p_123,E_lab] = max(abs(v_p2*n)/norm(n)); donde n es columna y v_p2 pues
            // v_p2
            // *******************************

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
            System.out.println("El mayor valor obtenido es : " + V_p_123 + "en indice: " + maxIndex);
            if (V_p_123 > 0.00001) {
                E[i] = maxIndex;
                Eaux[0] = maxIndex;
                ContadorE = ContadorE + 1;
                System.out.println("Contador : " + ContadorE);
            }
        }
        DMatrixRMaj mixedTransColumn = new DMatrixRMaj(matrixMixed.numCols, matrixMixed.numRows);
        CommonOps_DDRM.transpose(matrixMixed, mixedTransColumn);
        DMatrixRMaj matrixMixedTrans = mixedTransColumn;
        mixedTransColumn = getColumnas(mixedTransColumn, E);
        // System.out.println(mixedTrans);
        int Nx = 58;
        int Ny = 58;
        Nb = 188.0;

        //
        performance_metrics test = new performance_metrics();
        performance_metrics.performance_metrics_method(matrixA, mixedTransColumn, matrixMixedTrans, matrixABF, Nx, Ny,
                Nb, valueC);

    }
}
