package tfm;

import us.hebi.matlab.mat.ejml.Mat5Ejml;

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
import org.ejml.dense.row.decomposition.qr.QRDecompositionHouseholder_DDRM;
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

public class performance_metrics {

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

    public static void performance_metrics_method(DMatrixRMaj matrixA, DMatrixRMaj mixedTransColumn,
            DMatrixRMaj matrixMixed, DMatrixRMaj matrixABF, int Nx, int Ny, double Nb, double valueC) {
        DMatrixRMaj FirstMatrix = mixedTransColumn;
        DMatrixRMaj SecondMatrix = new DMatrixRMaj(1, mixedTransColumn.getNumCols());
        DMatrixRMaj AA = new DMatrixRMaj(mixedTransColumn.getNumRows(), mixedTransColumn.getNumCols());
        // function [E_aad,E_aid,E_sad,E_sid,E_rmse,Aest] =
        // performance_metrics(A,Aest,mixed,abf,M,N,D,c)
        // %%
        // warning off;
        // AA = [1e-5*Aest;ones(1,length(Aest(1,:)))];
        CommonOps_DDRM.fill(SecondMatrix, 1.0);
        CommonOps_DDRM.scale(1e-5, FirstMatrix);
        AA = CommonOps_DDRM.concatRowsMulti(FirstMatrix, SecondMatrix);
        // sest = zeros(length(Aest(1,:)),M*N);
        DMatrixRMaj sest = new DMatrixRMaj(mixedTransColumn.getNumCols(), Nx * Ny);
        CommonOps_DDRM.fill(sest, 0.0);
        // for j=1:M*N
        DMatrixRMaj r = new DMatrixRMaj(matrixMixed.getNumRows(), 1);
        for (int j = 0; j < (Nx * Ny) - 1; j++) {
            // r = [1e-5*mixed(:,j); 1];
            DMatrixRMaj rAux = null;
            rAux = getColumna(matrixMixed, j);
            DMatrixRMaj rAux2 = new DMatrixRMaj(1, 1);
            CommonOps_DDRM.fill(rAux2, 1.0);
            r = CommonOps_DDRM.concatRowsMulti(rAux, rAux2);
            DMatrixRMaj rP = null;
            // sest(:,j) = lsqnonneg(AA,r);
            
            // SimpleMatrix AA_sm = SimpleMatrix.wrap(AA); 
            // SimpleMatrix r_sm = SimpleMatrix.wrap(r);
            // SimpleMatrix sest_sm = AA_sm.solve
            // CommonOps_DDRM.extractColumn(sest_sm.getDDRM(), 0, sest, 0); 
                                               

        }
        // CRD = corrcoef([A Aest])
        // DD = abs(CRD(c+1:2*c,1:c))
        // perm_mtx = zeros(c,c);
        // aux=zeros(c,1);
        // for i=1:c
        // [ld cd]=find(max(DD(:))==DD);
        // ld=ld(1);cd=cd(1); 
        // perm_mtx(ld,cd)=1;
        // DD(:,cd)=aux; DD(ld,:)=aux';
        // end
        // Aest = Aest*perm_mtx;
        // sest = sest'*perm_mtx;
        // Sest = reshape(sest,[M,N,c]);
        // sest = sest';
        // E_rmse = sqrt(sum(sum(((abf-sest).*(abf-sest)).^2))/(M*N*c))

        // % the angle between abundances
        // nabf = diag(abf*abf');
        // nsest = diag(sest*sest');
        // ang_beta = 180/pi*acos( diag(abf*sest')./sqrt(nabf.*nsest));
        // E_aad = mean(ang_beta.^2)^.5

        // % cross entropy between abundance 丰度之间的交叉熵
        // E_entropy = sum(abf.*log((abf+1e-9)./(sest+1e-9))) +
        // sum(sest.*log((sest+1e-9)./(abf+1e-9)));
        // E_aid = mean(E_entropy.^2)^.5

        // % the angle between material signatures
        // nA = diag(A'*A);
        // nAest = diag(Aest'*Aest);
        // ang_theta = 180/pi*acos( diag(A'*Aest)./sqrt(nA.*nAest) );
        // E_sad = mean(ang_theta.^2)^.5

        // % the spectral information divergence
        // pA = A./(repmat(sum(A),[length(A(:,1)) 1]));
        // qA = Aest./(repmat(sum(Aest),[length(A(:,1)) 1]));
        // qA = abs(qA);
        // SID = sum(pA.*log((pA+1e-9)./(qA+1e-9))) +
        // sum(qA.*log((qA+1e-9)./(pA+1e-9)));
        // E_sid = mean(SID.^2)^.5

    }

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

    }
}
