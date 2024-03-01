#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic, *fundIN;
    int NRow, NCol, Ntic;
    double lambda1, lambda2, mu;

    int *Freq1, *Freq2, *fund;
    double *FreqOut1, *FreqOut2;
    double *Energy, *FVal;
    double minval, val;

    int multi1 = 1, multi2 = 2;
    int bw1, bw2;
    int i,j,k,j1,j2;

    /* Portal to matlab */
    InData = mxGetPr(prhs[0]);
    NRow = mxGetM(prhs[0]);
    NCol = mxGetN(prhs[0]);

    tic = mxGetPr(prhs[1]);
    Ntic = mxGetM(prhs[1]);

    fundIN = mxGetPr(prhs[2]);
    fund = (int *)mxMalloc(sizeof(int)*(NRow));
    for (i=0; i<NRow; ++i)
        // Transform to an integer array
        fund[i] = (int) fundIN[i];

    lambda1 = mxGetScalar(prhs[3]);
    lambda2 = mxGetScalar(prhs[4]);

    mu = mxGetScalar(prhs[5]);
    // multi1 = mxGetScalar(prhs[9]);
    // multi2 = mxGetScalar(prhs[10]);
    // multi3 = mxGetScalar(prhs[11]);
    bw1 = mxGetScalar(prhs[6]);
    bw2 = mxGetScalar(prhs[7]);
    
    plhs[0] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    // plhs[2] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut1 = mxGetPr(plhs[0]);
    FreqOut2 = mxGetPr(plhs[1]);
    printf("NRow = %d\n", NRow);
    printf("NCol = %d\n", NCol);

    /* Main operations start here */
    const double eps = 1e-8;
    double sum = 0, tmp;

    printf("BW1, BW2 = %d, %d\n", bw1, bw2);
    printf("Multiples are %d, %d\n", multi1, multi2);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(2*bw1)*(2*bw2));
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);

    int win = 1, start, end;
    double avg;

    // TFR matrix
    for (i=0; i<NRow; ++i) {
        for (j=0; j<NCol; ++j) {
            Energy[i*NCol+j] = InData[i+j*NRow];
            sum += Energy[i*NCol+j];
        }
    }

    for (i=0;i<NRow;++i) {
        for (j=0;j<NCol;++j)
            Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum+eps);
    }

    // Pre-allocation
    // for (i=0; i<NRow; ++i)
    int idxJ1, idxJ2;
    for (j1 = 0; j1 < 2*bw1; ++j1) {
        idxJ1 = j1+fund[0]-bw1;
        idxJ2 = multiIND(tic, Ntic, idxJ1, 2);
        for (j2 = 0; j2 < 2*bw2; ++j2)
            FVal[j1*(2*bw2) + j2] = Energy[idxJ1] + Energy[j2+idxJ2-bw2];
    }

    int idxK1, idxK2, k1, k2;
    for (i = 1; i < NRow; ++i) {
        for (j1 = 0; j1 < 2*bw1; ++j1) {
            idxJ1 = j1+fund[i]-bw1;
            idxJ2 = multiIND(tic, Ntic, idxJ1, 2);
            for (j2 = 0; j2 < 2*bw2; ++j2) {
                minval = 1e16; //FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*2*bw3+j3];
                for (k1 = 0; k1 < 2*bw1; ++k1) {
                    idxK1 = k1+fund[i-1]-bw1;
                    idxK2 = multiIND(tic, Ntic, idxK1, 2);
                    for (k2 = 0; k2 < 2*bw2; ++k2) {
                        tmp = FVal[(i-1)*(2*bw1*2*bw2) + k1*(2*bw2) + k2] + \
                        lambda1 * pow(idxK1-idxJ1, 2) + \
                        lambda2 * pow((k2+idxK2-bw2)-(j2+idxJ2-bw2), 2) + \
                        mu * pow(multi2*tic[idxK1-1] - tic[k2+idxK2-bw2-1], 2);
                        if(tmp < minval)
                            minval = tmp;
                    }
                }
                FVal[i*(2*bw1*2*bw2) + j1*(2*bw2) + j2] = minval + Energy[i*NCol+idxJ1] + Energy[i*NCol+(j2+idxJ2-bw2)];
            }
        }
    }

    minval = FVal[(NRow-1)*(2*bw1*2*bw2)];
    Freq1[NRow-1] = 0 + fund[NRow-1] - bw1;
    Freq2[NRow-1] = 0 + multiIND(tic, Ntic, Freq1[NRow-1], 2) - bw2;
    for (j1 = 0; j1 < 2*bw1; ++j1) {
        idxJ1 = j1 + fund[NRow-1] - bw1;
        idxJ2 = multiIND(tic, Ntic, idxJ1, 2);
        for (j2 = 0; j2 < 2*bw2; ++j2) {
            tmp = FVal[(NRow-1)*(2*bw1*2*bw2) + j1*(2*bw2) + j2];
            if (tmp < minval) {
                minval = tmp;
                Freq1[NRow-1] = idxJ1;
                Freq2[NRow-1] = j2 + idxJ2 - bw2;
            }
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(2*bw1*2*bw2) + (Freq1[i+1]-fund[i+1]+bw1)*(2*bw2) + (Freq2[i+1]-multiIND(tic,Ntic,Freq1[i+1],2)+bw2)] \
        - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]];

        for (j1 = 0; j1 < 2*bw1; ++j1) {
            idxJ1 = j1 + fund[i] - bw1;
            idxJ2 = multiIND(tic, Ntic, idxJ1, 2);
            for (j2 = 0; j2 < 2*bw2; ++j2) {
                tmp = FVal[i*(2*bw1*2*bw2) + j1*(2*bw2) + j2] + \
                lambda1 * pow(idxJ1-Freq1[i+1], 2) + \
                lambda2 * pow((j2+idxJ2-bw2) - Freq2[i+1], 2) + \
                mu * pow(multi2*tic[idxJ1-1] - tic[j2+idxJ2-bw2-1], 2);
                tmp = fabs(val-tmp);
                if(tmp < eps) {
                    Freq1[i] = idxJ1;
                    Freq2[i] = j2 + idxJ2 - bw2;
                    j1 = 2*bw1; j2 = 2*bw2; // get out of loops
                }
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut1[i] = (double)Freq1[i]+1.;
        FreqOut2[i] = (double)Freq2[i]+1.;
    }

    mxFree(Freq1); mxFree(Freq2);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}

int min(int a, int b) {
    return (a>b)?b:a;
}

int max(int a, int b) {
    return (a>b)?a:b;
}

int multiIND(double *tic, int Ntic, int ind, int multi) {
    for (int k = Ntic-1; k >= 1; --k) {
        // search for the true indices for multiple frequency in tfrtic
        if (tic[k]>=multi*tic[ind] && tic[k-1]<multi*tic[ind])
            return k;
    }
    return -1;
}