#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic;
    int NRow, NCol, Ntic;
    double lambda1, lambda2, lambda3, mu2, mu3;

    int *Freq1, *Freq2, *Freq3;
    double *FreqOut1, *FreqOut2, *FreqOut3;
    // double *FValOut;
    double *Energy, *FVal;
    double minval, val;

    double multi1, multi2;
    int bw2, bw3, bw1 = 75, LB = 25;
    int i,j,k;

    /* Portal to matlab */
    InData = mxGetPr(prhs[0]);
    NRow = mxGetM(prhs[0]);
    NCol = mxGetN(prhs[0]);

    tic = mxGetPr(prhs[1]);
    Ntic = mxGetM(prhs[1]);

    lambda1 = mxGetScalar(prhs[2]);
    lambda2 = mxGetScalar(prhs[3]);
    lambda3 = mxGetScalar(prhs[4]);
    mu2 = mxGetScalar(prhs[5]);
    mu3 = mxGetScalar(prhs[6]);

    multi1 = mxGetScalar(prhs[7]);
    multi2 = mxGetScalar(prhs[8]);

    bw2 = mxGetScalar(prhs[9]);
    bw3 = mxGetScalar(prhs[10]); 
    
    plhs[0] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(NRow, 1, mxDOUBLE_CLASS, mxREAL);

    FreqOut1 = mxGetPr(plhs[0]);
    FreqOut2 = mxGetPr(plhs[1]);
    FreqOut3 = mxGetPr(plhs[2]);
    printf("NRow = %d\n", NRow);
    printf("NCol = %d\n", NCol);

    /* Main operations start here */
    const double eps = 1e-8;
    double sum = 0, tmp;

    printf("BW1, BW2, BW3 = %d, %d, %d\n", bw1, bw2, bw3);
    int *baseLine;

    printf("Start here: %d\n", 3);
    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*bw1*(2*bw2)*(2*bw3));
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq3 = (int *)mxMalloc(sizeof(int)*NRow);
    baseLine = (int *)mxMalloc(sizeof(int)*2*bw1);


    // Build baseLine array
    // multi1 *= 1.01;
    // multi2 *= 1.01;
    int val1, val2;
    double freqVal1, freqVal2;
    for (int j1 = LB; j1 < bw1; ++j1) {
        freqVal1 = tic[j1]*multi1;
        freqVal2 = tic[j1]*multi2;
        for (int ii = Ntic-1; ii >= 1; --ii) {
            if (tic[ii]>=freqVal1 && tic[ii-1]<freqVal1)
                baseLine[j1] = ii;
            if (tic[ii]>=freqVal2 && tic[ii-1]<freqVal2)
                baseLine[bw1+j1] = ii;
        }
    }

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
    // for (i=0; i<NRow; ++i) {
    for (int j1=LB; j1<bw1; ++j1) {
        for (int j2=0; j2<2*bw2; ++j2) {
            for (int j3=0; j3<2*bw3; ++j3) {
                FVal[j1*(2*bw2*2*bw3)+j2*(2*bw3)+j3] = Energy[j1] + \
                Energy[j2+baseLine[j1]-bw2] + Energy[j3+baseLine[bw1+j1]-bw3];
            }
        }
    }

    for (i=1; i<NRow; ++i) {
        for (int j1=LB; j1<bw1; ++j1) {
            for (int j2=0; j2<2*bw2; ++j2) {
                for (int j3=0; j3<2*bw3; ++j3) {
                    minval = 1e16; //FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*2*bw3+j3];
                    for (int k1=LB; k1<bw1; ++k1) {
                        for (int k2=0; k2<2*bw2; ++k2) {
                            for (int k3=0; k3<2*bw3; ++k3) {
                                tmp = FVal[(i-1)*(bw1*2*bw2*2*bw3) + k1*(2*bw2*2*bw3) + k2*2*bw3 + k3] + \
                                lambda1*(k1-j1)*(k1-j1) + \
                                lambda2*((k2+baseLine[k1]-bw2)-(j2+baseLine[j1]-bw2))*((k2+baseLine[k1]-bw2)-(j2+baseLine[j1]-bw2)) + \
                                lambda3*((k3+baseLine[bw1+k1]-bw3)-(j3+baseLine[bw1+j1]-bw3))*((k3+baseLine[bw1+k1]-bw3)-(j3+baseLine[bw1+j1]-bw3)) + \
                                mu2*fabs(multi1*tic[k1]-tic[k2+baseLine[k1]-bw2]) + \
                                mu3*fabs(multi2*tic[k1]-tic[k3+baseLine[bw1+k1]-bw3]);
                                if(tmp < minval)
                                    minval = tmp;
                            }
                        }
                    }
                    FVal[i*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*2*bw3 + j3] = minval + Energy[i*NCol+j1] + \
                    Energy[i*NCol+j2+baseLine[j1]-bw2] + Energy[i*NCol+j3+baseLine[bw1+j1]-bw3];
                }
            }
        }
    }

    minval = FVal[(NRow-1)*(bw1*2*bw2*2*bw3) + LB*2*bw2*2*bw3];
    Freq1[NRow-1] = 0;
    Freq2[NRow-1] = 0;
    Freq3[NRow-1] = 0;
    for (int j1=LB; j1<bw1; ++j1) {
        for (int j2=0; j2<2*bw2; ++j2) {
            for (int j3=0; j3<2*bw3; ++j3) {
                tmp = FVal[(NRow-1)*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*2*bw3 + j3];
                if (tmp < minval) {
                    minval = tmp;
                    Freq1[NRow-1] = j1;
                    Freq2[NRow-1] = j2 + baseLine[j1] - bw2;
                    Freq3[NRow-1] = j3 + baseLine[bw1+j1] - bw3;
                }
            }
        }
    }

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(bw1*2*bw2*2*bw3) + Freq1[i+1]*(2*bw2*2*bw3) + \
        (Freq2[i+1]-baseLine[Freq1[i+1]]+bw2)*2*bw3 + Freq3[i+1]-baseLine[bw1+Freq1[i+1]]+bw3] \
        - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]] - Energy[(i+1)*NCol+Freq3[i+1]];
        minval = 100;
        for (int j1=LB; j1<bw1; ++j1) {
            for (int j2=0; j2<2*bw2; ++j2) {
                for (int j3=0; j3<2*bw3; ++j3) {
                    tmp = FVal[i*(bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*2*bw3 + j3] + \
                    lambda1*(j1-Freq1[i+1])*(j1-Freq1[i+1]) + \
                    lambda2*((j2+baseLine[j1]-bw2)-Freq2[i+1])*((j2+baseLine[j1]-bw2)-Freq2[i+1]) + \
                    lambda3*((j3+baseLine[bw1+j1]-bw3)-Freq3[i+1])*((j3+baseLine[bw1+j1]-bw3)-Freq3[i+1]) + \
                    mu2*fabs(multi1*tic[j1]-tic[j2+baseLine[j1]-bw2]) + \
                    mu3*fabs(multi2*tic[j1]-tic[j3+baseLine[bw1+j1]-bw3]);
                    tmp = fabs(val-tmp);
                    if(tmp < eps) {
                        Freq1[i] = j1;
                        Freq2[i] = j2 + baseLine[j1] - bw2;
                        Freq3[i] = j3 + baseLine[bw1+j1] - bw3;
                        j1 = bw1; j2 = 2*bw2; j3 = 2*bw3; // get out of loops
                    }
                }
            }
        }
    }

    for (i=0; i<NRow; ++i) {
        FreqOut1[i] = (double)Freq1[i]+1.;
        FreqOut2[i] = (double)Freq2[i]+1.;
        FreqOut3[i] = (double)Freq3[i]+1.;
    }

    mxFree(Freq1); mxFree(Freq2); mxFree(Freq3);
    mxFree(baseLine);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}

 
