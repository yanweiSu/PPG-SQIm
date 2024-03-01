#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *InData, *tic, *fundIN;
    int NRow, NCol, Ntic;
    double lambda1, lambda2, lambda3, mu1, mu2, mu3;

    int *Freq1, *Freq2, *Freq3, *fundM, *fund, *BW;
    double *FreqOut1, *FreqOut2, *FreqOut3;
    // double *FValOut;
    double *Energy, *FVal;
    double minval, val;

    double multi1, multi2, multi3;
    int bw2, bw3, bw1;
    int i,j,k,l;

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
    lambda3 = mxGetScalar(prhs[5]);

    mu1 = mxGetScalar(prhs[6]);
    mu2 = mxGetScalar(prhs[7]);
    mu3 = mxGetScalar(prhs[8]);

    multi1 = mxGetScalar(prhs[9]);
    multi2 = mxGetScalar(prhs[10]);
    multi3 = mxGetScalar(prhs[11]);

    bw1 = mxGetScalar(prhs[12]);
    bw2 = mxGetScalar(prhs[13]);
    bw3 = mxGetScalar(prhs[14]);
    
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
    printf("Multiples are %.1f, %.1f, %.1f\n", multi1, multi2, multi3);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(2*bw1)*(2*bw2)*(2*bw3));
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq3 = (int *)mxMalloc(sizeof(int)*NRow);
    fundM = (int *)mxMalloc(sizeof(int)*(3*NRow));

    int win = 5, start, end;
    double avg;

    // Build baseLine array
    double freqVal1, freqVal2, freqVal3;
    for (i = 0; i<NRow; ++i) {
        avg = 0.0;
        // Calculate moving average
        start = max(0, i-win);
        end = min(NRow-1, i+win);
        avg = tic[fund[start]-1];
        l = start;
        while (l < end) {
            avg += tic[fund[l+1]-1];
            l++;
        }
        avg = (double) (avg/(end-start+1));

        // avg = tic[fund[i]-1];
        freqVal1 = avg*multi1;
        freqVal2 = avg*multi2;
        freqVal3 = avg*multi3;

        for (k = Ntic-1; k >= 1; --k) {
            // search for the true indices for multiple frequency in tfrtic
            if (tic[k]>=freqVal1 && tic[k-1]<freqVal1)
                fundM[i] = k;
            if (tic[k]>=freqVal2 && tic[k-1]<freqVal2)
                fundM[NRow+i] = k;
            if (tic[k]>=freqVal3 && tic[k-1]<freqVal3)
                fundM[2*NRow+i] = k;
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
    // for (i=0; i<NRow; ++i)
    for (int j1=0; j1<2*bw1; ++j1) {
        for (int j2=0; j2<2*bw2; ++j2) {
            for (int j3=0; j3<2*bw3; ++j3) {
                FVal[j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3] = Energy[j1+fundM[0]-bw1] + \
                Energy[j2+fundM[NRow]-bw2] + Energy[j3+fundM[2*NRow]-bw3];
            }
        }
    }

    for (i = 1; i < NRow; ++i) {
        for (int j1=0; j1<2*bw1; ++j1) {
            for (int j2=0; j2<2*bw2; ++j2) {
                for (int j3=0; j3<2*bw3; ++j3) {
                    minval = 1e16; //FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*2*bw3+j3];
                    for (int k1=0; k1<2*bw1; ++k1) {
                        for (int k2=0; k2<2*bw2; ++k2) {
                            for (int k3=0; k3<2*bw3; ++k3) {
                                tmp = FVal[(i-1)*(2*bw1*2*bw2*2*bw3) + k1*(2*bw2*2*bw3) + k2*(2*bw3) + k3] + \
                                lambda1 * pow((k1+fundM[i-1])-(j1+fundM[i]), 2) + \
                                lambda2 * pow((k2+fundM[NRow+i-1])-(j2+fundM[NRow+i]), 2) + \
                                lambda3 * pow((k3+fundM[2*NRow+i-1])-(j3+fundM[2*NRow+i]), 2) + \
                                pow(mu1*fabs(multi1*tic[fund[i-1]-1]-tic[k1+fundM[i-1]-bw1-1]), 2) + \
                                pow(mu2*fabs(multi2*tic[fund[i-1]-1]-tic[k2+fundM[NRow+i-1]-bw2-1]), 2) + \
                                pow(mu3*fabs(multi3*tic[fund[i-1]-1]-tic[k3+fundM[2*NRow+i-1]-bw3-1]), 2);
                                if(tmp < minval)
                                    minval = tmp;
                            }
                        }
                    }
                    FVal[i*(2*bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3] = minval + Energy[i*NCol+(j1+fundM[i]-bw1)] + \
                    Energy[i*NCol+(j2+fundM[NRow+i]-bw2)] + Energy[i*NCol+(j3+fundM[2*NRow+i]-bw3)];
                }
            }
        }
    }

    minval = FVal[(NRow-1)*(2*bw1*2*bw2*2*bw3)];
    Freq1[NRow-1] = 0 + fundM[NRow-1] - bw1;
    Freq2[NRow-1] = 0 + fundM[2*NRow-1] - bw2;
    Freq3[NRow-1] = 0 + fundM[3*NRow-1] - bw3;
    for (int j1=0; j1<2*bw1; ++j1) {
        for (int j2=0; j2<2*bw2; ++j2) {
            for (int j3=0; j3<2*bw3; ++j3) {
                tmp = FVal[(NRow-1)*(2*bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3];
                // printf("j1,j2,j3,FVal = %d, %d, %d, %f\n", j1, j2, j3, tmp);
                if (tmp < minval) {
                    minval = tmp;
                    Freq1[NRow-1] = j1 + fundM[NRow-1] - bw1;
                    Freq2[NRow-1] = j2 + fundM[2*NRow-1] - bw2;
                    Freq3[NRow-1] = j3 + fundM[3*NRow-1] - bw3;
                }
            }
        }
    }

    // printf("check0 : %d, %d, %d\n", (Freq1[NRow-1]-fundM[NRow-1]+bw1), (Freq2[NRow-1]-fundM[2*NRow-1]+bw2), Freq3[NRow-1]-fundM[3*NRow-1]+bw3);

    for (i = NRow-2; i >= 0; --i) {
        val = FVal[(i+1)*(2*bw1*2*bw2*2*bw3) + (Freq1[i+1]-fundM[i+1]+bw1)*(2*bw2*2*bw3) + \
        (Freq2[i+1]-fundM[NRow+i+1]+bw2)*(2*bw3) + Freq3[i+1]-fundM[2*NRow+i+1]+bw3] \
        - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]] - Energy[(i+1)*NCol+Freq3[i+1]];
        // printf("check : %d, %d, %d\n", (Freq1[i+1]-fundM[i+1]+bw1), (Freq2[i+1]-fundM[NRow+i+1]+bw2), Freq3[i+1]-fundM[2*NRow+i+1]+bw3);

        for (int j1=0; j1<2*bw1; ++j1) {
            for (int j2=0; j2<2*bw2; ++j2) {
                for (int j3=0; j3<2*bw3; ++j3) {
                    tmp = FVal[i*(2*bw1*2*bw2*2*bw3) + j1*(2*bw2*2*bw3) + j2*(2*bw3) + j3] + \
                    lambda1 * pow((j1+fundM[i]-bw1) - Freq1[i+1], 2) + \
                    lambda2 * pow((j2+fundM[NRow+i]-bw2) - Freq2[i+1], 2) + \
                    lambda3 * pow((j3+fundM[2*NRow+i]-bw3) - Freq3[i+1], 2) + \
                    pow(mu1*fabs(multi1*tic[fund[i]-1]-tic[j1+fundM[i]-bw1-1]), 2) + \
                    pow(mu2*fabs(multi2*tic[fund[i]-1]-tic[j2+fundM[NRow+i]-bw2-1]), 2) + \
                    pow(mu3*fabs(multi3*tic[fund[i]-1]-tic[j3+fundM[2*NRow+i]-bw3-1]), 2);
                    tmp = fabs(val-tmp);
                    if(tmp < eps) {
                        Freq1[i] = j1 + fundM[i] - bw1;
                        Freq2[i] = j2 + fundM[NRow+i] - bw2;
                        Freq3[i] = j3 + fundM[2*NRow+i] - bw3;
                        j1 = 2*bw1; j2 = 2*bw2; j3 = 2*bw3; // get out of loops
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
    mxFree(fundM);
    mxFree(Energy);
    mxFree(FVal); 
    return;
}

void findIdx(int **array, int **band, int len, int value) {
    int temp = value;
    *array = mxMalloc(sizeof(int)*len);
    if (*array == NULL)
        return;
    // Get value's 'frequency location'
    for (int k=len-1; k>=0; --k) {
        (*array)[k] = temp%(2*(*band)[k]);
        temp = temp/(2*(*band)[k]);
    }
    // Now j1 = freqID[1], ..., jK = freqID[K].
}

int min(int a, int b) {
    return (a>b)?b:a;
}

int max(int a, int b) {
    return (a>b)?a:b;
}