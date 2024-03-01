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

    int multi1, multi2, multi3;
    int bw2, bw3, bw1;
    int i,j,k,l;
    int K = 3;

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

    int *lambda = (int *)mxMalloc(sizeof(int)*K);
    lambda1 = mxGetScalar(prhs[3]);
    lambda2 = mxGetScalar(prhs[4]);
    lambda3 = mxGetScalar(prhs[5]);
    lambda[0] = lambda3; lambda[1] = lambda2; lambda[2] = lambda1;

    int *mu = (int *)mxMalloc(sizeof(int)*K);
    mu1 = mxGetScalar(prhs[6]);
    mu2 = mxGetScalar(prhs[7]);
    mu3 = mxGetScalar(prhs[8]);
    mu[0] = mu3; mu[1] = mu2; mu[2] = mu1;

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
    printf("Multiples are %d, %d, %d\n", multi1, multi2, multi3);

    Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
    FVal = (double *)mxMalloc(sizeof(double)*NRow*(2*bw1)*(2*bw2)*(2*bw3));
    Freq1 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq2 = (int *)mxMalloc(sizeof(int)*NRow);
    Freq3 = (int *)mxMalloc(sizeof(int)*NRow);
    fundM = (int *)mxMalloc(sizeof(int)*(3*NRow));

    BW = (int *)mxMalloc(sizeof(int)*K);
    BW[0] = bw1; BW[1] = bw2; BW[2] = bw3;
    int* CW = (int *)mxMalloc(sizeof(int)*(K-1));
    CW[0] = (2*bw2)*(2*bw3); CW[1] = 2*bw3;
    int M = (2*bw1)*(2*bw2)*(2*bw3);

    int* idxJ = (int *)mxMalloc(sizeof(int)*K);
    int* idxK = (int *)mxMalloc(sizeof(int)*K);

    // Build baseLine array
    double freqVal1, freqVal2, freqVal3;
    for (i=0; i<NRow; ++i) {
        freqVal1 = tic[fund[i]]*multi1;
        freqVal2 = tic[fund[i]]*multi2;
        freqVal3 = tic[fund[i]]*multi3;
        for (k=Ntic-1; k>=1; --k) {
            // search for the true indices for multiple frequency in tfrtic
            if (tic[k]>=freqVal1 && tic[k-1]<freqVal1)
                fundM[i] = k;
            if (tic[k]>=freqVal2 && tic[k-1]<freqVal2)
                fundM[NRow+i] = k;
            if (tic[k]>=freqVal3 && tic[k-1]<freqVal3)
                fundM[2*NRow+i] = k;
        }
    }

    // // Test
    // for (i = 1; i < 20; ++i)
    //     printf("%f %f\n", 200*tic[i], 200*tic[multiIND(tic, Ntic, i, 3)]);

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
    // Initialization
    for (l = 0; l < K; ++l)
        idxJ[l] = 0;
    for (j = 0; j < M; ++j) {
        FVal[idxJ[0]*(2*bw2*2*bw3) + idxJ[1]*(2*bw3) + idxJ[2]] = Energy[idxJ[0]+fundM[0]-bw1] + \
        Energy[idxJ[1]+fundM[NRow]-bw2] + Energy[idxJ[2]+fundM[2*NRow]-bw3];

        idxJ[K-1] += 1;
        l = K-1;
        while ((l >= 1) && (idxJ[l] == 2*BW[l])) {
            idxJ[l] = 0;
            idxJ[l-1] += 1;
            l--;
        }
    }

    int countJ, countK;
    int j1,j2,j3,j4,j5,j6,j7,j8,j9;
    int k1,k2,k3,k4,k5,k6,k7,k8,k9;

    for (i = 1; i < NRow; ++i) {
        // for (l = 0; l < K; ++l) // Initialization
        //     idxJ[l] = 0;
        j1 = j2 = j3 = j4 = j5 = j6 = j7 = j8 = j9 = 0;
        countJ = K-1;
        for (j = 0; j < M; ++j) {
            minval = 1e16; //FVal[i*(bw1*2*bw2*2*bw3)+j1*(2*bw2*2*bw3)+j2*2*bw3+j3];
            // for (l = 0; l < K; ++l) // Initialization
            //     idxK[l] = 0;
            k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = k9 = 0;
            countK = K-1;
            for (k = 0; k < M; ++k) /*while(idxK[0] < 2*BW[0])*/ {
                tmp = FVal[(i-1)*(2*bw1*2*bw2*2*bw3) + k3*(2*bw2*2*bw3) + k2*(2*bw3) + k1] + \
                lambda3 * pow(((k1+fundM[2*NRow+i-1])-(j1+fundM[2*NRow+i])), 2) + \
                lambda2 * pow(((k2+fundM[NRow+i-1])-(j2+fundM[NRow+i])), 2) + \
                lambda1 * pow(((k3+fundM[i-1])-(j3+fundM[i])), 2) + \
                mu3 * fabs(multi3*tic[fund[i-1]] - tic[k1+fundM[2*NRow+i-1]-bw3]) + \
                mu2 * fabs(multi2*tic[fund[i-1]] - tic[k2+fundM[NRow+i-1]-bw2]) + \
                mu1 * fabs(multi1*tic[fund[i-1]] - tic[k3+fundM[i-1]-bw1]);
                if(tmp < minval)
                    minval = tmp;

                // Next frequency indicies
                // idxK[K-1]++;
                // if (idxK[K-1] == 2*BW[K-1]) {
                //     idxK[K-1] = 0;
                //     idxK[K-2]++;
                //     if (idxK[K-2] == 2*BW[K-2]) {
                //         idxK[K-2] = 0;
                //         temp = k;
                //         for (l = 0; l < K-2; ++l) {
                //             idxK[l] = temp/CW[l];
                //             temp = temp%CW[l];
                //         }
                //     }
                // }

                // l = K-1;
                // while (l >= 1) {
                //     if (BW[l]*2-idxK[l])) //Not equal zero
                //         break;
                //     idxK[l] = 0;
                //     idxK[l-1] += 1;
                //     l--;
                // }
                
                k1++;
                if (k1 == 2*BW[countK]) {
                    k1 = 0; k2++; countK--;
                    if (countK > 0) {
                        if (k2 == 2*BW[countK]) {
                            k2 = 0; k3++; countK--;
                            if (countK > 0) {
                                if (k3 == 2*BW[countK]) {
                                    k3 = 0; k4++; countK--;
                                    if (countK > 0) {
                                        if (k4 == 2*BW[countK]) {
                                            k4 = 0; k5++; countK--;
                                            if (countK > 0) {
                                                if (k5 == 2*BW[countK]) {
                                                    k5 = 0; k6++; countK--;
                                                    if (countK > 0) {
                                                        if (k6 == 2*BW[countK]) {
                                                            k6 = 0; k7++; countK--;
                                                            if (countK > 0) {
                                                                if (k7 == 2*BW[countK]) {
                                                                    k7 = 0; k8++; countK--;
                                                                    if (countK > 0) {
                                                                        if (k8 == 2*BW[countK]) {
                                                                            k8 = 0; k9++;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                countK = K-1;
            }
            FVal[i*(2*bw1*2*bw2*2*bw3) + j3*(2*bw2*2*bw3) + j2*(2*bw3) + j1] = minval + Energy[i*NCol+(j3+fundM[i]-bw1)] + \
            Energy[i*NCol+(j2+fundM[NRow+i]-bw2)] + Energy[i*NCol+(j1+fundM[2*NRow+i]-bw3)];

            // Next frequency indices
            // idxJ[K-1] += 1;
            // l = K-1;
            // while ((l >= 1) && (idxJ[l] == 2*BW[l])) {
            //     idxJ[l] = 0;
            //     idxJ[l-1] += 1;
            //     l--;
            // }
            j1++;
            if (j1 == 2*BW[countJ]) {
                j1 = 0; j2++; countJ--;
                if (countJ > 0) {
                    if (j2 == 2*BW[countJ]) {
                        j2 = 0; j3++; countJ--;
                        if (countJ > 0) {
                            if (j3 == 2*BW[countJ]) {
                                j3 = 0; j4++; countJ--;
                                if (countJ > 0) {
                                    if (j4 == 2*BW[countJ]) {
                                        j4 = 0; j5++; countJ--;
                                        if (countJ > 0) {
                                            if (j5 == 2*BW[countJ]) {
                                                j5 = 0; j6++; countJ--;
                                                if (countJ > 0) {
                                                    if (j6 == 2*BW[countJ]) {
                                                        j6 = 0; j7++; countJ--;
                                                        if (countJ > 0) {
                                                            if (j7 == 2*BW[countJ]) {
                                                                j7 = 0; j8++; countJ--;
                                                                if (countJ > 0) {
                                                                    if (j8 == 2*BW[countJ]) {
                                                                        j8 = 0; j9++;
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            countJ = K-1;
        }
    }

    // minval = FVal[(NRow-1)*(2*bw1*2*bw2*2*bw3)];
    // Freq1[NRow-1] = 0 + fundM[NRow-1] - bw1;
    // Freq2[NRow-1] = 0 + fundM[2*NRow-1] - bw2;
    // Freq3[NRow-1] = 0 + fundM[3*NRow-1] - bw3;
    // for (j = 0; j < M; ++j) {
    // // for (int j1=0; j1<2*bw1; ++j1) {
    // //     for (int j2=0; j2<2*bw2; ++j2) {
    // //         for (int j3=0; j3<2*bw3; ++j3) {
    //     tmp = FVal[(NRow-1)*(2*bw1*2*bw2*2*bw3) + freqID[j][0]*(2*bw2*2*bw3) + freqID[j][1]*(2*bw3) + freqID[j][2]];
    //     // printf("j1,j2,j3,FVal = %d, %d, %d, %f\n", j1, j2, j3, tmp);
    //     if (tmp < minval) {
    //         minval = tmp;
    //         Freq1[NRow-1] = freqID[j][0] + fundM[NRow-1] - bw1;
    //         Freq2[NRow-1] = freqID[j][1] + fundM[2*NRow-1] - bw2;
    //         Freq3[NRow-1] = freqID[j][2] + fundM[3*NRow-1] - bw3;
    //     }
    // //         }
    // //     }
    // // }
    // }
    // // printf("check0 : %d, %d, %d\n", (Freq1[NRow-1]-fundM[NRow-1]+bw1), (Freq2[NRow-1]-fundM[2*NRow-1]+bw2), Freq3[NRow-1]-fundM[3*NRow-1]+bw3);

    // for (i = NRow-2; i >= 0; --i) {
    //     val = FVal[(i+1)*(2*bw1*2*bw2*2*bw3) + (Freq1[i+1]-fundM[i+1]+bw1)*(2*bw2*2*bw3) + \
    //     (Freq2[i+1]-fundM[NRow+i+1]+bw2)*(2*bw3) + Freq3[i+1]-fundM[2*NRow+i+1]+bw3] \
    //     - Energy[(i+1)*NCol+Freq1[i+1]] - Energy[(i+1)*NCol+Freq2[i+1]] - Energy[(i+1)*NCol+Freq3[i+1]];
    //     // printf("check : %d, %d, %d\n", (Freq1[i+1]-fundM[i+1]+bw1), (Freq2[i+1]-fundM[NRow+i+1]+bw2), Freq3[i+1]-fundM[2*NRow+i+1]+bw3);

    //     for (j = 0; j < M; ++j) {
    //     // for (int j1=0; j1<2*bw1; ++j1) {
    //     //     for (int j2=0; j2<2*bw2; ++j2) {
    //     //         for (int j3=0; j3<2*bw3; ++j3) {
    //         tmp = FVal[i*(2*bw1*2*bw2*2*bw3) + freqID[j][0]*(2*bw2*2*bw3) + freqID[j][1]*(2*bw3) + freqID[j][2]] + \
    //         lambda1 * pow(((freqID[j][0]+fundM[i]-bw1) - Freq1[i+1]), 2) + \
    //         lambda2 * pow(((freqID[j][1]+fundM[NRow+i]-bw2) - Freq2[i+1]), 2) + \
    //         lambda3 * pow(((freqID[j][2]+fundM[2*NRow+i]-bw3) - Freq3[i+1]), 2) + \
    //         mu1 * fabs(multi1*tic[fund[i]] - tic[freqID[j][0]+fundM[i]-bw1]) + \
    //         mu2 * fabs(multi2*tic[fund[i]] - tic[freqID[j][1]+fundM[NRow+i]-bw2]) + \
    //         mu3 * fabs(multi3*tic[fund[i]] - tic[freqID[j][2]+fundM[2*NRow+i]-bw3]);
    //         tmp = fabs(val-tmp);
    //         if(tmp < eps) {
    //             Freq1[i] = freqID[j][0] + fundM[i] - bw1;
    //             Freq2[i] = freqID[j][1] + fundM[NRow+i] - bw2;
    //             Freq3[i] = freqID[j][2] + fundM[2*NRow+i] - bw3;
    //             j = M;
    //             // j1 = 2*bw1; j2 = 2*bw2; j3 = 2*bw3; // get out of loops
    //         }
    //     //         }
    //     //     }
    //     // }
    //     }
    // }

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

int multiIND(double *tic, int Ntic, int ind, int multi) {
    for (int k = Ntic-1; k >= 1; --k) {
        // search for the true indices for multiple frequency in tfrtic
        if (tic[k]>=multi*tic[ind] && tic[k-1]<multi*tic[ind])
            return k;
    }
    return -1;
}