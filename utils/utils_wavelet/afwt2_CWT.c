/* afwt2_CWT.c
 * Usage: u = fwt2_CWT(x, h0_r, h1_r, h0_c, h1_c, J, symr, symc);
 * Adjoint of 2D Complex wavelet transform.
 * Written by: M. Salman Asif, Georgia Tech December 2007
 *
 * ------------------------------
 * Copyright (c) 2011 Salman Asif
 * ------------------------------
 */

#include <stdlib.h>
#include "mex.h"

/* these are in afwt_level.c */
void adj_wavelet_one_level(double*, double*, int, double*, double*, int, int);
int checkPowerTwo(int, int);

void adj_wavelet_multi_level2D_CWT(double *u, double *w, int nR, int nC, double *h0_r, double *h1_r, double *h0_c, double *h1_c, int l, int J, int symr, int symc)
{
    int j, n, nRj, nCj, m1, m2;
    double *tmp, *tmpo;
    
    n = (nR > nC) ? nR: nC;    
    tmp = (double*) mxCalloc(n, sizeof(double));
    tmpo = (double*) mxCalloc(n, sizeof(double));
    
    nRj = nR;
    nCj = nC;
    for (m1=0; m1 < nR; m1++){
        for (m2=0; m2<nC; m2++){
            u[m2*nR+m1] = w[m2*nR+m1];
        }
    }
    
    for (j=J-1; j >= 0; j--) {
        nRj = nR>>j;
        nCj = nC>>j;
        
        /* rows */
        for (m1=0; m1 < nRj; m1++) {
            for (m2=0; m2 < nCj; m2++)
                tmp[m2] = u[nR*m2+m1];
            adj_wavelet_one_level(tmpo, tmp, nCj, h0_r, h1_r, l, symr);
            for (m2=0; m2 < nCj; m2++)
                u[nR*m2+m1] = tmpo[m2];
        }
        
        /* columns */
        for (m2=0; m2 < nCj; m2++) {
            for (m1=0; m1 < nRj; m1++)
                tmp[m1] = u[nR*m2+m1];
            adj_wavelet_one_level(tmpo, tmp, nRj, h0_c, h1_c, l, symc);
            for (m1=0; m1 < nRj; m1++)
                u[nR*m2+m1] = tmpo[m1];
        }
    }
    
    mxFree(tmp);
    mxFree(tmpo);
}


/* The gateway routine. */
/* u = fwt2_CWT(x, h0_r, h1_r, h0_c, h1_c, J, symr, symc); */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *u, *w, *h0_r, *h1_r, *h0_c, *h1_c;
    int nR, nC, l, l1, symr, symc, J, JmaxR, JmaxC;
    
  /* Check for the proper number of arguments. */
    if (nrhs != 8) {
        mexErrMsgTxt("Exactly five inputs required");
    }
    if (nlhs > 1) {
        mexErrMsgTxt("Too many output arguments");
    }
    
    nR = mxGetM(prhs[0]);
    nC = mxGetN(prhs[0]);
    /*
     * if (n != n1) {
        mexErrMsgTxt("Input w must be a square image.");
    }
     */
    if (mxIsComplex(prhs[0])) {
        mexErrMsgTxt("Input w must be real");
    }
    J = mxGetScalar(prhs[5]);
    JmaxR = checkPowerTwo(nR, J);
    JmaxC = checkPowerTwo(nC, J);
    if ((J < 0) || (J > JmaxR) || (J > JmaxC) ) {
        mxErrMsgTxt("Input J must be an integer between 0 and log2(n), and dyadic for ROW and COL---use smaller J.");
    }
    l = (mxGetM(prhs[1]) > mxGetN(prhs[1])) ? mxGetM(prhs[1]) : mxGetN(prhs[1]);
    l1 = (mxGetM(prhs[2]) > mxGetN(prhs[2])) ? mxGetM(prhs[2]) : mxGetN(prhs[2]);
    if (l != l1) {
        mexErrMsgTxt("Filters must be the same length");
    }
    
    symr = mxGetScalar(prhs[6]);
    symc = mxGetScalar(prhs[7]);
    if ((symr > 2) || (symc > 2)) {
        mexErrMsgTxt("Symmetry flag must be 0, 1, or 2");
    }
    
  /* Create matrix for the return argument. */
    plhs[0] = mxCreateDoubleMatrix(nR, nC, mxREAL);
    
  /* Assign pointers to each input and output. */
    w = mxGetPr(prhs[0]);
    h0_r = mxGetPr(prhs[1]);
    h1_r = mxGetPr(prhs[2]);
    h0_c = mxGetPr(prhs[3]);
    h1_c = mxGetPr(prhs[4]);
    
    u = mxGetPr(plhs[0]);
    
  /* Call the C subroutine. */
    adj_wavelet_multi_level2D_CWT(u, w, nR, nC, h0_r, h1_r, h0_c, h1_c, l, J, symr, symc);
    return;
}
