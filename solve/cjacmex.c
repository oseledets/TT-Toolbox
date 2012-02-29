/*=========================================================
 * matrixMultiply.c - Example for illustrating how to use
 * BLAS within a C MEX-file. matrixMultiply calls the
 * BLAS function dgemm.
 *
 * C = matrixMultiply(A,B) computes the product of A*B,
 *     where A, B, and C are matrices.
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2009 The MathWorks, Inc.
 *=======================================================*/
/* $Revision: 1.1.6.1 $ $Date: 2009/03/16 21:52:28 $ */

#if !defined(_WIN32)
#define dgemm dgemm_
#endif

#include "mex.h"
#include "blas.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C; /* pointers to input & output matrices*/
    mwSignedIndex m,n,p;      /* matrix dimensions */
    mwSize dimcount;
    mwSize *dims;
    const mwSize *rhsdims;
    /* form of op(A) & op(B) to use in matrix multiplication */
    char *chn = "N";
    /* scalar values to use in dgemm */
    double one = 1.0, zero = 0.0;

    A = mxGetPr(prhs[0]); /* first input matrix */
    B = mxGetPr(prhs[1]); /* second input matrix */

    dimcount=mxGetNumberOfDimensions(prhs[0]);
    dims = (mwSize *)malloc(sizeof(mwSize)*4);
    rhsdims = mxGetDimensions(prhs[0]);
    for (n=0; n<dimcount; n++) dims[n]=rhsdims[n];
    for (n=dimcount; n<4; n++) dims[n]=1;
//     for (n=0; n<4; n++) mexPrintf("dims[%d]=%d\n", n, dims[n]);
    /* dimensions of input matrices */
/*
     m = (mwSignedIndex)mxGetM(prhs[0]);
     p = (mwSignedIndex)mxGetN(prhs[0]);
     n = (mwSignedIndex)mxGetN(prhs[1]);

     if (p != mxGetM(prhs[1])) {
         mexErrMsgIdAndTxt("MATLAB:matrixMultiply:matchdims",
                 "Inner dimensions of matrix multiply do not match.");
     }
*/
    /* create output matrix C */
    /* n*r1*r2 */
    plhs[0] = mxCreateDoubleMatrix(dims[1]*dims[2]*dims[3], 1, mxREAL);
    C = mxGetPr(plhs[0]);

    /* Pass arguments to Fortran by reference */
    p = 1;
    for (m=0; m<dims[3]; m++) {
        for (n=0; n<dims[2]; n++) {
            dgemm(chn, chn, &dims[0], &p, &dims[1], &one, &A[n*dims[0]*dims[1] + m*dims[0]*dims[1]*dims[2]], &dims[0], &B[n*dims[1] + m*dims[1]*dims[2]], &dims[1], &zero, &C[n*dims[0] + m*dims[0]*dims[2]], &dims[0]);
        }
    }

    free(dims);
}
