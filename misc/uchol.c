#include "mex.h"
#include "lapack.h"
#include "blas.h"
/*
double dCompElem(double *dA, int i, int j, int m)
{
    double res=0.0;
    int k;
    for (k=0; k<m; k++) {
        res += dA[k+j*m]*dA[k+i*m];
    }
    return res;
}*/

long duchol(long m, long n, double *dA, long numvec, double *du, double *lam);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *dA, *scal;
    double *lam, *du;
    long numvec, i, j, m, n;
/*    double *dv, *du2, *ddiag;
    long m,n;
    long i,j,cnt;
    long maxi;
    double maxa, maxx;*/
//     double *dx;
//     char jobs;
//     char uplo = 'U';
//     double *work;
//     long lwork;

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    dA = mxGetPr(prhs[0]);

    scal = mxGetPr(prhs[1]);
    numvec = (long)(round(scal[0]));

    if (numvec>n) numvec = n;

    du = (double *)malloc(sizeof(double)*n*numvec);
    lam = (double *)malloc(sizeof(double)*numvec);

    numvec = duchol(m, n, dA, numvec, du, lam);

    plhs[0] = mxCreateDoubleMatrix(n, numvec, mxREAL);
    scal = mxGetPr(plhs[0]);
    for (j=0; j<numvec; j++) {
        for (i=0; i<n; i++) scal[i+j*n]=du[i+j*n];
    }

    free(du);

    if (nlhs>1) {
        plhs[1] = mxCreateDoubleMatrix(numvec,1, mxREAL);
        scal = mxGetPr(plhs[1]);
        for (j=0; j<numvec; j++) scal[j]=lam[j];
    }

    free(lam);
}
