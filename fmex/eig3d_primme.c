#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include "tt-fort/primme/primme.h"

int g_rx1, g_n, g_rx2, g_B, g_ra1, g_ra2;
double *g_Phi1, *g_A, *g_Phi2, *g_res1, *g_res2;
double done = 1.0;
double dzero = 0.0;

void dtransp(long m, long n, double *in, double *out)
{
  long i,j;
  j = 1;
  for (i=0; i<m; i++) {
    dcopy(&n, &in[i], &m, &out[i*n], &j);
  }
}


void PrimmeMatvec(void *x, void *y, int *blockSize, primme_params *primme) {
// bfun3
// sizes of res1, res2: max(rx1*m*ra2*ry2, rx1*ra1*n*ry2)
   double *xvec, *yvec;
   char op;
   long sz1, sz2, sz3;
   long i;

   xvec = (double *)x;
   yvec = (double *)y;

   op = 'N';

   //phi2(rx2,ra2,ry2)
   //phi1(ry1,rx1,ra1)

   for (i=0; i<*blockSize; i++) {
    sz1 = g_rx1*g_n;
    sz2 = g_ra2*g_rx2;
    sz3 = g_rx2;
    dgemm(&op, &op, &sz1, &sz2, &sz3, &done, &xvec[i*sz1*sz3], &sz1, g_Phi2, &sz3, &dzero, g_res2, &sz1);
    dtransp(g_rx1, g_n*g_ra2*g_rx2, g_res2, g_res1);
    sz1 = g_ra1*g_n;
    sz2 = g_rx2*g_rx1;
    sz3 = g_n*g_ra2;
    dgemm(&op, &op, &sz1, &sz2, &sz3, &done, g_A, &sz1, g_res1, &sz3, &dzero, g_res2, &sz1);
    dtransp(g_ra1*g_n*g_rx2, g_rx1, g_res2, g_res1);
    sz1 = g_rx1;
    sz2 = g_n*g_rx2;
    sz3 = g_rx1*g_ra1;
    dgemm(&op, &op, &sz1, &sz2, &sz3, &done, g_Phi1, &sz1, g_res1, &sz3, &dzero, &yvec[i*sz1*sz2], &sz1);
   }
}


// We needed to make int 8-byte in primme.
// However, Matlab wants it 4-byte. Undefine...
#undef int

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Input: Phi1 [r1',r1,ra1], A[ra1,n',n,ra2], Phi2[r2,ra2,r2'], tol, B, sol_prev, max_matvecs = 500, max_basis_size = 100
// Output: sol, lambda, num_matvecs
{

    double *scal, *dsol_prev;
    double tol;
    double *dsol, *dlambda, *dnum_matvecs, *drnorms;
    long max_matvecs, max_basis_size;
    long dimcount, *rhsdims;
    long dims[4];
    long rx1, rx2, ra1, ra2, n, i, j, B;

    primme_params primme;
    primme_preset_method method;


    if (nrhs<6) { mexPrintf("Specify at least Phi1,A,Phi2,tol,B,sol_prev\n"); return; }
    if (nrhs<7) max_matvecs=500; else { scal = mxGetPr(prhs[6]); max_matvecs = (long)round(scal[0]); }
    if (nrhs<8) max_basis_size=100; else { scal = mxGetPr(prhs[7]); max_basis_size = (long)round(scal[0]); }

    // Fetch the data
    g_Phi1 = mxGetPr(prhs[0]);
    dimcount=mxGetNumberOfDimensions(prhs[0]);
    rhsdims = mxGetDimensions(prhs[0]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<3; i++) dims[i]=1;
    rx1 = dims[0];
    if (dims[1]!=rx1) { mexPrintf("Phi1 is not square!\n"); return; }
    ra1 = dims[2];

    g_A = mxGetPr(prhs[1]);
    dimcount=mxGetNumberOfDimensions(prhs[1]);
    rhsdims = mxGetDimensions(prhs[1]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<4; i++) dims[i]=1;
    if (ra1 != dims[0]) { mexPrintf("ra1 in Phi1 and A are not consistent!\n"); return; }
    n = dims[1];
    if (n!=dims[2]) { mexPrintf("A is not square!\n"); return; }
    ra2 = dims[3];

    g_Phi2 = mxGetPr(prhs[2]);
    dimcount=mxGetNumberOfDimensions(prhs[2]);
    rhsdims = mxGetDimensions(prhs[2]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<3; i++) dims[i]=1;
    rx2 = dims[0];
    if (dims[2]!=rx2) { mexPrintf("Phi2 is not square!\n"); return; }
    if (ra2 != dims[1]) { mexPrintf("ra2 in Phi2 and A are not consistent!\n"); return; }

    scal = mxGetPr(prhs[3]);
    tol = scal[0];

    scal = mxGetPr(prhs[4]);
    B = (long)round(scal[0]);

    dsol_prev = mxGetPr(prhs[5]);
    if ((mxGetM(prhs[5])*mxGetN(prhs[5]))!=(rx1*n*rx2*B)) { mexPrintf("SOL_PREV size is not consistent!\n"); return; }


    // Allocate the output
    plhs[0] = mxCreateDoubleMatrix(rx1*n*rx2, B, mxREAL);
//     dsol = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(B, 1, mxREAL);
    dlambda = mxGetPr(plhs[1]);
    if (nlhs>2) {
      plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
      dnum_matvecs = mxGetPr(plhs[2]);
    }

    // Allocate the work arrays
    drnorms = (double *)malloc(sizeof(double)*B);
    dsol = (double *)malloc(sizeof(double)*rx1*n*rx2*(max_basis_size+B+1));

    g_rx1 = rx1;
    g_n = n;
    g_rx2 = rx2;
    g_ra1 = ra1;
    g_ra2 = ra2;
    g_B = B;

    i = ra1;
    if (ra2>ra1) { i=ra2; }
    i = rx1*n*rx2*i;
    g_res1 = (double *)malloc(sizeof(double)*i);
    g_res2 = (double *)malloc(sizeof(double)*i);

    // Call the solver
    //Initialization of the primme stuff;
    primme_initialize(&primme);
    method=DYNAMIC;
    primme_set_method(method, &primme);

    primme.n = rx1*n*rx2;
    primme.numEvals = B;

    primme.matrixMatvec = PrimmeMatvec;
    primme.printLevel = 1;
    primme.maxMatvecs = max_matvecs;
    primme.minRestartSize = B+1;
    if (B+1>rx1*n*rx2-1) { primme.minRestartSize = rx1*n*rx2-1; }
    primme.maxBasisSize = max_basis_size;

    primme.eps = tol;
//     primme.aNorm = 1.0; // It should be 1 to truncate in the residual norm. Default is fro

    primme.initSize = B;

    i = rx1*n*rx2*B;
    j = 1;
    dcopy(&i, dsol_prev, &j, dsol, &j);

// solution
    i = dprimme(dlambda, dsol, drnorms, &primme);

    if (nlhs>2) { dnum_matvecs[0] = (double)(primme.stats.numMatvecs); }

// finalize primme
    primme_Free(&primme);

    if (i < 0) {
      mexPrintf("Primme failed with err=%d",i);
      return;
    }

// output
    i = rx1*n*rx2*B;
    dcopy(&i, dsol, &j, mxGetPr(plhs[0]), &j);

// free the rest
    free(dsol);
    free(drnorms);
    free(g_res1);
    free(g_res2);
}
