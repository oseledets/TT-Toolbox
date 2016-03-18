#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <string.h>
#include <math.h>


/*
 * Cross approximation with full pivoting
 */

/* Fortran housekeeping */
mwIndex ione = 1;
char cN = 'N';
char cT = 'T';
char cR = 'R';
char cU = 'U';
double done = 1.0;
double dzero = 0.0;

/* Correctly determine three array dimensions */
void matlab_size3(const mxArray *prhs, mwIndex *dims, double **A)
{
  mwIndex dimcount, *rhsdims, i;

  *A = mxGetPr(prhs);
  dimcount=mxGetNumberOfDimensions(prhs);
  rhsdims = (mwIndex *)mxGetDimensions(prhs);
  for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
  for (i=dimcount; i<3; i++) dims[i]=1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/* RHS: Y, tol
   LHS: u,v,ind
*/
{
  double *Y, *u, *vt, *res, *ind, val_max, val, tol, *work;
  mwIndex dims[3], n, m, minsz, sz, r, piv1, piv2, i;

  if (nrhs<2) { mexPrintf("Give at least Y and tol\n"); return; }

  /* Fetch the data. We need to determine the sizes correctly */
  matlab_size3(prhs[0], dims, &Y);
  n = dims[0];
  m = dims[1]*dims[2];
  matlab_size3(prhs[1], dims, &u);
  tol = *u;

  minsz=n; if (m<minsz) minsz = m;

  u = (double *)malloc(sizeof(double)*n*minsz);
  vt = (double *)malloc(sizeof(double)*minsz*m);
  res = (double *)malloc(sizeof(double)*n*m);
  ind = (double *)malloc(sizeof(double)*minsz);

  sz = n*m;
  dcopy_(&sz, Y, &ione, res, &ione);
  /* i = idamax_(&sz, Y, &ione);
  val_max = fabs(Y[i]); */

  val_max = 0.0;
  for (i=0; i<n*m; i++) {
    if (fabs(Y[i])>val_max) {
      val_max = fabs(Y[i]);
    }
  }

  /* Main loop */
  for (r=0; r<minsz; r++) {
    /* Find the maximal element */
    val=0.0;
    for (i=0; i<n*m; i++) {
      if (fabs(res[i])>val) {
	val = fabs(res[i]);
	piv2 = i;
      }
    }

    if (val<=tol*val_max) break;
    piv1 = piv2 % n;
    piv2 = (int)(piv2/n);
    ind[r] = (double)(piv1+1);
    dcopy_(&n, &res[piv2*n], &ione, &u[r*n], &ione);
    dcopy_(&m, &res[piv1], &n, &vt[r*m], &ione);
    val = 1.0/res[piv1+piv2*n];
    dscal_(&m,&val,&vt[r*m],&ione);
    /* res = res - u(:,r)*v(r,:) */
    val = -1.0;
    dger_(&n, &m, &val, &u[r*n], &ione, &vt[r*m],&ione, res, &n);
  }
  if (r==0) {
    /* There was a zero matrix */
    r=1;
    memset(u, 0, sizeof(double)*n);
    memset(vt, 0, sizeof(double)*m);
  }

  /* QR u */
  sz = -1;
  dgeqrf_(&n, &r, u, &n, &res[1], &res[0], &sz, &i);
  sz = (mwIndex)res[0];
  work = (double *)malloc(sizeof(double)*sz);
  dgeqrf_(&n, &r, u, &n, res, work, &sz, &i);
  dtrmm_(&cR, &cU, &cT, &cN, &m, &r, &done, u, &n, vt, &m);
  dorgqr_(&n, &r, &r, u, &n, res, work, &sz, &i);
  
  free(work);
  free(res);

  /* Return the outputs */
  plhs[0] = mxCreateDoubleMatrix(n, r, mxREAL);
  res = mxGetPr(plhs[0]);
  sz = n*r;
  dcopy_(&sz, u, &ione, res, &ione);
  plhs[1] = mxCreateDoubleMatrix(r, m, mxREAL);
  res = mxGetPr(plhs[1]);
  /* vt should be transposed */
  for (i=0; i<r; i++) {
    dcopy_(&m, &vt[i*m], &ione, &res[i], &r);
  }
  /* return ind */
  if (nlhs>2) {
    plhs[2] = mxCreateDoubleMatrix(ione, r, mxREAL);
    res = mxGetPr(plhs[2]);
    dcopy_(&r, ind, &ione, res, &ione);
  }

  free(u);
  free(vt);
  free(ind);
}
