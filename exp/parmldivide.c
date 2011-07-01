#include "mex.h"
#include "matrix.h"
#include <string.h>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
/* Input: double matrices [A1, A2,...], double rhs [f1, f2, ...]
 Output: solution [u1, u2, ...] */
{
  int n, r;
  int i;
  int ret;
  char *funcname = "mldivide";
  char is_complex;
  double *dallM, *dallf, *dcurM, *dcurf, *dcuru, *dallu;
  mxArray *curMf[2],  *curu;
  double *callM, *callf, *ccurM, *ccurf, *ccuru, *callu;
  mwSize *dims;
  
  if (nrhs<2) {
    mexPrintf("Please specify matrices ([A2,A2,...]) and rhs ([f1,f2,...])\n");
    return;
  }
  
  dims = (mwSize *)malloc(sizeof(mwSize)*2);
  
  n = mxGetM(prhs[0]); 
  r = mxGetN(prhs[1]);
  
  if ((mxIsComplex(prhs[0]))||(mxIsComplex(prhs[1]))) {
    plhs[0] = mxCreateDoubleMatrix(n,r,mxCOMPLEX);
    dallu = mxGetPr(plhs[0]);
    callu = mxGetPi(plhs[0]);
    is_complex=1;
  } else {
    plhs[0] = mxCreateDoubleMatrix(n,r,mxREAL);
    dallu = mxGetPr(plhs[0]);
    is_complex=0;
  }
  
  if (is_complex==1) {
    dims[0]=n; dims[1]=n;
    curMf[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    dcurM = mxGetPr(curMf[0]);
    ccurM = mxGetPi(curMf[0]);
    dims[0]=n; dims[1]=1;
    curMf[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxCOMPLEX);
    dcurf = mxGetPr(curMf[1]);
    ccurf = mxGetPi(curMf[1]);    
  } else {
    dims[0]=n; dims[1]=n;
    curMf[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    dcurM = mxGetPr(curMf[0]);    
    dims[0]=n; dims[1]=1;
    curMf[1] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    dcurf = mxGetPr(curMf[1]);    
  }
  
  dallM = mxGetPr(prhs[0]);
  dallf = mxGetPr(prhs[1]);
  
  if (is_complex==1) {
    callM = mxGetPi(prhs[0]);
    callf = mxGetPi(prhs[1]);    
  }
  
  for (i=0; i<r; i++) {
    memcpy(dcurM, &(dallM[i*n*n]), sizeof(double)*n*n);
    memcpy(dcurf, &(dallf[i*n]), sizeof(double)*n);
    if (is_complex==1) {
      memcpy(ccurM, &(callM[i*n*n]), sizeof(double)*n*n);
      memcpy(ccurf, &(callf[i*n]), sizeof(double)*n);      
    }
    ret = mexCallMATLAB(1, &curu, 2, &curMf, funcname);
    if (ret!=0) {
      mexPrintf("mldivide returned with error %d on system %d\n", ret, i);
    }
    if (is_complex==1) {
      dcuru = mxGetPr(curu);
      ccuru = mxGetPi(curu);
      memcpy(&(dallu[i*n]), dcuru, sizeof(double)*n);
      memcpy(&(callu[i*n]), ccuru, sizeof(double)*n);
    } else {
      dcuru = mxGetPr(curu);
      memcpy(&(dallu[i*n]), dcuru, sizeof(double)*n);
    }
    mxDestroyArray(curu);
  }
  
  mxDestroyArray(curMf[0]);
  mxDestroyArray(curMf[1]);

}
