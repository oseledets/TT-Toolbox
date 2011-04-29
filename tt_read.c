#include "mex.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "complex.h"

// mex -v tt_write.c CC=icc CFLAGS="-fexceptions -fPIC -fno-omit-frame-pointer -pthread" COPTIMFLAGS="-O3" LDOPTIMFLAGS="-O3"

typedef struct {
  char txt[8];
  int ver[2];
  int inf[4];
  char comment[64];
  int i[8];
} tthead;

int tt_read(char *fname, tthead *hd, int *d, int **r, int **n, char **core, char *dtype, int *coresize)
{
  int fid;
  int mem=0;
  int curmem, lm[2], _coresize;

  fid = open(fname, O_RDONLY);
  if (fid<0) {
    mexPrintf("The file %s can not be opened\n", fname);
    return -1;
  }

  curmem = read(fid, hd, sizeof(tthead)); // read the head
  mem+=curmem;
  if (strncmp(hd[0].txt, "TT", 2)!=0) {
    mexPrintf("Not a SDV TT file %s\n", fname);
    return -1;
  }

  if (hd[0].ver[0]!=1) {
    mexPrintf("Wrong version of SDV TT file (current is 1), received %d\n", hd[0].ver[0]);
    return -1;
  }

  dtype[0]=hd[0].inf[1]; // we'll check the real/complex/double/float numbers

  curmem = read(fid, &lm, sizeof(int)*2); // read l and m
  mem+=curmem;

  d[0] = lm[1]-lm[0]+1;
  r[0] = (int *)malloc(sizeof(int)*(d[0]+1));
  n[0] = (int *)malloc(sizeof(int)*d[0]);

  curmem = read(fid, n[0], sizeof(int)*d[0]); // read n
  mem+=curmem;

  curmem = read(fid, r[0], sizeof(int)*(d[0]+1)); // read r
  mem+=curmem;

  // compute the size of the core
  _coresize=0;
  for (curmem=0; curmem<d[0]; curmem++) {
    _coresize+=r[0][curmem]*n[0][curmem]*r[0][curmem+1];
  }

  coresize[0]=_coresize;

  // Now, the most tricky part.
  // We have to fetch the core corresponding to the data type.
  switch (dtype[0]) {
    case 0: { // double, real
      core[0] = (char *)malloc(sizeof(double)*_coresize);
      curmem = read(fid, core[0], sizeof(double)*_coresize); // read core
      break;
    }
    case 1: { // double, complex
      core[0] = (char *)malloc(sizeof(double complex)*_coresize);
      curmem = read(fid, core[0], sizeof(double complex)*_coresize); // read core
      break;
    }
    case 2: { // float, real
      core[0] = (char *)malloc(sizeof(float)*_coresize);
      curmem = read(fid, core[0], sizeof(float)*_coresize); // read core
      break;
    }
    case 3: { // float, complex
      core[0] = (char *)malloc(sizeof(float complex)*_coresize);
      curmem = read(fid, core[0], sizeof(float complex)*_coresize); // read core
      break;
    }
    default: {
      mexPrintf("Wrong datatype %d\n", dtype[0]);
      return -1;
    }
  }

  mem+=curmem;


  close(fid);

  return mem;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
// // char *fname
// return: d, r, n, core
// if d<1, then smth is wrong
{
  tthead head;
  int d, *r, *n, mem, coresize, j;
  char dtype;
  double complex *zcore;
  double *dcore;
  float complex *ccore;
  float *score;
  double *num1, *num2;
  char *filename, *buf;

  if (nrhs<1) {
    mexPrintf("Please specify filename.\n");
    return;
  }

  if (sizeof(int)!=4) {
    mexPrintf("sizeof(int) is not 4. The header will be inconsistent. Try to play with compiler options.\n");
    return;
  }

// get filename
  coresize = mxGetM(prhs[0]);
  if (coresize==1) coresize = mxGetN(prhs[0]);
  filename = (char *)malloc(sizeof(char)*(coresize+1));
  memset(filename, 0, sizeof(char)*(coresize+1));
  mxGetString(prhs[0], filename, coresize+1);

  mem=tt_read(filename, &head, &d, &r, &n, &buf, &dtype, &coresize);

  // create and report d
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  num1 = mxGetPr(plhs[0]);
  num1[0]=(double)d;  // if evrthg is ok - d
  if (mem<0) {
    num1[0]=(double)mem; // if error - return it
    free(filename);
    return;
  }

  // create and report r
  plhs[1] = mxCreateDoubleMatrix(d+1,1,mxREAL);
  num1 = mxGetPr(plhs[1]);
  // create and report n
  plhs[2] = mxCreateDoubleMatrix(d,1,mxREAL);
  num2 = mxGetPr(plhs[2]);
  for (j=0; j<d; j++) {
    num1[j]=(double)r[j];
    num2[j]=(double)n[j];
  }
  num1[d]=(double)r[d];


  switch(dtype) {
    case 0: { // double, real
      dcore = (double *)buf;
      // create and report core
      plhs[3] = mxCreateDoubleMatrix(coresize,1,mxREAL);
      num1 = mxGetPr(plhs[3]);
      memcpy(num1, dcore, sizeof(double)*coresize);
      break;
    }
    case 1: { // double, complex
      zcore = (double complex *)buf;
      // create and report core
      plhs[3] = mxCreateDoubleMatrix(coresize,1,mxCOMPLEX);
      num1 = mxGetPr(plhs[3]);
      num2 = mxGetPi(plhs[3]);
      for (j=0; j<coresize; j++) {
	num1[j] = creal(zcore[j]);
	num2[j] = cimag(zcore[j]);
      }
      break;
    }
    case 2: { // float, real
      score = (float *)buf;
      // create and report core
      plhs[3] = mxCreateDoubleMatrix(coresize,1,mxREAL);
      num1 = mxGetPr(plhs[3]);
      for (j=0; j<coresize; j++) {
	num1[j] = score[j];
      }
      break;
    }
    case 3: { // float, complex
      ccore = (float complex *)buf;
      // create and report core
      plhs[3] = mxCreateDoubleMatrix(coresize,1,mxCOMPLEX);
      num1 = mxGetPr(plhs[3]);
      num2 = mxGetPi(plhs[3]);
      for (j=0; j<coresize; j++) {
	num1[j] = crealf(ccore[j]);
	num2[j] = cimagf(ccore[j]);
      }
      break;
    }
    default: {
      mexPrintf("Wrong datatype %d\n", dtype);
      return;
    }
  }

  free(r);
  free(n);
  free(filename);
  free(buf);
}
