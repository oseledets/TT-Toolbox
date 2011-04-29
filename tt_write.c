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

int tt_write(char *fname, tthead *hd, int d, int *r, int *n, char *core)
{
  int fid;
  int mem=0;
  int curmem, lm[2], coresize;

  fid = open(fname, O_WRONLY | O_CREAT | O_TRUNC, 00644);
  if (fid<0) {
    mexPrintf("The file %s can not be opened\n", fname);
    return -1;
  }

  // compute the size of the core
  coresize=0;
  for (curmem=0; curmem<d; curmem++) {
    coresize+=r[curmem]*n[curmem]*r[curmem+1];
  }

  curmem = write(fid, hd, sizeof(tthead)); // write out the head
  mem+=curmem;

  lm[0]=1; lm[1]=d;
  curmem = write(fid, &lm, sizeof(int)*2); // write out l and m
  mem+=curmem;

  curmem = write(fid, n, sizeof(int)*d); // write out n
  mem+=curmem;

  curmem = write(fid, r, sizeof(int)*(d+1)); // write out r
  mem+=curmem;

  switch (hd[0].inf[1]) {
    case 0: { // double, real
      curmem = write(fid, core, sizeof(double)*coresize); // write out core
      break;
    }
    case 1: { // double, complex
      curmem = write(fid, core, sizeof(double complex)*coresize); // write out core
      break;
    }
    case 2: { // float, real
      curmem = write(fid, core, sizeof(float)*coresize); // write out core
      break;
    }
    case 3: { // float, complex
      curmem = write(fid, core, sizeof(float complex)*coresize); // write out core
      break;
    }
    default: {
      mexPrintf("Wrong datatype %d\n", hd[0].inf[1]);
      return -1;
    }
  }
  mem+=curmem;

  close(fid);

  return mem;
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
// char *fname, double d, double r[d], double n[d], double [complex] core[sososo], double precision: 0,1 (double, float)
// return: bytes written
{
  tthead head;
  int d, *r, *n, j, coresize, mem;
  char precision;
  double complex *zcore;
  double *dcore;
  float complex *ccore;
  float *score;
  double *num1, *num2;
  char *filename;

  if (nrhs<6) {
    mexPrintf("Please specify filename, d, r, n, core and precision.\n");
    return;
  }

  if (sizeof(int)!=4) {
    mexPrintf("sizeof(int) is not 4. The header will be inconsistent.\n");
    return;
  }

// initialize head
  sprintf(head.txt, "TT      ");
  head.ver[0]=1; head.ver[1]=0;
  memset(head.inf, 0, sizeof(int)*4);
  head.inf[0]=128;
  memset(head.i, 0, sizeof(int)*8);
  memset(head.comment, 0, sizeof(char)*64);

// get filename
  coresize = mxGetM(prhs[0]);
  if (coresize==1) coresize = mxGetN(prhs[0]);
  filename = (char *)malloc(sizeof(char)*(coresize+1));
  memset(filename, 0, sizeof(char)*(coresize+1));
  mxGetString(prhs[0], filename, coresize+1);

// get d
  num1 = mxGetPr(prhs[1]);
  d = (int)(round(num1[0]));

  head.i[0]=1; head.i[2]=d;

// prepare ints for ranks and mode sizes
  r = (int *)malloc(sizeof(int)*(d+1));
  n = (int *)malloc(sizeof(int)*d);
// now fill'em
  num1 = mxGetPr(prhs[2]);
  num2 = mxGetPr(prhs[3]);
  for (j=0; j<d; j++) {
    r[j]=(int)(round(num1[j]));
    n[j]=(int)(round(num2[j]));
  }
  r[d]=(int)(round(num1[d]));

// get precision
  num1 = mxGetPr(prhs[5]);
  precision = (char)(round(num1[0]));

// now, the core
  coresize = mxGetM(prhs[4]);
  if (coresize==1) coresize = mxGetN(prhs[4]);

  if (mxIsComplex(prhs[4])) {
    num1 = mxGetPr(prhs[4]); // real part
    num2 = mxGetPi(prhs[4]); // imag part

    head.inf[1]=1+2*precision;
    if (precision==0) {
      zcore = (double complex *)malloc(sizeof(double complex)*coresize);
      for (j=0; j<coresize; j++) zcore[j]=num1[j]+I*num2[j];
      // write out the complex tt
      mem = tt_write(filename, &head, d, r, n, (char *)zcore);
      free(zcore);
    }
    else {
      ccore = (float complex *)malloc(sizeof(float complex)*coresize);
      for (j=0; j<coresize; j++) ccore[j]=num1[j]+I*num2[j];
      // write out the complex tt
      mem = tt_write(filename, &head, d, r, n, (char *)ccore);
      free(ccore);
    }
  }
  else {
    dcore = mxGetPr(prhs[4]);

    head.inf[1]=0+2*precision;
    if (precision==0) {
      // write out the real tt
      mem = tt_write(filename, &head, d, r, n, (char *)dcore);
    }
    else {
      score = (float *)malloc(sizeof(float)*coresize);
      for (j=0; j<coresize; j++) score[j]=dcore[j];
      // write out the complex tt
      mem = tt_write(filename, &head, d, r, n, (char *)score);
      free(score);
    }
  }

  // return the bytes written
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  num1 = mxGetPr(plhs[0]);
  num1[0]=(double)mem;

  free(r);
  free(n);
  free(filename);
}
