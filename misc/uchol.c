#include "mex.h"
#include "lapack.h"
#include "blas.h"

double dCompElem(double *dA, int i, int j, int m)
{
    double res=0.0;
    int k;
    for (k=0; k<m; k++) {
        res += dA[k+j*m]*dA[k+i*m];
    }
    return res;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *dA, *scal, *ddiag;
    double *lam, *du, *dv, *du2;
    long m,n, numvec;
    long i,j,cnt;
    long maxi;
    double maxa, maxx;
    double *dx;
    char jobs;
    char uplo = 'U';
    double *work;
    long lwork;

    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);
    dA = mxGetPr(prhs[0]);

    scal = mxGetPr(prhs[1]);
    numvec = (int)(round(scal[0]));

    if (numvec>n) numvec = n;

    dx = (double *)malloc(sizeof(double)*numvec);
    dv = (double *)malloc(sizeof(double)*numvec*numvec);

    lwork = 3*numvec;
    work = (double *)malloc(sizeof(double)*lwork);

    du = (double *)malloc(sizeof(double)*n*numvec);
    lam = (double *)malloc(sizeof(double)*numvec);

    for (j=0; j<numvec; j++) lam[j]=0.0;

    du2 = (double *)malloc(sizeof(double)*n*numvec);

    ddiag = (double *)malloc(sizeof(double)*n);
    /* Initial diagonal */
    for (i=0; i<n; i++) {
        ddiag[i] = dCompElem(dA, i, i, m);
//         if (ddiag[i]>maxa) {
//             maxi = i;
//             maxa = ddiag[i];
//         }
    }

    /* Main cycle */
    for (cnt=0; cnt<numvec; cnt++) {
        /* maxel */
        maxi = 0; maxa = 0.0;
        for (i=0; i<n; i++) {
            if (fabs(ddiag[i])>maxa) {
                maxi = i;
                maxa = fabs(ddiag[i]);
            }
        }

        /* residual */
        for (i=0; i<n; i++) {
            du[i+cnt*n] = dCompElem(dA, i, maxi, m);
            for (j=0; j<cnt; j++) {
                du[i+cnt*n] -= du[i+j*n]*lam[j]*du[maxi+j*n];
            }
        }
        maxa = sqrt(fabs(du[maxi+cnt*n]));
        if (maxa<1e-300) { // We've got the exact zero matrix
            break;
        }
        for (i=0; i<n; i++) {
            du[i+cnt*n] /= maxa;
            ddiag[i] = ddiag[i]-(du[i+cnt*n]*du[i+cnt*n]);
        }
        /* reorth */
// 	for (i=0; i<n; i++) du2[i]=du[i+cnt*n];
	j = cnt+1;
	dgeqrf_(&n, &j, du, &n, du2, work, &lwork, &i);
        for (j=0; j<=cnt; j++) dx[j] = du[j+cnt*n];
	dorgqr_(&n, &j, &j, du, &n, du2, work, &lwork, &i);

// 	for (j=0; j<=cnt; j++) {
// 	  dx[j]=0.0;
// 	  for (i=0; i<n; i++) {
// 	    dx[j] +=du[i+j*n]*du2[i];
// 	  }
// 	}


//         for (j=0; j<cnt; j++) {
//             dx[j] = 0.0;
//             for (i=0; i<n; i++) {
//                 dx[j]+=du[i+j*n]*du[i+cnt*n];
//             }
//         }
//         maxa = 0.0;
//         for (i=0; i<n; i++) {
//             for (j=0; j<cnt; j++) {
//                 du[i+cnt*n] -= du[i+j*n]*dx[j];
//             }
//             maxa += du[i+cnt*n]*du[i+cnt*n];
//         }
//         maxa = sqrt(maxa);
//         for (i=0; i<n; i++) {
//             du[i+cnt*n] /= maxa;
//         }
        /* new eig */
//         dx[cnt]=maxa;
        for (i=0; i<=cnt; i++) {
            for (j=0; j<=cnt; j++) {
                if (i==j)
                    dv[j+j*numvec]=lam[j];
                else
                    dv[i+j*numvec]=0.0;
            }
        }
//         dv[cnt+cnt*numvec]=0.0;
        for (i=0; i<=cnt; i++) {
            for (j=0; j<=cnt; j++) {
                dv[i+j*numvec]+=dx[i]*dx[j];
//                 dv[i+i*numvec]+=dx[cnt]*dx[j];
            }
        }
        j=cnt+1;
        jobs = 'V';
        dsyev_(&jobs, &uplo, &j, dv, &numvec, lam, work, &lwork, &i);
//         for (j=0; j<=cnt; j++) {
//             mexPrintf("cnt=%d, lam[%d]=%g\n", cnt, j, lam[j]);
//         }
        /* update u */
        jobs = 'N';
        j=cnt+1;
        work[0]=1.0; work[1]=0.0;
        dgemm_(&jobs,&jobs,&n,&j,&j,&work[0],du,&n,dv,&numvec,&work[1],du2,&n);
/*        for (i=0; i<n; i++) {
            for (j=0; j<=cnt; j++) {
                du2[i+j*n]=0.0;
                for (maxi=0; maxi<=cnt; maxi++) {
                    du2[i+j*n] += du[i+maxi*n]*dv[maxi+j*numvec];
                }
            }
        }*/
        for (i=0; i<n; i++) {
            for (j=0; j<=cnt; j++) du[i+j*n]=du2[i+j*n];
        }
    }


    free(ddiag);
    free(dx);
    free(dv);
    free(work);
    free(du2);

    plhs[0] = mxCreateDoubleMatrix(n, cnt, mxREAL);
    du2 = mxGetPr(plhs[0]);
    for (j=0; j<cnt; j++) {
        for (i=0; i<n; i++) du2[i+j*n]=du[i+j*n];
    }

    free(du);

    if (nlhs>1) {
        plhs[1] = mxCreateDoubleMatrix(cnt,1, mxREAL);
        dx = mxGetPr(plhs[1]);
        for (j=0; j<cnt; j++) dx[j]=lam[j];
    }

    free(lam);
}
