#include "mex.h"
#include "lapack.h"

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
    char jobs = 'V';
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

    plhs[0] = mxCreateDoubleMatrix(n, numvec, mxREAL);
    du = mxGetPr(plhs[0]);

    if (nlhs>1) {
        plhs[1] = mxCreateDoubleMatrix(numvec,1, mxREAL);
        lam = mxGetPr(plhs[1]);
    }
    else {
        lam = (double *)malloc(sizeof(double)*numvec);
    }

    du2 = (double *)malloc(sizeof(double)*n*numvec);

    ddiag = (double *)malloc(sizeof(double)*n);
    /* Initial diagonal */
    for (i=0; i<n; i++) {
        ddiag[i] = dCompElem(dA, i, i, m);
        if (ddiag[i]>maxa) {
            maxi = i;
            maxa = ddiag[i];
        }
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
        for (i=0; i<n; i++) {
            du[i+cnt*n] /= maxa;
            ddiag[i] = ddiag[i]-(du[i+cnt*n]*du[i+cnt*n]);
        }
        /* reorth */
        for (j=0; j<cnt; j++) {
            dx[j] = 0.0;
            for (i=0; i<n; i++) {
                dx[j]+=du[i+j*n]*du[i+cnt*n];
            }
        }
        maxa = 0.0;
        for (i=0; i<n; i++) {
            for (j=0; j<cnt; j++) {
                du[i+cnt*n] -= du[i+j*n]*dx[j];
            }
            maxa += du[i+cnt*n]*du[i+cnt*n];
        }
        maxa = sqrt(maxa);
        for (i=0; i<n; i++) {
            du[i+cnt*n] /= maxa;
        }
        /* new eig */
        dx[cnt]=maxa;
        for (i=0; i<cnt; i++) {
            for (j=0; j<cnt; j++) {
                if (i==j)
                    dv[j+j*numvec]=lam[j];
                else
                    dv[i+j*numvec]=0.0;
            }
        }
        for (j=0; j<=cnt; j++) {
            dv[j+cnt*numvec]=dx[j]*dx[cnt];
            dv[cnt+j*numvec]=dx[cnt]*dx[j];
        }
        j=cnt+1;
        dsyev_(&jobs, &uplo, &j, dv, &numvec, lam, work, &lwork, &i);
//         for (j=0; j<=cnt; j++) {
//             mexPrintf("cnt=%d, lam[%d]=%g\n", cnt, j, lam[j]);
//         }
        /* update u */
        for (i=0; i<n; i++) {
            for (j=0; j<=cnt; j++) {
                du2[i+j*n]=0.0;
                for (maxi=0; maxi<=cnt; maxi++) {
                    du2[i+j*n] += du[i+maxi*n]*dv[maxi+j*numvec];
                }
            }
        }
        for (i=0; i<n; i++) {
            for (j=0; j<=cnt; j++) du[i+j*n]=du2[i+j*n];
        }
    }


    free(ddiag);
    if (nlhs==0) free(lam);
    free(dx);
    free(dv);
    free(work);
    free(du2);
}
