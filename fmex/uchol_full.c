double dCompElem(double *dA, int i, int j, int m)
{
    double res=0.0;
    int k;
    for (k=0; k<m; k++) {
        res += dA[k+j*m]*dA[k+i*m];
    }
    return res;
}

long duchol(long m, long n, double *dA, long numvec, double *du, double *lam)
{
    double *ddiag, *du2, *work, *dx, *dv;
    long i,j,cnt;
//     long m = m0[0];
//     long n = n0[0];
//     long numvec = numvec0[0];
    long maxi;
    double maxa, maxx;
    char jobs;
    char uplo = 'U';
    long lwork;

//     printf("m: %d, n: %d, nv: %d\n", m, n, numvec);

    dx = (double *)malloc(sizeof(double)*numvec);
    dv = (double *)malloc(sizeof(double)*numvec*numvec);

    lwork = 3*numvec;
//     printf("%d\n", sizeof(lwork));
    work = (double *)malloc(sizeof(double)*lwork);

//     du = (double *)malloc(sizeof(double)*n*numvec);
//     lam = (double *)malloc(sizeof(double)*numvec);

    for (j=0; j<numvec; j++) lam[j]=0.0;

    du2 = (double *)malloc(sizeof(double)*n*numvec);

    ddiag = (double *)malloc(sizeof(double)*n);
    /* Initial diagonal */
    for (i=0; i<n; i++) {
        ddiag[i] = dCompElem(dA, i, i, m);
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
	j = cnt+1;
	dgeqrf_(&n, &j, du, &n, du2, work, &lwork, &i);
        for (j=0; j<=cnt; j++) dx[j] = du[j+cnt*n];
	dorgqr_(&n, &j, &j, du, &n, du2, work, &lwork, &i);

        /* new eig */
        for (i=0; i<=cnt; i++) {
            for (j=0; j<=cnt; j++) {
                if (i==j)
                    dv[j+j*numvec]=lam[j];
                else
                    dv[i+j*numvec]=0.0;
            }
        }
        for (i=0; i<=cnt; i++) {
            for (j=0; j<=cnt; j++) {
                dv[i+j*numvec]+=dx[i]*dx[j];
            }
        }
        j=cnt+1;
        jobs = 'V';
        dsyev_(&jobs, &uplo, &j, dv, &numvec, lam, work, &lwork, &i);

        /* update u */
        jobs = 'N';
        j=cnt+1;
        work[0]=1.0; work[1]=0.0;
        dgemm_(&jobs,&jobs,&n,&j,&j,&work[0],du,&n,dv,&numvec,&work[1],du2,&n);
        for (j=0; j<=cnt; j++) {
             for (i=0; i<n; i++) du[i+j*n]=du2[i+j*n];
        }
//         printf("vector %d done\n", cnt);
    }


    free(ddiag);
    free(dx);
    free(dv);
    free(work);
    free(du2);

//     if (cnt<numvec) cnt++;
    return cnt;
}


// dummy wrapper for fortran
void duchol_(long *m, long *n, double *dA, long *numvec, double *du, double *lam, long *real_numvec)
{
  real_numvec[0] = duchol(m[0], n[0], dA, numvec[0], du, lam);
}
