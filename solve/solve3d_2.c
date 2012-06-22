#include "mex.h"
#include "lapack.h"
#include "blas.h"
// #include "time.h"


//     typedef struct {
//                time_t   tv_sec;        // seconds
//                long     tv_nsec;       // nanoseconds
//            } tsc;
//
// char trans = 'N';
// char uplo = 'U';
// double done=1.0;
// double dzero=0.0;
long ione = 1;

extern void tt_adapt_als_mp_djac_gen_(char *ptype, long *rx1, long *n, long *rx2, long *ra1, long *ra2, double *Phi1, double *A, double *Phi2, double *jacs);
extern void tt_adapt_als_mp_djac_apply_(char *ptype, long *rx1, long *n, long *rx2, double *jacs, double *x, double *y, double *work1);
extern void tt_adapt_als_mp_bfun3_(long *rx1, long *m, long *rx2, long *ry1, long *n, long *ry2, long *ra1, long *ra2, double *phi1, double *A, double *phi2, double *x, double *y, double *res1, double *res2);
extern void tt_adapt_als_mp_transp_(long *n, long *m, double *A, double *B);

#define djac_apply tt_adapt_als_mp_djac_apply_
#define bfun3 tt_adapt_als_mp_bfun3_
#define djac_gen tt_adapt_als_mp_djac_gen_
#define transp tt_adapt_als_mp_transp_

extern char gmres_3d_printf_buf_[128];


void dcjacgen(double *Phi1,double *A,double *Phi2, long rx1, long n, long rx2, long ra1, long ra2,double *jacs)
{
    long i, k,m,cnt;
    double *diag1;
    double *diag2;
    double *Id, *B, *A2;
    long *ipiv, info;
    char trans = 'N';
    double done=1.0;
    double dzero=0.0;
    long ione = 1;

    diag1 = (double *)malloc(sizeof(double)*rx1*ra1);
    diag2 = (double *)malloc(sizeof(double)*rx2*ra2);

    for (k=0; k<ra1; k++) {
        for (m=0; m<rx1; m++) {
            diag1[m+k*rx1] = Phi1[m+m*rx1+k*rx1*rx1]; // rx1,ra1
        }
    }
    for (k=0; k<ra2; k++) {
        for (m=0; m<rx2; m++) {
            diag2[k+m*ra2] = Phi2[m+k*rx2+m*rx2*ra2]; // ra2,rx2
        }
    }
/*
    for (i=0; i<rx1*rx2*n*n; i++) jacs[i]=0.0;

    for (i=0; i<ra2; i++) {
        for (j=0; j<ra1; j++) {
            for (k=0; k<rx2; k++) {
                for (m=0; m<rx1; m++) {
                    for (cnt=0; cnt<n*n; cnt++) {
                        jacs[cnt+m*n*n+k*n*n*rx1] += A[j+cnt*ra1+i*n*n*ra1]*diag1[m+j*rx1]*diag2[k+i*rx2];
                    }
                }
            }
        }
    }
 */
// Vectorize this bydlocode
    A2 = (double *)malloc(sizeof(double)*rx1*n*n*ra2);
    // with phi1: sizes rx1,ra1 - ra1,n*n*ra2
    i = n*n*ra2;
    dgemm_(&trans,&trans,&rx1,&i,&ra1,&done,diag1,&rx1,A,&ra1,&dzero,A2,&rx1);
    // with phi2: sizes rx1*n*n,ra2 - ra2,rx2
    i = rx1*n*n;
    dgemm_(&trans,&trans,&i,&rx2,&ra2,&done,A2,&i,diag2,&ra2,&dzero,jacs,&i);
    free(A2);


    // permute so that n.n.rx1.rx2
    A2 = (double *)malloc(sizeof(double)*n*n*rx1*rx2);
    i = rx1; k = n*n*rx2;
    transp(&i, &k, jacs, A2);
//     dperm213(jacs,rx1, n*n, rx2, A2);

    Id = (double *)malloc(sizeof(double)*n*n);
    for (k=0; k<n; k++) {
        for (m=0; m<n; m++) {
            if (m==k) Id[m+k*n]=1.0;
            else      Id[m+k*n]=0.0;
        }
    }
    B = (double *)malloc(sizeof(double)*n*n);
    ipiv = (long*)malloc(sizeof(long)*n);

    for (k=0; k<rx2*rx1; k++) {
        i = n*n;
        dcopy_(&i,Id,&ione,B,&ione);
        dgesv_(&n, &n, &(A2[k*n*n]), &n, ipiv, B, &n, &info);
        if (info>0) mexPrintf("The component %d at ranks [%d,%d] is singular\n", info,m,k);
        dcopy_(&i,B,&ione,&(jacs[k*n*n]),&ione);
    }

    free(A2);
    free(diag1);
    free(diag2);
    free(Id);
    free(B);
    free(ipiv);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Phi1 [r1,r1',ra1], A[ra1,n,n',ra2], Phi2[r2,r2',ra2], rhs, tol, trunc_norm, sol_prev, prec, nrestart, niters, verb
{
    double *dPhi1, *dA, *dPhi2, *drhs, *scal, *dsol_prev, *dres;
    double tol, tol_prev;
    double *dsol, *djacs, *work1, *work2;
    double dbeta;
//     long *ipiv;
    long nrestart, niters;
    long dimcount, *rhsdims;
    long dims[4];
    long rx1, rx2, ra1, ra2, n, i;
    char prec, verb, trunc_norm;
    char ptype = 'n';
//     tsc ts0, ts1;

    if (nrhs<7) { mexPrintf("Specify at least Phi1,A,Phi2,rhs,tol,trunc_norm,sol_prev\n"); return; }
    if (nrhs<8) prec=0; else { scal = mxGetPr(prhs[7]); prec = (char)round(scal[0]); }
    if (nrhs<9) nrestart=40; else { scal = mxGetPr(prhs[8]); nrestart = (long)round(scal[0]); }
    if (nrhs<10) niters=2; else { scal = mxGetPr(prhs[9]); niters = (long)round(scal[0]); }
    if (nrhs<11) verb=0; else { scal = mxGetPr(prhs[10]); verb = (char)round(scal[0]); }

    // Fetch the data
    dPhi1 = mxGetPr(prhs[0]);
    dimcount=mxGetNumberOfDimensions(prhs[0]);
    rhsdims = mxGetDimensions(prhs[0]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<3; i++) dims[i]=1;
    rx1 = dims[0];
    if (dims[1]!=rx1) { mexPrintf("Phi1 is not square!\n"); return; }
    ra1 = dims[2];

    dA = mxGetPr(prhs[1]);
    dimcount=mxGetNumberOfDimensions(prhs[1]);
    rhsdims = mxGetDimensions(prhs[1]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<4; i++) dims[i]=1;
    if (ra1 != dims[0]) { mexPrintf("ra1 in Phi1 and A are not consistent!\n"); return; }
    n = dims[1];
    if (n!=dims[2]) { mexPrintf("A is not square!\n"); return; }
    ra2 = dims[3];

    dPhi2 = mxGetPr(prhs[2]);
    dimcount=mxGetNumberOfDimensions(prhs[2]);
    rhsdims = mxGetDimensions(prhs[2]);
    for (i=0; i<dimcount; i++) dims[i]=rhsdims[i];
    for (i=dimcount; i<3; i++) dims[i]=1;
    rx2 = dims[0];
    if (dims[2]!=rx2) { mexPrintf("Phi2 is not square!\n"); return; }
    if (ra2 != dims[1]) { mexPrintf("ra2 in Phi2 and A are not consistent!\n"); return; }

    drhs = mxGetPr(prhs[3]);
    if ((mxGetM(prhs[3])*mxGetN(prhs[3]))!=(rx1*n*rx2)) { mexPrintf("RHS size is not consistent!\n"); return; }

    scal = mxGetPr(prhs[4]);
    tol = scal[0];

    scal = mxGetPr(prhs[5]);
    trunc_norm = (char)round(scal[0]);

    dsol_prev = mxGetPr(prhs[6]);
    if ((mxGetM(prhs[6])*mxGetN(prhs[6]))!=(rx1*n*rx2)) { mexPrintf("SOL_PREV size is not consistent!\n"); return; }

    if (nrestart>rx1*n*rx2) nrestart=rx1*n*rx2;

//     mexPrintf("Sizes: rx1: %d, n: %d, rx2: %d, ra1: %d, ra2: %d\n", rx1, n, rx2, ra1, ra2);

//     dsol = (double *)malloc(sizeof(double)*rx1*n*rx2);

    plhs[0] = mxCreateDoubleMatrix(rx1*n*rx2, 1, mxREAL);
    dsol = mxGetPr(plhs[0]);

//     plhs[1] = mxCreateDoubleMatrix(rx1*ra1, 1, mxREAL);

//     prec = 1;

    if (prec==1) {
        ptype = 'c';
        i = n*n*rx1*rx2;
        if (n*n*rx1*ra2>i) i = n*n*rx1*ra2;
    }
    if (prec==2) {
        ptype = 'l';
        i = rx1*n*rx1*n*rx2;
        if (rx1*rx1*n*n*ra2>i) i = rx1*rx1*n*n*ra2;
    }
    if (prec==3) {
        ptype = 'r';
        i = n*n*rx2*rx2*rx1;
        if (rx1*n*n*ra2>i) i = rx1*n*n*ra2;
    }
    djacs=NULL;
    if (prec>0) {
         djacs = (double *)malloc(sizeof(double)*i);
         djac_gen(&ptype, &rx1, &n, &rx2, &ra1, &ra2, dPhi1, dA, dPhi2, djacs);
    }
//     if (prec==1) {
//         ptype = 'c';
// 	i = n*n*rx1*rx2;
// 	if (n*n*rx1*ra2>i) i = n*n*rx1*ra2;
//         /*if (rx1*n*ra2*rx2>i) i = rx1*n*ra2*rx2;
//         if (rx1*ra1*n*rx2>i) i = rx1*ra1*n*rx2;
//         if (ra1*rx1>i) i=ra1*rx1;
//         if (ra2*rx2>i) i=ra2*rx2;
//         if (ra2*rx2*rx2>i) i=ra2*rx2*rx2;
//         work1 = (double *)malloc(sizeof(double)*i);
//         work2 = (double *)malloc(sizeof(double)*i);*/
//         djacs = (double *)malloc(sizeof(double)*i);
// //         ipiv = (long *)malloc(sizeof(long)*n);
//
// //         printf("djac_gen started...\n");
// 	djac_gen(&ptype, &rx1, &n, &rx2, &ra1, &ra2, dPhi1, dA, dPhi2, djacs);
// //         clock_gettime(CLOCK_REALTIME, &ts0);
// //         dcjacgen(dPhi1,dA,dPhi2, rx1, n, rx2, ra1, ra2, djacs);
// //         clock_gettime(CLOCK_REALTIME, &ts1);
// //         if (verb>0) mexPrintf("JacGen time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);
//
// //         dims[0]=n; dims[1]=n; dims[2]=rx1; dims[3]=rx2;
// //         plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
// //         dcjacapply(djacs, n, rx1, rx2, drhs, dsol);
// //         return;
// //         memcpy(dsol, djacs, sizeof(double)*n*n*rx1*rx2);
//     }
//     else {
    i = rx1*n*ra2*rx2;
    if (ra1>ra2) i = rx1*ra1*n*rx2;
    work1 = (double *)malloc(sizeof(double)*i);
    work2 = (double *)malloc(sizeof(double)*i);
//     }


//     verb = 2;

    if (trunc_norm==1) { // residual
        // prepare initial residual - for right prec
        dres = (double *)malloc(sizeof(double)*rx1*n*rx2);
	bfun3(&rx1,&n,&rx2, &rx1,&n,&rx2, &ra1, &ra2, dPhi1,dA,dPhi2, dsol_prev, dres, work1, work2);
//         bfun3(dPhi1,dA,dPhi2,rx1,n,rx2,ra1,ra2,dsol_prev,dres);
        dbeta = -1.0;
        i = rx1*n*rx2;
        daxpy_(&i,&dbeta,drhs,&ione,dres,&ione);
        dscal_(&i,&dbeta,dres,&ione);
        tol_prev = dnrm2_(&i, dres, &ione) / dnrm2_(&i, drhs, &ione);
        if (tol_prev<tol) tol_prev = tol;
//         mexPrintf("tol0: %3.5e, tol: %3.5e\n", tol_prev, tol);
//         clock_gettime(CLOCK_REALTIME, &ts0);
        dgmresr_hh(dPhi1, dA, dPhi2, dres, rx1, n, rx2, ra1, ra2, nrestart, tol/tol_prev, niters, ptype, djacs, dsol, verb);
	if (verb>0) mexPrintf("%s", gmres_3d_printf_buf_);
//         clock_gettime(CLOCK_REALTIME, &ts1);
//         if (verb>0) mexPrintf("gmres time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);

        if (prec>0) {
	    djac_apply(&ptype, &rx1, &n, &rx2, djacs, dsol, dsol, work1);
//             dcjacapply(djacs, n, rx1, rx2, dsol, dres);
//             dcopy_(&i,dres,&ione,dsol,&ione);
        }
        dbeta = 1.0;
        daxpy_(&i,&dbeta,dsol_prev,&ione,dsol,&ione);
        free(dres);
    }
    else { // fro
        // left-prec gmres until the new correction is less than tol
        i = rx1*n*rx2;
        dcopy_(&i, dsol_prev, &ione, dsol, &ione);
//         for (i=0; i<rx1*n*rx2; i++) mexPrintf("%g\n", dsol[i]);
//         clock_gettime(CLOCK_REALTIME, &ts0);
        dgmresl_hh(dPhi1, dA, dPhi2, drhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, ptype, djacs, dsol, verb);
	if (verb>0) mexPrintf("%s", gmres_3d_printf_buf_);
//         clock_gettime(CLOCK_REALTIME, &ts1);
//         if (verb>0) mexPrintf("gmres time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);

    }

//     i = rx1*ra1;
//     if (prec==1) { dcopy_(&i, scal, &ione, mxGetPr(plhs[1]), &ione);  /*free(scal);*/ }

    if (prec>0) { free(djacs); }

    free(work1); free(work2);
}
