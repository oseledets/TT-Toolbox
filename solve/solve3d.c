#include "mex.h"
#include "lapack.h"
#include "blas.h"
// #include "time.h"


//     typedef struct {
//                time_t   tv_sec;        // seconds
//                long     tv_nsec;       // nanoseconds
//            } tsc;
//
char trans = 'N';
char uplo = 'U';
double done=1.0;
double dzero=0.0;
long ione = 1;

void dperm213(double *in, long n1, long n2, long n3, double *out)
{
    long j,k;
    for (k=0; k<n3; k++) {
        for (j=0; j<n2; j++) {
            dcopy_(&n1, &in[j*n1+k*n1*n2], &ione, &out[j+k*n1*n2], &n2);
//             for (i=0; i<n1; i++) {
//                 out[j+i*n2+k*n1*n2] = in[i+j*n1+k*n1*n2];
//             }
        }
    }
}

void dperm312(double *in, long n1, long n2, long n3, double *out)
{
    long i,j; //,k;
    j = n1*n2;
    for (i=0; i<n3; i++) {
        dcopy_(&j, &in[i*j], &ione, &out[i], &n3);
    }
//     for (k=0; k<n3; k++) {
//         for (j=0; j<n2; j++) {
//             for (i=0; i<n1; i++) {
//                 out[k+i*n3+j*n3*n1] = in[i+j*n1+k*n1*n2];
//             }
//         }
//     }
}

void dcjacgen(double *Phi1,double *A,double *Phi2, long rx1, long n, long rx2, long ra1, long ra2,double *jacs)
{
    long i, k,m,cnt;
    double *diag1;
    double *diag2;
    double *Id, *B, *A2;
    long *ipiv, info;

    diag1 = (double *)malloc(sizeof(double)*rx1*ra1);
    diag2 = (double *)malloc(sizeof(double)*rx2*ra2);

    for (k=0; k<ra1; k++) {
        for (m=0; m<rx1; m++) {
            diag1[m+k*rx1] = Phi1[m+m*rx1+k*rx1*rx1]; // rx1,ra1
        }
    }
    for (k=0; k<ra2; k++) {
        for (m=0; m<rx2; m++) {
            diag2[k+m*ra2] = Phi2[m+m*rx2+k*rx2*rx2]; // ra2,rx2
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
    dperm213(jacs,rx1, n*n, rx2, A2);

    Id = (double *)malloc(sizeof(double)*n*n);
    for (k=0; k<n; k++) {
        for (m=0; m<n; m++) {
            if (m==k) Id[m+k*n]=1.0;
            else      Id[m+k*n]=0.0;
        }
    }
    B = (double *)malloc(sizeof(double)*n*n);
    ipiv = (long*)malloc(sizeof(long)*n);

    for (k=0; k<rx2; k++) {
        for (m=0; m<rx1; m++) {
            i = n*n;
            dcopy_(&i,Id,&ione,B,&ione);
            dgesv_(&n, &n, &(A2[m*n*n+k*n*n*rx1]), &n, ipiv, B, &n, &info);
            if (info>0) mexPrintf("The component %d at ranks [%d,%d] is singular\n", info,m,k);
            dcopy_(&i,B,&ione,&(jacs[m*n*n+k*n*n*rx1]),&ione);
        }
    }

    free(A2);
    free(diag1);
    free(diag2);
    free(Id);
    free(B);
    free(ipiv);
}


void dcjacapply(double *jacs, long n, long rx1, long rx2, double *x, double *y)
// x will be broken, &x must be != &y!
{
    long i,j;
    dperm213(x, rx1, n, rx2, y);
    for (j=0; j<rx2; j++) {
        for (i=0; i<rx1; i++) {
            dgemv_(&trans,&n,&n,&done,&jacs[i*n*n+j*n*n*rx1],&n,&y[i*n+j*n*rx1],&ione, &dzero, &x[i*n+j*n*rx1], &ione);
        }
    }
    dperm213(x, n, rx1, rx2, y);
}


void bfun3(double *Phi1, double *A, double *Phi2, long rx1, long n, long rx2, long ra1, long ra2, double *x, double *y)
// x is preserved, &x==&y is possible
{
    // we need A in form ra1,n,n',ra2
    long i,j,sz1, sz2, sz3;
    double *z1, *z2;

// bydlocoded max. size calculus
    sz1 = rx1*n*rx2;
    sz2 = rx2*rx1*n*ra2;
    if (sz2>sz1) sz1 = sz2;
    sz2 = ra1*n*rx1*rx2;
    if (sz2>sz1) sz1 = sz2;

    z1 = (double *)malloc(sizeof(double)*sz1);
    z2 = (double *)malloc(sizeof(double)*sz1);

//     z1 = (double *)malloc(sizeof(double)*rx1*n*rx2);
    dperm213(x, rx1*n, rx2, 1, z1); // rx2', rx1', n

    // Phi2*x
//     z2 = (double *)malloc(sizeof(double)*rx2*rx1*n*ra2);
    for (j=0; j<ra2; j++) {
        for (i=0; i<rx1*n; i++) {
            dgemv_(&trans,&rx2,&rx2,&done,&Phi2[j*rx2*rx2],&rx2,&z1[i*rx2],&ione, &dzero, &z2[i*rx2+j*rx2*rx1*n], &ione); // indices rx2, rx1', n', ra2
        }
    }
//     free(z1);
//     z1 = (double *)malloc(sizeof(double)*n*ra2*rx2*rx1);
    dperm213(z2, rx2*rx1, n*ra2, 1, z1); // now n'*ra2, rx2*rx1'
//     free(z2);

    // A*x
//     z2 = (double *)malloc(sizeof(double)*ra1*n*rx1*rx2);
    sz1 = ra1*n; sz2 = n*ra2; sz3 = rx2*rx1;
    dgemm_(&trans, &trans, &sz1, &sz3, &sz2, &done, A, &sz1, z1, &sz2, &dzero, z2, &sz1);
/*    for (i=0; i<rx2*rx1; i++) {
        dgemv_(&trans,&sz1,&sz2,&done,A,&sz1,&z1[i*n*ra2],&ione, &dzero, &z2[i*ra1*n], &ione);
        // indices ra1,n,rx2,rx1'
    }*/
//     free(z1);
//     z1 = (double *)malloc(sizeof(double)*rx1*ra1*n*rx2);
    dperm312(z2, ra1, n*rx2, rx1, z1); // now rx1'*ra1*n*rx2
//     free(z2);

    // Phi1*x
    sz1 = rx1; sz2 = rx1*ra1; sz3 = n*rx2;
    dgemm_(&trans, &trans, &sz1, &sz3, &sz2, &done, Phi1, &sz1, z1, &sz2, &dzero, y, &sz1);
//     for (i=0; i<n*rx2; i++) {
//         dgemv_(&trans,&sz1,&sz2,&done,Phi1,&sz1,&z1[i*rx1*ra1],&ione, &dzero, &y[i*rx1], &ione); // indices rx1, n, rx2
//     }

    free(z1);
    free(z2);
}


void dgmresr(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, int nrestart, double tol, int niters, double *jacs, double *sol, char verb)
// right preconditioned - for residual tolerance
{
    long i,j,it, sz, sz2, lda, lwork;
    double *v, *w, *w2;
    double *H,*Q,*R, *work;
    double nrmr, curres, nrmrhs;
    double dalpha, dbeta;
    char last_iter = 0;

    sz = rx1*n*rx2;

    v = (double *)malloc(sizeof(double)*sz*nrestart);
    w = (double *)malloc(sizeof(double)*sz);
    w2 = (double *)malloc(sizeof(double)*sz);
    H = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart+1)*(nrestart+1));

    lda = nrestart+1;

    for (i=0; i<nrestart*lda; i++) H[i]=0.0;

    lwork = 3*(nrestart+1);
    work = (double *)malloc(sizeof(double)*lwork);

    for (j=0; j<sz; j++) sol[j]=0.0;

    for (it=0; it<niters; it++) {
        // r0
        if (jacs!=NULL) {
            dcopy_(&sz,sol,&ione,w,&ione);
            dcjacapply(jacs, n, rx1, rx2, w, w2);
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w2,&v[0]);
        }
        else {
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,&v[0]);
        }
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,&v[0],&ione);
        nrmr = dnrm2_(&sz,&v[0],&ione);
        if (verb>1) mexPrintf("restart %d, res: %3.5e\n", it, nrmr);
        if (it==0) nrmrhs = nrmr;
        dbeta = -1.0/nrmr;
        dscal_(&sz,&dbeta,&v[0],&ione);

        for (j=0; j<nrestart; j++) {
            // precvec, matvec
            if (jacs!=NULL) {
                dcopy_(&sz,&v[j*sz],&ione,w,&ione);
                dcjacapply(jacs, n, rx1, rx2, w, w2);
                bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w2,w);
            }
            else {
                bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,&v[j*sz],w);
//                 dgemv_(&trans,&sz,&sz,&dalpha,A,&sz,w,&ione, &dbeta, w2, &ione);
//                 dcopy_(&sz,w2,&ione,w,&ione);
            }

            // Arnoldi
            for (i=0; i<=j; i++) {
                H[i+j*lda] = ddot_(&sz,&v[i*sz],&ione,w,&ione);
                dbeta = -H[i+j*lda];
//                 mexPrintf("iter [%d,%d], v[%d]'*w: %3.5e\n", it, j, i, dbeta);
                daxpy_(&sz,&dbeta,&v[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
//             mexPrintf("iter [%d,%d], v[%d]'*w: %3.5e\n", it, j, j+1, dbeta);
            H[j+1+j*lda] = dbeta;
            if (j<nrestart-1) {
                dbeta = 1.0/dbeta;
                dscal_(&sz,&dbeta,w,&ione);
                dcopy_(&sz,w,&ione,&v[(j+1)*sz],&ione);
            }

            // QR
            /*for (sz=0; sz<j+2; sz++) {
                for (i=0; i<j+1; i++)  mexPrintf("%3.10e\t", H[sz+lda*i]);
                mexPrintf("\n");
            }*/
            sz = nrestart*lda;
            dcopy_(&sz,H,&ione,Q,&ione);
            sz = j+2; sz2 = j+1;
            dgeqrf_(&sz, &sz2, Q, &lda, &R[nrestart*lda], work, &lwork, &i);
            for (sz=0; sz<=j; sz++) {
                for (i=0; i<=sz; i++) {
                    R[i+lda*sz]=Q[i+lda*sz];
                }
            }
            sz = j+2; sz2 = j+1;
            dorgqr_(&sz, &sz2, &sz2, Q, &lda, &R[nrestart*lda], work, &lwork, &i);
//             for (sz=0; sz<j+2; sz++) {
//                 for (i=0; i<j+1; i++)  mexPrintf("%3.10e\t", Q[sz+lda*i]);
//                 mexPrintf("\n");
//             }
//             for (sz=0; sz<j+2; sz++) {
//                 for (i=0; i<j+1; i++) mexPrintf("%3.10e\t", R[sz+lda*i]);
//                 mexPrintf("\n");
//             }
//             sz = j+2;
            // solution
            for (i=0; i<=j; i++) Q[i] = Q[i*lda]*nrmr;
            dtrsv_(&uplo,&trans,&trans,&sz2,R,&lda,Q,&ione);

            // residual
            dalpha = -1.0; dbeta = 1.0;
            for (i=0; i<=(j+1); i++) R[i] = 0.0;
            R[0]=nrmr;
            dgemv_(&trans,&sz,&sz2,&dalpha,H,&lda,Q,&ione, &dbeta, R, &ione);
            curres = dnrm2_(&sz2,R,&ione);
            curres = curres/nrmrhs;
            if (verb>1) mexPrintf("iter [%d,%d], res: %3.5e\n", it, j, curres);

            sz = rx1*n*rx2;

            if (last_iter==1) break;
            if ((curres<tol)) last_iter=1;
        }
        if (j==nrestart) j=nrestart-1;
        // Correction
        for (i=0; i<=j; i++) {
//             if (jacs!=NULL) dcjacapply(jacs, n, rx1, rx2, &v[i*sz], w);
//             else            dcopy_(&sz,&v[i*sz],&ione,w,&ione);
//             daxpy_(&sz,&Q[i],w,&ione,sol,&ione);
            daxpy_(&sz,&Q[i],&v[i*sz],&ione,sol,&ione);
        }
        if (last_iter==1) break;
        if ((curres<tol)) last_iter=1;
    }
    if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(v);
    free(w);
    free(w2);
    free(H);
    free(Q);
    free(R);
    free(work);
}


void dgmresr_hh(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, long nrestart, double tol, long niters, double *jacs, double *sol, char verb)
// right preconditioned - for residual tolerance
// This one with Householder tranforms
{
    long i,j,it, sz;
    double *U, *w, *w2;
    double *R, *J, *tau;
    double nrmr, curres, nrmrhs;
    double dalpha, dbeta;
    char last_iter = 0;

    sz = rx1*n*rx2;

    U = (double *)malloc(sizeof(double)*sz*(nrestart+1));
    w = (double *)malloc(sizeof(double)*sz);
    w2 = (double *)malloc(sizeof(double)*sz);
    tau = (double *)malloc(sizeof(double)*(nrestart+1));
    J = (double *)malloc(sizeof(double)*2*(nrestart));
//     Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart)*(nrestart));

//     lda = nrestart+1;


//     lwork = 3*(nrestart+1);
//     work = (double *)malloc(sizeof(double)*lwork);

    for (j=0; j<sz; j++) sol[j]=0.0;

    for (it=0; it<niters; it++) {
        // r0
        if (jacs!=NULL) {
            dcopy_(&sz,sol,&ione,w,&ione);
            dcjacapply(jacs, n, rx1, rx2, w, w2);
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w2,w);
        }
        else {
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,w);
        }
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,w,&ione);
        dscal_(&sz,&dbeta,w,&ione);
        nrmr = dnrm2_(&sz,w,&ione);
        if (verb>1) mexPrintf("restart %d, res: %3.5e\n", it, nrmr);
        if (it==0) nrmrhs = nrmr;
        // initial HHT
        dbeta = nrmr;
        if (w[0]<0.0) dbeta = -dbeta;
        w[0]+=dbeta;
        tau[0]=-dbeta;
        nrmr = dnrm2_(&sz,w,&ione);
        dbeta = 1.0/nrmr;
        dscal_(&sz,&dbeta,w,&ione);
        dcopy_(&sz,w,&ione,&U[0],&ione);

        for (j=0; j<nrestart; j++) {
            // HHT on last U
            dcopy_(&sz,&U[sz*j],&ione,w,&ione);
            dbeta = -2.0*U[j+sz*j];
            dscal_(&sz,&dbeta,w,&ione);
            w[j]+=1.0;
            for (i=j-1; i>=0; i--) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
            dbeta = 1.0/dbeta;
            dscal_(&sz,&dbeta,w,&ione); // w=w/norm(w);

//             mexPrintf("norm_v0=%3.7e\n", 1.0/dbeta);

            // precvec, matvec
            if (jacs!=NULL) {
                dcopy_(&sz,w,&ione,w2,&ione);
                dcjacapply(jacs, n, rx1, rx2, w2, w);
                bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            }
            else {
                bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            }

            // Orthog w to j projectors
            for (i=0; i<=j; i++) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }

            // new P_{j+1}
            if (j<sz-1) {
                for (i=0; i<=j; i++) U[i+(j+1)*sz]=0.0;
                i = sz-j-1;
                dcopy_(&i, &w[j+1], &ione, &U[j+1+(j+1)*sz], &ione);
                dalpha = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                if (dalpha!=0.0) {
                    if (w[j+1]<0.0) dalpha = -dalpha;
                    U[j+1+(j+1)*sz]+=dalpha;
                    dbeta = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                    dbeta = 1.0/dbeta;
                    dscal_(&i,&dbeta,&U[j+1+(j+1)*sz],&ione);

                    w[j+1]=-dalpha;
                    for (i=j+2; i<sz; i++) w[i]=0.0;
                }
            }

            // Givens rotators to the top of w
            for (i=0; i<=j-1; i++) {
                dbeta = w[i];
                w[i] = J[0+i*2]*w[i] + J[1+i*2]*w[i+1];
                w[i+1] = -J[1+i*2]*dbeta + J[0+i*2]*w[i+1];
            }

//            mexPrintf("w[%d]=%3.7e\tw[%d]=%3.7e\n", i, w[j], j+1, w[j+1]);

            // New rotator
            if (j<sz-1) {
                dalpha = sqrt((w[j]*w[j])+(w[j+1]*w[j+1]));
//                 mexPrintf("rho=%3.7e\n", dalpha);
                J[0+j*2] = w[j]/dalpha;
                J[1+j*2] = w[j+1]/dalpha;
//                 mexPrintf("Jnew=%3.7e\t%3.7e\n", J[j*2], J[j*2+1]);
                tau[j+1] = -J[1+j*2]*tau[j];
                tau[j] = J[0+j*2]*tau[j];
//                 mexPrintf("tau=%3.7e\t%3.7e\n", tau[j], tau[j+1]);
                w[j] = dalpha;
                w[j+1] = 0.0;
            }

//             for (i=0; i<=j+1; i++) mexPrintf("tau[%d]=%3.7e\n", i, tau[i]);

            dcopy_(&nrestart, w, &ione, &R[j*nrestart], &ione);


            // residual
            curres = fabs(tau[j+1])/nrmrhs;
            if (verb>1) mexPrintf("iter [%d,%d], res: %3.5e\n", it, j, curres);

            if ((curres<tol)) break;
        }

//         for (i=0; i<=j; i++) mexPrintf("%g\t", tau[i]);
//         mexPrintf("\n");

        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;
        dtrsv_(&uplo,&trans,&trans,&i,R,&nrestart,tau,&ione);

/*        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;*/

        // Correction
        dcopy_(&sz, &U[j*sz], &ione, w, &ione);
        dbeta = -2.0*U[j+j*sz]*tau[j];
        dscal_(&sz, &dbeta, w, &ione);
        w[j]+=tau[j];
        for (i=j-1; i>=0; i--) {
            w[i]+=tau[i];
            dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
            daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
        }
        dalpha=1.0;
//         dcopy_(&sz, w, &ione, sol, &ione);
        daxpy_(&sz,&dalpha,w,&ione,sol,&ione);
        if ((curres<tol)) break;
    }
    if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(U);
    free(w);
    free(w2);
    free(tau);
    free(J);
    free(R);
}


void dgmresl(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, int nrestart, double tol, int niters, double *jacs, double *sol, char verb)
// left preconditioned - for fro tolerance
{
    long i,j,it, sz, sz2, lda, lwork;
    double *v, *w, *w2;
    double *H,*Q,*R, *work;
    double nrmr, curres, nrmrhs, nrmsol;
    double dalpha, dbeta;

    sz = rx1*n*rx2;

    v = (double *)malloc(sizeof(double)*sz*nrestart);
    w = (double *)malloc(sizeof(double)*sz);
    w2 = (double *)malloc(sizeof(double)*sz);
    H = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart+1)*(nrestart+1));

    lda = nrestart+1;

    for (i=0; i<nrestart*lda; i++) H[i]=0.0;

    lwork = 3*(nrestart+1);
    work = (double *)malloc(sizeof(double)*lwork);

    nrmsol = dnrm2_(&sz, sol, &ione);
    if (jacs!=NULL) {
        dcopy_(&sz,rhs,&ione,w,&ione);
        dcjacapply(jacs, n, rx1, rx2, w, w2);
        nrmrhs = dnrm2_(&sz, w2, &ione);
    } else
        nrmrhs = dnrm2_(&sz, rhs, &ione);
//     for (j=0; j<sz; j++) sol[j]=0.0;

    for (it=0; it<niters; it++) {
        // r0
        bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,&v[0]);
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,&v[0],&ione);
        if (jacs!=NULL) {
            dcjacapply(jacs, n, rx1, rx2, &v[0], w);
            dcopy_(&sz,w,&ione,&v[0],&ione);
        }
        nrmr = dnrm2_(&sz,&v[0],&ione);
        if (verb>1) mexPrintf("restart %d, res: %3.5e\n", it, nrmr);
        dbeta = -1.0/nrmr;
        dscal_(&sz,&dbeta,&v[0],&ione);

        for (j=0; j<nrestart; j++) {
            // matvec, precvec
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,&v[j*sz],w);
            if (jacs!=NULL) {
                dcjacapply(jacs, n, rx1, rx2, w, w2);
                dcopy_(&sz,w2,&ione,w,&ione);
            }

            // Arnoldi
            for (i=0; i<=j; i++) {
                H[i+j*lda] = ddot_(&sz,&v[i*sz],&ione,w,&ione);
                dbeta = -H[i+j*lda];
//                 mexPrintf("iter [%d,%d], v[%d]'*w: %3.5e\n", it, j, i, dbeta);
                daxpy_(&sz,&dbeta,&v[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
//             mexPrintf("iter [%d,%d], v[%d]'*w: %3.5e\n", it, j, j+1, dbeta);
            H[j+1+j*lda] = dbeta;
            if (j<nrestart-1) {
                dbeta = 1.0/dbeta;
                dscal_(&sz,&dbeta,w,&ione);
                dcopy_(&sz,w,&ione,&v[(j+1)*sz],&ione);
            }

            // QR
            /*for (sz=0; sz<j+2; sz++) {
                for (i=0; i<j+1; i++)  mexPrintf("%3.10e\t", H[sz+lda*i]);
                mexPrintf("\n");
            }*/
            sz = nrestart*lda;
            dcopy_(&sz,H,&ione,Q,&ione);
            sz = j+2; sz2 = j+1;
            dgeqrf_(&sz, &sz2, Q, &lda, &R[nrestart*lda], work, &lwork, &i);
            for (sz=0; sz<=j; sz++) {
                for (i=0; i<=sz; i++) {
                    R[i+lda*sz]=Q[i+lda*sz];
                }
            }
            sz = j+2; sz2 = j+1;
            dorgqr_(&sz, &sz2, &sz2, Q, &lda, &R[nrestart*lda], work, &lwork, &i);
//             for (sz=0; sz<j+2; sz++) {
//                 for (i=0; i<j+1; i++)  mexPrintf("%3.10e\t", Q[sz+lda*i]);
//                 mexPrintf("\n");
//             }
//             for (sz=0; sz<j+2; sz++) {
//                 for (i=0; i<j+1; i++) mexPrintf("%3.10e\t", R[sz+lda*i]);
//                 mexPrintf("\n");
//             }
//             sz = j+2;
            // solution
            for (i=0; i<=j; i++) Q[i] = Q[i*lda]*nrmr;
            dtrsv_(&uplo,&trans,&trans,&sz2,R,&lda,Q,&ione);

            // residual
            dalpha = -1.0; dbeta = 1.0;
            for (i=0; i<=(j+1); i++) R[i] = 0.0;
            R[0]=nrmr;
            dgemv_(&trans,&sz,&sz2,&dalpha,H,&lda,Q,&ione, &dbeta, R, &ione);
            curres = dnrm2_(&sz2,R,&ione);
            curres = curres/nrmrhs;
            if (verb>1) mexPrintf("iter [%d,%d], res: %3.5e\n", it, j, curres);

            sz = rx1*n*rx2;

            if ((curres<tol)&&(fabs(Q[j])/nrmsol<tol)) break;
        }
        if (j==nrestart) j=nrestart-1;
        // Correction
        for (i=0; i<=j; i++) {
/*            if (jacs!=NULL) dcjacapply(jacs, n, rx1, rx2, &v[i*sz], w);
            else            dcopy_(&sz,&v[i*sz],&ione,w,&ione);
            daxpy_(&sz,&Q[i],w,&ione,sol,&ione);*/
            daxpy_(&sz,&Q[i],&v[i*sz],&ione,sol,&ione);
        }
//         j=j+1;
        if ((curres<tol)&&(fabs(Q[j])/nrmsol<tol)) break;
//         if ((curres<tol)&&(dnrm2_(&sz,Q,&ione)/nrmsol<tol)) break;
//         j=j-1;
    }
    if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(v);
    free(w);
    free(w2);
    free(H);
    free(Q);
    free(R);
    free(work);
}


void dgmresl_hh(double *Phi1, double *A, double *Phi2, double *rhs, long rx1, long n, long rx2, long ra1, long ra2, long nrestart, double tol, long niters, double *jacs, double *sol, char verb)
// left preconditioned - for fro tolerance
// This one with Householder tranforms
{
    long i,j,it, sz;
    double *U, *w, *w2;
    double *R, *J, *tau;
    double nrmr, curres, nrmrhs, nrmsol;
    double dalpha, dbeta;
    char last_iter = 0;

    sz = rx1*n*rx2;

    U = (double *)malloc(sizeof(double)*sz*(nrestart+1));
    w = (double *)malloc(sizeof(double)*sz);
    w2 = (double *)malloc(sizeof(double)*sz);
    tau = (double *)malloc(sizeof(double)*(nrestart+1));
    J = (double *)malloc(sizeof(double)*2*(nrestart));
//     Q = (double *)malloc(sizeof(double)*nrestart*(nrestart+1));
    R = (double *)malloc(sizeof(double)*(nrestart)*(nrestart));

    if (jacs!=NULL) {
        dcopy_(&sz,rhs,&ione,w,&ione);
        dcjacapply(jacs, n, rx1, rx2, w, w2);
        nrmrhs = dnrm2_(&sz, w2, &ione);
    } else
        nrmrhs = dnrm2_(&sz, rhs, &ione);

    for (it=0; it<niters; it++) {
        // r0
        nrmsol = dnrm2_(&sz, sol, &ione);
        bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,sol,w);
        dbeta = -1.0;
        daxpy_(&sz,&dbeta,rhs,&ione,w,&ione);
        dscal_(&sz,&dbeta,w,&ione);
        if (jacs!=NULL) {
            dcopy_(&sz,w,&ione,w2,&ione);
            dcjacapply(jacs, n, rx1, rx2, w2, w);
        }
        nrmr = dnrm2_(&sz,w,&ione);
        if (verb>1) mexPrintf("restart %d, res: %3.5e\n", it, nrmr);
        if (it==0) nrmrhs = nrmr;
        // initial HHT
        dbeta = nrmr;
        if (w[0]<0.0) dbeta = -dbeta;
        w[0]+=dbeta;
        tau[0]=-dbeta;
        nrmr = dnrm2_(&sz,w,&ione);
        dbeta = 1.0/nrmr;
        dscal_(&sz,&dbeta,w,&ione);
        dcopy_(&sz,w,&ione,&U[0],&ione);

        for (j=0; j<nrestart; j++) {
            // HHT on last U
            dcopy_(&sz,&U[sz*j],&ione,w,&ione);
            dbeta = -2.0*U[j+sz*j];
            dscal_(&sz,&dbeta,w,&ione);
            w[j]+=1.0;
            for (i=j-1; i>=0; i--) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }
            dbeta = dnrm2_(&sz,w,&ione);
            dbeta = 1.0/dbeta;
            dscal_(&sz,&dbeta,w,&ione); // w=w/norm(w);

//             mexPrintf("norm_v0=%3.7e\n", 1.0/dbeta);

            // precvec, matvec
            bfun3(Phi1,A,Phi2,rx1,n,rx2,ra1,ra2,w,w);
            if (jacs!=NULL) {
                dcopy_(&sz,w,&ione,w2,&ione);
                dcjacapply(jacs, n, rx1, rx2, w2, w);
            }

            // Orthog w to j projectors
            for (i=0; i<=j; i++) {
                dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
                daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
            }

            // new P_{j+1}
            if (j<sz-1) {
                for (i=0; i<=j; i++) U[i+(j+1)*sz]=0.0;
                i = sz-j-1;
                dcopy_(&i, &w[j+1], &ione, &U[j+1+(j+1)*sz], &ione);
                dalpha = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                if (dalpha!=0.0) {
                    if (w[j+1]<0.0) dalpha = -dalpha;
                    U[j+1+(j+1)*sz]+=dalpha;
                    dbeta = dnrm2_(&i, &U[j+1+(j+1)*sz], &ione);
                    dbeta = 1.0/dbeta;
                    dscal_(&i,&dbeta,&U[j+1+(j+1)*sz],&ione);

                    w[j+1]=-dalpha;
                    for (i=j+2; i<sz; i++) w[i]=0.0;
                }
            }

            // Givens rotators to the top of w
            for (i=0; i<=j-1; i++) {
                dbeta = w[i];
                w[i] = J[0+i*2]*w[i] + J[1+i*2]*w[i+1];
                w[i+1] = -J[1+i*2]*dbeta + J[0+i*2]*w[i+1];
            }

//            mexPrintf("w[%d]=%3.7e\tw[%d]=%3.7e\n", i, w[j], j+1, w[j+1]);

            // New rotator
            if (j<sz-1) {
                dalpha = sqrt((w[j]*w[j])+(w[j+1]*w[j+1]));
//                 mexPrintf("rho=%3.7e\n", dalpha);
                J[0+j*2] = w[j]/dalpha;
                J[1+j*2] = w[j+1]/dalpha;
//                 mexPrintf("Jnew=%3.7e\t%3.7e\n", J[j*2], J[j*2+1]);
                tau[j+1] = -J[1+j*2]*tau[j];
                tau[j] = J[0+j*2]*tau[j];
//                 mexPrintf("tau=%3.7e\t%3.7e\n", tau[j], tau[j+1]);
                w[j] = dalpha;
                w[j+1] = 0.0;
            }

//             for (i=0; i<=j+1; i++) mexPrintf("tau[%d]=%3.7e\n", i, tau[i]);

            dcopy_(&nrestart, w, &ione, &R[j*nrestart], &ione);


            // residual
            curres = fabs(tau[j+1])/nrmrhs;
            if (verb>1) mexPrintf("iter [%d,%d], res: %3.5e\n", it, j, curres);

            if ((curres<tol)) break;
        }

//         for (i=0; i<=j; i++) mexPrintf("%g\t", tau[i]);
//         mexPrintf("\n");

        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;
        dtrsv_(&uplo,&trans,&trans,&i,R,&nrestart,tau,&ione);

/*        if (j==nrestart) {
            j=nrestart-1;
            i=nrestart;
        }
        else i=j+1;*/

        // Correction
        dcopy_(&sz, &U[j*sz], &ione, w, &ione);
        dbeta = -2.0*U[j+j*sz]*tau[j];
        dscal_(&sz, &dbeta, w, &ione);
        w[j]+=tau[j];
        for (i=j-1; i>=0; i--) {
            w[i]+=tau[i];
            dbeta = -2.0*ddot_(&sz,&U[i*sz],&ione,w,&ione);
            daxpy_(&sz,&dbeta,&U[i*sz],&ione,w,&ione);
        }
        dalpha=1.0;
//         dcopy_(&sz, w, &ione, sol, &ione);
        daxpy_(&sz,&dalpha,w,&ione,sol,&ione);
        if ((curres<tol)&&(fabs(tau[j])/nrmsol<tol)) break;
    }
    if (verb>0) mexPrintf("gmres conducted [%d,%d] iters to relres %3.3e\n", it, j, curres);

    free(U);
    free(w);
    free(w2);
    free(tau);
    free(J);
    free(R);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
// Phi1 [r1,r1',ra1], A[ra1,n,n',ra2], Phi2[r2,r2',ra2], rhs, tol, trunc_norm, sol_prev, prec, nrestart, niters, verb
{
    double *dPhi1, *dA, *dPhi2, *drhs, *scal, *dsol_prev, *dres;
    double tol, tol_prev;
    double *dsol, *djacs;
    double dbeta;
    long nrestart, niters;
    long dimcount, *rhsdims;
    long dims[4];
    long rx1, rx2, ra1, ra2, n, i;
    char prec, verb, trunc_norm;
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
    if (dims[1]!=rx2) { mexPrintf("Phi2 is not square!\n"); return; }
    if (ra2 != dims[2]) { mexPrintf("ra2 in Phi2 and A are not consistent!\n"); return; }

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


    djacs=NULL;
    if (prec==1) {
        djacs = (double *)malloc(sizeof(double)*n*n*rx1*rx2);
//         clock_gettime(CLOCK_REALTIME, &ts0);
        dcjacgen(dPhi1,dA,dPhi2, rx1, n, rx2, ra1, ra2, djacs);
//         clock_gettime(CLOCK_REALTIME, &ts1);
//         if (verb>0) mexPrintf("JacGen time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);

//         dims[0]=n; dims[1]=n; dims[2]=rx1; dims[3]=rx2;
//         plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
//         dcjacapply(djacs, n, rx1, rx2, drhs, dsol);
//         return;
//         memcpy(dsol, djacs, sizeof(double)*n*n*rx1*rx2);
    }

    if (trunc_norm==1) { // residual
        // prepare initial residual - for right prec
        dres = (double *)malloc(sizeof(double)*rx1*n*rx2);
        bfun3(dPhi1,dA,dPhi2,rx1,n,rx2,ra1,ra2,dsol_prev,dres);
        dbeta = -1.0;
        i = rx1*n*rx2;
        daxpy_(&i,&dbeta,drhs,&ione,dres,&ione);
        dscal_(&i,&dbeta,dres,&ione);
        tol_prev = dnrm2_(&i, dres, &ione) / dnrm2_(&i, drhs, &ione);
        if (tol_prev<tol) tol_prev = tol;
//         mexPrintf("tol0: %3.5e, tol: %3.5e\n", tol_prev, tol);
//         clock_gettime(CLOCK_REALTIME, &ts0);
        dgmresr_hh(dPhi1, dA, dPhi2, dres, rx1, n, rx2, ra1, ra2, nrestart, tol/tol_prev, niters, djacs, dsol, verb);
//         clock_gettime(CLOCK_REALTIME, &ts1);
//         if (verb>0) mexPrintf("gmres time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);

        if (prec==1) {
            dcjacapply(djacs, n, rx1, rx2, dsol, dres);
            dcopy_(&i,dres,&ione,dsol,&ione);
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
        dgmresl_hh(dPhi1, dA, dPhi2, drhs, rx1, n, rx2, ra1, ra2, nrestart, tol, niters, djacs, dsol, verb);
//         clock_gettime(CLOCK_REALTIME, &ts1);
//         if (verb>0) mexPrintf("gmres time: %g\n", difftime(ts1.tv_sec, ts0.tv_sec) + ((double)(ts1.tv_nsec-ts0.tv_nsec))*1e-9);

    }

//     i = rx1*ra1;
//     if (prec==1) { dcopy_(&i, scal, &ione, mxGetPr(plhs[1]), &ione);  /*free(scal);*/ }

    if (prec==1) free(djacs);

}
