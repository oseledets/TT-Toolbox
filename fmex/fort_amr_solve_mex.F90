#include "fintrf.h"

# define mwpointer integer(8)
# define mwPointer integer(8)
# define MWPOINTER INTEGER(8)
# define mwsize  mwpointer
# define mwSize  mwpointer
# define MWSIZE  MWPOINTER
# define mwindex mwpointer
# define mwIndex mwpointer
# define MWINDEX MWPOINTER
# define mwsignedindex mwpointer
# define mwSignedIndex mwpointer
# define MWSIGNEDINDEX MWPOINTER





      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
! RHS format: d, n, m, ry, ra, crA, crY, crX, rx, eps, params, prec
! PARAMS: nswp, kickrank, nrestart, niters, verb
      use tt_adapt_als

!     Declarations
      implicit none

!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer(4) nlhs, nrhs

!     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr
      integer(4) mxIsNumeric
      mwPointer mxCreateNumericMatrix
      integer(4) mxClassIDFromClassName
      integer(4) dumme_kuh
!       mwPointer mxGetM, mxGetN
!       real(8) dummy
!       mwPointer ptr1
      mwSize sz, sz2

      integer d
      character prec
!       integer,pointer :: n(:), r(:), ps(:)
      real(8), allocatable :: crA(:), crX(:), crY(:)
!       real(8) :: scal(*)
      integer,allocatable :: n(:),rx(:),ra(:),ry(:),m(:),params(:)
      real(8) eps

! subroutine tt_amr_solve(d,n,m,ry,ra,crA, crY, crX0, rx, eps, kickrank, nswp, verb, prec, nrestart, niters)

      sz = 1
!       d = 5
!       eps = 1e-6
      call mxCopyPtrToInteger8(mxGetPr(prhs(1)), d, sz)
!        call get_d(, d)
!        call get_eps(eps, %VAL(mxGetPr(prhs(6))))
!        print *, eps
!       d = %VAL(%VAL(mxGetPr(prhs(1))));

        allocate(rx(d+1), ry(d+1), ra(d+1), n(d), m(d), params(5))

        sz = d;
        call mxCopyPtrToInteger8(mxGetPr(prhs(2)), n, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(3)), m, sz)

        sz = d+1;
        call mxCopyPtrToInteger8(mxGetPr(prhs(4)), ry, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(5)), ra, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(9)), rx, sz)

!         print *,  n(1:d), ra(1:d+1), ry(1:d+1), rx(1:d+1)

        sz = sum(ra(2:d+1)*ra(1:d)*n(1:d)*m(1:d))
        allocate (crA(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(6)),crA,sz)
        sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
        allocate (crY(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(7)),crY,sz)


        sz = sum(rx(2:d+1)*rx(1:d)*m(1:d))
        allocate (crX(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(8)),crX,sz)

        sz = 1
        call mxCopyPtrToReal8(mxGetPr(prhs(10)), eps, sz)
        call mxCopyPtrToCharacter(mxGetPr(prhs(12)), prec, sz)
        sz = 5
        call mxCopyPtrToInteger8(mxGetPr(prhs(11)), params, sz)

!         print *, params, n(1:d), ra(1:d+1), ry(1:d+1), rx(1:d+1), prec

! PARAMS: nswp, kickrank, nrestart, niters, prec, verb

	call tt_amr_solve(d,n,m,ry,ra,crA, crY, crX, rx, eps, kickrank=params(2), nswp=params(1), verb=params(5), nrestart=params(3), niters=params(4), prec=prec)


! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10));


!     Create matrix for the return argument.
      sz = d+1; sz2=1
      dumme_kuh=0
      plhs(1) = mxCreateNumericMatrix(sz,sz2, mxClassIDFromClassName('int64'), dumme_kuh)

      sz = sum(rx(2:d+1)*rx(1:d)*n(1:d))
      plhs(2) = mxCreateDoubleMatrix(sz,sz2,dumme_kuh)

      sz = d+1
      call mxCopyInteger8ToPtr(rx, mxGetPr(plhs(1)), sz)
      sz = sum(rx(2:d+1)*rx(1:d)*n(1:d))
      call mxCopyReal8ToPtr(result_core, mxGetPr(plhs(2)), sz)

      deallocate(result_core)
      deallocate(ry, rx,ra, n, m, crA, crX, crY, params)

      return
      end

