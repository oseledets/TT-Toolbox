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


!       subroutine copy_int(pt,n,sz)
!       integer :: pt(*)
!       integer,allocatable :: n(:)
!       integer :: sz
!       integer :: i
!
!       allocate(n(sz))
!       do i=1,sz
!         n(i)=pt(i)
!       end do
!       end
!
!       subroutine copy_real(pt,n,sz)
!       real(8) :: pt(*)
!       real(8),allocatable :: n(:)
!       integer :: sz
!       integer :: i
!
!       allocate(n(sz))
!       do i=1,sz
!         n(i)=pt(i)
!       end do
!       end
!
!       subroutine get_int(pt,d)
!       integer :: pt(*)
!       integer :: d
!
! !       d = pt(1)
!       d = pt(1)
!       end
!
!       subroutine get_real(pt,eps)
!       real(8) :: eps
!       real(8) :: pt(*)
!
!       eps = pt(1)
!       !       d = int(pt(1))
!       end
!


      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

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

      integer d, rmax, nswp, kickrank, verb
!       integer,pointer :: n(:), r(:), ps(:)
      real(8), allocatable :: crA(:), crX(:), crY(:)
!       real(8) :: scal(*)
      integer,allocatable :: n(:),rx(:),ra(:),ry(:),m(:)
      real(8) eps

! subroutine tt_mvk4(d,n,m,rx,ra,crA, crX, crY0, ry, eps)



      sz = 1
!       d = 5
!       eps = 1e-6
      call mxCopyPtrToInteger8(mxGetPr(prhs(1)), d, sz)
!        call get_d(, d)
!        print *, d
!        call get_eps(eps, %VAL(mxGetPr(prhs(6))))
!        print *, eps
!       d = %VAL(%VAL(mxGetPr(prhs(1))));

        allocate(rx(d+1), ry(d+1), ra(d+1), n(d), m(d))

        sz = d;
        call mxCopyPtrToInteger8(mxGetPr(prhs(2)), n, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(3)), m, sz)

        sz = d+1;
        call mxCopyPtrToInteger8(mxGetPr(prhs(4)), rx, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(5)), ra, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(9)), ry, sz)

        sz = sum(ra(2:d+1)*ra(1:d)*n(1:d)*m(1:d))
        allocate (crA(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(6)),crA,sz)
        sz = sum(rx(2:d+1)*rx(1:d)*m(1:d))
        allocate (crX(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(7)),crX,sz)


        sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
        allocate (crY(sz))
        call mxCopyPtrToReal8(mxGetPr(prhs(8)),crY,sz)

        sz = 1
        call mxCopyPtrToReal8(mxGetPr(prhs(10)), eps, sz)
        call mxCopyPtrToInteger8(mxGetPr(prhs(11)), rmax, sz)
	call mxCopyPtrToInteger8(mxGetPr(prhs(12)), nswp, sz)
	call mxCopyPtrToInteger8(mxGetPr(prhs(13)), kickrank, sz)
	call mxCopyPtrToInteger8(mxGetPr(prhs(14)), verb, sz)

!        call copy_int(%VAL(mxGetPr(prhs(3))), r, d+1);
!        print *, ry
!        print *, rmax
!        call copy_int(%VAL(mxGetPr(prhs(4))), ps, d+1);
!        print *, ps
!        call copy_int(%VAL(mxGetPr(prhs(1))), n, d);
!        print *, n
!        call copy_real(%VAL(mxGetPr(prhs(5))), cr, ps(d+1)-1);
!        print *, cr

	call tt_mvk4(d,n,m,rx,ra,crA,crX,crY,ry,eps,rmax,nswp=nswp,kickrank=kickrank,verb=verb)

! 	print *, 'mvk4 done'

! 	call mexPrintf(matlab_ist_dumme_kuh//achar(10));
!         call tt_mvk4(d=d,n=n,m=m,rx=rx,ra=ra,crA=crA, crX=crX, crY0=crY, ry=ry, eps=eps, rmax=rmax, nswp=nswp, kickrank=kickrank)

!         print *, ry
!         print *, result_core

!       call mxCopyPtrToInteger4(mxGetPr(prhs(0)), d, 1)
!       ptr1 = mxGetPr(prhs(1))
!       call mxCopyPtrToReal8(ptr1,scal1,sz)
!       scal = mxGetPr(prhs(1))
!       d = int(scal1)
!       print *, d

!       allocate(n(d), r(d+1), ps(d+1))

!       core_mex = mxGetPr(prhs(4));




!     Create matrix for the return argument.
      sz = d+1; sz2=1
      dumme_kuh=0
      plhs(1) = mxCreateNumericMatrix(sz,sz2, mxClassIDFromClassName('int64'), dumme_kuh)

      sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
      plhs(2) = mxCreateDoubleMatrix(sz,sz2,dumme_kuh)
!
!       scal = mxGetPr(plhs(1))
!       scal(1)=real(d)
!       x_pr = mxGetPr(prhs(1))
!       y_pr = mxGetPr(plhs(1))

!       call mxCopyPtrToReal8(x_pr,x,size)

      sz = d+1
      call mxCopyInteger8ToPtr(ry, mxGetPr(plhs(1)), sz)
      sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
      call mxCopyReal8ToPtr(result_core, mxGetPr(plhs(2)), sz)

!       print *, sz

      deallocate(result_core)
      deallocate(ry, rx,ra, n, m, crA, crX, crY)

      return
      end

