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
! PARAMS: nswp, kickrank, nrestart, niters, max_full_size, trunc_norm, verb
      use ttamen_lib
      use python_conv_lib
      use tt_lib

!     Declarations
      implicit none

!     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer(4) nlhs, nrhs

!     Function declarations:
      mwPointer mxCreateDoubleMatrix, mxGetPr, mxGetPi
      integer(4) mxIsNumeric
      integer(4) mxIsComplex
      mwPointer mxCreateNumericMatrix
      integer(4) mxClassIDFromClassName
      integer(4) complex_data
      mwSize sz, sz2
      real(8) dznrm2

      integer d, i
      character prec
      real(8), allocatable :: dcrA(:), dcrX(:), dcrY(:)
      type(dtt) :: dA, dY, dX
      complex(8), allocatable :: zcrA(:), zcrX(:), zcrY(:)
      type(ztt) :: zA, zY, zX
      integer,allocatable :: n(:),rx(:),ra(:),ry(:),m(:),params(:), psa(:), psy(:), psx(:)
      real(8) eps

      sz = 1
      call mxCopyPtrToInteger8(mxGetPr(prhs(1)), d, sz)
      allocate(rx(d+1), ry(d+1), ra(d+1), n(d), m(d), params(7), psa(d+1), psy(d+1), psx(d+1))

      sz = d;
      call mxCopyPtrToInteger8(mxGetPr(prhs(2)), n, sz)
      call mxCopyPtrToInteger8(mxGetPr(prhs(3)), m, sz)

      sz = d+1;
      call mxCopyPtrToInteger8(mxGetPr(prhs(4)), ry, sz)
      call mxCopyPtrToInteger8(mxGetPr(prhs(5)), ra, sz)
      call mxCopyPtrToInteger8(mxGetPr(prhs(9)), rx, sz)

      sz = 1
      call mxCopyPtrToReal8(mxGetPr(prhs(10)), eps, sz)
      call mxCopyPtrToCharacter(mxGetPr(prhs(12)), prec, sz)
      sz = 7
      call mxCopyPtrToInteger8(mxGetPr(prhs(11)), params, sz)

      ! We need to determine whether we have complex data
      if ((mxIsComplex(prhs(6))).or.(mxIsComplex(prhs(7))).or.(mxIsComplex(prhs(8)))) then

	sz = sum(ra(2:d+1)*ra(1:d)*n(1:d)*m(1:d))
	allocate (zcrA(sz))
	call mxCopyPtrToComplex16(mxGetPr(prhs(6)),mxGetPi(prhs(6)),zcrA,sz)
	sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
	allocate (zcrY(sz))
	call mxCopyPtrToComplex16(mxGetPr(prhs(7)),mxGetPi(prhs(7)),zcrY,sz)
	sz = sum(rx(2:d+1)*rx(1:d)*m(1:d))
	allocate (zcrX(sz))
	call mxCopyPtrToComplex16(mxGetPr(prhs(8)),mxGetPi(prhs(8)),zcrX,sz)
	! Hey Mathworks, this naming sucks, because Complex16 is the quadruple precision!

	! Compute positions of the cores
	psa(1)=1
	psy(1)=1
	psx(1)=1
	do i=1,d
	  psa(i+1) = psa(i)+ra(i)*n(i)*m(i)*ra(i+1)
	  psx(i+1) = psx(i)+rx(i)*m(i)*rx(i+1)
	  psy(i+1) = psy(i)+ry(i)*n(i)*ry(i+1)
	end do
	! Copy to SDV format
	call zarrays_to_sdv(n*m,ra,d,psa,zcrA,zA)
	call zarrays_to_sdv(n,ry,d,psy,zcrY,zY)
	call zarrays_to_sdv(m,rx,d,psx,zcrX,zX)

	call ztt_amen_solve(zA, zY, eps, zX, kickrank=params(2), nswp=params(1), local_prec=prec, local_iters=params(4), local_restart=params(3), trunc_norm=params(6), max_full_size=params(5), verb=params(7))

	deallocate(zcrX) ! next line will allocate it accordingly
        call zsdv_to_arrays(m,rx,d,psx,zcrX,zX)

        ! Create matrix for the return argument.
	sz = d+1; sz2=1
	complex_data=0
	! This is rx
	plhs(1) = mxCreateNumericMatrix(sz,sz2, mxClassIDFromClassName('int64'), complex_data)
	call mxCopyInteger8ToPtr(rx, mxGetPr(plhs(1)), sz)

	! This is core
	complex_data=1
	sz = sum(rx(2:d+1)*rx(1:d)*n(1:d))
	plhs(2) = mxCreateDoubleMatrix(sz,sz2,complex_data)
	call mxCopyComplex16ToPtr(zcrX, mxGetPr(plhs(2)), mxGetPi(plhs(2)), sz)

	deallocate(zcrA, zcrX, zcrY)
	call dealloc(zA)
	call dealloc(zX)
	call dealloc(zY)

      else ! Real data

	sz = sum(ra(2:d+1)*ra(1:d)*n(1:d)*m(1:d))
	allocate (dcrA(sz))
	call mxCopyPtrToReal8(mxGetPr(prhs(6)),dcrA,sz)
	sz = sum(ry(2:d+1)*ry(1:d)*n(1:d))
	allocate (dcrY(sz))
	call mxCopyPtrToReal8(mxGetPr(prhs(7)),dcrY,sz)
	sz = sum(rx(2:d+1)*rx(1:d)*m(1:d))
	allocate (dcrX(sz))
	call mxCopyPtrToReal8(mxGetPr(prhs(8)),dcrX,sz)

	! Compute positions of the cores
	psa(1)=1
	psy(1)=1
	psx(1)=1
	do i=1,d
	  psa(i+1) = psa(i)+ra(i)*n(i)*m(i)*ra(i+1)
	  psx(i+1) = psx(i)+rx(i)*m(i)*rx(i+1)
	  psy(i+1) = psy(i)+ry(i)*n(i)*ry(i+1)
	end do
	! Copy to SDV format
	call darrays_to_sdv(n*m,ra,d,psa,dcrA,dA)
	call darrays_to_sdv(n,ry,d,psy,dcrY,dY)
	call darrays_to_sdv(m,rx,d,psx,dcrX,dX)

        call dtt_amen_solve(dA, dY, eps, dX, kickrank=params(2), nswp=params(1), local_prec=prec, local_iters=params(4), local_restart=params(3), trunc_norm=params(6), max_full_size=params(5), verb=params(7))

	deallocate(dcrX) ! next line will allocate it accordingly
        call dsdv_to_arrays(m,rx,d,psx,dcrX,dX)

        ! Create matrix for the return argument.
	sz = d+1; sz2=1
	complex_data=0
	! This is rx
	plhs(1) = mxCreateNumericMatrix(sz,sz2, mxClassIDFromClassName('int64'), complex_data)
	call mxCopyInteger8ToPtr(rx, mxGetPr(plhs(1)), sz)

	! This is core
	sz = sum(rx(2:d+1)*rx(1:d)*n(1:d))
	plhs(2) = mxCreateDoubleMatrix(sz,sz2,complex_data)
	call mxCopyReal8ToPtr(dcrX, mxGetPr(plhs(2)), sz)

	deallocate(dcrA, dcrX, dcrY)
	call dealloc(dA)
	call dealloc(dX)
	call dealloc(dY)

      end if

      deallocate(ry, rx,ra, n, m, params, psa, psx, psy)

      end subroutine

