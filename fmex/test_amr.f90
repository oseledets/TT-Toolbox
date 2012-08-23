program test_amr

use tt_lib
use ttop_lib
use ttio_lib
use python_conv_lib
use tt_adapt_als

! d,n,m,ry,ra,crA, crY, crX0, rx, eps, kickrank, nswp, verb, prec, nrestart, niters
integer :: d
integer, allocatable :: n(:),m(:),ry(:),ra(:), rx(:), ps(:), na(:), psa(:)
integer :: kickrank, nswp, verb, nrestart, niters
! verb: 0: off; 1: matlab; 2: console
character :: prec ! 'n', 'l', 'r', 'c'
real(8),allocatable :: crA(:), crb(:)

type(dtt) :: tt_A, tt_b

call dtt_read(tt_A,'A.sdv')
call dtt_read(tt_b,'b.sdv')

d = tt_A%m-tt_A%l+1

!call svd(tt_A,1d-6)
allocate(n(d),na(d),m(d),ry(d+1),ra(d+1),rx(d+1),ps(d+1),psa(d+1))

 call dsdv_to_arrays(na,ra,d,psa,crA,tt_A)
!print *,na(1:d)
call dsdv_to_arrays(n,rx,d,ps,crb,tt_b)

call dealloc(tt_A)
call dealloc(tt_b)
ry(1:d+1) = rx(1:d+1)
call tt_amr_solve(d,n,n,ry,ra,crA, crb, crb, rx, 1d-6)

! deallocate (crA, crb)

end program
