program test_amr

use tt_lib
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

print *, tt_A%m, tt_A%l

! allocate(n(d),na(d),m(d),ry(d+1),ra(d+1),rx(d+1),ps(d+1),psa(d+1))

! call dsdv_to_arrays(na,ra,d,psa,crA,tt_A)
! call dsdv_to_arrays(n,rx,d,ps,crb,tt_b)



! deallocate (crA, crb)

end program