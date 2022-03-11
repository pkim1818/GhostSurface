














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine bicubic( li, lj, lx, ly, wgt )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  use globals, only : ounit, one, zero, small, dx, dy, BC, bcf, FD, fdf
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  implicit none  

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, intent(in)  :: li, lj ! lower, left adjacent grid point;
  real   , intent(in)  :: lx, ly ! local interpolation variables;
  real   , intent(out) :: wgt(1:16) ! weights;
  integer              :: ierr, ii, jj, kk, nn(1:16), mm(1:16), iin(1:16), jin(1:16) ! cubic-powers;
  real                 :: lxn(0:3), lym(0:3)  ! local interpolation;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  nn(1:16) = (/  0,  1,  2,  3,  0,  1,  2,  3,  0,  1,  2,  3,  0,  1,  2,  3 /)
  mm(1:16) = (/  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3 /)
  
  iin(1:16) = (/ -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2 /) ! this must be consistent with matrix.h; probably should make global ; 12/17/20;
  jin(1:16) = (/ -1,-1,-1,-1,  0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2 /)
  
  lxn(0) = one ; lxn(1) = lx * lxn(0) ; lxn(2) = lx * lxn(1) ; lxn(3) = lx * lxn(2) ! powers of xx, local interpolation ; 12/17/20;
  lym(0) = one ; lym(1) = ly * lym(0) ; lym(2) = ly * lym(1) ; lym(3) = ly * lym(2) ! powers of yy, local interpolation ; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  wgt(1:16) = zero
  
  do kk = 1, 16
   
   do jj = 1, 16
    
    if( FD(jj,kk).eq.0 ) cycle ! matrices are sparse;
    
    do ii = 1, 16
     
     if( BC(ii,jj).eq.0 ) cycle ! matrices are sparse;
     
     wgt(kk) = wgt(kk) + lxn(nn(ii)) * lym(mm(ii)) * BC(ii,jj) * bcf(jj) * FD(jj,kk) * fdf(jj)
     
    enddo
    
   enddo
   
  enddo
  
  if(abs(sum(wgt)-one).gt.small ) then
          write(ounit,'("bicubic  :    fatal    : abs(sum(wgt)-one).gt.small  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine bicubic

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
