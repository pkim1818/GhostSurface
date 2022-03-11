














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine sumpolynomial( mn, Np, kmn, yy, kk ) 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  use globals, only : zero, one

  implicit none

  include "mpif.h"

  integer, intent(in ) :: mn, Np
  real   , intent(in ) :: kmn(1:mn,0:Np), yy
  real   , intent(out) :: kk(1:mn,0:2) 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  integer :: imn, nn
  real    :: yn(0:Np)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  ;             ; yn( 0) = one
  do nn = 1, Np ; yn(nn) = yn(nn-1) * yy
  enddo
  
  kk(1:mn,0:2) = zero
  
  do imn = 1, mn
   
   do nn = 0, Np

                kk(imn,0) = kk(imn,0) + kmn(imn,nn) * yn(nn  )
    if(nn.gt.0) kk(imn,1) = kk(imn,1) + kmn(imn,nn) * yn(nn-1) * nn
    if(nn.gt.1) kk(imn,2) = kk(imn,2) + kmn(imn,nn) * yn(nn-2) * nn * (nn-1)
   !if(nn.gt.2) kk(imn,3) = kk(imn,3) + kmn(imn,nn) * yn(nn-3) * nn * (nn-1) * (nn-2)

   enddo

  enddo

  return

end subroutine sumpolynomial

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
