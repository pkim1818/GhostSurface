














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine bfield( zeta, xy, Bxy )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  use globals, only : zero, one, chi, Bzeta
  
  implicit none  

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, parameter :: Node = 6
  real, intent(in)   :: zeta, xy(1:Node)
  real, intent(out)  :: Bxy(1:Node)
  
  integer            :: imn
  real               :: lrho, teta, arg, carg, sarg, TM(1:2,1:2), DB(1:2,1:2)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  Bxy(1:Node) = zero ! summation variable;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
! lrho = min( chi%upp, max( xy(2), chi%low ) )
  lrho = xy(2) ; teta = xy(1)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  call sumpolynomial( chi%mn, chi%Np, chi%kmn(1:chi%mn,0:chi%Np), lrho, chi%k(1:chi%mn,0:2) ) 
  
  do imn = 1, chi%mn

   if( imn.gt.1 .and. ( lrho.lt.chi%low .or. lrho.gt.chi%upp ) ) chi%k(imn,0:2) = zero ! force perturbation to zero outside domain ; 02/19/2021 ;
   
   arg = chi%im(imn) * teta - chi%in(imn) * zeta ; sarg = sin(arg) ; carg = cos(arg)
   
   Bxy(1) = Bxy(1) + chi%k(imn,1) * carg                                     !      \dot xx
   Bxy(2) = Bxy(2) + chi%k(imn,0) * sarg * ( chi%im(imn) )                   !      \dot yy
   
!  Bxy(3) = Bxy(3) + chi%k(imn,1) * sarg * (-chi%im(imn) )                   ! d_xx \dot xx
!  Bxy(4) = Bxy(4) + chi%k(imn,2) * carg                                     ! d_yy \dot xx
   
!  Bxy(5) = Bxy(5) + chi%k(imn,0) * carg * ( chi%im(imn) ) * ( chi%im(imn) ) ! d_xx \dot yy
!  Bxy(6) = Bxy(6) + chi%k(imn,1) * sarg * ( chi%im(imn) )                   ! d_yy \dot yy
   
  enddo
  
  Bzeta = one
  
! TM(1,1) = Bxy(3) ; TM(1,2) = Bxy(4) ; TM(2,1) = Bxy(5) ; TM(2,2) = Bxy(6)
! DB(1,1) =  xy(3) ; DB(1,2) =  xy(4) ; DB(2,1) =  xy(5) ; DB(2,2) =  xy(6)
  
! DB = matmul( TM , DB )

! Bxy(3) = DB(1,1) ; Bxy(4) = DB(1,2) ; Bxy(5) = DB(2,1) ; Bxy(6) = DB(2,2)
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine bfield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 


! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine afield( zeta, xy, Bxy )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  use globals, only : zero, one, chi, Bzeta
  
  implicit none  

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, parameter :: Node = 6
  real, intent(in)   :: zeta, xy(1:Node)
  real, intent(out)  :: Bxy(1:Node)
  
  integer            :: imn
  real               :: lrho, teta, arg, carg, sarg, TM(1:2,1:2), DB(1:2,1:2)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  Bxy(1:Node) = zero ! summation variable;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  lrho = xy(2) ; teta = xy(1)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  call sumpolynomial( chi%mn, chi%Np, chi%kmn(1:chi%mn,0:chi%Np), lrho, chi%k(1:chi%mn,0:2) ) 
  
  do imn = 1, chi%mn
   
   arg = chi%im(imn) * teta - chi%in(imn) * zeta ; sarg = sin(arg) ; carg = cos(arg)
   
   Bxy(1) = Bxy(1) + chi%k(imn,1) * carg                                     !      \dot xx
   Bxy(2) = Bxy(2) + chi%k(imn,0) * sarg * ( chi%im(imn) )                   !      \dot yy
   
  enddo
  
  Bzeta = one
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine afield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

