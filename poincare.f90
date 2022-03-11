














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine poincare
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : ounit, tunit, zero, pi2, one, myid, tstart, &
                      hmnfile, chi, &
                      NPtrj, NPpts, odetol, &
                      Nqfms, pp, qq, Pxyz, &
                      Linterpol
  
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, parameter :: Node = 6
  
  integer            :: ifail, jj, ipq, ii, ipoincare, lNode, ierr
  real               :: cput, stz(1:3), xy(1:Node), Bxy(1:Node), zst, zed, d02bjfwk(1:20*Node), lxy(1:2,0:NPpts), xyz(1:3,0:3), lodetol
  
  EXTERNAL           :: bfield, pfield, D02BJX, D02BJW
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  lodetol = odetol
  
  do ipoincare = 0, 1 ! construct the Poincare plot in chaotic coordinates and then in background coordinates ; 20/10/21 ;
   
   cput = MPI_WTIME()
   
   write(ounit,'("poincare  : ",f10.1,"s : writing poincare ; NPpts =",i6," ; NPtrj =",i4," ; odetol =",es8.1," ;")') cput-tstart, NPpts, NPtrj, odetol
   
   if( ipoincare.eq.0 ) open( tunit, file = "."//trim(hmnfile)//".poincare.s", status='unknown', form='unformatted')
   if( ipoincare.eq.1 ) open( tunit, file = "."//trim(hmnfile)//".poincare"  , status='unknown', form='unformatted')
   
   write(tunit) NPpts
   
   do ipq = 2, Nqfms-1 + NPtrj-1
    
    do ii = 0, 2 ! the starting points are the X, midpoint and O points of the islands, which are located in action.h ; 20/10/21 ;
     
     if( ipq.lt.Nqfms ) then
      if(.not.allocated(qq) ) then
          write(ounit,'("poincare  :    fatal    : .not.allocated(qq)  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
      if( qq(ipq).gt.8 ) cycle ! don't include high-order islands ; 02/19/2021 ;
      stz(1:3) = Pxyz(1:3,ii,ipq)
     else
      if( ii.gt.0 ) cycle ! only one field line per radial value ; 02/19/2021 ;
      stz(1:3) = (/ chi%low + ( ipq - Nqfms + 1 ) * ( chi%upp - chi%low ) / NPtrj, zero, zero /)
     endif
     
     call stzxyz( stz(1:3), xyz(1:3,0:3), Linterpol ) ! the starting points are known in chaotic coordinates, and need to be mapped ; 20/10/21 ;
     
     if( ipoincare.eq.0 ) xy(1:2) = (/ stz(2)  , stz(1)   /)
     if( ipoincare.eq.1 ) xy(1:2) = (/ xyz(2,0), xyz(1,0) /)
     
     lxy(1:2,0) = xy(1:2)

     do jj = 1, NPpts
      
      zst = zero ; zed = pi2 ; xy(1:Node) = (/ xy(1), xy(2), one, zero, zero, one /)
      
      ifail = 0 ; lNode = Node
      
      if( ipoincare.eq.0 ) call D02BJF( zst, zed, lNode, xy, pfield, lodetol, 'D', D02BJX, D02BJW, d02bjfwk, ifail )
      if( ipoincare.eq.1 ) call D02BJF( zst, zed, lNode, xy, bfield, lodetol, 'D', D02BJX, D02BJW, d02bjfwk, ifail )
      
      xy(1) = mod( xy(1), pi2 )
      
      if( xy(1).lt.zero ) xy(1) = xy(1) + pi2
      
      lxy(1:2,jj) = xy(1:2)
      
     enddo
     
     write(tunit) lxy(1:2,0:NPpts)
     
    enddo ! end of do ii ; 12/17/20;
    
   enddo ! end of do ipq ; 12/17/20;
   
   close(tunit)
   
  enddo ! end of do ipoincare ; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine poincare

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine pfield( zeta, xys, Bxys )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  use globals, only : zero, one, tstart, ounit, &
                      Bzeta, &
                      Linterpol
  
  implicit none  
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, parameter   :: Node = 6
  real   , intent(in)  :: zeta, xys(1:Node)
  real   , intent(out) :: Bxys(1:Node)
  
  integer              :: lLinterpol
  real                 :: xy(1:Node), Bxy(1:Node), Delta, AA(1:3,1:3), Brtp(1:3), stz(1:3), xyz(1:3,0:3)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  stz(1:3) = (/ xys(2), xys(1), zeta /)

  lLinterpol = Linterpol
! lLinterpol = 2         ! using C1 coordinate transformation ; 02/19/2021 ;

  call stzxyz( stz(1:3), xyz(1:3,0:3), lLinterpol ) ! coordinate transformation ;

  Delta = xyz(2,2) * xyz(1,1) - xyz(2,1) * xyz(1,2)
  
  AA(1,1) = + xyz(2,2) / Delta  ; AA(1,2) = - xyz(1,2) / Delta ; AA(1,3) = (   xyz(2,3) * xyz(1,2) - xyz(2,2) * xyz(1,3) ) / Delta ! vector transformation ; 
  AA(2,1) = - xyz(2,1) / Delta  ; AA(2,2) = + xyz(1,1) / Delta ; AA(2,3) = ( - xyz(2,3) * xyz(1,1) + xyz(2,1) * xyz(1,3) ) / Delta
  AA(3,1) =   zero              ; AA(3,2) =   zero             ; AA(3,3) =     one

  xy(1:Node) = (/ xyz(2,0), xyz(1,0), one, zero, zero, one /)

  call bfield( zeta, xy(1:Node), Bxy(1:Node) )

  Brtp(1:3) = matmul( AA(1:3,1:3), (/ Bxy(2), Bxy(1), Bzeta /) )
    
  Bxys(1:Node) = (/ Brtp(2), Brtp(1), zero, zero, zero, zero /)

  Bzeta = one

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine pfield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
