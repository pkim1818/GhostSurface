














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine volume
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, half, pi2, &
                      ounit, &
                      myid, tstart, &
                      chi, tkperp, &
                      Ni, Nj, Nk, dx, dy, dz, &
                      odetol, &
                      ijk
  
  implicit none

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  integer, parameter :: Node = 9, MM = 1
  integer            :: ii, jj, ifail, ie02def, iintegration
  real               :: cput, cpuo, xyv(1:Node), xx, yy, zst, zed, d02bjfwk(1:20*Node), lvolume, evolume, pvolume, dvolume, verr
  real               :: lx(1:MM), ly(1:MM), ff(1:MM)
  
  EXTERNAL           :: vfield, D02BJX, D02BJW
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cpuo = MPI_WTIME()
  
  write(ounit,1000) cpuo-tstart, "computing volume integrals"
  
1000 format("volume    : ",f10.2,"s : ",a," ;",:," volume =",f19.15," ; err =",es8.1," ; pvolume =",es12.5," ; dvolume =",es12.5," ; time =",f10.2,"s ;")
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  lvolume = zero
  pvolume = zero
  dvolume = zero
  
  do ii = 0, Ni-1

   do jj = 0, Nj-1 ; ijk(1:3) = (/ ii, jj,  0 /)
    
    xx = ( ii + half ) * dx ; yy = chi%low + ( jj + half ) * dy
    
    xyv(1:Node) = (/ xx , yy , one , zero , zero , one, zero, zero, zero /)
    
    zst = zero ; zed = pi2
    
    ifail = 0 ; call D02BJF( zst, zed, Node, xyv(1:Node), vfield, odetol, 'D', D02BJX, D02BJW, d02bjfwk, ifail )
    
    lvolume = lvolume + xyv(7)
    pvolume = pvolume + xyv(8)
    dvolume = dvolume + xyv(9)
    
   enddo ! end of do jj
   
  enddo ! end of do ii
  
  lvolume = lvolume * dx * dy
  pvolume = pvolume * dx * dy
  dvolume = dvolume * dx * dy * tkperp
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  evolume = pi2 * pi2 * ( chi%upp - chi%low ) ; verr = lvolume-evolume ! this is the exact volume ; 12/17/20;
  
  cput = MPI_WTIME()
  
  write(ounit,1000) cput-tstart, "computed  volume integrals", lvolume, abs(verr), pvolume, dvolume, cput-cpuo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine volume

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine vfield( zeta, xyv, Bxyv )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  use globals, only : zero, one, pi2, &
                     chi, &
                     Bzeta, &
                     C0, D0, E0, C1, D1, E1, &
                     ijk
  
  implicit none  

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, parameter   :: Node = 9
  real   , intent(in)  :: zeta, xyv(1:Node)
  real   , intent(out) :: Bxyv(1:Node)
  
  integer              :: imn, ii, jj
  real                 :: teta, arg, carg, sarg, TM(1:2,1:2), DB(1:2,1:2), Tp, Ta, Tz, Bsqd, AA(1:3,1:3), g11, g12, g13, g22, g23, g33, Delta
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  Bxyv(1:Node) = zero ; teta = xyv(1) ; ii = ijk(1) ; jj = ijk(2) ! initialize summation variable; shorthand counters;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  call sumpolynomial( chi%mn, chi%Np, chi%kmn(1:chi%mn,0:chi%Np), xyv(2), chi%k(1:chi%mn,0:2) ) 
  
  do imn = 1, chi%mn
   
   arg = chi%im(imn)*teta - chi%in(imn)*zeta ; sarg = sin(arg) ; carg = cos(arg)
   
   Bxyv(1) = Bxyv(1) + chi%k(imn,1) * carg                                     !      \dot xx
   Bxyv(2) = Bxyv(2) + chi%k(imn,0) * sarg * ( chi%im(imn) )                   !      \dot yy
   
   Bxyv(3) = Bxyv(3) + chi%k(imn,1) * sarg * (-chi%im(imn) )                   ! d_xx \dot xx
   Bxyv(4) = Bxyv(4) + chi%k(imn,2) * carg                                     ! d_yy \dot xx
   
   Bxyv(5) = Bxyv(5) + chi%k(imn,0) * carg * ( chi%im(imn) ) * ( chi%im(imn) ) ! d_xx \dot yy
   Bxyv(6) = Bxyv(6) + chi%k(imn,1) * sarg * ( chi%im(imn) )                   ! d_yy \dot yy
   
  enddo
  
  Bzeta = one
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  TM(1,1) = Bxyv(3) ; TM(1,2) = Bxyv(4) ; TM(2,1) = Bxyv(5) ; TM(2,2) = Bxyv(6) ! packing local      tangent map into matrix ; 20/10/21 ;
  DB(1,1) =  xyv(3) ; DB(1,2) =  xyv(4) ; DB(2,1) =  xyv(5) ; DB(2,2) =  xyv(6) ! packing integrated tangent map into matrix ; 20/10/21 ;
  
  DB(1:2,1:2) = matmul( TM(1:2,1:2) , DB(1:2,1:2) )                             ! product tangent map ; 20/10/21 ;

  Bxyv(3) = DB(1,1) ; Bxyv(4) = DB(1,2) ; Bxyv(5) = DB(2,1) ; Bxyv(6) = DB(2,2) ! product tangent map; unpacking ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

! Bxyv(7) = ( one + Bxyv(6) ) * ( one + Bxyv(3) ) - Bxyv(4) * Bxyv(5) ! Jacobian ; 12/17/20; ! this seems to have worked, but it must be an error ; 20/10/21 ;
  Bxyv(7) =          xyv(6)   *          xyv(3)   -  xyv(4) *  xyv(5) ! Jacobian ; 12/17/20;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  Tp = E0(ii,jj) + ( E1(ii,jj) - E0(ii,jj) ) * zeta / pi2 ! radial   derivative ; 20/10/21 ;
  Ta = D0(ii,jj) + ( D1(ii,jj) - D0(ii,jj) ) * zeta / pi2 ! poloidal derivative ; 20/10/21 ;
  Tz = C0(ii,jj) + ( C1(ii,jj) - C0(ii,jj) )        / pi2 ! parallel derivative ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  Bsqd = Bxyv(1)**2 + Bxyv(2)**2 + Bzeta**2

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  Bxyv(8) = Bxyv(7) * ( Bzeta * Tz ) * ( Bzeta * Tz ) / Bsqd ! parallel derivative ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  Delta = Bxyv(7)
  
  AA(1,1) = + xyv(3)   / Delta  ; AA(1,2) = - xyv(5)   / Delta ; AA(1,3) = (   xyv(5)   * Bxyv(1)  - Bxyv(2)  * xyv(3)   ) / Delta
  AA(2,1) = - xyv(4)   / Delta  ; AA(2,2) = + xyv(6)   / Delta ; AA(2,3) = (   Bxyv(2)  * xyv(4)   - Bxyv(1)  * xyv(6)   ) / Delta
  AA(3,1) =   zero              ; AA(3,2) =   zero             ; AA(3,3) =     one
  
  g11 = AA(1,1) * AA(1,1) + AA(1,2) * AA(1,2) + AA(1,3) * AA(1,3)
  g12 = AA(1,1) * AA(2,1) + AA(1,2) * AA(2,2) + AA(1,3) * AA(2,3)
  g13 = AA(1,1) * AA(3,1) + AA(1,2) * AA(3,2) + AA(1,3) * AA(3,3)
  g22 = AA(2,1) * AA(2,1) + AA(2,2) * AA(2,2) + AA(2,3) * AA(2,3)
  g23 = AA(2,1) * AA(3,1) + AA(2,2) * AA(3,2) + AA(2,3) * AA(3,3)
  g33 = AA(3,1) * AA(3,1) + AA(3,2) * AA(3,2) + AA(3,3) * AA(3,3)
  
  Bxyv(9) = Tp * Tp * g11 + 2 * Tp * Ta * g12 + 2 * Tp * Tz * g13 &
                          + 1 * Ta * Ta * g22 + 2 * Ta * Tz * g23 &
                                              + 1 * Tz * Tz * g33

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine vfield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

