














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine diffusion( Ni, Nj, Nk, TT, diff, dT )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, half, ounit, tstart, &
                      tkperp, Tlow, Tupp, &
                      Ndof, dx, dy, dz, dv, &
                      Bx, By, Bz, BB
  
  implicit none
  
  include "mpif.h"
  
  integer            :: Ni, Nj, Nk
  real               :: TT(-1:Ni,-1:Nj,0:Nk), diff(0:4), dT(0:Ni-1,0:Nj-1,1:Nk-1)

  real               :: Tx(-1:Ni-1,-1:Nj-1,0:Nk-1), Ty(-1:Ni-1,-1:Nj-1,0:Nk-1), Tz(-1:Ni-1,-1:Nj-1,0:Nk-1)
  real               :: Bd(-1:Ni-1,-1:Nj-1,0:Nk-1), Td(-1:Ni-1,-1:Nj-1,0:Nk-1)

  real               :: dx4, dy4, dz4
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  dx4 = 4 * dx ; dy4 = 4 * dy ; dz4 = 4 * dz

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  TT(Ni,:,:) = TT(   0,:,:) ! periodicity in x = teta ; 20/10/21 ;
  TT(-1,:,:) = TT(Ni-1,:,:) ! "ghost" cell for convenience ; 20/10/21 ;  
  TT(:,Nj,:) = TT(:,   0,:) ! periodicity in y = zeta ; 20/10/21 ;
  TT(:,-1,:) = TT(:,Nj-1,:) ! "ghost" cell for convenience ; 20/10/21 ;
  TT(:,:, 0) = Tlow         ! not required ; 20/10/21 ;
  TT(:,:,Nk) = Tupp         ! not required ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) = ( TT( 0:Ni  ,-1:Nj-1,0:Nk-1) + TT( 0:Ni  ,-1:Nj-1,1:Nk  ) + TT( 0:Ni  , 0:Nj  ,0:Nk-1) + TT( 0:Ni  , 0:Nj  ,1:Nk  ) &
                               - TT(-1:Ni-1,-1:Nj-1,0:Nk-1) - TT(-1:Ni-1,-1:Nj-1,1:Nk  ) - TT(-1:Ni-1, 0:Nj  ,0:Nk-1) - TT(-1:Ni-1, 0:Nj  ,1:Nk  ) ) / dx4
    
  Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) = ( TT(-1:Ni-1, 0:Nj  ,0:Nk-1) + TT(-1:Ni-1, 0:Nj  ,1:Nk  ) + TT( 0:Ni  , 0:Nj  ,0:Nk-1) + TT( 0:Ni  , 0:Nj  ,1:Nk  ) &
                               - TT(-1:Ni-1,-1:Nj-1,0:Nk-1) - TT(-1:Ni-1,-1:Nj-1,1:Nk  ) - TT( 0:Ni  ,-1:Nj-1,0:Nk-1) - TT( 0:Ni  ,-1:Nj-1,1:Nk  ) ) / dy4
 
  Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) = ( TT(-1:Ni-1,-1:Nj-1,1:Nk  ) + TT(-1:Ni-1, 0:Nj  ,1:Nk  ) + TT( 0:Ni  ,-1:Nj-1,1:Nk  ) + TT( 0:Ni  , 0:Nj  ,1:Nk  ) &
                               - TT(-1:Ni-1,-1:Nj-1,0:Nk-1) - TT(-1:Ni-1, 0:Nj  ,0:Nk-1) - TT( 0:Ni  ,-1:Nj-1,0:Nk-1) - TT( 0:Ni  , 0:Nj  ,0:Nk-1) ) / dz4
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) = Bx(-1:Ni-1,-1:Nj-1,0:Nk-1) * Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) &
                             + By(-1:Ni-1,-1:Nj-1,0:Nk-1) * Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) &
                             + Bz(-1:Ni-1,-1:Nj-1,0:Nk-1) * Tz(-1:Ni-1,-1:Nj-1,0:Nk-1)

  Td(-1:Ni-1,-1:Nj-1,0:Nk-1) = Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) * Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) &
                             + Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) * Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) &
                             + Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) * Tz(-1:Ni-1,-1:Nj-1,0:Nk-1)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
! diff(1) = half *          sum( Bd(0:Ni-1,0:Nj-1,0:Nk-1) * Bd(0:Ni-1,0:Nj-1,0:Nk-1) / BB(0:Ni-1,0:Nj-1,0:Nk-1) ) * dv
! diff(2) = half * tkperp * sum( Td(0:Ni-1,0:Nj-1,0:Nk-1)                                                       ) * dv
  diff(1) = half *          sum( Bd(0:Ni-1,0:Nj-1,0:Nk-1) * Bd(0:Ni-1,0:Nj-1,0:Nk-1) / BB(0:Ni-1,0:Nj-1,0:Nk-1) ) * dv / tkperp
  diff(2) = half *          sum( Td(0:Ni-1,0:Nj-1,0:Nk-1)                                                       ) * dv
  
  diff(0) = diff(1) + diff(2)

  diff(4) = zero ! this should be the heatflux ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) = Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) / BB(-1:Ni-1,-1:Nj-1,0:Nk-1)

! Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) = tkperp * Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) + Bx(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1)
! Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) = tkperp * Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) + By(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1)
! Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) = tkperp * Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) + Bz(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1)
  Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) =          Tx(-1:Ni-1,-1:Nj-1,0:Nk-1) + Bx(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) / tkperp
  Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) =          Ty(-1:Ni-1,-1:Nj-1,0:Nk-1) + By(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) / tkperp
  Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) =          Tz(-1:Ni-1,-1:Nj-1,0:Nk-1) + Bz(-1:Ni-1,-1:Nj-1,0:Nk-1) * Bd(-1:Ni-1,-1:Nj-1,0:Nk-1) / tkperp

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  dT( 0:Ni-1, 0:Nj-1,1:Nk-1) = ( Tx( 0:Ni-1,-1:Nj-2,0:Nk-2) + Tx( 0:Ni-1,-1:Nj-2,1:Nk-1) + Tx( 0:Ni-1, 0:Nj-1,0:Nk-2) + Tx( 0:Ni-1, 0:Nj-1,1:Nk-1) &
                               - Tx(-1:Ni-2,-1:Nj-2,0:Nk-2) - Tx(-1:Ni-2,-1:Nj-2,1:Nk-1) - Tx(-1:Ni-2, 0:Nj-1,0:Nk-2) - Tx(-1:Ni-2, 0:Nj-1,1:Nk-1) ) / dx4 &
                             + ( Ty(-1:Ni-2, 0:Nj-1,0:Nk-2) + Ty(-1:Ni-2, 0:Nj-1,1:Nk-1) + Ty( 0:Ni-1, 0:Nj-1,0:Nk-2) + Ty( 0:Ni-1, 0:Nj-1,1:Nk-1) &
                               - Ty(-1:Ni-2,-1:Nj-2,0:Nk-2) - Ty(-1:Ni-2,-1:Nj-2,1:Nk-1) - Ty( 0:Ni-1,-1:Nj-2,0:Nk-2) - Ty( 0:Ni-1,-1:Nj-2,1:Nk-1) ) / dy4 &
                             + ( Tz(-1:Ni-2,-1:Nj-2,1:Nk-1) + Tz(-1:Ni-2, 0:Nj-1,1:Nk-1) + Tz( 0:Ni-1,-1:Nj-2,1:Nk-1) + Tz( 0:Ni-1, 0:Nj-1,1:Nk-1) &
                               - Tz(-1:Ni-2,-1:Nj-2,0:Nk-2) - Tz(-1:Ni-2, 0:Nj-1,0:Nk-2) - Tz( 0:Ni-1,-1:Nj-2,0:Nk-2) - Tz( 0:Ni-1, 0:Nj-1,0:Nk-2) ) / dz4

  dT = - dT * dv

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  diff(3) = sqrt( sum( dT(0:Ni-1,0:Nj-1,1:Nk-1)*dT(0:Ni-1,0:Nj-1,1:Nk-1) ) / Ndof ) / diff(0) ! gradient magnitude ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine diffusion

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
