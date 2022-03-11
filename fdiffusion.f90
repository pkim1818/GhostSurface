














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine fdiffusion( mn, Ni, fTT, diff, dfT )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, half, one, pi2, small, ounit, tstart, &
                      tkperp, Tlow, Tupp, &
                      Ndof, dx, dy, dz, dv, &
                      fBs, fBt, fBz, fBB, &
                      im, in, Nteta, Nzeta, Ntz
  
  implicit none
  
  include "mpif.h"
  
  integer            :: mn, Ni
  real               :: fTT(1:mn,0:Ni), diff(0:4), dfT(1:mn,1:Ni-1)

  integer            :: ii, jj, kk, jk
  real               :: teta, zeta, heatflux
  real               :: ff(1:Ntz,0:Ni-1), fs(1:Ntz,0:Ni-1), ft(1:Ntz,0:Ni-1), fz(1:Ntz,0:Ni-1)
  real               :: BT(1:Ntz,0:Ni-1)
  real               :: efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), ijreal(1:Ntz), ijimag(1:Ntz), jireal(1:Ntz), jiimag(1:Ntz)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!  ijreal(1:Ntz) = zero
!  ijimag(1:Ntz) = zero
!  do jj = 0, Nteta-1 ; teta = jj * pi2 / Nteta
!   do kk = 0, Nzeta-1 ; zeta = kk * pi2 / Nzeta ; jk = 1 + jj + kk * Nteta
!    ijreal(jk) = one + cos( 1*teta- 1*zeta) + sin( 2*teta- 2*zeta)
!    ijimag(jk) = zero
!   enddo
!  enddo
!  
!  call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn) )
!
!  do ii = 1, mn
!   if( sqrt(efmn(ii)**2+ofmn(ii)**2+cfmn(ii)**2+sfmn(ii)**2).gt.small ) write(ounit,'(2i,4es13.5)') im(ii), in(ii), efmn(ii), ofmn(ii), cfmn(ii), sfmn(ii)
!  enddo
!
!  call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nteta, Nzeta, jireal(1:Ntz), jiimag(1:Ntz) )
!
!  write(ounit,'(2es13.5)') sqrt(sum((ijreal-jireal)**2)/Ntz), sqrt(sum((ijimag-jiimag)**2)/Ntz)
!
!  pause

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  fTT(1:mn, 0) = zero ; fTT(1, 0) = Tlow ! lower boundary condition ; 20/10/21 ;
  fTT(1:mn,Ni) = zero ; fTT(1,Ni) = Tupp ! upper boundary condition ; 20/10/21 ;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  diff(1) = zero ! summation ; parallel diffusion; 20/10/21 ;
  diff(2) = zero ! summation ; perpendi diffusion; 20/10/21 ;

  heatflux = zero

  do ii = 0, Ni-1
   
   efmn(1:mn) =   ( fTT(1:mn,ii+1) + fTT(1:mn,ii) ) * half ! temperature on radial half-grid; not actually required; 20/10/21 ;
   ofmn(1:mn) = zero
   cfmn(1:mn) =   ( fTT(1:mn,ii+1) - fTT(1:mn,ii) ) / dz   ! d Temp / ds on radial half-grid;
   sfmn(1:mn) = zero
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nteta, Nzeta, ff(1:Ntz,ii), fs(1:Ntz,ii) )
   
   efmn(1:mn) = zero
   ofmn(1:mn) = - ( fTT(1:mn,ii+1) + fTT(1:mn,ii) ) * half * ( + im(1:mn) ) ! d Temp / d teta on radial half-grid;
   cfmn(1:mn) = zero
   sfmn(1:mn) = - ( fTT(1:mn,ii+1) + fTT(1:mn,ii) ) * half * ( - in(1:mn) ) ! d Temp / d zeta on radial half-grid;
   
   call invfft( mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn), Nteta, Nzeta, ft(1:Ntz,ii), fz(1:Ntz,ii) )
   
   BT(1:Ntz,ii) = fBs(1:Ntz,ii) * fs(1:Ntz,ii) + fBt(1:Ntz,ii) * ft(1:Ntz,ii) + fBz(1:Ntz,ii) * fz(1:Ntz,ii) ! {\bf B}  \cdot \nabla T ; 20/10/21 ;

   ijreal(1:Ntz) = BT(1:Ntz,ii) * BT(1:Ntz,ii) * fBB(1:Ntz,ii)
   ijimag(1:Ntz) = fs(1:Ntz,ii) * fs(1:Ntz,ii) + ft(1:Ntz,ii) * ft(1:Ntz,ii) + fz(1:Ntz,ii) * fz(1:Ntz,ii) ! \nabla T \cdot \nabla T ; 20/10/21 ;

   heatflux = heatflux + sum( fs(1:Ntz,ii) + ( BT(1:Ntz,ii) * fBB(1:Ntz,ii) ) * fBs(1:Ntz,ii) / tkperp )

   call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn) )

   diff(1) = diff(1) + efmn(1)
   diff(2) = diff(2) + cfmn(1)
   
!  write(ounit,'("fdiffusion : ", 10x,"s : heatflux =",es23.15," ;")') sum( fs(1:Ntz,ii) + ( BT(1:Ntz,ii) * fBB(1:Ntz,ii) ) * fBs(1:Ntz,ii) / tkperp )

  enddo ! end of do ii ; 20/10/21 ;
  
! pause

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  heatflux = heatflux / ( Ntz * Ni )

  diff(1) = half * diff(1) * ( pi2 * pi2 * dz )
  diff(2) = half * diff(2) * ( pi2 * pi2 * dz ) * tkperp

  diff(0) = ( diff(1) + diff(2) ) / tkperp

  diff(4) = heatflux

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  BT(1:Ntz,0:Ni-1) = ( BT(1:Ntz,0:Ni-1) * fBB(1:Ntz,0:Ni-1) ) / tkperp
  
  do ii = 1, Ni-1
   
   ijreal = zero
   ijimag = ( ( fs(1:Ntz,ii  ) + BT(1:Ntz,ii  ) * fBs(1:Ntz,ii  ) ) &
            - ( fs(1:Ntz,ii-1) + BT(1:Ntz,ii-1) * fBs(1:Ntz,ii-1) ) ) / dz
   
   call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn) )

   dfT(1:mn,ii) = cfmn(1:mn) ! Fourier harmonics of radial derivative of d Temp / d s ; 20/10/21 ;
   
   ijreal = ( ( ft(1:Ntz,ii  ) + BT(1:Ntz,ii  ) * fBt(1:Ntz,ii  ) ) &
            + ( ft(1:Ntz,ii-1) + BT(1:Ntz,ii-1) * fBt(1:Ntz,ii-1) ) ) * half
   ijimag = ( ( fz(1:Ntz,ii  ) + BT(1:Ntz,ii  ) * fBz(1:Ntz,ii  ) ) &
            + ( fz(1:Ntz,ii-1) + BT(1:Ntz,ii-1) * fBz(1:Ntz,ii-1) ) ) * half
   
   call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn) )

   dfT(1:mn,ii) = dfT(1:mn,ii) + ofmn(1:mn) * im(1:mn) - sfmn(1:mn) * in(1:mn)

   dfT(2:mn,ii) = dfT(2:mn,ii) * half ! THIS IS A FUDGE ; THIS NEEDS TO BE CHECKED ; 20/10/21 ;

  enddo ! end of do ii ; 20/10/21 ;

  dfT(1:mn,1:Ni-1) = - dfT(1:mn,1:Ni-1) * ( pi2 * pi2 * dz )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 
  diff(3) = sqrt( sum( dfT(1:mn,1:Ni-1)*dfT(1:mn,1:Ni-1) ) / Ndof ) / diff(0) ! gradient magnitude ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine fdiffusion

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
