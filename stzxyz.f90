














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine stzxyz( stz, xyz, Linterpol )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  use globals, only : zero, one, machprec, tstart, ounit, &
                      oNqfms, oqfmsmn, oqfmsim, oqfmsin, oiota, oflux, orcmn, orsmn, otcmn, otsmn, drcmn, drsmn, dtcmn, dtsmn
  
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  real, intent(in)  :: stz(1:3)
  real, intent(out) :: xyz(1:3,0:3)
  integer           :: Linterpol
  
  integer           :: ii, mm, ie01bgf, ierr, ifail, iqfms, imn
  real              :: drho(0:3), dthe(0:3), lorcmn, lorsmn, lotcmn, lotsmn, ldrcmn, ldrsmn, ldtcmn, ldtsmn, lflux, carg, sarg, teta, zeta
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  lflux = stz(1) ; teta = stz(2) ; zeta = stz(3)
  
  mm = 1 ! number of interpolation points ; 01/06/21 ; 
  
  select case( Linterpol )
   
  case( 1 ) ! Linterpol = 1 ; 02/19/2021 ;
   
   lflux = min( max(lflux,oflux(1)+machprec), oflux(oNqfms)-machprec )

   drho(0:3) = (/ zero, zero, zero, zero /) ! initialize summation ; 01/06/21 ; 
   dthe(0:3) = (/ teta, zero,  one, zero /)
   
   ifail = 1
   do iqfms = 1, oNqfms-1
    if( lflux.ge.oflux(iqfms) .and. lflux.le.oflux(iqfms+1) ) then
     do imn = 1, oqfmsmn
      lorcmn = orcmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      ldrcmn =                                               ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      lotsmn = otsmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      ldtsmn =                                               ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      
      carg = cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      sarg = sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      
!     drho(0) = drho(0) + (   lorcmn * carg + lorsmn * sarg )
!     dthe(0) = dthe(0) + (   lotcmn * carg + lotsmn * sarg )   
      drho(0) = drho(0) + (   lorcmn * carg                 )
      dthe(0) = dthe(0) + (                   lotsmn * sarg )
      
!     drho(1) = drho(1) + (   ldrcmn * carg + ldrsmn * sarg )
!     dthe(1) = dthe(1) + (   ldtcmn * carg + ldtsmn * sarg )
      drho(1) = drho(1) + (   ldrcmn * carg                 )
      dthe(1) = dthe(1) + (                   ldtsmn * sarg )
      
!     drho(2) = drho(2) + ( - lorcmn * sarg + lorsmn * carg ) * ( + oqfmsim(imn) )
!     dthe(2) = dthe(2) + ( - lotcmn * sarg + lotsmn * carg ) * ( + oqfmsim(imn) )
      drho(2) = drho(2) + ( - lorcmn * sarg                 ) * ( + oqfmsim(imn) )
      dthe(2) = dthe(2) + (                   lotsmn * carg ) * ( + oqfmsim(imn) )
      
!     drho(3) = drho(3) + ( - lorcmn * sarg + lorsmn * carg ) * ( - oqfmsin(imn) )
!     dthe(3) = dthe(3) + ( - lotcmn * sarg + lotsmn * carg ) * ( - oqfmsin(imn) )
      drho(3) = drho(3) + ( - lorcmn * sarg                 ) * ( - oqfmsin(imn) )
      dthe(3) = dthe(3) + (                   lotsmn * carg ) * ( - oqfmsin(imn) )
      
     enddo ! end of do imn ; 02/19/2021 ;
     
     ifail = 0 
     
     exit
     
    endif

    if( ifail.eq.0 ) exit
    
   enddo ! end of do iqfms ; 02/19/2021 ;

   if(ifail.ne.0 ) then
          write(ounit,'("stzxyz :    fatal    : ifail.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
   
  case( 2 ) ! Linterpol = 2 ; 02/19/2021 ;

   drho(0:3) = (/ zero, zero, zero, zero /) ! initialize summation ; 01/06/21 ; 
   dthe(0:3) = (/ teta, zero,  one, zero /)
   
   do imn = 1, oqfmsmn
    
    ie01bgf = 1 ; call E01BGF( oNqfms, oflux(1:oNqfms), orcmn(1:oNqfms,imn), drcmn(1:oNqfms,imn), mm, lflux, lorcmn, ldrcmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), orsmn(1:oNqfms,imn), drsmn(1:oNqfms,imn), mm, lflux, lorsmn, ldrsmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), otcmn(1:oNqfms,imn), dtcmn(1:oNqfms,imn), mm, lflux, lotcmn, ldtcmn, ie01bgf )
    ie01bgf = 1 ; call E01BGF( oNqfms, oflux(1:oNqfms), otsmn(1:oNqfms,imn), dtsmn(1:oNqfms,imn), mm, lflux, lotsmn, ldtsmn, ie01bgf )
    
    carg = cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    sarg = sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    
!   drho(0) = drho(0) + (   lorcmn * carg + lorsmn * sarg )
!   dthe(0) = dthe(0) + (   lotcmn * carg + lotsmn * sarg )   
    drho(0) = drho(0) + (   lorcmn * carg                 )
    dthe(0) = dthe(0) + (                   lotsmn * sarg )
    
!   drho(1) = drho(1) + (   ldrcmn * carg + ldrsmn * sarg )
!   dthe(1) = dthe(1) + (   ldtcmn * carg + ldtsmn * sarg )
    drho(1) = drho(1) + (   ldrcmn * carg                 )
    dthe(1) = dthe(1) + (                   ldtsmn * sarg )
    
!   drho(2) = drho(2) + ( - lorcmn * sarg + lorsmn * carg ) * ( + oqfmsim(imn) )
!   dthe(2) = dthe(2) + ( - lotcmn * sarg + lotsmn * carg ) * ( + oqfmsim(imn) )
    drho(2) = drho(2) + ( - lorcmn * sarg                 ) * ( + oqfmsim(imn) )
    dthe(2) = dthe(2) + (                   lotsmn * carg ) * ( + oqfmsim(imn) )
    
!   drho(3) = drho(3) + ( - lorcmn * sarg + lorsmn * carg ) * ( - oqfmsin(imn) )
!   dthe(3) = dthe(3) + ( - lotcmn * sarg + lotsmn * carg ) * ( - oqfmsin(imn) )
    drho(3) = drho(3) + ( - lorcmn * sarg                 ) * ( - oqfmsin(imn) )
    dthe(3) = dthe(3) + (                   lotsmn * carg ) * ( - oqfmsin(imn) )
    
   enddo
   
!  write(ounit,'(8f15.10)') drho(0:3), dthe(0:3)

  case default ! Linterpol ; 

   if(.true. ) then
          write(ounit,'("stzxyz :    fatal    : .true.  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
   
  end select
  
  xyz(1:3,0) = (/ drho(0), dthe(0), zeta /)
  xyz(1:3,1) = (/ drho(1), dthe(1), zero /)
  xyz(1:3,2) = (/ drho(2), dthe(2), zero /)
  xyz(1:3,3) = (/ drho(3), dthe(3),  one /)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine stzxyz

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
