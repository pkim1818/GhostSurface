














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine viscous
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, two, goldmean, ounit, myid, tstart, &
                      Lmoser, kpert, npert, Mpol, NFarey, viscosity, &
                      munit, &
                      Nqfms, pp, qq
                                
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer              :: ifail, astat, ipert, ipq
  real                 :: cput
  
  integer              :: Ndof, ibound, liwork, lrwork, iuser(1:1), mm
  integer, allocatable :: iwork(:)
  real                 :: action, daction(1:Mpol), ruser(0:32*Mpol+1), laction, uaction, fdiff = 1.0E-04, lx(1:Mpol), xx(1:Mpol), bl(1:Mpol), bu(1:Mpol), error
  real   , allocatable :: rwork(:), iota(:)

  EXTERNAL             :: Fmoser, Dmoser
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  Ndof = Mpol ; ibound = 1 ; bl(1:Mpol) = - one ; bl(1:Mpol) = + one
  
  liwork = Ndof + 2 ; lrwork = maxval( (/ Ndof * ( Ndof - 1 ) / 2 + 12 * Ndof + 13, 10 * Ndof + Ndof * ( Ndof - 1 ) / 2, 11 /) )  
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( iwork ) ) then
    write(0,'("macros    : 0123456789 : iwork already allocated ;")') 
    stop      'macros    : 0123456789 : iwork already allocated ;'
   endif
!#endif

   allocate( iwork(1:liwork), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating iwork ;")') 
    stop      'macros   : 0123456789 : error allocating iwork ;'
   endif
!#endif

   iwork(1:liwork) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rwork ) ) then
    write(0,'("macros    : 0123456789 : rwork already allocated ;")') 
    stop      'macros    : 0123456789 : rwork already allocated ;'
   endif
!#endif

   allocate( rwork(1:lrwork), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rwork ;")') 
    stop      'macros   : 0123456789 : error allocating rwork ;'
   endif
!#endif

   rwork(1:lrwork) = zero 

 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!  cput = MPI_WTIME()
!  
!  Nqfms = 2
!  SALLOCATE( pp, (0:Nqfms-1), (/ 0, 1/) )
!  SALLOCATE( qq, (0:NPq-1), (/ 1, 1/) )
!  
!  do iFarey = 1, NFarey
!   
!   lNqfms = 2*Nqfms-1
!
!   SALLOCATE( lp, (0:lNqfms-1), 0 )
!   SALLOCATE( lq, (0:lNqfms-1), 0 )
!
!   lp(0:lNqfms-1:2) = pp(0:Nqfms-1) ; lp(1:lNqfms-2:2) = pp(0:Nqfms-2) + pp(1:Nqfms-1)
!   lq(0:lNqfms-1:2) = qq(0:Nqfms-1) ; lq(1:lNqfms-2:2) = qq(0:Nqfms-2) + qq(1:Nqfms-1)
!   
!   DALLOCATE( pp )
!   DALLOCATE( qq )
!   
!   Nqfms = lNqfms
!
!   SALLOCATE( pp, (0:Nqfms-1), lp(0:Nqfms-1) )
!   SALLOCATE( qq, (0:Nqfms-1), lq(0:Nqfms-1) )
!   
!   DALLOCATE( lp )
!   DALLOCATE( lq )
!   
!  enddo ! end of do iFarey ; 12/17/20;
!  
!  SALLOCATE( iota, (1:2*Nqfms-1), zero )
!
!  do ipq = 0, Nqfms-2
!   
!   iota(2*ipq+1) = ( pp(ipq+1) + goldmean * pp(ipq  ) ) / ( qq(ipq+1) + goldmean * qq(ipq  ) )
!   iota(2*ipq+2) = ( pp(ipq  ) + goldmean * pp(ipq+1) ) / ( qq(ipq  ) + goldmean * qq(ipq+1) )
!   
!  enddo
!
!1000 format("viscous   : "f10.2"s : "i4" : ("i3","i3" ) + ("i3","i3" ) = "f20.16" ; ")

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
!  do ipq = 1, 2*Nqfms-2
!   
!   xx = zero ; ruser(32*Mpol) = iota(ipq) ; ruser(32*Mpol+1) = kpert
!   
!   do mm = 1, Mpol ; lx = xx
!    
!    lx(mm) = xx(mm) - fdiff ; call Fmoser( Ndof, lx, laction,          iuser, ruser )
!    lx(mm) = xx(mm) + fdiff ; call Fmoser( Ndof, lx, uaction,          iuser, ruser )
!    ;                       ; call Dmoser( Ndof, xx,  action, daction, iuser, ruser )
!    
!    write(ounit,'("gradient =",f20.15," =",f20.15," ; err =",es10.2," ;")') daction(mm),  ( uaction - laction ) / ( two * fdiff ), &
!                                                                            daction(mm) - ( uaction - laction ) / ( two * fdiff )
!   enddo
!
!   pause
!
!  enddo

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  open( munit, file = ".moser", status='unknown', form='unformatted' )   

  write(munit) kpert, Mpol, 2*Nqfms-2, viscosity
  
  do ipq = 1, 2*Nqfms-2
   
   xx(1:Mpol) = zero ; ruser(32*Mpol) = iota(ipq)

  !cput = MPI_WTIME()
  !write(ounit,1100) cput-tstart, ipq, 2*Nqfms-2, iota(ipq)
   
   do ipert = 1, npert

    ruser(32*Mpol+1) = ipert * kpert / npert

    ifail = 1
    
    select case( Lmoser )    
    case( 1 ) ; call E04JYF( Ndof, ibound, Fmoser, bl, bu, xx, action,          iwork, Liwork, rwork, Lrwork, iuser, ruser, ifail )
    case( 2 ) ; call E04KYF( Ndof, ibound, Dmoser, bl, bu, xx, action, daction, iwork, Liwork, rwork, Lrwork, iuser, ruser, ifail )
    end select
    
    error = sqrt(  sum( daction(1:Mpol) * daction(1:Mpol) ) / Mpol ) / action

   !cput = MPI_WTIME()
   !write(ounit,1100) cput-tstart, ipq, 2*Nqfms-2, iota(ipq), ruser(32*Mpol+1), ifail

1100 format("viscous   : "f10.2"s : "i4" /"i4" : iota =",f19.16," ; ",:,"kpert =",f5.2," ; ",:,"|DA| =",es7.0," ; ",:,"ifail =",i3," ;")

   enddo ! end of do ipert ; 12/17/20;

   cput = MPI_WTIME()
   write(ounit,1100) cput-tstart, ipq, 2*Nqfms-2, iota(ipq), ruser(32*Mpol+1), error, ifail
   
   write(munit) ruser(      0: 8*Mpol-1)
   write(munit) ruser( 8*Mpol:16*Mpol-1)
   write(munit) ruser(16*Mpol:24*Mpol-1)
   write(munit) ruser(24*Mpol:32*Mpol-1)

!  do mm = 1, Mpol ; lx = xx
    
!   lx(mm) = xx(mm) - fdiff ; call Fmoser( Ndof, lx, laction,          iuser, ruser )
!   lx(mm) = xx(mm) + fdiff ; call Fmoser( Ndof, lx, uaction,          iuser, ruser )
!   ;                       ; call Dmoser( Ndof, xx,  action, daction, iuser, ruser )
    
!   write(ounit,'("gradient =",f20.15," =",f20.15," ; err =",es10.2," ;")') daction(mm),  ( uaction - laction ) / ( two * fdiff ), &
!                                                                           daction(mm) - ( uaction - laction ) / ( two * fdiff )
!  enddo

  enddo ! end of do ipq ; 12/17/20;

  close(munit)

   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( iwork  ) ) then
    write(0,'("macros   : 0123456789 : iwork  not already allocated ;")') 
    stop      'macros   : 0123456789 : iwork  not already allocated ;'
   endif
!#endif

   deallocate( iwork , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating iwork  ;")') 
    stop      'macros   : 0123456789 : error de-allocating iwork  ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rwork  ) ) then
    write(0,'("macros   : 0123456789 : rwork  not already allocated ;")') 
    stop      'macros   : 0123456789 : rwork  not already allocated ;'
   endif
!#endif

   deallocate( rwork , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rwork  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rwork  ;'
   endif
!#endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine viscous

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine Fmoser( Ndof, xx, action, iuser, ruser )
  
  use globals, only : ounit, zero, one, half, pi2, &
                      myid, tstart, &
                      Mpol, viscosity
  
  implicit none
  
  include "mpif.h"
  
  integer :: Ndof, iuser(1:1)
  real    :: xx(1:Ndof), action, ruser(0:32*Mpol+1), iota, kpert
  
  integer :: ii, mm
  real    :: teta, uu, up, ut
  
  iota = ruser(32*Mpol) ; kpert = ruser(32*Mpol+1)
  
  action = zero
  
  do ii = 0, 8 * Mpol-1
   
   teta = ii * pi2 / ( 8 * Mpol )
   
   uu = teta
   ut =  one
   up = teta + iota * pi2
   
   do mm = 1, Mpol
    
    uu = uu + xx(mm) * sin( mm *   teta                )
    ut = ut + xx(mm) * cos( mm *   teta                ) * mm
    up = up + xx(mm) * sin( mm * ( teta + iota * pi2 ) )
    
   enddo ! end of do mm ; 12/17/20;
   
   ruser(ii        ) =        uu
   ruser(ii+ 8*Mpol) = ( up - uu ) + kpert * sin( uu )
   ruser(ii+16*Mpol) =   up       
   ruser(ii+24*Mpol) = ( up - uu )
   
   action = action + ( - viscosity * log(ut) + half * ( up - uu ) * ( up - uu ) + kpert * cos( uu ) ) 
   
  enddo ! end of do ii ; 12/17/20;
  
  action = action * pi2 / ( 8 * Mpol )
  
  return
  
end subroutine Fmoser

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine Dmoser( Ndof, xx, action, daction, iuser, ruser )
  
  use globals, only : ounit, zero, one, half, pi2, &
                      myid, tstart, &
                      Mpol, viscosity
  
  implicit none
  
  include "mpif.h"
  
  integer :: Ndof, iuser(1:1)
  real    :: xx(1:Ndof), action, daction(1:Mpol), ruser(0:32*Mpol+1), iota, kpert
  
  integer :: ii, mm
  real    :: teta, uu, up, ut
  
  iota = ruser(32*Mpol) ; kpert = ruser(32*Mpol+1)
  
  action          = zero
  daction(1:Mpol) = zero
  
  do ii = 0, 8 * Mpol-1
   
   teta = ii * pi2 / ( 8 * Mpol )
   
   uu = teta
   ut =  one
   up = teta + iota * pi2
   
   do mm = 1, Mpol
    
    uu = uu + xx(mm) * sin( mm *   teta                )
    ut = ut + xx(mm) * cos( mm *   teta                ) * mm
    up = up + xx(mm) * sin( mm * ( teta + iota * pi2 ) )
    
   enddo ! end of 12/17/20;
   
   ruser(ii        ) =        uu
   ruser(ii+ 8*Mpol) = ( up - uu ) + kpert * sin( uu )
   ruser(ii+16*Mpol) =   up       
   ruser(ii+24*Mpol) = ( up - uu )
   
   action = action + ( - viscosity * log(ut) + half * ( up - uu ) * ( up - uu ) + kpert * cos( uu ) ) 
   
   do mm = 1, Mpol
    daction(mm) = daction(mm) + ( &
                                  - viscosity * cos( mm * teta ) * mm / ut &
                                  + ( up - uu ) * ( sin( mm * ( teta + iota * pi2) ) - sin( mm * teta ) ) &
                                  - kpert * sin( uu ) * sin( mm * teta ) &
                                  )

   enddo ! end of do mm ; 12/17/20;

  enddo ! end of do ii ; 12/17/20;
  
  action          = action          * pi2 / ( 8 * Mpol )

  daction(1:Mpol) = daction(1:Mpol) * pi2 / ( 8 * Mpol )
  
  return
  
end subroutine Dmoser

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
