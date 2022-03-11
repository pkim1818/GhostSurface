














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine savedata( lNi, lNj, lNk, lTijk, err, rnorm )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only: zero, half, pi2, ounit, tunit, tstart, myid, kperp, chi, imethod, &
                     infile, hmnfile, Tlow, Tupp, iNi, iNj, iNk, metrixsuff, ksuff, lsuff, &
                     Ni, Nj, dx, dy, &
                     px, py, lamda, mu, cc, e02defrwrk, e02defiwrk, &
                     dxydabf, T0, C0, D0, E0, T1, C1, D1, E1
  
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, intent(in)  :: lNi, lNj, lNk
  real   , intent(in)  :: lTijk(0:lNi-1,1:lNj-1,0:lNk-1), err, rnorm

  LOGICAL              :: Te_exist
  integer              :: ii, jj, astat, lwrk, liwrk, MX, MY, nx, ny, ie01daf, ie02def, MM, ierr
  real                 :: cput, lx(1:1), ly(1:1), wTijk, lerr
  real   , allocatable :: xx(:), yy(:), ff(:)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( imethod )
   
  case( 0 )
   
  case( 1 )
   
   open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Te", status='replace', form='unformatted' )
   write(tunit) iNi, iNj, iNk, kperp, err, rnorm, Tlow, Tupp, chi%low, chi%upp, lNi, lNj, lNk
   write(tunit) lTijk(0:lNi-1,1:lNj-1,0:lNk-1)
   close(tunit)
   
   inquire(file="."//trim(hmnfile)//"."//metrixsuff//"."//lsuff//".Te", exist=Te_exist )
   
   if( .not.Te_exist ) then ! construct initial guess for lower kperp calculation ; 12/17/20;
    
    open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//lsuff//".Te", status='replace', form='unformatted' )
    write(tunit) iNi, iNj, iNk, kperp, err, rnorm, Tlow, Tupp, chi%low, chi%upp, lNi, lNj, lNk
    write(tunit) lTijk(0:lNi-1,1:lNj-1,0:lNk-1)
    close(tunit)
    
   endif
   
  case( 2 )
   
  end select
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  select case( imethod )
   
  case( 0 )
   
  case( 1 )
   
   MX = Ni + 1 ; MY = Nj + 1 ; lwrk = ( MX + 6 ) * ( MY + 6 )
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( xx ) ) then
    write(0,'("macros    : 0123456789 : xx already allocated ;")') 
    stop      'macros    : 0123456789 : xx already allocated ;'
   endif
!#endif

   allocate( xx(1:MX), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating xx ;")') 
    stop      'macros   : 0123456789 : error allocating xx ;'
   endif
!#endif

   xx(1:MX) = (/ ( ii, ii = 0, Ni ) /) * dx           


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( yy ) ) then
    write(0,'("macros    : 0123456789 : yy already allocated ;")') 
    stop      'macros    : 0123456789 : yy already allocated ;'
   endif
!#endif

   allocate( yy(1:MY), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating yy ;")') 
    stop      'macros   : 0123456789 : error allocating yy ;'
   endif
!#endif

   yy(1:MY) = (/ ( jj, jj = 0, Nj ) /) * dy + chi%low 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ff ) ) then
    write(0,'("macros    : 0123456789 : ff already allocated ;")') 
    stop      'macros    : 0123456789 : ff already allocated ;'
   endif
!#endif

   allocate( ff(1:MX*MY), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ff ;")') 
    stop      'macros   : 0123456789 : error allocating ff ;'
   endif
!#endif

   ff(1:MX*MY) = zero 


   
   do ii = 0, Ni-1 ;    jj = 0       ; ff(1+jj+ii*(Nj+1)) = Tlow
    ;              ; do jj = 1, Nj-1 ; ff(1+jj+ii*(Nj+1)) = lTijk(ii,jj,0)
    ;              ; enddo    
    ;              ;    jj =    Nj   ; ff(1+jj+ii*(Nj+1)) = Tupp
   enddo
   ;  ii =    Ni   ;    jj = 0       ; ff(1+jj+ii*(Nj+1)) = Tlow
   ;               ; do jj = 1, Nj-1 ; ff(1+jj+ii*(Nj+1)) = lTijk( 0,jj,0)
   ;               ; enddo
   ;               ;    jj =    Nj   ; ff(1+jj+ii*(Nj+1)) = Tupp
   
  case( 2 ) 
   
  case default
   
  end select
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( lamda      ) ) then
    write(0,'("macros    : 0123456789 : lamda      already allocated ;")') 
    stop      'macros    : 0123456789 : lamda      already allocated ;'
   endif
!#endif

   allocate( lamda     (1:MX+4 ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating lamda      ;")') 
    stop      'macros   : 0123456789 : error allocating lamda      ;'
   endif
!#endif

   lamda     (1:MX+4 ) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( mu         ) ) then
    write(0,'("macros    : 0123456789 : mu         already allocated ;")') 
    stop      'macros    : 0123456789 : mu         already allocated ;'
   endif
!#endif

   allocate( mu        (1:MY+4 ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating mu         ;")') 
    stop      'macros   : 0123456789 : error allocating mu         ;'
   endif
!#endif

   mu        (1:MY+4 ) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( cc         ) ) then
    write(0,'("macros    : 0123456789 : cc         already allocated ;")') 
    stop      'macros    : 0123456789 : cc         already allocated ;'
   endif
!#endif

   allocate( cc        (1:MX*MY), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating cc         ;")') 
    stop      'macros   : 0123456789 : error allocating cc         ;'
   endif
!#endif

   cc        (1:MX*MY) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( e02defrwrk ) ) then
    write(0,'("macros    : 0123456789 : e02defrwrk already allocated ;")') 
    stop      'macros    : 0123456789 : e02defrwrk already allocated ;'
   endif
!#endif

   allocate( e02defrwrk(1:lwrk ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating e02defrwrk ;")') 
    stop      'macros   : 0123456789 : error allocating e02defrwrk ;'
   endif
!#endif

   e02defrwrk(1:lwrk ) = zero 


  
  ie01daf = 0
  
  call E01DAF( MX, MY, xx, yy, ff, px, py, lamda, mu, cc, e02defrwrk, ie01daf ) ! construct bi-cubic interpolation ; 20/10/21 ;
  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( xx  ) ) then
    write(0,'("macros   : 0123456789 : xx  not already allocated ;")') 
    stop      'macros   : 0123456789 : xx  not already allocated ;'
   endif
!#endif

   deallocate( xx , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating xx  ;")') 
    stop      'macros   : 0123456789 : error de-allocating xx  ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( yy  ) ) then
    write(0,'("macros   : 0123456789 : yy  not already allocated ;")') 
    stop      'macros   : 0123456789 : yy  not already allocated ;'
   endif
!#endif

   deallocate( yy , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating yy  ;")') 
    stop      'macros   : 0123456789 : error de-allocating yy  ;'
   endif
!#endif


  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( e02defrwrk  ) ) then
    write(0,'("macros   : 0123456789 : e02defrwrk  not already allocated ;")') 
    stop      'macros   : 0123456789 : e02defrwrk  not already allocated ;'
   endif
!#endif

   deallocate( e02defrwrk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating e02defrwrk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating e02defrwrk  ;'
   endif
!#endif


  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( e02defrwrk ) ) then
    write(0,'("macros    : 0123456789 : e02defrwrk already allocated ;")') 
    stop      'macros    : 0123456789 : e02defrwrk already allocated ;'
   endif
!#endif

   allocate( e02defrwrk(1:py-4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating e02defrwrk ;")') 
    stop      'macros   : 0123456789 : error allocating e02defrwrk ;'
   endif
!#endif

   e02defrwrk(1:py-4) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( e02defiwrk ) ) then
    write(0,'("macros    : 0123456789 : e02defiwrk already allocated ;")') 
    stop      'macros    : 0123456789 : e02defiwrk already allocated ;'
   endif
!#endif

   allocate( e02defiwrk(1:py-4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating e02defiwrk ;")') 
    stop      'macros   : 0123456789 : error allocating e02defiwrk ;'
   endif
!#endif

   e02defiwrk(1:py-4) = zero 


  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
! MM = 1 ; lx = zero ; ly = zero
  
! lerr = zero
  
! do ii = 0, Ni-1 ; lx = ii * dx
!  do jj = 0, Nj   ; ly = jj * dy + chi%low
    
!   ie02def = 0 
    
!   call E02DEF( MM, px, py, lx(1:MM), ly(1:MM), lamda(1:px), mu(1:py), cc(1:(px-4)*(py-4)), ff(1:MM), e02defrwrk(1:py-4), e02defiwrk(1:py-4), ie02def )
    
!   if( ii.lt.Ni ) then
!    if    ( jj.eq.0                ) then ; wTijk = Tlow
!    elseif( jj.gt.0 .and. jj.lt.Nj ) then ; wTijk = lTijk(ii,jj,0)
!    elseif(               jj.eq.Nj ) then ; wTijk = Tupp
!    endif
!   else
!    if    ( jj.eq.0                ) then ; wTijk = Tlow
!    elseif( jj.gt.0 .and. jj.lt.Nj ) then ; wTijk = lTijk( 0,jj,0)
!    elseif(               jj.eq.Nj ) then ; wTijk = Tupp
!    endif
!   endif
    
!   lerr = lerr + abs( wTijk-ff(1) )

!  enddo
! enddo
  
! cput = MPI_WTIME()
  
! write(ounit,'("savedata  : "f10.2"s : average interpolation error =",es08.1," ;")') cput-tstart, lerr / ((Ni+1)*(Nj+1))
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( imethod )
   
  case( 0 )
   
  case( 1 )
   
   MM = 1
   
   do ii = 0, Ni
    
    do jj = 0, Nj
     
     lx(1) = ii * dx ; ly(1) = chi%low + jj * dy
     
     ie02def = 0
     
     call E02DEF( MM, px, py, lx(1:MM), ly(1:MM), lamda(1:px), mu(1:py), cc(1:(px-4)*(py-4)), ff(1:MM), e02defrwrk(1:py-4), e02defiwrk(1:py-4), ie02def )
     
     T0(ii,jj) = ff(1) ! this should be exactly the original grid data ; 20/10/21 ;
     
     lx(1) = dxydabf(1,2,ii,jj, 0) ; lx(1) = dxydabf(2,2,ii,jj, 0)
     
     if( lx(1).gt.pi2 ) lx(1) = lx(1) - pi2
     
     ie02def = 0
     
     call E02DEF( MM, px, py, lx(1:MM), ly(1:MM), lamda(1:px), mu(1:py), cc(1:(px-4)*(py-4)), ff(1:MM), e02defrwrk(1:py-4), e02defiwrk(1:py-4), ie02def )
     
     T1(ii,jj) = ff(1)
     
    enddo ! end of do jj ; 20/10/21 ;
    
   enddo ! end of do ii ; 20/10/21 ;
   
   do ii = 0, Ni-1
    
    do jj = 0, Nj-1
     
     C0(ii,jj) = ( + T0(ii,jj) + T0(ii,jj+1) + T0(ii+1,jj) + T0(ii+1,jj+1) ) / ( 4      )
     D0(ii,jj) = ( - T0(ii,jj) - T0(ii,jj+1) + T0(ii+1,jj) + T0(ii+1,jj+1) ) / ( 2 * dx )
     E0(ii,jj) = ( - T0(ii,jj) + T0(ii,jj+1) - T0(ii+1,jj) + T0(ii+1,jj+1) ) / ( 2 * dy )
     
     C1(ii,jj) = ( + T1(ii,jj) + T1(ii,jj+1) + T1(ii+1,jj) + T1(ii+1,jj+1) ) / ( 4      )
     D1(ii,jj) = ( - T1(ii,jj) - T1(ii,jj+1) + T1(ii+1,jj) + T1(ii+1,jj+1) ) / ( 2 * dx )
     E1(ii,jj) = ( - T1(ii,jj) + T1(ii,jj+1) - T1(ii+1,jj) + T1(ii+1,jj+1) ) / ( 2 * dy )
     
    enddo ! end of do jj ; 20/10/21 ;
    
   enddo ! end of do ii ; 20/10/21 ;
   
  case( 2 )
   
  case default
   
  end select
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

end subroutine savedata

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
