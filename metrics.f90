














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine metrics
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, &
                              ounit, munit, &
                              myid, tstart, icheck, &
                              iNi, iNj, iNk, Ni, Nj, Nk, dx, dy, dz, metrixsuff, &
                              Bzeta, Bfh, dxydabf, &
                              hmnfile, chi, &
                              odetol, &
                              Lbck, Lfwd, Lfull, Lhalf, ijk, kbf
  
  implicit none

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  LOGICAL            :: metrics_exist
  integer, parameter :: Node = 6
  integer            :: ierr, astat, ios, liNi, liNj, liNk, lNi, lNj, lNk, ii, jj, kk, ifail, rfail
  real               :: cpui, cpuo, cput, lodetol, xx, yy, zz, xy(1:Node), Bxy(1:Node), zst, zed, d02bjfwk(1:20*Node)
  
  EXTERNAL           :: bfield, bfieldout, bfieldend
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cpui = MPI_WTIME() ; cpuo = cpui

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! SALLOCATE( Bfh    , (    -2:2,0:Ni-1,0:Nj,0:Nk-1), zero )
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Bfh     ) ) then
    write(0,'("macros    : 0123456789 : Bfh     already allocated ;")') 
    stop      'macros    : 0123456789 : Bfh     already allocated ;'
   endif
!#endif

   allocate( Bfh    (    -2:2,0:Ni  ,0:Nj,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Bfh     ;")') 
    stop      'macros   : 0123456789 : error allocating Bfh     ;'
   endif
!#endif

   Bfh    (    -2:2,0:Ni  ,0:Nj,0:Nk-1) = zero 


! SALLOCATE( dxydabf, (1:2,-4:4,0:Ni-1,0:Nj,0:Nk-1), zero ) ! field line mapped points d(x,y) / d(a,b) forward/backward on full grid;  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dxydabf ) ) then
    write(0,'("macros    : 0123456789 : dxydabf already allocated ;")') 
    stop      'macros    : 0123456789 : dxydabf already allocated ;'
   endif
!#endif

   allocate( dxydabf(1:2,-4:4,0:Ni  ,0:Nj,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dxydabf ;")') 
    stop      'macros   : 0123456789 : error allocating dxydabf ;'
   endif
!#endif

   dxydabf(1:2,-4:4,0:Ni  ,0:Nj,0:Nk-1) = zero 

 ! field line mapped points d(x,y) / d(a,b) forward/backward on full grid;  
 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  inquire( file = "."//trim(hmnfile)//"."//metrixsuff//".metrics", exist = metrics_exist )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if( metrics_exist ) then
   
   rfail = 0
   
   if( myid.eq.0 ) then

    cput = MPI_WTIME()

    open( munit, file = "."//trim(hmnfile)//"."//metrixsuff//".metrics", status='old', form='unformatted' )
    
    if( rfail.eq.0 ) then

     read(munit,iostat=ios) liNi, liNj, liNk, lodetol

     if( ios.ne.0 .or. liNi.ne.iNi .or. liNj.ne.iNj .or. liNk.ne.iNk ) rfail = 1

     if( lodetol.gt.odetol ) rfail = -1

    endif
    
    if( rfail.eq.0 ) then

     read( munit, iostat=ios ) Bfh

     if( ios.ne.0 ) rfail = 2

    endif
    
    if( rfail.eq.0 ) then

     read( munit, iostat=ios ) dxydabf

     if( ios.ne.0 ) rfail = 3

    endif
    
    close(munit)
    
   endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   call MPI_BCAST(rfail,1 ,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   cput = MPI_WTIME()
   
   if( myid.eq.0 ) write(ounit,1007) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//".metrics", rfail, lodetol, odetol, cput-cpuo

   cpuo = cput

1007 format("metrics   : "f10.2"s : "a" ; reading ; rfail ="i3" ; old odetol ="es8.1" ; odetol ="es8.1" ; time ="f10.2" ;")
   
   if( rfail.eq.0 ) then
    
!   RBCAST( Bfh    ,     5 * Ni * (Nj+1) * Nk )
    call MPI_BCAST(Bfh    ,5 * (Ni+1) * (Nj+1) * Nk ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
!   RBCAST( dxydabf, 2 * 9 * Ni * (Nj+1) * Nk )
    call MPI_BCAST(dxydabf,2 * 9 * (Ni+1) * (Nj+1) * Nk ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    
   else
    
    metrics_exist = .false.    
    
   endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
  endif ! end of if( metrics_exist ) ; 
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if( .not.metrics_exist ) then
   
   do kk = 0, Nk-1 ; ijk(3) = kk
    
!   do ii = 0, Ni-1 ; ijk(1) = ii
    do ii = 0, Ni   ; ijk(1) = ii
     do jj = 0, Nj   ; ijk(2) = jj

     !if( ind(ii,jj,kk).eq.0 ) cycle ! don't follow points not ind computational domain;
      
      if( jj.ge.0 ) then
       
       Lfull = .true. ; Lhalf = .not.Lfull
       
       xx = ii * dx ; yy = chi%low + jj * dy
       
       Lbck = .true.  ; Lfwd = .not.Lbck ; zst = kk*dz ; zed = zst - 2*dz ; kbf = 0

       xy(1:Node) = (/ xx , yy , one , zero , zero , one /)

       call bfield( zst, xy(1:Node), Bxy(1:Node) ) ! xy(3) = one / Bzeta ! ensures $\sqrt g = 1 / B^\phi$;
       
       dxydabf( 1:2, kbf, ii, jj, kk ) = xy(1:2)
       
       ifail = 0 ; call D02BJF( zst, zed, Node, xy(1:Node), bfield, odetol, 'D', bfieldout, bfieldend, d02bjfwk(1:20*Node), ifail )
       
       dxydabf( 1:2,  -4, ii, jj, kk ) = xy(1:2)
       
       Lbck = .false. ; Lfwd = .not.Lbck ; zst = kk*dz ; zed = zst + 2*dz ; kbf = 0

       xy(1:Node) = (/ xx , yy , one , zero , zero , one /)

       call bfield( zst, xy(1:Node), Bxy(1:Node) ) ! ; xy(3) = one / Bzeta ! ensures $\sqrt g = 1 / B^\phi$;

       ifail = 0 ; call D02BJF( zst, zed, Node, xy(1:Node), bfield, odetol, 'D', bfieldout, bfieldend, d02bjfwk(1:20*Node), ifail )
       
       dxydabf( 1:2,   4, ii, jj, kk ) = xy(1:2)
       
      endif ! end of if( jj.ge. 0 ) then ; Lfull=.true. ; Lhalf=.not.Lfull ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

     !if( jj.lt.Nj ) then

     ! Lfull=.false. ; Lhalf=.not.Lfull

     ! xx = dx*half + ii * dx ; yy = dy*half + chi%yylu(1) + jj * dy

     ! Lbck= .true. ; Lfwd=.not.Lbck ; zst = kk*dz ; zed = zst - dz*half ; kbf=0 ; xy = (/ xx , yy , one , zero , zero , one /)
       
     ! call bfield( zst , xy , Bxy ) ; xy(3) = one / Bzeta ! ensures $\sqrt g = 1 / B^\phi$;

     ! dxydabh( 1:2 ,  0  , kbf , ii , jj , kk ) = xy(1:2)
     ! dxydabh(  1  , 1:2 , kbf , ii , jj , kk ) = xy(3:4)
     ! dxydabh(  2  , 1:2 , kbf , ii , jj , kk ) = xy(5:6)
       
     ! ifail=0 ; call D02BJF( zst , zed , Node , xy(1:6) , bfield , odetol , 'D' , bfieldout , bfieldend , d02bjfwk , ifail )
       
     ! dxydabh( 1:2 ,  0  , -1    , ii , jj , kk ) = xy(1:2)
     ! dxydabh(  1  , 1:2 , -1    , ii , jj , kk ) = xy(3:4)
     ! dxydabh(  2  , 1:2 , -1    , ii , jj , kk ) = xy(5:6)

     ! Lbck=.false. ; Lfwd=.not.Lbck ; zst = kk*dz ; zed = zst + dz*half ; kbf=0 ; xy = (/ xx , yy , one , zero , zero , one /)
       
     ! call bfield( zst , xy , Bxy ) ; xy(3) = one / Bzeta ! ensures $\sqrt g = 1 / B^\phi$;

     ! ifail=0 ; call D02BJF( zst , zed , Node , xy(1:6) , bfield , odetol , 'D' , bfieldout , bfieldend , d02bjfwk , ifail )
        
     ! dxydabh( 1:2 ,  0  ,  1    , ii , jj , kk ) = xy(1:2)
     ! dxydabh(  1  , 1:2 ,  1    , ii , jj , kk ) = xy(3:4)
     ! dxydabh(  2  , 1:2 ,  1    , ii , jj , kk ) = xy(5:6)
       
     !endif ! end of if( jj.lt.Nj ) then ; Lfull=.true. ; Lhalf=.not.Lfull ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      
     enddo ! end of do jj
    enddo ! end of do ii

   enddo ! end of do kk

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   if( myid.eq.0 ) then

    cput = MPI_WTIME()

    open( munit, file = "."//trim(hmnfile)//"."//metrixsuff//".metrics", status='unknown', form='unformatted' )

    write(munit) iNi, iNj, iNk, odetol
    write(munit) Bfh ! Bfh was constructed in bfieldout ;
    write(munit) dxydabf
    close(munit)

   endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   cput = MPI_WTIME()
   
   write(ounit,1001) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//".metrics", odetol, cput-cpuo

1001 format("metrics   : "f10.2"s : "a" ; followed field lines ; odetol ="es9.2" ; time ="f10.2"s ;")
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
  endif
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine metrics

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

real function bfieldend( zz , xy )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  use globals, only : two, half, Lbck, Lfwd, Lfull, Lhalf, dz, ijk

  implicit none

  include "mpif.h"

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  integer, parameter  :: Node = 6

  real   , intent(in) :: zz , xy(1:Node)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  if( Lfull .and. Lbck ) bfieldend = zz - ( ijk(3) * dz - dz * two  )
  if( Lfull .and. Lfwd ) bfieldend = zz - ( ijk(3) * dz + dz * two  )

  if( Lhalf .and. Lbck ) bfieldend = zz - ( ijk(3) * dz - dz * half )
  if( Lhalf .and. Lfwd ) bfieldend = zz - ( ijk(3) * dz + dz * half )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  return

end function bfieldend

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine bfieldout( zz , xy )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  use globals, only : one, half, &
                      ounit, icheck, Lbck, Lfwd, Lfull, Lhalf, ijk, kbf, dz, Bfh, dxydabf, Bzeta

  implicit none

  include "mpif.h"

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  integer, parameter  :: Node = 6

  real, intent(inout) :: zz
  real, intent(in)    :: xy(1:Node)

  integer             :: ierr
  real                :: Bxy(1:Node)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  call bfield( zz, xy(1:Node), Bxy(1:Node) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  if( Lfull ) then
   
   dxydabf( 1:2, kbf, ijk(1), ijk(2), ijk(3) ) = xy(1:2) ! 12/17/20;
   
   if( abs(kbf).le.2 ) Bfh( kbf, ijk(1), ijk(2), ijk(3) ) = one / ( Bxy(1)**2 + Bxy(2)**2 + Bzeta**2 ) ! on toroidal half-grid;
   
  else

   if(.true. ) then
          write(ounit,'("bfieldout :    fatal    : .true.  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
   
  endif
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  if( Lbck ) then ; zz = zz - dz * half ; kbf = kbf - 1 ! set next intermediate output location;
  else            ; zz = zz + dz * half ; kbf = kbf + 1
  endif
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  return

end subroutine bfieldout

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
