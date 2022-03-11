














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine tsolve
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, pi2, ten, half, ounit, tunit, myid, numprocs, readinput, infile, hmnfile, &
                      skperp, siNi, siNj, siNk, sMpol, sNtor, &
                      ekperp, eiNi,             eMpol, eNtor, &
                      tstart, imethod, &
                      kperp, Tlow, Tupp, chi, &
                      Ni, Nj, Nk, Ndof, &
                      oNqfms, oqfmsmn, oqfmsim, oqfmsin, oqfmsNteta, oqfmsNzeta, oqfmsNtz, &
                      oflux, orcmn, orsmn, otcmn, otsmn, &
                      le04dgf, dTtol, ld02bjf, lrkutta, Ntau, Mtau, dtau, difforig, dodetol, dend, &
                      NPpts, Laction, &
                      PX, PY, lamda, mu, cc, e02defrwrk, e02defiwrk, &
                      dx, dy, dz, Mpol, Ntor, mn, im, in, Nteta, &
                      Linterpol, dcomp, LTprofile, itlimit, linesearch
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include "mpif.h"
  
  include 'fftw3.f03'
  
  LOGICAL              :: T_exist
  integer              :: iok, ierr, itau, jtau, ii, jj, kk, astat, ifd, liw, lrw, mm, nn, iuser(1:5), iter, ie04dgf, id02bjf, oNi, oNj, oNk, omn, il, lNi
  integer              :: okperp, MX, MY, lwrk, ie01daf, ie01bef, ie01bgf
  integer, allocatable :: iw(:)
  real                 :: cput, cpul, diff(0:4), diffl(0:4), diffu(0:4), fdiff, fdest, totaldiffusion, ruser(1:1), timestart, timeend, teta
  real                 :: stz(1:3), xyz(1:3,0:3), srho, lTmn, lDmn
  real                 :: oTlow, oTupp, olow, oupp, lodetol
  real, allocatable    :: TT(:,:,:), dT(:,:,:), TL(:,:,:), TU(:,:,:), rw(:), Tijk(:), dTijk(:), d02bjfwk(:)
  real, allocatable    :: fTT(:,:), dfT(:,:), fTL(:,:), fTU(:,:), oTT(:)
  real, allocatable    :: Mijk(:,:,:)
  real                 :: rho, lorcmn, lorsmn, lotcmn, lotsmn
  real, allocatable    :: ijreal(:), ijimag(:), efmn(:), ofmn(:), cfmn(:), sfmn(:)
  
  real   , allocatable :: xx(:), yy(:), ff(:)
  real   , allocatable :: Smn(:), Tmn(:,:),  Dmn(:,:)

  character            :: stext*28

  EXTERNAL             :: fdiffuse, gdiffuse, dfield, dfieldout, D02BJX, D02BJW
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  select case( imethod )
  
  case( 1 ) ! imethod = 1; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   cput = MPI_WTIME()

   write(ounit,'("tsolve    : "f10.1"s ; "a" ;")') cput-tstart, "calling metrics ;"
   
   call metrics ! performs fieldline tracing to construct metric elements of local field-alligned coordinates ; 12/17/20;
   
   call readdata( iok ) ! this might be able to go before reading metrics ; 
   
   call matrix ! construct the sparse matrix that represents the anisotropic diffusion operator ; 12/17/20;
   
   call relax ! given the matrices, iteratively solves sparse linear system ; 12/17/20;
   
   call volume
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
  case( 2 ) ! imethod = 2; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( TT ) ) then
    write(0,'("macros    : 0123456789 : TT already allocated ;")') 
    stop      'macros    : 0123456789 : TT already allocated ;'
   endif
!#endif

   allocate( TT(-1:Ni  ,-1:Nj  ,0:Nk  ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating TT ;")') 
    stop      'macros   : 0123456789 : error allocating TT ;'
   endif
!#endif

   TT(-1:Ni  ,-1:Nj  ,0:Nk  ) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dT ) ) then
    write(0,'("macros    : 0123456789 : dT already allocated ;")') 
    stop      'macros    : 0123456789 : dT already allocated ;'
   endif
!#endif

   allocate( dT( 0:Ni-1, 0:Nj-1,1:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dT ;")') 
    stop      'macros   : 0123456789 : error allocating dT ;'
   endif
!#endif

   dT( 0:Ni-1, 0:Nj-1,1:Nk-1) = zero 


   
   do kk = 0, Nk ; TT(0:Ni-1,0:Nj-1,kk) = Tlow + kk * ( Tupp - Tlow ) / Nk
   enddo
   
   inquire(file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//siNj//":"//siNk//".T", exist=T_exist)

   if( T_exist ) then
    
    cput = MPI_WTIME()
    write(ounit,'("tsolve    : "f10.1"s : "a" ;")') cput-tstart, "reading "//"."//trim(hmnfile)//"."//skperp//"."//siNi//":"//siNj//":"//siNk//".T"
    cpul = cput
    
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//siNj//":"//siNk//".T", status='old', form='unformatted' )
    read(tunit) okperp, oTlow, oTupp, olow, oupp, oNi, oNj, oNk
    read(tunit) TT(0:Ni-1,0:Nj-1,1:Nk-1)
    close(tunit)
    
   endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   cpul = MPI_WTIME()
   
   call diffusion( Ni, Nj, Nk, TT, diff, dT )
   
   difforig(0:2) = diff(0:2)

   cput = MPI_WTIME()

   write(ounit,2000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul

   cpul = cput
   
2000 format("tsolve    : "f10.1"s : total ="es13.05" ("es10.2" ) ; para ="es13.05" ("es10.2" ) ; perp ="es13.05" ("es10.2" ) ; |dD|/D ="es9.2" ; Q =",es13.5, &
            " ; "f10.1"s ;",:," ie04dgf ="i3" ; iter ="i6" ;")

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
!  SALLOCATE( TL, (-1:Ni,-1:Nj,0:Nk), zero )
!  SALLOCATE( TU, (-1:Ni,-1:Nj,0:Nk), zero )
   
!  do ii = 0, Ni-1
!   do jj = 0, Nj-1
!    do kk = 1, Nk-1
      
!     do ifd = -1, -6, -1
       
!      fdiff = ten**ifd
       
!      TL = TT ; TL(ii,jj,kk) = TT(ii,jj,kk) - fdiff * half ; call diffusion( Ni, Nj, Nk, TL, diffl, dT )
!      TU = TT ; TU(ii,jj,kk) = TT(ii,jj,kk) + fdiff * half ; call diffusion( Ni, Nj, Nk, TU, diffu, dT )
!      ;       ;                                            ; call diffusion( Ni, Nj, Nk, TT, diff , dT )
       
!      fdest = ( diffu(0) - diffl(0) ) / fdiff
       
!      write(ounit,'("(",3i6," ) : fdiff =",es8.1," ; deriv =",2es23.15," ; ratio =",f15.10," ;")') ii, jj, kk, fdiff, fdest, dT(ii,jj,kk), dT(ii,jj,kk)/fdest
       
!     enddo
      
!     pause
      
!    enddo
!   enddo
!  enddo
   
!  DALLOCATE( TL )
!  DALLOCATE( TU )
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   select case( le04dgf )
    
   case( 0 )
    
   case( 1 )
    
    iter = 0
    
    iuser(1:5) = (/ Ni, Nj, Nk, mn, 0 /) ; ruser(1:1) = one
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( iw ) ) then
    write(0,'("macros    : 0123456789 : iw already allocated ;")') 
    stop      'macros    : 0123456789 : iw already allocated ;'
   endif
!#endif

   allocate( iw(1:Ndof+1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating iw ;")') 
    stop      'macros   : 0123456789 : error allocating iw ;'
   endif
!#endif

   iw(1:Ndof+1) = 0 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rw ) ) then
    write(0,'("macros    : 0123456789 : rw already allocated ;")') 
    stop      'macros    : 0123456789 : rw already allocated ;'
   endif
!#endif

   allocate( rw(1:13*Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rw ;")') 
    stop      'macros   : 0123456789 : error allocating rw ;'
   endif
!#endif

   rw(1:13*Ndof) = zero 


    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Tijk ) ) then
    write(0,'("macros    : 0123456789 : Tijk already allocated ;")') 
    stop      'macros    : 0123456789 : Tijk already allocated ;'
   endif
!#endif

   allocate( Tijk(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Tijk ;")') 
    stop      'macros   : 0123456789 : error allocating Tijk ;'
   endif
!#endif

   Tijk(1:Ndof) = zero 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dTijk ) ) then
    write(0,'("macros    : 0123456789 : dTijk already allocated ;")') 
    stop      'macros    : 0123456789 : dTijk already allocated ;'
   endif
!#endif

   allocate( dTijk(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dTijk ;")') 
    stop      'macros   : 0123456789 : error allocating dTijk ;'
   endif
!#endif

   dTijk(1:Ndof) = zero 


    
    nn = 0
    do ii = 0, Ni-1
     do jj = 0, Nj-1
      do kk = 1, Nk-1 ; nn = nn + 1 ; Tijk(nn) = TT(ii,jj,kk)
      enddo
     enddo
    enddo
    
    call E04DKF( 'Nolist' )
    call E04DKF( 'Print Level = 0' )
    call E04DKF( 'Verify Objective Gradients = 0' )
    call E04DKF( 'Iteration Limit = 500' )
    call E04DKF( 'Function Precision = 1.0E-14' )
    call E04DKF( 'Linesearch Tolerance = 0.10' )
    call E04DKF( 'Optimality Tolerance = 1.0E-14' )
    
    ie04dgf = 1
    
    call E04DGF( Ndof, fdiffuse, iter, totaldiffusion, dTijk(1:Ndof), Tijk(1:Ndof), iw(1:Ndof+1), rw(1:13*Ndof), iuser(1:5), ruser(1:1), ie04dgf )
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( iw  ) ) then
    write(0,'("macros   : 0123456789 : iw  not already allocated ;")') 
    stop      'macros   : 0123456789 : iw  not already allocated ;'
   endif
!#endif

   deallocate( iw , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating iw  ;")') 
    stop      'macros   : 0123456789 : error de-allocating iw  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rw  ) ) then
    write(0,'("macros   : 0123456789 : rw  not already allocated ;")') 
    stop      'macros   : 0123456789 : rw  not already allocated ;'
   endif
!#endif

   deallocate( rw , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rw  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rw  ;'
   endif
!#endif


    
    nn = 0   
    do ii = 0, Ni-1
     do jj = 0, Nj-1
      do kk = 1, Nk-1 ; nn = nn + 1 ; TT(ii,jj,kk) = Tijk(nn)
      enddo
     enddo
    enddo
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( Tijk  ) ) then
    write(0,'("macros   : 0123456789 : Tijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : Tijk  not already allocated ;'
   endif
!#endif

   deallocate( Tijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating Tijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating Tijk  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( dTijk  ) ) then
    write(0,'("macros   : 0123456789 : dTijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : dTijk  not already allocated ;'
   endif
!#endif

   deallocate( dTijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating dTijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating dTijk  ;'
   endif
!#endif


    
    call diffusion( Ni, Nj, Nk, TT, diff, dT )
    
    cput = MPI_WTIME()
    write(ounit,2000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul, &
                      ie04dgf, iter
    cpul = cput
    
    difforig(0:2) = diff(0:2)
    
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//siNj//":"//siNk//".T", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, Nj, Nk
    write(tunit) TT(0:Ni-1,0:Nj-1,1:Nk-1)
    close(tunit)
    
   case default
    
   end select ! end of select case( le04dgf ) ; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   select case( ld02bjf )
    
   case( 0 )
    
   case( 1 )
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Tijk ) ) then
    write(0,'("macros    : 0123456789 : Tijk already allocated ;")') 
    stop      'macros    : 0123456789 : Tijk already allocated ;'
   endif
!#endif

   allocate( Tijk(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Tijk ;")') 
    stop      'macros   : 0123456789 : error allocating Tijk ;'
   endif
!#endif

   Tijk(1:Ndof) = zero 


    
    nn = 0
    do ii = 0, Ni-1
     do jj = 0, Nj-1
      do kk = 1, Nk-1 ; nn = nn + 1 ; Tijk(nn) = TT(ii,jj,kk)
      enddo
     enddo
    enddo
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( d02bjfwk ) ) then
    write(0,'("macros    : 0123456789 : d02bjfwk already allocated ;")') 
    stop      'macros    : 0123456789 : d02bjfwk already allocated ;'
   endif
!#endif

   allocate( d02bjfwk(1:20*Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating d02bjfwk ;")') 
    stop      'macros   : 0123456789 : error allocating d02bjfwk ;'
   endif
!#endif

   d02bjfwk(1:20*Ndof) = zero 


    
    lodetol = dodetol ; timestart = zero ; timeend = dend
    
    id02bjf = 0
    
    call D02BJF( timestart, timeend, Ndof, Tijk(1:Ndof), dfield, lodetol, 'D', dfieldout, D02BJW, d02bjfwk(1:20*Ndof), id02bjf )
    
    nn = 0   
    do ii = 0, Ni-1
     do jj = 0, Nj-1
      do kk = 1, Nk-1 ; nn = nn + 1 ; TT(ii,jj,kk) = Tijk(nn)
      enddo
     enddo
    enddo
    
    call diffusion( Ni, Nj, Nk, TT, diff, dT )
    
    cput = MPI_WTIME()
    write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul
    cpul = cput
    
    difforig(0:2) = diff(0:2)
    
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//siNj//":"//siNk//".T", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, Nj, Nk
    write(tunit) TT(0:Ni-1,0:Nj-1,1:Nk-1)
    close(tunit)
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( Tijk  ) ) then
    write(0,'("macros   : 0123456789 : Tijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : Tijk  not already allocated ;'
   endif
!#endif

   deallocate( Tijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating Tijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating Tijk  ;'
   endif
!#endif


    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( d02bjfwk  ) ) then
    write(0,'("macros   : 0123456789 : d02bjfwk  not already allocated ;")') 
    stop      'macros   : 0123456789 : d02bjfwk  not already allocated ;'
   endif
!#endif

   deallocate( d02bjfwk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating d02bjfwk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating d02bjfwk  ;'
   endif
!#endif


    
   case default
    
   end select ! end of select case( ld02bjf ) ; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
! the following is a hand-made descent algorithm and can be deleted ; 04/1/21 ;

!   cput = MPI_WTIME()
!   
!   cpul = cput
!   
!   do itau = 0, Ntau
!    
!    do jtau = 0, Mtau
!     
!     call diffusion( Ni, Nj, Nk, TT, diff, dT )
!     
!!    gradienterror = sqrt( sum( dT(0:Ni-1,0:Nj-1,1:Nk-1)*dT(0:Ni-1,0:Nj-1,1:Nk-1) ) / Ndof ) / diff(0)
!     
!     TT(0:Ni-1,0:Nj-1,1:Nk-1) = TT(0:Ni-1,0:Nj-1,1:Nk-1) - dtau * dT(0:Ni-1,0:Nj-1,1:Nk-1)
!     
!    enddo ! end of do jtau ; 20/10/21 ;
!    
!    cput = MPI_WTIME()
!
!    write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul
!
!    cpul = cput
!
!   enddo ! end of do itau ; 20/10/21 ;
   
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
   
!  select case( lrkutta )
!   
!  case( 0 )
!   
!  case( 1 )
!   
!   call gradient(        Ndof, xdof(0:Ndof), fdof(0:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)
!   
!   do irksave = 0, NRKsave-1 ! allows intermediate output;
!    
!    do irkstep = 1, NRKstep ! allows intermediate output;
!     
!     ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,1) * RKstep * half ; call gradient(        Ndof, ydof(0:Ndof), fdof(0:Ndof) ) ; fk(1:Ndof,2) = - fdof(1:Ndof)
!     ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,2) * RKstep * half ; call gradient(        Ndof, ydof(0:Ndof), fdof(0:Ndof) ) ; fk(1:Ndof,3) = - fdof(1:Ndof)
!     ydof(1:Ndof) = xdof(1:Ndof) + fk(1:Ndof,3) * RKstep        ; call gradient(        Ndof, ydof(0:Ndof), fdof(0:Ndof) ) ; fk(1:Ndof,4) = - fdof(1:Ndof)
!     
!     xdof(1:Ndof) = xdof(1:Ndof) + ( fk(1:Ndof,1) + 2 * fk(1:Ndof,2) + 2 * fk(1:Ndof,3) + fk(1:Ndof,4) ) * RKstep / six ! Runge-Kutta step;
!     
!     call gradient(        Ndof, xdof(0:Ndof), fdof(0:Ndof) ) ; fk(1:Ndof,1) = - fdof(1:Ndof)
!     
!     ferr = sqrt( sum(fdof(1:Ndof)*fdof(1:Ndof)) / Ndof ) / fdof(0)
!     
!    enddo ! end of do irkstep ;
!    
!   enddo ! end of do irksave ;
!   
!  case default
!   
!  end select
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
!  DALLOCATE( TT ) ! will probably need later for interpolations etc. ; 20/10/21 ;
    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( dT  ) ) then
    write(0,'("macros   : 0123456789 : dT  not already allocated ;")') 
    stop      'macros   : 0123456789 : dT  not already allocated ;'
   endif
!#endif

   deallocate( dT , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating dT  ;")') 
    stop      'macros   : 0123456789 : error de-allocating dT  ;'
   endif
!#endif


   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   MX = Ni + 1 ; MY = Nk + 1 ; lwrk = ( MX + 6 ) * ( MY + 6 )
   
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

   yy(1:MY) = (/ ( kk, kk = 0, Nk ) /) * dz + chi%low 


   
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


   
   do ii = 0, Ni-1 ;    kk = 0       ; ff(1+kk+ii*(Nk+1)) = Tlow
    ;              ; do kk = 1, Nk-1 ; ff(1+kk+ii*(Nk+1)) = TT(ii,0,kk)
    ;              ; enddo    
    ;              ;    kk =    Nk   ; ff(1+kk+ii*(Nk+1)) = Tupp
   enddo
   ;  ii =    Ni   ;    kk = 0       ; ff(1+kk+ii*(Nk+1)) = Tlow
   ;               ; do kk = 1, Nk-1 ; ff(1+kk+ii*(Nk+1)) = TT( 0,0,kk)
   ;               ; enddo
   ;               ;    kk =    Nk   ; ff(1+kk+ii*(Nk+1)) = Tupp
   
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
   
   call E01DAF( MX, MY, xx, yy, ff, PX, PY, lamda, mu, cc, e02defrwrk, ie01daf ) ! construct bi-cubic interpolation ; 20/10/21 ;
   
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

   allocate( e02defrwrk(1:PY-4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating e02defrwrk ;")') 
    stop      'macros   : 0123456789 : error allocating e02defrwrk ;'
   endif
!#endif

   e02defrwrk(1:PY-4) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( e02defiwrk ) ) then
    write(0,'("macros    : 0123456789 : e02defiwrk already allocated ;")') 
    stop      'macros    : 0123456789 : e02defiwrk already allocated ;'
   endif
!#endif

   allocate( e02defiwrk(1:PY-4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating e02defiwrk ;")') 
    stop      'macros   : 0123456789 : error allocating e02defiwrk ;'
   endif
!#endif

   e02defiwrk(1:PY-4) = zero 


   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
  case( 3 ) ! imethod = 3; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( fTT ) ) then
    write(0,'("macros    : 0123456789 : fTT already allocated ;")') 
    stop      'macros    : 0123456789 : fTT already allocated ;'
   endif
!#endif

   allocate( fTT(1:mn,0:Ni  ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating fTT ;")') 
    stop      'macros   : 0123456789 : error allocating fTT ;'
   endif
!#endif

   fTT(1:mn,0:Ni  ) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dfT ) ) then
    write(0,'("macros    : 0123456789 : dfT already allocated ;")') 
    stop      'macros    : 0123456789 : dfT already allocated ;'
   endif
!#endif

   allocate( dfT(1:mn,1:Ni-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dfT ;")') 
    stop      'macros   : 0123456789 : error allocating dfT ;'
   endif
!#endif

   dfT(1:mn,1:Ni-1) = zero 


   
   do ii = 0, Ni ; fTT(1,ii) = Tlow + ii * ( Tupp - Tlow ) / Ni
   enddo
   
   inquire(file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".F", exist=T_exist )
   
   if( T_exist ) then
    
    cput = MPI_WTIME()
    write(ounit,'("tsolve    : "f10.1"s : "a" ;")') cput-tstart, "reading "//"."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".F"
    cpul = cput
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oTT ) ) then
    write(0,'("macros    : 0123456789 : oTT already allocated ;")') 
    stop      'macros    : 0123456789 : oTT already allocated ;'
   endif
!#endif

   allocate( oTT(0:Ni), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oTT ;")') 
    stop      'macros   : 0123456789 : error allocating oTT ;'
   endif
!#endif

   oTT(0:Ni) = zero 


    
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".F", status='old', form='unformatted' )
    read(tunit) okperp, oTlow, oTupp, olow, oupp, oNi, omn
    if(oNi.ne.Ni ) then
          write(ounit,'("tsolve :    fatal    : oNi.ne.Ni  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
    do jj = 1, omn
     read(tunit) mm, nn, oTT(0:oNi)
     do ii = 1, mn
      if( im(ii).eq.mm .and. in(ii).eq.nn ) fTT(ii,0:Ni) = oTT(0:Ni)
     enddo
    enddo
    close(tunit)
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( oTT  ) ) then
    write(0,'("macros   : 0123456789 : oTT  not already allocated ;")') 
    stop      'macros   : 0123456789 : oTT  not already allocated ;'
   endif
!#endif

   deallocate( oTT , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating oTT  ;")') 
    stop      'macros   : 0123456789 : error de-allocating oTT  ;'
   endif
!#endif



   else

    cput = MPI_WTIME()
    write(ounit,'("tsolve    : "f10.1"s : "a" ;")') cput-tstart, "does not exist "//"."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".F"
    cpul = cput
        
   endif ! end of if( T_exist ) ; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   cpul = MPI_WTIME()
   
   call fdiffusion( mn, Ni, fTT(1:mn,0:Ni), diff, dfT(1:mn,1:Ni-1) )
   
   difforig(0:2) = diff(0:2)
   
   cput = MPI_WTIME()

!  write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul
   write(ounit,1001) cput-tstart, diff(0), diff(0)-difforig(0),                                                                      diff(4), cput-cpul
   
1000 format("tsolve    : "f10.1"s : total ="es13.05" ("es10.2" ) ; para ="es13.05" ("es10.2" ) ; perp ="es13.05" ("es10.2" ) ; |dD|/D ="es9.2" ; Q =",es13.5, &
            " ; "f10.1"s ;",:," ie04dgf ="i3" ; iter ="i6" ;")
1001 format("tsolve    : "f10.1"s : total ="es13.05" ("es10.2" ) ;       "  13x  "  "  10x "           "  13x  "  "  10x "   ; |dD|/D ="es9.2" ; Q =",es13.5, &
            " ; "f10.1"s ;",:," ie04dgf ="i3" ; iter ="i6" ;")

   dcomp(0:2,2) = diff(0:2)

   cpul = cput

   if( diff(3).lt.dTtol ) le04dgf = 0 ! initial guess (from previous calculation) is sufficiently converged ; 04/1/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
! compare analytic derivatives to finite-difference estimate ; 04/1/21 ;

!   SALLOCATE( fTU, (1:mn,0:Ni  ), zero )
!   SALLOCATE( fTL, (1:mn,0:Ni  ), zero )
!   
!   do ii = 1, Ni-1
!    do jj = 1, mn
!     
!     do ifd = -2, -4, -2
!      
!      fdiff = ten**ifd
!      
!      fTL = fTT ; fTL(jj,ii) = fTT(jj,ii) - fdiff * half ; call fdiffusion( mn, Ni, fTL(1:mn,0:Ni), diffl, dfT(1:mn,1:Ni-1) )
!      fTU = fTT ; fTU(jj,ii) = fTT(jj,ii) + fdiff * half ; call fdiffusion( mn, Ni, fTU(1:mn,0:Ni), diffu, dfT(1:mn,1:Ni-1) )
!      ;                                                  ; call fdiffusion( mn, Ni, fTT(1:mn,0:Ni), diff , dfT(1:mn,1:Ni-1) )
!      
!      fdest = ( diffu(0) - diffl(0) ) / fdiff
!      
!      write(ounit,'("(",2i6," ) : fdiff =",es8.1," ; est, analytic =",2es23.15," ; ratio =",f15.10," ;")') ii, jj, fdiff, fdest, dfT(jj,ii), fdest/dfT(jj,ii)
!      
!     enddo
!     
!     pause
!     
!    enddo
!   enddo
!   
!   DALLOCATE( fTL )
!   DALLOCATE( fTU )
!   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   select case( le04dgf )
    
   case( 0 ) ! le04dgf = 0 ; 20/10/21 ;
    
   case( 1 ) ! le04dgf = 1 ; 20/10/21 ;
    
    iter = 0
    
    iuser(1:5) = (/ Ni, Nj, Nk, mn, 0 /) ; ruser(1:1) = one
        
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( iw ) ) then
    write(0,'("macros    : 0123456789 : iw already allocated ;")') 
    stop      'macros    : 0123456789 : iw already allocated ;'
   endif
!#endif

   allocate( iw(1:  Ndof+1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating iw ;")') 
    stop      'macros   : 0123456789 : error allocating iw ;'
   endif
!#endif

   iw(1:  Ndof+1) = 0 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rw ) ) then
    write(0,'("macros    : 0123456789 : rw already allocated ;")') 
    stop      'macros    : 0123456789 : rw already allocated ;'
   endif
!#endif

   allocate( rw(1:13*Ndof  ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rw ;")') 
    stop      'macros   : 0123456789 : error allocating rw ;'
   endif
!#endif

   rw(1:13*Ndof  ) = zero 


    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Tijk ) ) then
    write(0,'("macros    : 0123456789 : Tijk already allocated ;")') 
    stop      'macros    : 0123456789 : Tijk already allocated ;'
   endif
!#endif

   allocate( Tijk(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Tijk ;")') 
    stop      'macros   : 0123456789 : error allocating Tijk ;'
   endif
!#endif

   Tijk(1:Ndof) = zero 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dTijk ) ) then
    write(0,'("macros    : 0123456789 : dTijk already allocated ;")') 
    stop      'macros    : 0123456789 : dTijk already allocated ;'
   endif
!#endif

   allocate( dTijk(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dTijk ;")') 
    stop      'macros   : 0123456789 : error allocating dTijk ;'
   endif
!#endif

   dTijk(1:Ndof) = zero 


    
    nn = 0
    do ii = 1, Ni-1
     do jj = 1, mn ; nn = nn + 1 ; Tijk(nn) = fTT(jj,ii) ! this "packing" seems clumsy ; 20/10/21 ;
     enddo
    enddo

    call E04DKF( 'Nolist')
    call E04DKF( 'Print Level = 0')
    call E04DKF( 'Verify Objective Gradients = 0')
!   call E04DKF( 'Iteration Limit = 100')
    write(stext,'("Iteration Limit      ="i6  )') itlimit    ; call E04DKF(stext)
    call E04DKF( 'Function Precision = 1.0E-14')
!   call E04DKF( 'Linesearch Tolerance = 0.50')
    write(stext,'("Linesearch Tolerance ="f6.2)') linesearch ; call E04DKF(stext)
    call E04DKF( 'Optimality Tolerance = 1.0E-14')
    
    ie04dgf = 1
    
    call E04DGF( Ndof, gdiffuse, iter, totaldiffusion, dTijk(1:Ndof), Tijk(1:Ndof), iw(1:Ndof+1), rw(1:13*Ndof), iuser(1:5), ruser(1:1), ie04dgf )
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( iw  ) ) then
    write(0,'("macros   : 0123456789 : iw  not already allocated ;")') 
    stop      'macros   : 0123456789 : iw  not already allocated ;'
   endif
!#endif

   deallocate( iw , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating iw  ;")') 
    stop      'macros   : 0123456789 : error de-allocating iw  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rw  ) ) then
    write(0,'("macros   : 0123456789 : rw  not already allocated ;")') 
    stop      'macros   : 0123456789 : rw  not already allocated ;'
   endif
!#endif

   deallocate( rw , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rw  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rw  ;'
   endif
!#endif


    
    nn = 0   
    do ii = 1, Ni-1
     do jj = 1, mn ; nn = nn + 1 ; fTT(jj,ii) = Tijk(nn) ! this "unpacking" seems clumsy ; 20/10/21 ;
     enddo
    enddo
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( Tijk  ) ) then
    write(0,'("macros   : 0123456789 : Tijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : Tijk  not already allocated ;'
   endif
!#endif

   deallocate( Tijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating Tijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating Tijk  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( dTijk  ) ) then
    write(0,'("macros   : 0123456789 : dTijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : dTijk  not already allocated ;'
   endif
!#endif

   deallocate( dTijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating dTijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating dTijk  ;'
   endif
!#endif


    
    call fdiffusion( mn, Ni, fTT, diff, dfT )
    
    cput = MPI_WTIME()

    write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul, &
                      ie04dgf, iter

    cpul = cput
    
    difforig(0:2) = diff(0:2)
    
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".F", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, mn, Mpol, Ntor
    do ii = 1, mn ; write(tunit) im(ii), in(ii), fTT(ii,0:Ni)
    enddo
    close(tunit)
  
    dcomp(0:2,2) = diff(0:2)

   case default ! le04dgf ;
    
   end select ! end of select case( le04dgf ) ; 20/10/21 ;
   
   write(ounit,'("tsolve    : ",11x," : reduction =",f10.6," ;")') (dcomp(0,0)-dcomp(0,1))/(dcomp(0,0)-dcomp(0,2))
   write(ounit,'("tsolve    : ",11x," : reduction =",f10.6," ;")') (dcomp(1,0)-dcomp(1,1))/(dcomp(1,0)-dcomp(1,2))
   write(ounit,'("tsolve    : ",11x," : reduction =",f10.6," ;")') (dcomp(2,0)-dcomp(2,1))/(dcomp(2,0)-dcomp(2,2))

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!  DALLOCATE( fTT ) ! will probably need later for interpolations etc. ; 20/10/21 ;
    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( dfT  ) ) then
    write(0,'("macros   : 0123456789 : dfT  not already allocated ;")') 
    stop      'macros   : 0123456789 : dfT  not already allocated ;'
   endif
!#endif

   deallocate( dfT , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating dfT  ;")') 
    stop      'macros   : 0123456789 : error de-allocating dfT  ;'
   endif
!#endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! save output to file ; 04/1/21 ;

   inquire(file="."//trim(hmnfile)//"."//ekperp//"."//siNi//":"//sMpol//":"//sNtor//".F", exist=T_exist )
   
   if( .not.T_exist ) then    
    open( tunit, file="."//trim(hmnfile)//"."//ekperp//"."//siNi//":"//sMpol//":"//sNtor//".F", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, mn, Mpol, Ntor
    do ii = 1, mn ; write(tunit) im(ii), in(ii), fTT(ii,0:Ni)
    enddo
    close(tunit)
   endif

   inquire(file="."//trim(hmnfile)//"."//skperp//"."//eiNi//":"//sMpol//":"//sNtor//".F", exist=T_exist )
   
   if( .not.T_exist ) then
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//eiNi//":"//sMpol//":"//sNtor//".F", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, 2*Ni, mn, Mpol, Ntor
    do ii = 1, mn ; write(tunit) im(ii), in(ii), ( fTT(ii,jj), (fTT(ii,jj)+fTT(ii,jj+1))*half, jj = 0, Ni-1 ), fTT(ii,Ni)
    enddo
    close(tunit)
   endif

   inquire(file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//eMpol//":"//eNtor//".F", exist=T_exist )
   
   if( .not.T_exist ) then
    open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//eMpol//":"//eNtor//".F", status='replace', form='unformatted' )
    write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, mn, Mpol, Ntor
    do ii = 1, mn ; write(tunit) im(ii), in(ii), fTT(ii,0:Ni)
    enddo
    close(tunit)
   endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   cput = MPI_WTIME()
   
   write(ounit,'("tsolve    : "f10.1"s : radial interpolation of Fourier harmonics ;")') cput-tstart
   
   cpul = cput
   
   lNi = Ni+1

    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Smn ) ) then
    write(0,'("macros    : 0123456789 : Smn already allocated ;")') 
    stop      'macros    : 0123456789 : Smn already allocated ;'
   endif
!#endif

   allocate( Smn(0:Ni     ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Smn ;")') 
    stop      'macros   : 0123456789 : error allocating Smn ;'
   endif
!#endif

   Smn(0:Ni     ) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Tmn ) ) then
    write(0,'("macros    : 0123456789 : Tmn already allocated ;")') 
    stop      'macros    : 0123456789 : Tmn already allocated ;'
   endif
!#endif

   allocate( Tmn(0:Ni,1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Tmn ;")') 
    stop      'macros   : 0123456789 : error allocating Tmn ;'
   endif
!#endif

   Tmn(0:Ni,1:mn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Dmn ) ) then
    write(0,'("macros    : 0123456789 : Dmn already allocated ;")') 
    stop      'macros    : 0123456789 : Dmn already allocated ;'
   endif
!#endif

   allocate( Dmn(0:Ni,1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Dmn ;")') 
    stop      'macros   : 0123456789 : error allocating Dmn ;'
   endif
!#endif

   Dmn(0:Ni,1:mn) = zero 


   
   do ii = 0, Ni ; Smn(ii) = chi%low + ii * dz
    do jj = 1, mn ; Tmn(ii,jj) = fTT(jj,ii)
    enddo
   enddo
   
   do ii = 1, mn ; ie01bef = 0 ; call E01BEF( lNi, Smn(0:Ni), Tmn(0:Ni,ii), Dmn(0:Ni,ii), ie01bef )
   enddo
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
   cput = MPI_WTIME()
   
   cpul = cput
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Mijk ) ) then
    write(0,'("macros    : 0123456789 : Mijk already allocated ;")') 
    stop      'macros    : 0123456789 : Mijk already allocated ;'
   endif
!#endif

   allocate( Mijk(0:Nteta,0:Ni,0:4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Mijk ;")') 
    stop      'macros   : 0123456789 : error allocating Mijk ;'
   endif
!#endif

   Mijk(0:Nteta,0:Ni,0:4) = zero 


   
   do jj = 0, Ni    ; srho = chi%low + jj * ( chi%upp - chi%low ) / Ni
    
    do ii = 0, Nteta ; teta = ii * pi2 / Nteta
     
     stz(1:3) = (/ srho, teta, zero /) ; call stzxyz( stz(1:3), xyz(1:3,0:3), Linterpol )
     
     Mijk(ii,jj,0:4) = (/ zero, teta, srho, xyz(1,0), xyz(2,0) /)

     do kk = 1, mn

      ie01bgf = 1 ; mm = 1 ; call E01BGF( lNi, Smn(0:Ni), Tmn(0:Ni,kk), Dmn(0:Ni,kk), mm, xyz(1,0), lTmn, lDmn, ie01bgf )

      Mijk(ii,jj,0) = Mijk(ii,jj,0) + lTmn * cos( im(kk) * xyz(2,0) - in(kk) * xyz(3,0) )

     enddo ! end of do kk Fourier summation ; 2020/02/16 ;

    enddo ! end of do ii ; 20/10/21 ;
    
   enddo ! end of do jj ; 20/10/21 ;
   
   open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".M", status='replace', form='unformatted' )
   write(tunit) Nteta, Ni
   write(tunit) Mijk(0:Nteta,0:Ni,0:4)
   close(tunit)
   
    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( Mijk  ) ) then
    write(0,'("macros   : 0123456789 : Mijk  not already allocated ;")') 
    stop      'macros   : 0123456789 : Mijk  not already allocated ;'
   endif
!#endif

   deallocate( Mijk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating Mijk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating Mijk  ;'
   endif
!#endif



   cput = MPI_WTIME()
   
   write(ounit,'("tsolve    : ",f10.1,"s : mapped Temp to chaotic coordinates ; ",f10.2,"s ;")') cput-tstart, cput-cpul
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  case default ! imethod ; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
  end select ! end select case( imethod ) ; 20/10/21 ;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
9998 continue
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine  

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine dfield( time, Tijk, dTijk )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : ounit, tstart, &
                      Ni, Nj, Nk, Ndof
  
  implicit none
  
  include "mpif.h"
  
  real               :: time
  real               :: Tijk(1:Ndof), dTijk(1:Ndof)
  
  integer            :: ii, jj, kk, nn
  real               :: TT(-1:Ni,-1:Nj,0:Nk), diff(0:4), dT(0:Ni-1,0:Nj-1,1:Nk-1)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nn = 0
  do ii = 0, Ni-1
   do jj = 0, Nj-1
    do kk = 1, Nk-1 ; nn = nn + 1 ; TT(ii,jj,kk) = Tijk(nn)
    enddo
   enddo
  enddo

  call diffusion( Ni, Nj, Nk, TT, diff, dT )
  
  nn = 0
  do ii = 0, Ni-1
   do jj = 0, Nj-1
    do kk = 1, Nk-1 ; nn = nn + 1 ; dTijk(nn) = - dT(ii,jj,kk)
    enddo
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine dfield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine dfieldout( time, Tijk )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : ounit, tunit, hmnfile, metrixsuff, msuff, tstart, &
                      Ni, Nj, Nk, Ndof, &
                      dend, difforig, Ntau, &
                      kperp, Tlow, Tupp, chi
  
  implicit none
  
  include "mpif.h"

  real               :: time
  real               :: Tijk(1:Ndof), dTijk(1:Ndof)
 
  integer            :: ii, jj, kk, nn
  real               :: cput, cpul, TT(-1:Ni,-1:Nj,0:Nk), diff(0:4), dT(0:Ni-1,0:Nj-1,1:Nk-1)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nn = 0
  do ii = 0, Ni-1
   do jj = 0, Nj-1
    do kk = 1, Nk-1 ; nn = nn + 1 ; TT(ii,jj,kk) = Tijk(nn)
    enddo
   enddo
  enddo

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  call diffusion( Ni, Nj, Nk, TT, diff, dT )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cput = MPI_WTIME()
  write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4)
  cpul = cput

1000 format("dfieldout : "f10.1"s : total ="es13.05" ("es10.2" ) ; para ="es13.05" ("es10.2" ) ; perp ="es13.05" ("es10.2" ) ; |dD|/D ="es9.2" ; Q =",es13.5, &
            " ; "f10.1"s ;")
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//msuff//".T", status='replace', form='unformatted' )
  write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, Nj, Nk
  write(tunit) TT(0:Ni-1,0:Nj-1,1:Nk-1)
  close(tunit)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  time = time + dend / Ntau

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine dfieldout

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine fdiffuse( mode, Ndof, Tijk, totaldiffusion, dTijk, nstate, iuser, ruser )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only :
  
  implicit none
  
  include "mpif.h"
  
  integer            :: mode, Ndof, nstate, iuser(1:5)
  real               :: Tijk(1:Ndof), totaldiffusion, dTijk(1:Ndof), ruser(1:1)
  
  integer            :: Ni, Nj, Nk, ii, jj, kk, nn, mn
  real               :: TT(-1:iuser(1),-1:iuser(2),0:iuser(3)), diff(0:4), dT(0:iuser(1)-1,0:iuser(2)-1,1:iuser(3)-1)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Ni = iuser(1) ; Nj = iuser(2) ; Nk = iuser(3) ; mn = iuser(4) ; iuser(5) = iuser(5) + 1
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nn = 0
  do ii = 0, Ni-1
   do jj = 0, Nj-1
    do kk = 1, Nk-1 ; nn = nn + 1 ; TT(ii,jj,kk) = Tijk(nn)
    enddo
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call diffusion( Ni, Nj, Nk, TT, diff, dT )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  totaldiffusion = diff(0)
  
  if( mode.le.0 ) return
  
  nn = 0
  do ii = 0, Ni-1
   do jj = 0, Nj-1
    do kk = 1, Nk-1 ; nn = nn + 1 ; dTijk(nn) = dT(ii,jj,kk)
    enddo
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine fdiffuse

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine gdiffuse( mode, Ndof, Tijk, totaldiffusion, dTijk, nstate, iuser, ruser )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only :
  
  implicit none
  
  include "mpif.h"
  
  integer            :: mode, Ndof, nstate, iuser(1:5)
  real               :: Tijk(1:Ndof), totaldiffusion, dTijk(1:Ndof), ruser(1:1)
  
  integer            :: Ni, Nj, Nk, ii, jj, kk, nn, mn
  real               :: fTT(1:iuser(4),0:iuser(1)), diff(0:4), dfT(1:iuser(4),1:iuser(1)-1)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Ni = iuser(1) ; Nj = iuser(2) ; Nk = iuser(3) ; mn = iuser(4) ; iuser(5) = iuser(5) + 1
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nn = 0
  do ii = 1, Ni-1
   do jj = 1, mn ; nn = nn + 1 ; fTT(jj,ii) = Tijk(nn)
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call fdiffusion( mn, Ni, fTT, diff, dfT )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  totaldiffusion = diff(0)
  
  if( mode.le.0 ) return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nn = 0
  do ii = 1, Ni-1
   do jj = 1, mn ; nn = nn + 1 ; dTijk(nn) = dfT(jj,ii)
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine gdiffuse

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
