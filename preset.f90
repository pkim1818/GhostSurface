














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine preset
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals
  
  implicit none
  
  include "mpif.h"

  LOGICAL              :: qfms_exist, pq_exist
  integer, parameter   :: Node = 6, NFT = 2
  integer              :: astat, ii, jj, kk, nn, mm, ie01bef, jk, lpp, lqq, ipq, iFarey, lNqfms, lNode, ifail, iostat, ierr, ie01bgf, MX, MY, lwrk, ie01daf
  integer              :: iuser(1:3), ic05rbf, Nii, id02bjf, iflag, lorder
  integer, allocatable :: lp(:), lq(:)
  real                 :: cput, cpul, xy(1:Node), Bxy(1:Node), teta, zeta, srho, stz(1:3), xyz(1:3,0:3), norm(1:3), sqrtg, lodetol
  real                 :: lorcmn, ldrcmn, lotsmn, ldtsmn
  real                 :: grads(1:3), bdotn, Temp(1:1), rhost, rhoend, d02bjfwk(1:20), rho, dTdz(1:1), lflux, ldtemp(0:1)
  real                 :: fluxteta(1:NFT,0:3), fvec(1:NFT), fjac(1:NFT,1:NFT), xtol, ruser(1:1), ltemp, dltemp, diff(0:4)
  real   , allocatable :: xx(:), yy(:), ff(:)
  real   , allocatable :: efmn(:), ofmn(:), cfmn(:), sfmn(:), ijreal(:), ijimag(:), fTT(:,:), dfT(:,:)

  EXTERNAL             :: dTemp, dTempout, D02BJX, D02BJW, invcoords

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cput = MPI_WTIME()
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Ni = 2**iNi ; Nj = 2**iNj ; Nk = 2**iNk

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  mn = 1 + Ntor + Mpol * ( 2 * Ntor + 1 ) ! Fourier resolution of temperature ; 02/19/2021 ;
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( im ) ) then
    write(0,'("macros    : 0123456789 : im already allocated ;")') 
    stop      'macros    : 0123456789 : im already allocated ;'
   endif
!#endif

   allocate( im(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating im ;")') 
    stop      'macros   : 0123456789 : error allocating im ;'
   endif
!#endif

   im(1:mn) = 0 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( in ) ) then
    write(0,'("macros    : 0123456789 : in already allocated ;")') 
    stop      'macros    : 0123456789 : in already allocated ;'
   endif
!#endif

   allocate( in(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating in ;")') 
    stop      'macros   : 0123456789 : error allocating in ;'
   endif
!#endif

   in(1:mn) = 0 


  
  ii = 0
  
  ;  mm = 0  
  ;do nn = 0, Ntor
  ; ii = ii+1 ; im(ii) = mm ; in(ii) = nn * nfp
  ;enddo
  ;
  
  do mm = 1, Mpol
   do nn = -Ntor, Ntor
    ii = ii+1 ; im(ii) = mm ; in(ii) = nn * nfp
   enddo
  enddo
  
  Nteta = 8 * max(Mpol,1) ; Nzeta = 8 * max(Ntor,1) ; Ntz = Nteta * Nzeta
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( cplxin  ) ) then
    write(0,'("macros    : 0123456789 : cplxin  already allocated ;")') 
    stop      'macros    : 0123456789 : cplxin  already allocated ;'
   endif
!#endif

   allocate( cplxin (1:Nteta,1:Nzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating cplxin  ;")') 
    stop      'macros   : 0123456789 : error allocating cplxin  ;'
   endif
!#endif

   cplxin (1:Nteta,1:Nzeta) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( cplxout ) ) then
    write(0,'("macros    : 0123456789 : cplxout already allocated ;")') 
    stop      'macros    : 0123456789 : cplxout already allocated ;'
   endif
!#endif

   allocate( cplxout(1:Nteta,1:Nzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating cplxout ;")') 
    stop      'macros   : 0123456789 : error allocating cplxout ;'
   endif
!#endif

   cplxout(1:Nteta,1:Nzeta) = zero 


  
  planb = fftw_plan_dft_2d( Nzeta, Nteta, cplxin, cplxout, fftw_backward, fftw_measure + fftw_destroy_input )
  planf = fftw_plan_dft_2d( Nzeta, Nteta, cplxin, cplxout, fftw_forward , fftw_measure + fftw_destroy_input )
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  qfmsmn = 1 + qfmsNtor + qfmsMpol * ( 2 * qfmsNtor + 1 ) ! Fourier resolution of qfm surfaces ; 02/19/2021 ;
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qfmsim ) ) then
    write(0,'("macros    : 0123456789 : qfmsim already allocated ;")') 
    stop      'macros    : 0123456789 : qfmsim already allocated ;'
   endif
!#endif

   allocate( qfmsim(1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qfmsim ;")') 
    stop      'macros   : 0123456789 : error allocating qfmsim ;'
   endif
!#endif

   qfmsim(1:qfmsmn) = 0 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qfmsin ) ) then
    write(0,'("macros    : 0123456789 : qfmsin already allocated ;")') 
    stop      'macros    : 0123456789 : qfmsin already allocated ;'
   endif
!#endif

   allocate( qfmsin(1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qfmsin ;")') 
    stop      'macros   : 0123456789 : error allocating qfmsin ;'
   endif
!#endif

   qfmsin(1:qfmsmn) = 0 


  
  ii = 0
  
  ;  mm = 0  
  ;do nn = 0, qfmsNtor
  ; ii = ii+1 ; qfmsim(ii) = mm ; qfmsin(ii) = nn * nfp
  ;enddo
  ;
  
  do mm = 1, qfmsMpol
   do nn = -qfmsNtor, qfmsNtor
    ii = ii+1 ; qfmsim(ii) = mm ; qfmsin(ii) = nn * nfp
   enddo
  enddo
  
  qfmsNteta = 8 * max(qfmsMpol,1) ; qfmsNzeta = 8 * max(qfmsNtor,1) ; qfmsNtz = qfmsNteta * qfmsNzeta
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qfmscplxin  ) ) then
    write(0,'("macros    : 0123456789 : qfmscplxin  already allocated ;")') 
    stop      'macros    : 0123456789 : qfmscplxin  already allocated ;'
   endif
!#endif

   allocate( qfmscplxin (1:qfmsNteta,1:qfmsNzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qfmscplxin  ;")') 
    stop      'macros   : 0123456789 : error allocating qfmscplxin  ;'
   endif
!#endif

   qfmscplxin (1:qfmsNteta,1:qfmsNzeta) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qfmscplxout ) ) then
    write(0,'("macros    : 0123456789 : qfmscplxout already allocated ;")') 
    stop      'macros    : 0123456789 : qfmscplxout already allocated ;'
   endif
!#endif

   allocate( qfmscplxout(1:qfmsNteta,1:qfmsNzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qfmscplxout ;")') 
    stop      'macros   : 0123456789 : error allocating qfmscplxout ;'
   endif
!#endif

   qfmscplxout(1:qfmsNteta,1:qfmsNzeta) = zero 


  
  qfmsplanb = fftw_plan_dft_2d( qfmsNzeta, qfmsNteta, qfmscplxin, qfmscplxout, fftw_backward, fftw_measure + fftw_destroy_input )
  qfmsplanf = fftw_plan_dft_2d( qfmsNzeta, qfmsNteta, qfmscplxin, qfmscplxout, fftw_forward , fftw_measure + fftw_destroy_input )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  select case( imethod )
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  case( 1 ) ! imethod = 1 ; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
   
   dx = pi2 / Ni ; dy = ( chi%upp - chi%low ) / Nj ; dz = pi2 / Nk ; dv = dx * dy * dz

   Ndof = Ni * (Nj-1) * Nk ! #degrees-of-freedom;

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
   if( allocated( T0 ) ) then
    write(0,'("macros    : 0123456789 : T0 already allocated ;")') 
    stop      'macros    : 0123456789 : T0 already allocated ;'
   endif
!#endif

   allocate( T0(0:Ni  ,0:Nj  ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating T0 ;")') 
    stop      'macros   : 0123456789 : error allocating T0 ;'
   endif
!#endif

   T0(0:Ni  ,0:Nj  ) = zero 

 ! full grid; regular plane; 20/10/21 ;
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( C0 ) ) then
    write(0,'("macros    : 0123456789 : C0 already allocated ;")') 
    stop      'macros    : 0123456789 : C0 already allocated ;'
   endif
!#endif

   allocate( C0(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating C0 ;")') 
    stop      'macros   : 0123456789 : error allocating C0 ;'
   endif
!#endif

   C0(0:Ni-1,0:Nj-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( D0 ) ) then
    write(0,'("macros    : 0123456789 : D0 already allocated ;")') 
    stop      'macros    : 0123456789 : D0 already allocated ;'
   endif
!#endif

   allocate( D0(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating D0 ;")') 
    stop      'macros   : 0123456789 : error allocating D0 ;'
   endif
!#endif

   D0(0:Ni-1,0:Nj-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( E0 ) ) then
    write(0,'("macros    : 0123456789 : E0 already allocated ;")') 
    stop      'macros    : 0123456789 : E0 already allocated ;'
   endif
!#endif

   allocate( E0(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating E0 ;")') 
    stop      'macros   : 0123456789 : error allocating E0 ;'
   endif
!#endif

   E0(0:Ni-1,0:Nj-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( T1 ) ) then
    write(0,'("macros    : 0123456789 : T1 already allocated ;")') 
    stop      'macros    : 0123456789 : T1 already allocated ;'
   endif
!#endif

   allocate( T1(0:Ni  ,0:Nj  ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating T1 ;")') 
    stop      'macros   : 0123456789 : error allocating T1 ;'
   endif
!#endif

   T1(0:Ni  ,0:Nj  ) = zero 

 ! full grid; mapped  plane; 20/10/21 ;
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( C1 ) ) then
    write(0,'("macros    : 0123456789 : C1 already allocated ;")') 
    stop      'macros    : 0123456789 : C1 already allocated ;'
   endif
!#endif

   allocate( C1(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating C1 ;")') 
    stop      'macros   : 0123456789 : error allocating C1 ;'
   endif
!#endif

   C1(0:Ni-1,0:Nj-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( D1 ) ) then
    write(0,'("macros    : 0123456789 : D1 already allocated ;")') 
    stop      'macros    : 0123456789 : D1 already allocated ;'
   endif
!#endif

   allocate( D1(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating D1 ;")') 
    stop      'macros   : 0123456789 : error allocating D1 ;'
   endif
!#endif

   D1(0:Ni-1,0:Nj-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( E1 ) ) then
    write(0,'("macros    : 0123456789 : E1 already allocated ;")') 
    stop      'macros    : 0123456789 : E1 already allocated ;'
   endif
!#endif

   allocate( E1(0:Ni-1,0:Nj-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating E1 ;")') 
    stop      'macros   : 0123456789 : error allocating E1 ;'
   endif
!#endif

   E1(0:Ni-1,0:Nj-1) = zero 


   
   xxn(0) = one ; xxn(1) = dx*xxn(0) ; xxn(2) = dx*xxn(1) ; xxn(3) = dx*xxn(2)
   yym(0) = one ; yym(1) = dy*yym(0) ; yym(2) = dy*yym(1) ; yym(3) = dy*yym(2)
   
   bcf( 1: 4) = one ; bcf( 5: 8) =               dx ; bcf( 9:12) =               dy ; bcf(13:16) =               dx*dy
   fdf( 1: 4) = one ; fdf( 5: 8) = one / eight / dx ; fdf( 9:12) = one / eight / dy ; fdf(13:16) = one / four / (dx*dy)
   
   BC(  1 , 1:16 ) = (/ +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   BC(  2 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   BC(  3 , 1:16 ) = (/ -3 , +3 , +0 , +0 , -2 , -1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   BC(  4 , 1:16 ) = (/ +2 , -2 , +0 , +0 , +1 , +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   BC(  5 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   BC(  6 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 /)
   BC(  7 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , -3 , +3 , +0 , +0 , -2 , -1 , +0 , +0 /)
   BC(  8 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +2 , -2 , +0 , +0 , +1 , +1 , +0 , +0 /)
   BC(  9 , 1:16 ) = (/ -3 , +0 , +3 , +0 , +0 , +0 , +0 , +0 , -2 , +0 , -1 , +0 , +0 , +0 , +0 , +0 /)
   BC( 10 , 1:16 ) = (/ +0 , +0 , +0 , +0 , -3 , +0 , +3 , +0 , +0 , +0 , +0 , +0 , -2 , +0 , -1 , +0 /)
   BC( 11 , 1:16 ) = (/ +9 , -9 , -9 , +9 , +6 , +3 , -6 , -3 , +6 , -6 , +3 , -3 , +4 , +2 , +2 , +1 /)
   BC( 12 , 1:16 ) = (/ -6 , +6 , +6 , -6 , -3 , -3 , +3 , +3 , -4 , +4 , -2 , +2 , -2 , -2 , -1 , -1 /)
   BC( 13 , 1:16 ) = (/ +2 , +0 , -2 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +1 , +0 , +0 , +0 , +0 , +0 /)
   BC( 14 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +2 , +0 , -2 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +1 , +0 /)
   BC( 15 , 1:16 ) = (/ -6 , +6 , +6 , -6 , -4 , -2 , +4 , +2 , -3 , +3 , -3 , +3 , -2 , -1 , -2 , -1 /)
   BC( 16 , 1:16 ) = (/ +4 , -4 , -4 , +4 , +2 , +2 , -2 , -2 , +2 , -2 , +2 , -2 , +1 , +1 , +1 , +1 /)
!                        1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
   FD(  1 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   FD(  2 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 /)
   FD(  3 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 , +0 /)
   FD(  4 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +0 , +1 , +0 , +0 , +0 , +0 , +0 /)
   FD(  5 , 1:16 ) = (/ -1 , +0 , +1 , +0 , -2 , +0 , +2 , +0 , -1 , +0 , +1 , +0 , +0 , +0 , +0 , +0 /)
   FD(  6 , 1:16 ) = (/ +0 , -1 , +0 , +1 , +0 , -2 , +0 , +2 , +0 , -1 , +0 , +1 , +0 , +0 , +0 , +0 /)
   FD(  7 , 1:16 ) = (/ +0 , +0 , +0 , +0 , -1 , +0 , +1 , +0 , -2 , +0 , +2 , +0 , -1 , +0 , +1 , +0 /)
   FD(  8 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , -1 , +0 , +1 , +0 , -2 , +0 , +2 , +0 , -1 , +0 , +1 /)
   FD(  9 , 1:16 ) = (/ -1 , -2 , -1 , +0 , +0 , +0 , +0 , +0 , +1 , +2 , +1 , +0 , +0 , +0 , +0 , +0 /)
   FD( 10 , 1:16 ) = (/ +0 , -1 , -2 , -1 , +0 , +0 , +0 , +0 , +0 , +1 , +2 , +1 , +0 , +0 , +0 , +0 /)
   FD( 11 , 1:16 ) = (/ +0 , +0 , +0 , +0 , -1 , -2 , -1 , +0 , +0 , +0 , +0 , +0 , +1 , +2 , +1 , +0 /)
   FD( 12 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , -1 , -2 , -1 , +0 , +0 , +0 , +0 , +0 , +1 , +2 , +1 /)
   FD( 13 , 1:16 ) = (/ +1 , +0 , -1 , +0 , +0 , +0 , +0 , +0 , -1 , +0 , +1 , +0 , +0 , +0 , +0 , +0 /)
   FD( 14 , 1:16 ) = (/ +0 , +1 , +0 , -1 , +0 , +0 , +0 , +0 , +0 , -1 , +0 , +1 , +0 , +0 , +0 , +0 /)
   FD( 15 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +1 , +0 , -1 , +0 , +0 , +0 , +0 , +0 , -1 , +0 , +1 , +0 /)
   FD( 16 , 1:16 ) = (/ +0 , +0 , +0 , +0 , +0 , +1 , +0 , -1 , +0 , +0 , +0 , +0 , +0 , -1 , +0 , +1 /)
   
   nperp(0:2) = (/ 5, 9, 9 /) ! #points used in perpendicular finite differencing for idiff = 0, 1, 2;
   
   npint(1:2) = (/ 4, 16 /) ! #points used in perpendicular interpolation for opint = 1, 2;
   
   ifd(0,1:5) = (/   0, - 1,   0, + 1,   0 /)
   jfd(0,1:5) = (/ - 1,   0,   0,   0,   1 /)
  
   wfd(0,1:5) = (/   0,   1, - 2,   1,   0 /) / dx / dx &
              + (/   1,   0, - 2,   0,   1 /) / dy / dy  

   
   ifd(1,1:9) = (/ - 1,   0,   1, - 1,   0,   1, - 1,   0,   1 /)
   jfd(1,1:9) = (/ - 1, - 1, - 1,   0,   0,   0,   1,   1,   1 /)
   
   wfd(1,1:9) = (/   1, - 2,   1,   2, - 4,   2,   1, - 2,   1 /) / dx / dx / four &
              + (/   1,   2,   1, - 2, - 4, - 2,   1,   2,   1 /) / dy / dy / four
   
   ifd(2,1:9) = (/   0,   0, - 2, - 1,   0,   1,   2,   0,   0 /)
   jfd(2,1:9) = (/ - 2, - 1,   0,   0,   0,   0,   0,   1,   2 /)
   
   wfd(2,1:9) = (/   0,   0, - 1,  16, -30,  16, - 1,   0,   0 /) / dx / dx / twelve &
              + (/ - 1,  16,   0,   0, -30,   0,   0,  16, - 1 /) / dy / dy / twelve
    
   iin(1,1: 4) = (/  0, 1, 0, 1 /) ! linear interpolation? ; 12/17/20;
   jin(1,1: 4) = (/  0, 0, 1, 1 /)
   
   iin(2,1:16) = (/ -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2 /) ! cubic interpolation? ; 12/17/20;
   jin(2,1:16) = (/ -1,-1,-1,-1,  0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2 /)
   
   wpd(1,-2,-2:2) = (/   0,   0,   0,   0,   0 /) * kpara / dz**2
   wpd(1,-1,-2:2) = (/   0,   1,   0,   0,   0 /) * kpara / dz**2
   wpd(1, 0,-2:2) = (/   0, - 1,   0, - 1,   0 /) * kpara / dz**2
   wpd(1, 1,-2:2) = (/   0,   0,   0,   1,   0 /) * kpara / dz**2
   wpd(1, 2,-2:2) = (/   0,   0,   0,   0,   0 /) * kpara / dz**2
   
   wpd(2,-2,-2:2) = (/ - 1,   0,   0,   0,   0 /) * kpara / dz**2 / 12
   wpd(2,-1,-2:2) = (/   0,  16,   0,   0,   0 /) * kpara / dz**2 / 12
   wpd(2, 0,-2:2) = (/   1, -16,   0, -16,   1 /) * kpara / dz**2 / 12
   wpd(2, 1,-2:2) = (/   0,   0,   0,  16,   0 /) * kpara / dz**2 / 12
   wpd(2, 2,-2:2) = (/   0,   0,   0,   0, - 1 /) * kpara / dz**2 / 12

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
   
  case( 2 ) ! imethod = 2; 20/10/21 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

   Ndof = Ni * Nj * (Nk-1) ! #degrees-of-freedom;

   dx = pi2 / Ni ; dy = pi2 / Nj ; dz = ( chi%upp - chi%low ) / Nk ; dv = dx * dy * dz
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Bx ) ) then
    write(0,'("macros    : 0123456789 : Bx already allocated ;")') 
    stop      'macros    : 0123456789 : Bx already allocated ;'
   endif
!#endif

   allocate( Bx(-1:Ni-1,-1:Nj-1,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Bx ;")') 
    stop      'macros   : 0123456789 : error allocating Bx ;'
   endif
!#endif

   Bx(-1:Ni-1,-1:Nj-1,0:Nk-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( By ) ) then
    write(0,'("macros    : 0123456789 : By already allocated ;")') 
    stop      'macros    : 0123456789 : By already allocated ;'
   endif
!#endif

   allocate( By(-1:Ni-1,-1:Nj-1,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating By ;")') 
    stop      'macros   : 0123456789 : error allocating By ;'
   endif
!#endif

   By(-1:Ni-1,-1:Nj-1,0:Nk-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Bz ) ) then
    write(0,'("macros    : 0123456789 : Bz already allocated ;")') 
    stop      'macros    : 0123456789 : Bz already allocated ;'
   endif
!#endif

   allocate( Bz(-1:Ni-1,-1:Nj-1,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Bz ;")') 
    stop      'macros   : 0123456789 : error allocating Bz ;'
   endif
!#endif

   Bz(-1:Ni-1,-1:Nj-1,0:Nk-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( BB ) ) then
    write(0,'("macros    : 0123456789 : BB already allocated ;")') 
    stop      'macros    : 0123456789 : BB already allocated ;'
   endif
!#endif

   allocate( BB(-1:Ni-1,-1:Nj-1,0:Nk-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating BB ;")') 
    stop      'macros   : 0123456789 : error allocating BB ;'
   endif
!#endif

   BB(-1:Ni-1,-1:Nj-1,0:Nk-1) = zero 



   do ii = -1, Ni-1 ; teta = (ii+half) * pi2                  / Ni
    do jj = -1, Nj-1 ; zeta = (jj+half) * pi2                  / Nj
     do kk =  0, Nk-1 ; srho = (kk+half) * (chi%upp - chi%low ) / Nk + chi%low
      
      xy(1:Node) = (/ teta, srho, zero, zero, zero, zero /)

      call bfield( zeta, xy(1:Node), Bxy(1:Node) )

      Bx(ii,jj,kk) = Bxy(1)
      By(ii,jj,kk) = Bzeta 
      Bz(ii,jj,kk) = Bxy(2)

      BB(ii,jj,kk) = Bx(ii,jj,kk)*Bx(ii,jj,kk) + By(ii,jj,kk)*By(ii,jj,kk) + Bz(ii,jj,kk)*Bz(ii,jj,kk)
   
     enddo
    enddo
   enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  case( 3 ) ! imethod = 3; 20/10/21 ;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
   
   Ndof = mn * ( Ni-1 ) ! #degrees-of-freedom ;
   
   dz = ( chi%upp - chi%low ) / Ni
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( fBs ) ) then
    write(0,'("macros    : 0123456789 : fBs already allocated ;")') 
    stop      'macros    : 0123456789 : fBs already allocated ;'
   endif
!#endif

   allocate( fBs(1:Ntz, 0:Ni-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating fBs ;")') 
    stop      'macros   : 0123456789 : error allocating fBs ;'
   endif
!#endif

   fBs(1:Ntz, 0:Ni-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( fBt ) ) then
    write(0,'("macros    : 0123456789 : fBt already allocated ;")') 
    stop      'macros    : 0123456789 : fBt already allocated ;'
   endif
!#endif

   allocate( fBt(1:Ntz, 0:Ni-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating fBt ;")') 
    stop      'macros   : 0123456789 : error allocating fBt ;'
   endif
!#endif

   fBt(1:Ntz, 0:Ni-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( fBz ) ) then
    write(0,'("macros    : 0123456789 : fBz already allocated ;")') 
    stop      'macros    : 0123456789 : fBz already allocated ;'
   endif
!#endif

   allocate( fBz(1:Ntz, 0:Ni-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating fBz ;")') 
    stop      'macros   : 0123456789 : error allocating fBz ;'
   endif
!#endif

   fBz(1:Ntz, 0:Ni-1) = zero 



    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( fBB ) ) then
    write(0,'("macros    : 0123456789 : fBB already allocated ;")') 
    stop      'macros    : 0123456789 : fBB already allocated ;'
   endif
!#endif

   allocate( fBB(1:Ntz, 0:Ni-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating fBB ;")') 
    stop      'macros   : 0123456789 : error allocating fBB ;'
   endif
!#endif

   fBB(1:Ntz, 0:Ni-1) = zero 


   
   do ii =  0, Ni-1    ; srho = ( ii + half ) * ( chi%upp - chi%low ) / Ni + chi%low
    do kk =  0, Nzeta-1 ; zeta =  kk           * pi2                   / Nzeta
     do jj =  0, Nteta-1 ; teta =  jj           * pi2                   / Nteta
      
      jk = 1 + jj + kk * Nteta

      xy(1:Node) = (/ teta, srho, zero, zero, zero, zero /)

      call bfield( zeta, xy(1:Node), Bxy(1:Node) )

      fBs(jk,ii) = Bxy(2)
      fBt(jk,ii) = Bxy(1)
      fBz(jk,ii) = Bzeta 
      
      fBB(jk,ii) = fBs(jk,ii) * fBs(jk,ii) + fBt(jk,ii) * fBt(jk,ii) + fBz(jk,ii) * fBz(jk,ii)
      
     enddo ! end of do jj ; 02/19/2021 ;
    enddo ! end of do kk ; 02/19/2021 ;
   enddo ! end of do ii ; 02/19/2021 ;
   
   fBB = one / fBB

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  case default

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  end select ! end select( imethod ) ; 0/10/21 ;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine preset

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
