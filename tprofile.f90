














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine tprofile( itprofile ) 
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals
  
  implicit none
  
  include "mpif.h"

  LOGICAL              :: qfms_exist, pq_exist
  integer, parameter   :: Node = 6, NFT = 2
  integer              :: astat, ii, jj, kk, nn, mm, ie01bef, jk, lpp, lqq, ipq, iFarey, lNqfms, lNode, ifail, iostat, ierr, ie01bgf, MX, MY, lwrk, ie01daf
  integer              :: iuser(1:3), ic05rbf, Nii, id02bjf, iflag, lorder, itprofile
  integer, allocatable :: lp(:), lq(:)
  real                 :: cput, cpul, xy(1:Node), Bxy(1:Node), teta, zeta, srho, stz(1:3), xyz(1:3,0:3), norm(1:3), sqrtg, lodetol
  real                 :: lorcmn, ldrcmn, lotsmn, ldtsmn
  real                 :: grads(1:3), bdotn, Temp(1:1), rhost, rhoend, d02bjfwk(1:20), rho, dTdz(1:1), lflux, ldtemp(0:1)
  real                 :: fluxteta(1:NFT,0:3), fvec(1:NFT), fjac(1:NFT,1:NFT), xtol, ruser(1:1), ltemp, dltemp, diff(0:4)
  real   , allocatable :: xx(:), yy(:), ff(:)
  real   , allocatable :: efmn(:), ofmn(:), cfmn(:), sfmn(:), ijreal(:), ijimag(:), fTT(:,:), dfT(:,:)

  EXTERNAL             :: dTemp, dTempout, D02BJX, D02BJW, invcoords

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cpul = MPI_WTIME()
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Tempprofile ) ) then
    write(0,'("macros    : 0123456789 : Tempprofile already allocated ;")') 
    stop      'macros    : 0123456789 : Tempprofile already allocated ;'
   endif
!#endif

   allocate( Tempprofile(0:Ni,0:4), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Tempprofile ;")') 
    stop      'macros   : 0123456789 : error allocating Tempprofile ;'
   endif
!#endif

   Tempprofile(0:Ni,0:4) = zero 


  
  rhost = chi%low ; rhoend = chi%upp ; lNode = 1 ; lodetol = odetol ; Temp(1) = Tlow ; ijk(1:3) = (/ 0, 0, 0 /)
  
  id02bjf = 0 ; call D02BJF( rhost, rhoend, lNode, Temp, dTemp, lodetol, 'D', dTempout, D02BJW, d02bjfwk, id02bjf ) ! need to record intermediate output;
  
  Tempprofile(0:Ni,1:2) = Tempprofile(0:Ni,1:2) / Tempprofile(Ni,1) ! rescaling temperature profile ; 02/19/2021 ;
  
  cput = MPI_WTIME()

  select case( id02bjf )
  case( 0 ) 
  case default ; write(ounit,'("tprofile : "f10.1"s : error integrating temperature profile ;"f10.2"s ;")') cput-tstart, cput-cpul
  end select

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  open(tunit, file = "."//trim(hmnfile)//".profile"  , status='unknown', form='unformatted')
  write(tunit) Ni
  write(tunit) Tempprofile
  close(tunit)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  Nii = Ni + 1 ; ie01bef = 0 ; call E01BEF( Nii, Tempprofile(0:Ni,0), Tempprofile(0:Ni,1), Tempprofile(0:Ni,2), ie01bef ) ! interpolate temperature ; 04/1/21 ;
  
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


  
  do ii = 0, Ni ; fTT(1,ii) = Tlow + ii * ( Tupp - Tlow ) / Ni ! 04/1/21 ;
  enddo
  
  call fdiffusion( mn, Ni, fTT, diff, dfT ) ! this is the "only-perpendicular-diffusion limit" ; 02/19/2021 ;
  
 !diff(0:2) = half * tkperp * pi2 * pi2 / tkperp

  difforig(0:2) = diff(0:2)
  
  cput = MPI_WTIME()
  
  if( itprofile.eq.1 ) then
!  write(ounit,1000) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul
   write(ounit,1001) cput-tstart, diff(0), diff(0)-difforig(0)                                                                              , cput-cpul
  endif

1000 format("tprofile  : "f10.1"s : total ="es13.05" ("es10.2" ) ; para ="es13.05" ("es10.2" ) ; perp ="es13.05" ("es10.2" ) ; |dD|/D ="es9.2" ; o =",es13.5, &
            " ; "f10.1"s ;")

1001 format("tprofile  : "f10.1"s : total ="es13.05" ("es10.2" ) ; "f10.1"s ;")

  dcomp(0:2,0) = diff(0:2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( efmn ) ) then
    write(0,'("macros    : 0123456789 : efmn already allocated ;")') 
    stop      'macros    : 0123456789 : efmn already allocated ;'
   endif
!#endif

   allocate( efmn(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating efmn ;")') 
    stop      'macros   : 0123456789 : error allocating efmn ;'
   endif
!#endif

   efmn(1:mn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ofmn ) ) then
    write(0,'("macros    : 0123456789 : ofmn already allocated ;")') 
    stop      'macros    : 0123456789 : ofmn already allocated ;'
   endif
!#endif

   allocate( ofmn(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ofmn ;")') 
    stop      'macros   : 0123456789 : error allocating ofmn ;'
   endif
!#endif

   ofmn(1:mn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( cfmn ) ) then
    write(0,'("macros    : 0123456789 : cfmn already allocated ;")') 
    stop      'macros    : 0123456789 : cfmn already allocated ;'
   endif
!#endif

   allocate( cfmn(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating cfmn ;")') 
    stop      'macros   : 0123456789 : error allocating cfmn ;'
   endif
!#endif

   cfmn(1:mn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( sfmn ) ) then
    write(0,'("macros    : 0123456789 : sfmn already allocated ;")') 
    stop      'macros    : 0123456789 : sfmn already allocated ;'
   endif
!#endif

   allocate( sfmn(1:mn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating sfmn ;")') 
    stop      'macros   : 0123456789 : error allocating sfmn ;'
   endif
!#endif

   sfmn(1:mn) = zero 


  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ijreal ) ) then
    write(0,'("macros    : 0123456789 : ijreal already allocated ;")') 
    stop      'macros    : 0123456789 : ijreal already allocated ;'
   endif
!#endif

   allocate( ijreal(1:Ntz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ijreal ;")') 
    stop      'macros   : 0123456789 : error allocating ijreal ;'
   endif
!#endif

   ijreal(1:Ntz) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ijimag ) ) then
    write(0,'("macros    : 0123456789 : ijimag already allocated ;")') 
    stop      'macros    : 0123456789 : ijimag already allocated ;'
   endif
!#endif

   allocate( ijimag(1:Ntz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ijimag ;")') 
    stop      'macros   : 0123456789 : error allocating ijimag ;'
   endif
!#endif

   ijimag(1:Ntz) = zero 


  
  xtol = odetol ; mm = 1 ; ltemp = zero
  
  do ii = 0, Ni ! map to background coordinates; need to invert ; 04/1/21 ;
   
   ijreal = zero ; ijimag = zero
   
   do kk = 0, Nzeta-1
    
    fluxteta(1:NFT,0:3) = zero
    
    do jj = 0, Nteta-1 ; jk = 1 + jj + kk * Nteta ; iuser(1:3) = (/ ii, jj, kk /)
     
     if( jj.eq.0 ) fluxteta(1:NFT,0) = (/ oflux(1) + ii * ( oflux(oNqfms)-oflux(1) ) / Ni, jj * pi2 / Nteta /) ! set initial guess ; 02/19/2021 ;
     
     if( jj.eq.1 ) fluxteta(1:NFT,0) =     fluxteta(1:NFT,1) + (/ 0, 1 /) * pi2 / Nteta
     if( jj.eq.2 ) fluxteta(1:NFT,0) = 2 * fluxteta(1:NFT,1) - 1 * fluxteta(1:NFT,2)
     if( jj.ge.3 ) fluxteta(1:NFT,0) = 3 * fluxteta(1:NFT,1) - 3 * fluxteta(1:NFT,2) + 1 * fluxteta(1:NFT,3)
     
     ic05rbf = 1
     
     call C05RBF( invcoords, NFT, fluxteta(1:NFT,0), fvec(1:NFT), fjac(1:NFT,1:NFT), xtol, iuser(1:3), ruser(1:1), ic05rbf )
     
     fluxteta(1:NFT,3) = fluxteta(1:NFT,2)
     fluxteta(1:NFT,2) = fluxteta(1:NFT,1)
     fluxteta(1:NFT,1) = fluxteta(1:NFT,0)
     
1100 format("tprofile  : ", 10x,"s : ",3i6," (",f20.16,",",f20.16," ) : fvec =",2es23.15," ;")
     
     lflux = fluxteta(1,0)
     
     ie01bgf = 1
     
     call E01BGF( Nii, Tempprofile(0:Ni,0), Tempprofile(0:Ni,1), Tempprofile(0:Ni,2), mm, lflux, ltemp, dltemp, ie01bgf )
     
     ijreal(jk) = ltemp
     ijimag(jk) =  zero
     
    enddo ! end of do jj ; 02/19/2021 ;
    
   enddo ! end of do kk ;02/19/2021 ;
   
   call tfft( Nteta, Nzeta, ijreal(1:Ntz), ijimag(1:Ntz), mn, im(1:mn), in(1:mn), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn) )
   
   fTT(1:mn,ii) = efmn(1:mn)
   
  enddo ! end of do ii ; 02/19/2021 ;
  
  cput = MPI_WTIME()
  
! write(ounit,'("tprofile  : ",f10.1,"s : inverted  coordinates ;",f10.2,"s ;")') cput-tstart, cput-cpul
  
  open( tunit, file="."//trim(hmnfile)//"."//skperp//"."//siNi//":"//sMpol//":"//sNtor//".E", status='replace', form='unformatted' )
  write(tunit) kperp, Tlow, Tupp, chi%low, chi%upp, Ni, mn, Mpol, Ntor
  do ii = 1, mn ; write(tunit) im(ii), in(ii), fTT(ii,0:Ni)
  enddo
  close(tunit)
  
  call fdiffusion( mn, Ni, fTT, diff, dfT )
  
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


  
  cput = MPI_WTIME()
  
! write(ounit,1002) cput-tstart, diff(0), diff(0)-difforig(0), diff(1), diff(1)-difforig(1), diff(2), diff(2)-difforig(2), diff(3), diff(4), cput-cpul
  write(ounit,1003) cput-tstart, diff(0), diff(0)-difforig(0),                                                                               cput-cpul
  
1002 format("tprofile  : "f10.1"s : total ="es13.05" ("es10.2" ) ; para ="es13.05" ("es10.2" ) ; perp ="es13.05" ("es10.2" ) ; |dD|/D ="es9.2" ; q =",es13.5, &
            " ; "f10.1"s ;")
1003 format("tprofile  : "f10.1"s : total ="es13.05" ("es10.2" ) ; "f10.1"s ;")
  
  dcomp(0:2,1) = diff(0:2)
  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( fTT ) ) then
    write(0,'("macros   : 0123456789 : fTT not already allocated ;")') 
    stop      'macros   : 0123456789 : fTT not already allocated ;'
   endif
!#endif

   deallocate( fTT, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating fTT ;")') 
    stop      'macros   : 0123456789 : error de-allocating fTT ;'
   endif
!#endif


  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( efmn ) ) then
    write(0,'("macros   : 0123456789 : efmn not already allocated ;")') 
    stop      'macros   : 0123456789 : efmn not already allocated ;'
   endif
!#endif

   deallocate( efmn, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating efmn ;")') 
    stop      'macros   : 0123456789 : error de-allocating efmn ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( ofmn ) ) then
    write(0,'("macros   : 0123456789 : ofmn not already allocated ;")') 
    stop      'macros   : 0123456789 : ofmn not already allocated ;'
   endif
!#endif

   deallocate( ofmn, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating ofmn ;")') 
    stop      'macros   : 0123456789 : error de-allocating ofmn ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( cfmn ) ) then
    write(0,'("macros   : 0123456789 : cfmn not already allocated ;")') 
    stop      'macros   : 0123456789 : cfmn not already allocated ;'
   endif
!#endif

   deallocate( cfmn, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating cfmn ;")') 
    stop      'macros   : 0123456789 : error de-allocating cfmn ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( sfmn ) ) then
    write(0,'("macros   : 0123456789 : sfmn not already allocated ;")') 
    stop      'macros   : 0123456789 : sfmn not already allocated ;'
   endif
!#endif

   deallocate( sfmn, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating sfmn ;")') 
    stop      'macros   : 0123456789 : error de-allocating sfmn ;'
   endif
!#endif


  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( ijreal ) ) then
    write(0,'("macros   : 0123456789 : ijreal not already allocated ;")') 
    stop      'macros   : 0123456789 : ijreal not already allocated ;'
   endif
!#endif

   deallocate( ijreal, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating ijreal ;")') 
    stop      'macros   : 0123456789 : error de-allocating ijreal ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( ijimag ) ) then
    write(0,'("macros   : 0123456789 : ijimag not already allocated ;")') 
    stop      'macros   : 0123456789 : ijimag not already allocated ;'
   endif
!#endif

   deallocate( ijimag, stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating ijimag ;")') 
    stop      'macros   : 0123456789 : error de-allocating ijimag ;'
   endif
!#endif


  
  if(.not.allocated(Tempprofile)) stop 'Tempprofile not allocated'
	deallocate(Tempprofile,stat=astat)
	if(astat.ne.0) stop 'error deallocating Tempprofile'

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine tprofile

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine dTemp( lflux, Temp, dTdz )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, pi2, ounit, tkperp, tstart, &
                      Bzeta, &
                      ijoteta, ijozeta, &
                      Linterpol, &
                      oNqfms, oqfmsmn, oqfmsim, oqfmsin, oqfmsNteta, oqfmsNzeta, oqfmsNtz, oflux, orcmn, drcmn, otsmn, dtsmn, &
                      qflux, gflux
  
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, parameter   :: lNode = 1, Node = 6
  real   , intent(in)  :: lflux, Temp(1:lNode)
  real   , intent(out) :: dTdz(1:lNode)
  
  integer              :: jj, kk, ii, ie01bgf, mm, jk, ifail, imn, ierr, iqfms
  real                 :: cput, teta, zeta, stz(1:3), xyz(1:3,0:3), xy(1:Node), Bxy(1:Node), norm(1:3), grads(1:3), bdotn, sqrtg
  real                 :: lorcmn, ldrcmn, lotsmn, ldtsmn
  real                 :: efmn(1:oqfmsmn), ofmn(1:oqfmsmn), cfmn(1:oqfmsmn), sfmn(1:oqfmsmn)
  real                 :: efnm(1:oqfmsmn), ofnm(1:oqfmsmn), cfnm(1:oqfmsmn), sfnm(1:oqfmsmn)
  real                 :: femn(1:oqfmsmn), fomn(1:oqfmsmn), fcmn(1:oqfmsmn), fsmn(1:oqfmsmn)
  real                 :: fenm(1:oqfmsmn), fonm(1:oqfmsmn), fcnm(1:oqfmsmn), fsnm(1:oqfmsmn)
  real                 :: ijreal(1:oqfmsNtz), ijimag(1:oqfmsNtz)
  real                 :: jireal(1:oqfmsNtz), jiimag(1:oqfmsNtz)
  real                 :: jkreal(1:oqfmsNtz), jkimag(1:oqfmsNtz)
  real                 :: kjreal(1:oqfmsNtz), kjimag(1:oqfmsNtz)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  qflux = zero
  gflux = zero
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( Linterpol )
   
  case( 1 ) ! Linterpol = 1 ; 02/19/2021 ;
   
   ifail = 1
   do iqfms = 1, oNqfms-1
    if( lflux.ge.oflux(iqfms) .and. lflux.le.oflux(iqfms+1) ) then
     do imn = 1, oqfmsmn
      efmn(imn) = orcmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      efnm(imn) =                                               ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      ofmn(imn) =  zero
      ofnm(imn) =  zero
      cfmn(imn) =  zero
      cfnm(imn) =  zero
      sfmn(imn) = otsmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      sfnm(imn) =                                               ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
     enddo ! end of do imn ; 02/19/2021 ;
     ifail = 0
     exit
    endif
    if( ifail.eq.0 ) exit
   enddo ! end of do iqfms ; 02/19/2021 ;
  
   if(ifail.ne.0 ) then
          write(ounit,'("initializ :    fatal    : ifail.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif

  case( 2 ) ! Linterpol = 2 ; 02/19/2021 ;
   
   do imn = 1, oqfmsmn ; mm = 1
    
    ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), orcmn(1:oNqfms,imn), drcmn(1:oNqfms,imn), mm, lflux, lorcmn, ldrcmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), orsmn(1:oNqfms,imn), drsmn(1:oNqfms,imn), mm, lflux, lorsmn, ldrsmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), otcmn(1:oNqfms,imn), dtcmn(1:oNqfms,imn), mm, lflux, lotcmn, ldtcmn, ie01bgf )
    ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), otsmn(1:oNqfms,imn), dtsmn(1:oNqfms,imn), mm, lflux, lotsmn, ldtsmn, ie01bgf )
    
    efmn(imn) = lorcmn ; ofmn(imn) = zero ; cfmn(imn) = zero ; sfmn(imn) = lotsmn
    efnm(imn) = ldrcmn ; ofnm(imn) = zero ; cfnm(imn) = zero ; sfnm(imn) = ldtsmn
    
   enddo
   
  case default ! Linterpol ; 02/19/2021 ;

   if(.true. ) then
          write(ounit,'("initialize :    fatal    : .true.  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif

  end select ! end select( Linterpol ) ; 02/19/2021 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call oqfmsinvfft( oqfmsmn, oqfmsim(1:oqfmsmn), oqfmsin(1:oqfmsmn), efmn(1:oqfmsmn), ofmn(1:oqfmsmn), cfmn(1:oqfmsmn), sfmn(1:oqfmsmn), &
                    oqfmsNteta, oqfmsNzeta, ijreal(1:oqfmsNtz), ijimag(1:oqfmsNtz) )

  ijimag = ijoteta + ijimag

  call oqfmsinvfft( oqfmsmn, oqfmsim(1:oqfmsmn), oqfmsin(1:oqfmsmn), efnm(1:oqfmsmn), ofnm(1:oqfmsmn), cfnm(1:oqfmsmn), sfnm(1:oqfmsmn), &
                    oqfmsNteta, oqfmsNzeta, jireal(1:oqfmsNtz), jiimag(1:oqfmsNtz) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  fomn = - oqfmsim * efmn ; femn = zero
  fcmn = + oqfmsim * sfmn ; fsmn = zero

  call oqfmsinvfft( oqfmsmn, oqfmsim(1:oqfmsmn), oqfmsin(1:oqfmsmn), femn(1:oqfmsmn), fomn(1:oqfmsmn), fcmn(1:oqfmsmn), fsmn(1:oqfmsmn), &
                    oqfmsNteta, oqfmsNzeta, jkreal(1:oqfmsNtz), jkimag(1:oqfmsNtz) )

  jkimag = one + jkimag

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  fonm = + oqfmsin * efmn ; fenm = zero
  fcnm = - oqfmsin * sfmn ; fsnm = zero

  call oqfmsinvfft( oqfmsmn, oqfmsim(1:oqfmsmn), oqfmsin(1:oqfmsmn), fenm(1:oqfmsmn), fonm(1:oqfmsmn), fcnm(1:oqfmsmn), fsnm(1:oqfmsmn), &
                    oqfmsNteta, oqfmsNzeta, kjreal(1:oqfmsNtz), kjimag(1:oqfmsNtz) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

!  call oqfmstfft( oqfmsNteta, oqfmsNzeta, rjk(1:oqfmsNtz), tjk(1:oqfmsNtz), &
!                 oqfmsmn, oqfmsim, oqfmsin, efmn(1:oqfmsmn), ofmn(1:oqfmsmn), cfmn(1:oqfmsmn), sfmn(1:oqfmsmn) )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  do jj = 0, oqfmsNteta-1 ; teta = jj * pi2 / oqfmsNteta
   
   do kk = 0, oqfmsNzeta-1 ; zeta = kk * pi2 / oqfmsNzeta ; jk = 1 + jj + kk * oqfmsNteta
    
    xyz(1:3,0) = (/ ijreal(jk), ijimag(jk), ijozeta(jk) /)
    xyz(1:3,1) = (/ jireal(jk), jiimag(jk),    zero     /)
    xyz(1:3,2) = (/ jkreal(jk), jkimag(jk),    zero     /)
    xyz(1:3,3) = (/ kjreal(jk), kjimag(jk),     one     /)

    xy(1:Node) = (/ xyz(2,0), xyz(1,0), one, zero, zero, one /)
    
    call bfield( zeta, xy(1:Node), Bxy(1:Node) )
    
    norm(1:3) = (/ xyz(2,2) * xyz(3,3) - xyz(3,2) * xyz(2,3), &
                   xyz(3,2) * xyz(1,3) - xyz(1,2) * xyz(3,3), &
                   xyz(1,2) * xyz(2,3) - xyz(2,2) * xyz(1,3) /)
    
    sqrtg = xyz(1,1) * norm(1) + xyz(2,1) * norm(2) + xyz(3,1) * norm(3)
    
    grads(1:3) = norm(1:3) / sqrtg
    
    norm(1:3) = norm(1:3) / sqrt( norm(1)*norm(1) + norm(2)*norm(2) + norm(3)*norm(3) )
    
    bdotn = ( Bxy(2)*norm(1) + Bxy(1)*norm(2) + Bzeta*norm(3) ) / sqrt( Bxy(2)*Bxy(2) + Bxy(1)*Bxy(1) + Bzeta*Bzeta )
    
    qflux = qflux + ( bdotn * bdotn                                             ) * sqrtg
    gflux = gflux + ( grads(1)*grads(1) + grads(2)*grads(2) + grads(3)*grads(3) ) * sqrtg

    dTdz(1) = dTdz(1) + ( bdotn * bdotn / tkperp + grads(1)*grads(1) + grads(2)*grads(2) + grads(3)*grads(3) ) * sqrtg
    
!   if( jj.eq. oqfmsNteta/4 .and. kk.eq. oqfmsNzeta/4 ) then
!    write(ounit,'("dTemp     : "10x"  : ("f8.5","f8.5","f8.5" ) ("f8.5","f8.5","f8.5" ) : B(r,t,z) ="3f14.10" ; norm ="3es13.5" ; bdotn ="es13.5" ;")') &
! lflux, teta, zeta, xy(2), xy(1), zeta, Bxy(2), Bxy(1), Bzeta, norm(1:3), bdotn
!   endif

   enddo ! end of do kk ; 2020/02/16 ;
   
  enddo ! end of do jj ; 2020/02/16 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! probably better to use FFT to perform summation ; ( probably faster and will also get Fourier harmonics of parallel/perpendicular terms ) ; 02/19/2021 ;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! write(ounit,'("dTemp     : ",10x,"s : qflux ="f18.10" ; gflux ="f18.10" ;")') qflux, gflux
  
  qflux = qflux * (pi2/oqfmsNteta) * (pi2/oqfmsNzeta) / tkperp
  gflux = gflux * (pi2/oqfmsNteta) * (pi2/oqfmsNzeta)

  dTdz(1) = one / ( qflux + gflux )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine dTemp

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine dTempout( lflux, Temp )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  use globals, only : chi, Ni, Tempprofile, ijk, qflux, gflux
  
  implicit none

  include "mpif.h"

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, parameter   :: lNode = 1
  real                 :: lflux, Temp(1:lNode), dTdz(1:lNode)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call dTemp( lflux, Temp, dTdz )

  Tempprofile(ijk(1),0) = lflux
  Tempprofile(ijk(1),1) = Temp(1)
  Tempprofile(ijk(1),2) = dTdz(1)
  Tempprofile(ijk(1),3) = qflux
  Tempprofile(ijk(1),4) = gflux
  
  ijk(1) = ijk(1) + 1 ; lflux = chi%low + ijk(1) * ( chi%upp - chi%low ) / Ni  

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  return

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

end subroutine dTempout

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine invcoords( NFT, fluxteta, fvec, fjac, iuser, ruser, iflag )
  
  use globals, only : zero, one, pi2, machprec, small, sqrtmachprec, ounit, &
                      chi, &
                      Ni, Nteta, Nzeta, &
                      Linterpol, &
                      oNqfms, oqfmsmn, oqfmsim, oqfmsin, oflux, orcmn, drcmn, otsmn, dtsmn
  
  implicit none
  
  include "mpif.h"
  
  integer :: NFT, iuser(1:3), iflag
  real    :: fluxteta(1:NFT), fvec(1:NFT), fjac(1:NFT,1:NFT), ruser(1:1)
  
  integer :: ii, jj, kk, imn, mm, ie01bgf, ifail, ierr, iqfms
  real    :: lflux, teta, zeta, lorcmn, ldrcmn, lotsmn, ldtsmn, xx(0:2), zz(0:2), zzi, xxj, ff(0:2), gg(0:2)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  ii = iuser(1) ; jj = iuser(2) ; kk = iuser(3)

  lflux = min( max( fluxteta(1), oflux(1) + machprec ), oflux(oNqfms) - machprec ) ; teta = fluxteta(2) ; zeta = kk * pi2 / Nzeta
  
  zzi = chi%low + ii * ( chi%upp - chi%low ) / Ni

  xxj = jj * pi2 / Nteta

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( Linterpol )
   
  case( 1 ) ! Linterpol = 1 ; 02/19/2021 ;
   
   zz(0:2) = (/ zero, zero, zero /)
   xx(0:2) = (/ teta, zero,  one /)
   
   ifail = 1

   do iqfms = 1, oNqfms-1
    
    if( lflux.ge.oflux(iqfms) .and. lflux.le.oflux(iqfms+1) ) then
     
     do imn = 1, oqfmsmn
      
      lorcmn = orcmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      ldrcmn =                                               ( orcmn(iqfms+1,imn) - orcmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      lotsmn = otsmn(iqfms,imn) + ( lflux - oflux(iqfms) ) * ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      ldtsmn =                                               ( otsmn(iqfms+1,imn) - otsmn(iqfms,imn) ) / ( oflux(iqfms+1) - oflux(iqfms) )
      
      zz(0) = zz(0) + lorcmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      zz(1) = zz(1) + ldrcmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      zz(2) = zz(2) + lorcmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta ) * ( - oqfmsim(imn) )
      
      xx(0) = xx(0) + lotsmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      xx(1) = xx(1) + ldtsmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
      xx(2) = xx(2) + lotsmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta ) * ( + oqfmsim(imn) )
      
     enddo ! end of do imn ; 02/19/2021 ;
     
     ifail = 0 

     exit

    endif

   enddo ! end of do iqfms ; 02/19/2021 ;
     
   if( ifail.ne.0 ) write(ounit,'(3es23.15)') oflux(1), lflux-one, oflux(oNqfms)

   if(ifail.ne.0 ) then
          write(ounit,'("invcoords :    fatal    : ifail.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif

  case( 2 ) ! Linterpol = 2 ; 02/19/2021 ;

   zz(0:2) = (/ zero, zero, zero /)
   xx(0:2) = (/ teta, zero,  one /)
   
   do imn = 1, oqfmsmn ; mm = 1
    
    ie01bgf = 1 ; call E01BGF( oNqfms, oflux(1:oNqfms), orcmn(1:oNqfms,imn), drcmn(1:oNqfms,imn), mm, lflux, lorcmn, ldrcmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), orsmn(1:oNqfms,imn), drsmn(1:oNqfms,imn), mm, lflux, lorsmn, ldrsmn, ie01bgf )
!   ie01bgf = 0 ; call E01BGF( oNqfms, oflux(1:oNqfms), otcmn(1:oNqfms,imn), dtcmn(1:oNqfms,imn), mm, lflux, lotcmn, ldtcmn, ie01bgf )
    ie01bgf = 1 ; call E01BGF( oNqfms, oflux(1:oNqfms), otsmn(1:oNqfms,imn), dtsmn(1:oNqfms,imn), mm, lflux, lotsmn, ldtsmn, ie01bgf )
    
    zz(0) = zz(0) + lorcmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    zz(1) = zz(1) + ldrcmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    zz(2) = zz(2) - lorcmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta ) * oqfmsim(imn)
    
    xx(0) = xx(0) + lotsmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    xx(1) = xx(1) + ldtsmn * sin( oqfmsim(imn) * teta - oqfmsin(imn) * zeta )
    xx(2) = xx(2) + lotsmn * cos( oqfmsim(imn) * teta - oqfmsin(imn) * zeta ) * oqfmsim(imn)
    
   enddo ! end of do imn ; 02/19/2021 ;
  
  case default

  end select ! end select( Linterpol ) ; 02/19/2021 ;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  ff(0) = zz(0) - zzi ; ff(1) = zz(1) ; ff(2) = zz(2)
  gg(0) = xx(0) - xxj ; gg(1) = xx(1) ; gg(2) = xx(2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( iflag )
   
  case( 1 )
   
   fvec(1:2) = (/ ff(0), gg(0) /)
   
  case( 2 )
   
   fjac(1,1) = ff(1) ; fjac(1,2) = ff(2)
   fjac(2,1) = gg(1) ; fjac(2,2) = gg(2)
   
  case( 3 )
   
   fvec(1:2) = (/ ff(0), gg(0) /)
   
   fjac(1,1) = ff(1) ; fjac(1,2) = ff(2)
   fjac(2,1) = gg(1) ; fjac(2,2) = gg(2)
   
  case default
   
  end select

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  if( sqrt(ff(0)*ff(0)+gg(0)*gg(0)).lt.sqrtmachprec ) iflag = -1 ! it might not be the case that fluxteta is returned ; 02/19/2021 ;

! write(ounit,1100) ii, kk, jj, fluxteta, ff(0), gg(0), ff(1)*gg(2)-ff(2)*gg(1), sqrt(ff(0)*ff(0)+gg(0)*gg(0)), iflag

1100 format("initializ : ", 10x,"s : ",3i6," (",f20.16,",",f20.16," ) : fvec =",2es23.15," ;",2es13.05,i3)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  return

end subroutine invcoords

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
