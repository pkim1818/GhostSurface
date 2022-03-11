














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

module globals
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include 'fftw3.f03'
  
  real, parameter :: zero    =   0.0
  real, parameter :: one     =   1.0 
  real, parameter :: two     =   2.0 
  real, parameter :: three   =   3.0 
  real, parameter :: four    =   4.0 
  real, parameter :: five    =   5.0 
  real, parameter :: six     =   6.0 
  real, parameter :: seven   =   7.0 
  real, parameter :: eight   =   8.0 
  real, parameter :: nine    =   9.0 
  real, parameter :: ten     =  10.0
  real, parameter :: twelve  =  12.0
  real, parameter :: hundred = 100.0

  real, parameter :: half   =  one / two
  real, parameter :: pi2    = 6.283185307179586 
  real, parameter :: pi     = pi2 / two

  real, parameter :: goldmean = ( one + sqrt(five) ) / two
  
  real            :: machprec, sqrtmachprec, small

  integer         :: ounit =  0
  integer         :: iunit = 10
  integer         :: munit = 11
  integer         :: tunit = 12
  
  TYPE(C_PTR)                            :: planf, planb
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: cplxin(:,:), cplxout(:,:)

  TYPE(C_PTR)                            :: qfmsplanf, qfmsplanb
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: qfmscplxin(:,:), qfmscplxout(:,:)

  TYPE(C_PTR)                            :: oqfmsplanf, oqfmsplanb
  COMPLEX(C_DOUBLE_COMPLEX), allocatable :: oqfmscplxin(:,:), oqfmscplxout(:,:)

  integer         :: icheck     =      0
  integer         :: Imethod    =      3
  integer         :: itlimit    =    100
  real            :: linesearch =     0.5
  integer         :: Nfp        =      1

  real            :: Tlow       =      0.0 
  real            :: Tupp       =      1.0 
  real            :: kpara      =      1.0 
  integer         :: kperp      =      4      

  integer         :: iNi        =      7 
  integer         :: iNj        =      7 
  integer         :: iNk        =      7
  integer         :: Mpol       =      4
  integer         :: Ntor       =      2

  real            :: odetol     =      1.0e-08
  integer         :: idiff      =      2 
  integer         :: kdiff      =      2
  integer         :: opint      =      2
  integer         :: iupdate    =      0
  integer         :: Mits       =    100
  real            :: lintol     =      1.0e-12
  integer         :: lbasis     =     10
  real            :: omega      =      1.0
  character       :: method*8   =     'BICGSTAB'
  character       :: precon     =     'S'

  integer         :: le04dgf    =      0
  real            :: dTtol      =      1.0e-06
  integer         :: ld02bjf    =      0
  integer         :: lrkutta    =      0
  real            :: dend       =      1.0e+00
  real            :: dodetol    =      1.0e-06
  integer         :: Ntau       =     10
  integer         :: Mtau       =     10
  real            :: dtau       =      1.0e-06

  integer         :: NPtrj      =     50
  integer         :: NPpts      =      0

  integer         :: LTprofile  =      0 ! integration of temperature profile ; LTprofile = 1 : D02BJF ; LTprofile = 2 : 1st order ; 04/1/21 ;
  integer         :: Linterpol  =      1 ! interpolation between qfm surfaces ; Linterpol = 1 : piecewise linear ; Linterpol = 2 : E01BGF ; 04/1/21 ;

  integer         :: Laction    =      1
  integer         :: Naction    =      3
  integer         :: NFarey     =      4
  integer         :: pqMpol     =      2
  integer         :: pqNtor     =      2
  integer         :: qfmsMpol   =      4
  integer         :: qfmsNtor   =      2

  integer         :: Lmoser     =      0
  real            :: kpert      =      0.0e-00
  integer         :: npert      =    100
  real            :: viscosity  =      0.0e-00

  namelist/heatlist/    &
           icheck      ,&
           imethod     ,&
           itlimit     ,&
           linesearch  ,&
           Nfp         ,&
           Tlow        ,&
           Tupp        ,&
           kpara       ,&
           kperp       ,&
           NPtrj       ,&
           NPpts       ,&
           iNi         ,&
           iNj         ,&
           iNk         ,&
           Mpol        ,&
           Ntor        ,&
           odetol      ,&
           idiff       ,&
           kdiff       ,&
           opint       ,&
           Mits        ,&
           lintol      ,&
           lbasis      ,&
           omega       ,&
           method      ,&
           precon      ,&
           le04dgf     ,&
           dTtol       ,&
           ld02bjf     ,&
           lrkutta     ,&
           Ntau        ,&
           Mtau        ,&
           dtau        ,&
           dend        ,&
           dodetol     ,&
           LTprofile   ,&
           Linterpol   ,&
           Laction     ,&
           Naction     ,&
           NFarey      ,&
           pqMpol      ,&
           pqNtor      ,&
           qfmsMpol    ,&
           qfmsNtor    ,&
!          maxiter     ,&
           Lmoser      ,&
           kpert       ,&
           npert       ,&
           viscosity   ,&
           iupdate    

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  type magneticfieldlinehamiltonian
     integer              :: mn, Np                 ! #Fourier harmonics, degree of polynomial approximation;
     integer, allocatable :: im(:), in(:)           ! mode identification;
     real   , allocatable :: kmn(:,:)                  ! magnetic-field-line harmonics;
     real                 :: low, upp, scale         ! lower and upper boundary;
     real   , allocatable :: k(:,:)                    ! harmonic and derivatives;
  end type magneticfieldlinehamiltonian
  
  type(magneticfieldlinehamiltonian) :: chi ! this structure contains all information regarding the magnetic field;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  real   , allocatable :: zfcmn(:,:), zfsmn(:,:), zhcmn(:,:), zhsmn(:,:)

  integer              :: mn, qfmsmn, Nteta, Nzeta, Ntz, qfmsNteta, qfmsNzeta, qfmsNtz, ipqfms, oqfmsNteta, oqfmsNzeta, oqfmsNtz
  integer, allocatable :: im(:), in(:), qfmsim(:), qfmsin(:)
  real   , allocatable :: zbc(:,:), zbs(:,:), cosij(:,:,:), sinij(:,:,:)

  real                 :: Bzeta ! zz field returned through common;

  integer              :: myid, numprocs
  
  real                 :: pi2nfp

  real                 :: tstart
  
  integer              :: Ni, Nj, Nk, Ndof ! grid resolution in x,y, and z; #degrees-of-freedom;

  integer              :: ijk(1:3), kbf ! counters used in field-line integration;

  real                 :: dx, dy, dz, dv ! grid differences;

  real                 :: xxn(0:3), yym(0:3), bcf(1:16), fdf(1:16)
  integer              :: BC(1:16,1:16), FD(1:16,1:16)

  real   , allocatable :: Tijk(:) ! solution on full interior (which may include inner boundary);
  real   , allocatable :: T0(:,:), C0(:,:), D0(:,:), E0(:,:), T1(:,:), C1(:,:), D1(:,:), E1(:,:)

  real   , allocatable :: Bx(:,:,:), By(:,:,:), Bz(:,:,:), BB(:,:,:)
  real   , allocatable :: fBs(:,:), fBt(:,:), fBz(:,:), fBB(:,:)
  real                 :: difforig(0:2)

  LOGICAL              :: matrix_exist

  integer              :: Mnz, nz ! used for coordinate storage representation of linear operator;
  integer, allocatable :: irow(:), jcol(:)
  real   , allocatable :: Acs(:), rhs(:) ! matrix representation of anisotropic diffusion coordinate storage format;
  
  real   , allocatable :: Bfh(:,:,:,:) ! B^phi / B^2 on forward/backward full (i,j) half (k) grid;

  real   , allocatable :: dxydabf(:,:,:,:,:) ! coordinate derivatives dR/da, dR/db, dZ/da, dZ/db on forward/backward full-grid;
    
  LOGICAL              :: Lbck, Lfwd, Lfull, Lhalf ! controls field line integration construction of coordinates;
  
  character            :: infile*100, hmnfile*100

  character            :: skperp*2, siNi*2, siNj*2, siNk*2, sMpol*2, sNtor*2
  character            :: ekperp*2, eiNi*2,                 eMpol*2, eNtor*2

  character            :: metrixsuff*8, fouriersuff*2
  character            :: ksuff*9, lsuff*9, csuff*14, msuff*2, nsuff*2

  integer              :: Rdof
  integer, allocatable :: ijkn(:,:), nijk(:,:,:) ! inverts (i,j,k) identification of grid point;
  
  integer, parameter   :: NDpq = 4
  integer              :: Nqfms
  integer, allocatable :: pp(:), qq(:)
! REAL   , allocatable :: iota(:), avetemp(:), flux(:)
  real   , allocatable :: iota(:), flux(:)

  integer              :: oNqfms, oqfmsmn, oqfmsMpol, oqfmsNtor
  integer, allocatable :: oqfmsim(:), oqfmsin(:)
  integer, allocatable :: opp(:),     oqq(:)     

  real   , allocatable :: oiota(:), oflux(:), dflux(:)
  real   , allocatable :: orcmn(:,:), orsmn(:,:), otcmn(:,:), otsmn(:,:)
  real   , allocatable :: drcmn(:,:), drsmn(:,:), dtcmn(:,:), dtsmn(:,:)

  real   , allocatable :: ijoteta(:), ijozeta(:)

  real   , allocatable ::  Pxyz(:,:,:)
  real   , allocatable :: oPxyz(:,:,:)

  real                 :: tkperp

  integer              :: nperp(0:2), ifd(0:2,1:9), jfd(0:2,1:9), npint(1:2), iin(1:2,1:16), jin(1:2,1:16)
  real                 :: wfd(0:2,1:9), wpd(1:2,-2:2,-2:2)

  integer              :: PX, PY
  integer, allocatable :: e02defiwrk(:)
  real   , allocatable :: lamda(:), mu(:), cc(:), e02defrwrk(:)

  real   , allocatable :: Tempprofile(:,:)

  real                 :: qflux, gflux
  
  real                 :: dcomp(0:2,0:2)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

contains
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  subroutine readinput
    
    implicit none
    
    include "mpif.h"
    
    integer              :: ierr, ios, getarg, astat, imn, iov, iargc
    real                 :: cput, X02AJF
    
    LOGICAL              :: qfms_exist
    integer              :: ii, jj, kk, mm, nn, lNqfms, iFarey, ifail, ipq

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    machprec = X02AJF() ; small = 1.0e+02 * machprec ; sqrtmachprec=sqrt(machprec)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    if( myid.eq.0 ) then
     
     iov = iargc() ! determine number of command line arguments;

     if(iov.lt.2 ) then
          write(ounit,'("readinput :             : two command line arguments required :    fatal    : iov.lt.2  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
     
     ios = getarg( 1, infile ) ! first command line argument is input file extension;
     ios = getarg( 2, hmnfile ) ! second is magnetic field line Hamiltonian;

     open( iunit, file=trim(infile)//".hr", status="old", iostat=ios )

     if( ios.ne.0 ) close(iunit)

     if(ios.ne.0 ) then
          write(ounit,'("readinput : error opening ext.hr :    fatal    : ios.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif

     read( iunit, heatlist )

     close( iunit )

     cput = MPI_WTIME()

     tkperp = one / (ten**kperp)
     
     write(ounit,1000) cput-tstart, trim(infile), trim(hmnfile), icheck, imethod, itlimit, linesearch
     write(ounit,1001) cput-tstart, Tlow, Tupp, kpara, kperp
     write(ounit,1002) cput-tstart, iNi, iNj, iNk, 2**iNi, 2**iNj, 2**iNk, Mpol, Ntor
!    write(ounit,1003) cput-tstart, odetol, idiff, kdiff, opint, Mits, lintol
!    write(ounit,1004) cput-tstart, lbasis, omega, method, precon, iupdate
     write(ounit,1005) cput-tstart, le04dgf, dTtol, ld02bjf, lrkutta, dend, dodetol
!    write(ounit,1006) cput-tstart, Ntau, Mtau, dtau
     write(ounit,1007) cput-tstart, NPtrj, NPpts, LTprofile, Linterpol
     write(ounit,1008) cput-tstart, Laction, Naction, NFarey, pqMpol, pqNtor, qfmsMpol, qfmsNtor
!    write(ounit,1009) cput-tstart, Lmoser, kpert, npert, viscosity

1000 format("readinput : "f10.1"s : ",a," ; ",a," ; icheck =",i2," ; imethod =",i2," ; itlimit ="i6" ; linesearch ="e12.5" ;")
1001 format("readinput : "f10.1"s : Tlow =",f7.3," ; Tupp =",f7.3," ; kpara =",es8.1," ; kperp = 1/ten^(",i2.2,") ; ")
1002 format("readinput : "f10.1"s : iNi =",i3," ; iNj =",i3," ; iNk =",i3," ; Ni =",i5," ; Nj =",i5," ; Nk =",i5," ; Mpol =",i4," ; Ntor =",i3," ;")
1003 format("readinput : "f10.1"s : odetol =",es8.1," ; idiff =",i2," ; kdiff =",i2," ; opint =",i2," ; Mits =",i6," ; lintol =",es8.1," ;")
1004 format("readinput : "f10.1"s : lbasis =",i3," ; omega =",f6.3," ; method = ",a," ; precon = ",a," ; iupdate =",i3," ;")
1005 format("readinput : "f10.1"s : le04dgf =",i2," ; dTtol =",es9.1," ; ld02bjf =",i2," ; lrukkta =",i2," ; dend =",es10.3," ; dodetol =",es10.3," ;")
1006 format("readinput : "f10.1"s : Ntau =",i9," ; Mtau =",i9," ; dtau =",es10.3," ;")
1007 format("readinput : "f10.1"s : NPtrj =",i6," ; NPpts =",i6," ; LTprofile ="i2" ; Linterpol ="i2" ;")
1008 format("readinput : "f10.1"s : Laction =",i2," ; Naction ="i3" ; NFarey =",i3," ; pqMpol =",i3," ; pqNtor =",i3," ; qfmsMpol =",i4," ; qfmsNtor =",i3," ;")
1009 format("readinput : "f10.1"s : Lmoser =",i3," ; kpert =",f8.4," ; npert =",i4," ; viscosity =",f8.4," ;")
     
    endif
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    write( skperp,'(i2.2)') kperp ; write( ekperp,'(i2.2)') kperp+1
    write( siNi  ,'(i2.2)') iNi   ; write( eiNi  ,'(i2.2)') iNi+1
    write( siNj  ,'(i2.2)') iNj   ; 
    write( siNk  ,'(i2.2)') iNk   ; 
    write( sMpol ,'(i2.2)') Mpol  ; write( eMpol ,'(i2.2)') Mpol+2
    write( sNtor ,'(i2.2)') Ntor  ; write( eNtor ,'(i2.2)') Ntor+2

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

    write( metrixsuff, '(i2.2":"i2.2":"i2.2)' ) iNi, iNj, iNk ! suffix for metrics file;
    write( fouriersuff, '(              i2.2)' )             iNi ! suffix for metrics file;
    write( msuff, '(i2.2)') kperp
    write( nsuff, '(i2.2)') kperp+1
    write( ksuff, '(i3.2":"i1":"i1":"i1)' ) kperp  , idiff, kdiff, opint ! suffix for matrix file;
    write( lsuff, '(i3.2":"i1":"i1":"i1)' ) kperp+1, idiff, kdiff, opint ! suffix for matrix file;
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    if( myid.eq.0 ) then
     
     cput = MPI_WTIME()
     
     open( iunit, file=trim(hmnfile)//".hmn", status="old", iostat=ios )
     
     if( ios.ne.0 ) close(iunit)
     
     if(ios.ne.0 ) then
          write(ounit,'("readinput :             : error opening ext.hmn :    fatal    : ios.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
     
     read( iunit, *, iostat=ios ) chi%mn, chi%Np, Nfp, chi%low, chi%upp, chi%scale
     
     if(ios.ne.0 ) then
          write(ounit,'("readinput :             : error reading mn Np Nfp yyl yyu scale :    fatal    : ios.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
     
     write(ounit,1100) cput-tstart, chi%mn, chi%Np, Nfp, chi%low, chi%upp, chi%scale
     
1100 format("readinput : "f10.1"s : chi%[ mn ="i3", Np ="i3", Nfp ="i3", low ="f9.5", upp ="f9.5", scale ="f9.5" ] ;")
     
     if(chi%mn.lt.1 ) then
          write(ounit,'("readinput :             : :    fatal    : chi%mn.lt.1  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif ! need at least one harmonic;
     if(chi%Np.lt.2 ) then
          write(ounit,'("readinput :             : :    fatal    : chi%Np.lt.2  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif ! need at least a quadratic to give shear;
     
    endif
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    if(allocated(chi%im     )) stop 'chi%im      allocated'
	allocate(chi%im     (1:chi%mn)          ,stat=astat)
	if(astat.ne.0) stop 'error allocating chi%im     ' ! poloidal mode identification;
    if(allocated(chi%in     )) stop 'chi%in      allocated'
	allocate(chi%in     (1:chi%mn)          ,stat=astat)
	if(astat.ne.0) stop 'error allocating chi%in     ' ! toroidal mode identification;
    if(allocated(chi%kmn    )) stop 'chi%kmn     allocated'
	allocate(chi%kmn    (1:chi%mn,0:chi%Np) ,stat=astat)
	if(astat.ne.0) stop 'error allocating chi%kmn    ' ! polynomial coefficients;
    if(allocated(chi%k      )) stop 'chi%k       allocated'
	allocate(chi%k      (1:chi%mn,0:3)      ,stat=astat)
	if(astat.ne.0) stop 'error allocating chi%k      ' ! perturbation harmonics and derivatives;
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    if( myid.eq.0 ) then
     
     do imn = 1, chi%mn
      
      read( iunit, *, iostat=ios ) chi%im(imn), chi%in(imn), chi%kmn(imn,0:chi%Np)
      
      if(ios.ne.0 ) then
          write(ounit,'("readinput : error reading m n kmn(0:Np) :    fatal    : ios.ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
      
     enddo
     
     if(chi%im(1).ne.0 ) then
          write(ounit,'("readinput :    fatal    : chi%im(1).ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif ! ensure first harmonic is integrable component;
     if(chi%in(1).ne.0 ) then
          write(ounit,'("readinput :    fatal    : chi%in(1).ne.0  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
     
     close( iunit )
     
    endif
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    chi%kmn(2:chi%mn,0:chi%Np) = chi%kmn(2:chi%mn,0:chi%Np) * chi%scale
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    pi2nfp = pi2 / Nfp
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
    
    return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  end subroutine readinput
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end module globals

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine tfft( Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn )
  
  use globals, only : half, zero, pi2, nfp, cplxin, cplxout, planf
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include 'fftw3.f03'
  
  intrinsic aimag
  
  integer :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, ii, mm, nn
  real    :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  
  COMPLEX(C_DOUBLE_COMPLEX) :: z1, z2, z3
  
  do ii = 1, Nz ; cplxin(:,ii) = cmplx( ijreal((ii-1)*Nt+1:ii*Nt), ijimag((ii-1)*Nt+1:ii*Nt), kind=c_double_complex )
  enddo
  
  call fftw_execute_dft( planf, cplxin, cplxout ) ! forward transform ;
  
  Ntz = Nt * Nz
  
  cplxout = cplxout / Ntz ; cplxout(1,1) = half * cplxout(1,1)
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   
   z1 = cplxout( 1 + mod(Nt-mm,Nt), 1 + mod(Nz+nn,Nz) )
   z2 = cplxout( 1 +        mm,     1 + mod(Nz-nn,Nz) )
   
   z3 = z1 + z2 ; efmn(ii) =  real(z3) ; cfmn(ii) = aimag(z3)
   z3 = z1 - z2 ; ofmn(ii) = aimag(z3) ; sfmn(ii) = -real(z3)

  enddo
  
  return

end subroutine tfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine invfft( mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag )
  
  use globals, only : zero, two, half, nfp, cplxin, cplxout, planb

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  integer, intent(in)  :: mn, im(mn), in(mn)
  real   , intent(in)  :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn) ! input Fourier harmonics;
  integer, intent(in)  :: Nt, Nz
  real   , intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  
  integer              :: ii, mm, nn
  
  ijreal = zero ; ijimag = zero
  
  cplxin = zero
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   cplxin( 1+mod(Nt-mm,Nt), 1+mod(Nz+nn,Nz) ) = half * cmplx( efmn(ii) - sfmn(ii), cfmn(ii) + ofmn(ii), kind=c_double_complex )
   cplxin( 1+       mm,     1+mod(Nz-nn,Nz) ) = half * cmplx( efmn(ii) + sfmn(ii), cfmn(ii) - ofmn(ii), kind=c_double_complex )
  enddo

  cplxin(1,1) = two * cplxin(1,1) ! 29 Oct 19;
  
  call fftw_execute_dft( planb, cplxin, cplxout ) ! inverse transform ;
  
  do ii = 1, Nz
   ijreal( (ii-1)*Nt+1 : ii*Nt ) =  real( cplxout(:,ii) )
   ijimag( (ii-1)*Nt+1 : ii*Nt ) = aimag( cplxout(:,ii) )
  enddo
  
  return

end subroutine invfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine qfmstfft( Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn )
  
  use globals, only : half, zero, pi2, nfp, qfmscplxin, qfmscplxout, qfmsplanf
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include 'fftw3.f03'
  
  intrinsic aimag
  
  integer :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, ii, mm, nn
  real    :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  
  COMPLEX(C_DOUBLE_COMPLEX) :: z1, z2, z3
  
  do ii = 1, Nz ; qfmscplxin(:,ii) = cmplx( ijreal((ii-1)*Nt+1:ii*Nt), ijimag((ii-1)*Nt+1:ii*Nt), kind=c_double_complex )
  enddo
  
  call fftw_execute_dft( qfmsplanf, qfmscplxin, qfmscplxout ) ! forward transform ;
  
  Ntz = Nt * Nz
  
  qfmscplxout = qfmscplxout / Ntz ; qfmscplxout(1,1) = half * qfmscplxout(1,1)
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   
   z1 = qfmscplxout( 1 + mod(Nt-mm,Nt), 1 + mod(Nz+nn,Nz) )
   z2 = qfmscplxout( 1 +        mm,     1 + mod(Nz-nn,Nz) )
   
   z3 = z1 + z2 ; efmn(ii) =  real(z3) ; cfmn(ii) = aimag(z3)
   z3 = z1 - z2 ; ofmn(ii) = aimag(z3) ; sfmn(ii) = -real(z3)

  enddo
  
  return

end subroutine qfmstfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine qfmsinvfft( mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag )
  
  use globals, only : zero, two, half, nfp, qfmscplxin, qfmscplxout, qfmsplanb

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  integer, intent(in)  :: mn, im(mn), in(mn)
  real   , intent(in)  :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn) ! input Fourier harmonics;
  integer, intent(in)  :: Nt, Nz
  real   , intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  
  integer              :: ii, mm, nn
  
  ijreal = zero ; ijimag = zero
  
  qfmscplxin = zero
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   qfmscplxin( 1+mod(Nt-mm,Nt), 1+mod(Nz+nn,Nz) ) = half * cmplx( efmn(ii) - sfmn(ii), cfmn(ii) + ofmn(ii), kind=c_double_complex )
   qfmscplxin( 1+       mm,     1+mod(Nz-nn,Nz) ) = half * cmplx( efmn(ii) + sfmn(ii), cfmn(ii) - ofmn(ii), kind=c_double_complex )
  enddo

  qfmscplxin(1,1) = two * qfmscplxin(1,1) ! 29 Oct 19;
  
  call fftw_execute_dft( qfmsplanb, qfmscplxin, qfmscplxout ) ! inverse transform ;
  
  do ii = 1, Nz
   ijreal( (ii-1)*Nt+1 : ii*Nt ) =  real( qfmscplxout(:,ii) )
   ijimag( (ii-1)*Nt+1 : ii*Nt ) = aimag( qfmscplxout(:,ii) )
  enddo
  
  return

end subroutine qfmsinvfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine oqfmstfft( Nt, Nz, ijreal, ijimag, mn, im, in, efmn, ofmn, cfmn, sfmn )
  
  use globals, only : half, zero, pi2, nfp, oqfmscplxin, oqfmscplxout, oqfmsplanf
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include 'fftw3.f03'
  
  intrinsic aimag
  
  integer :: Nt, Nz, mn, im(1:mn), in(1:mn), Ntz, ii, mm, nn
  real    :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  
  COMPLEX(C_DOUBLE_COMPLEX) :: z1, z2, z3
  
  do ii = 1, Nz ; oqfmscplxin(:,ii) = cmplx( ijreal((ii-1)*Nt+1:ii*Nt), ijimag((ii-1)*Nt+1:ii*Nt), kind=c_double_complex )
  enddo
  
  call fftw_execute_dft( oqfmsplanf, oqfmscplxin, oqfmscplxout ) ! forward transform ;
  
  Ntz = Nt * Nz
  
  oqfmscplxout = oqfmscplxout / Ntz ; oqfmscplxout(1,1) = half * oqfmscplxout(1,1)
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   
   z1 = oqfmscplxout( 1 + mod(Nt-mm,Nt), 1 + mod(Nz+nn,Nz) )
   z2 = oqfmscplxout( 1 +        mm,     1 + mod(Nz-nn,Nz) )
   
   z3 = z1 + z2 ; efmn(ii) =  real(z3) ; cfmn(ii) = aimag(z3)
   z3 = z1 - z2 ; ofmn(ii) = aimag(z3) ; sfmn(ii) = -real(z3)

  enddo
  
  return

end subroutine oqfmstfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine oqfmsinvfft( mn, im, in, efmn, ofmn, cfmn, sfmn, Nt, Nz, ijreal, ijimag )
  
  use globals, only : zero, two, half, nfp, oqfmscplxin, oqfmscplxout, oqfmsplanb

  use, intrinsic :: iso_c_binding
  
  implicit none

  include 'fftw3.f03'
  
  integer, intent(in)  :: mn, im(mn), in(mn)
  real   , intent(in)  :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn) ! input Fourier harmonics;
  integer, intent(in)  :: Nt, Nz
  real   , intent(out) :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  
  integer              :: ii, mm, nn
  
  ijreal = zero ; ijimag = zero
  
  oqfmscplxin = zero
  
  do ii = 1, mn ; mm = im(ii) ; nn = in(ii) / nfp
   oqfmscplxin( 1+mod(Nt-mm,Nt), 1+mod(Nz+nn,Nz) ) = half * cmplx( efmn(ii) - sfmn(ii), cfmn(ii) + ofmn(ii), kind=c_double_complex )
   oqfmscplxin( 1+       mm,     1+mod(Nz-nn,Nz) ) = half * cmplx( efmn(ii) + sfmn(ii), cfmn(ii) - ofmn(ii), kind=c_double_complex )
  enddo

  oqfmscplxin(1,1) = two * oqfmscplxin(1,1) ! 29 Oct 19;
  
  call fftw_execute_dft( oqfmsplanb, oqfmscplxin, oqfmscplxout ) ! inverse transform ;
  
  do ii = 1, Nz
   ijreal( (ii-1)*Nt+1 : ii*Nt ) =  real( oqfmscplxout(:,ii) )
   ijimag( (ii-1)*Nt+1 : ii*Nt ) = aimag( oqfmscplxout(:,ii) )
  enddo
  
  return

end subroutine oqfmsinvfft

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
