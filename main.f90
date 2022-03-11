














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

program heat
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, pi2, ten, half, ounit, tunit, myid, numprocs, readinput, infile, hmnfile, &
                      skperp, siNi, siNj, siNk, sMpol, sNtor, &
                      ekperp, eiNi,             eMpol, eNtor, &
                      tstart, imethod, &
                      kperp, Tlow, Tupp, chi, &
                      Ni, Nj, Nk, Ndof, &
                      oNqfms, oqfmsmn, oqfmsim, oqfmsin, oqfmsNteta, oqfmsNzeta, oqfmsNtz, &
                      Nqfms, &
                      oflux, orcmn, orsmn, otcmn, otsmn, &
                      le04dgf, dTtol, ld02bjf, lrkutta, Ntau, Mtau, dtau, difforig, dodetol, dend, &
                      NPpts, Laction, &
                      PX, PY, lamda, mu, cc, e02defrwrk, e02defiwrk, &
                      dx, dy, dz, Mpol, Ntor, mn, im, in, Nteta, &
                      Linterpol, dcomp, LTprofile, &
                      opp, oqq, oiota, dflux, drcmn, drsmn, dtcmn, dtsmn, oqfmscplxin, oqfmscplxout, ijoteta, ijozeta, oPxyz, &
                      pp, qq, iota, flux, Pxyz, &
                      Naction
  
  use, intrinsic :: iso_c_binding
  
  implicit none
  
  include "mpif.h"
  
  include 'fftw3.f03'
  
  LOGICAL              :: T_exist
  integer              :: iok, ierr, itau, jtau, ii, jj, kk, astat, ifd, liw, lrw, mm, nn, iuser(1:5), iter, ie04dgf, id02bjf, oNi, oNj, oNk, omn, il, lNi
  integer              :: okperp, MX, MY, lwrk, ie01daf, ie01bef, ie01bgf, iaction, ipqadd, itprofile, ireadqfms, iterations
  integer, allocatable :: iw(:)
  real                 :: cput, cpul, diff(0:4), diffl(0:4), diffu(0:4), fdiff, fdest, totaldiffusion, ruser(1:1), timestart, timeend, teta
  real                 :: stz(1:3), xyz(1:3,0:3), srho, lTmn, lDmn, dcomparison(0:2)
  real                 :: oTlow, oTupp, olow, oupp, lodetol
  real, allocatable    :: TT(:,:,:), dT(:,:,:), TL(:,:,:), TU(:,:,:), rw(:), Tijk(:), dTijk(:), d02bjfwk(:)
  real, allocatable    :: fTT(:,:), dfT(:,:), fTL(:,:), fTU(:,:), oTT(:)
  real, allocatable    :: Mijk(:,:,:)
  real                 :: rho, lorcmn, lorsmn, lotcmn, lotsmn
  real, allocatable    :: ijreal(:), ijimag(:), efmn(:), ofmn(:), cfmn(:), sfmn(:)
  
  real   , allocatable :: xx(:), yy(:), ff(:)
  real   , allocatable :: Smn(:), Tmn(:,:),  Dmn(:,:)

  EXTERNAL             :: fdiffuse, gdiffuse, dfield, dfieldout, D02BJX, D02BJW
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call mpi_init( ierr )
  call mpi_comm_rank( mpi_comm_world, myid, ierr )
  call mpi_comm_size( mpi_comm_world, numprocs, ierr )
  
! c mpilation ! replaced with a write statement identifying details of compilation; see Makefile;
  
! call A00AAF() ! determine mark of NAG library;
  
  tstart = MPI_WTIME()
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  call readinput ! reads .hr, .hmn ; 05/21/21;
  
  call preset ! 05/21/21;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  ipqadd = 0

  ireadqfms = 1 ; iaction = 1 ; itprofile = 1
  
  do iterations = 1, Naction
   
   call readqfms( ipqadd, ireadqfms ) ; ireadqfms = 0 ! reads .qfms, interpolates; reads .pq, adds ipqadd mediant ; 05/21/21;
   
   if( Laction.ge.1 ) call action( ipqadd, iaction ) ; iaction = 0 ! constructs QFM surfaces; writes .qfms ; 05/21/21;
   
   call readqfms( ipqadd, ireadqfms )

   if( LTprofile.ge.1 ) call tprofile( itprofile ) ; itprofile = 0 ! integrates o.d.e. for temperature; maps to background; 05/21/21;
   
   dcomparison(0:2) = dcomp(0:2,1)

  enddo

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  if( Laction.ge.2 ) then
   
   do ipqadd = 1, Nqfms-1
    
    ireadqfms = 1 ; iaction = 0 ; itprofile = 0
    
    do iterations = 1, Naction
     
     call readqfms( ipqadd, ireadqfms ) ; ireadqfms = 0
     
     if( Laction.ge.1 ) call action( ipqadd, iaction ) ; iaction = 0
     
     call readqfms( ipqadd, ireadqfms )
     
     if( LTprofile.ge.1 ) call tprofile( itprofile ) ; itprofile = 0
     
    enddo ! end of do iterations ; 04/20/21 ;
    
    write(ounit,'("main      : "10x"  :        "es13.05"  "10x"           "es13.05"  "10x"           "es13.05"  " )') dcomp(0:2,1) - dcomparison(0:2)

   enddo ! end of do ipqadd ; 04/20/21 ;
   
  endif ! end of if( Laction.ge.2 ) ; 04/20/21 ;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  ipqadd = 0

  ireadqfms = 1 ; iaction = 1 ; itprofile = 1
  
  do iterations = 1, Naction
   
   call readqfms( ipqadd, ireadqfms ) ; ireadqfms = 0
   
   if( Laction.ge.1 ) call action( ipqadd, iaction ) ; iaction = 0
   
   call readqfms( ipqadd, ireadqfms )

   if( LTprofile.ge.1 ) call tprofile( itprofile ) ; itprofile = 0
   
  enddo

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  call tsolve
 
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
! call maptemp
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if( NPpts.gt.0 ) call poincare ! this comes after qfms surfaces and temperature is constructed ; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
9998 continue
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cput = MPI_WTIME()
  write(ounit,9999) cput-tstart, (cput-tstart) / (/ 1, 60, 60**2, 24*60**2 /), trim(infile), trim(hmnfile)
  cpul = cput
  
9999 format("main      : "f10.1"s : completion ; time =",f10.1,"s = ",f8.2,"m = ",f6.2,"h = ",f5.2,"d ; ",a,":",a," ;")
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call MPI_FINALIZE(ierr)
  
  stop
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end program heat

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

