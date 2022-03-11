














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine readqfms( ipqadd, ireadqfms )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals
  
  implicit none
  
  include "mpif.h"

  LOGICAL              :: qfms_exist, pq_exist
  integer, parameter   :: Node = 6, NFT = 2
  integer              :: astat, ii, jj, kk, nn, mm, ie01bef, jk, lpp, lqq, ipq, iFarey, lNqfms, lNode, ifail, iostat, ierr, ie01bgf, MX, MY, lwrk, ie01daf
  integer              :: iuser(1:3), ic05rbf, Nii, id02bjf, iflag, lorder, ipqadd, ireadqfms
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

  if( allocated( oqfmsim ) ) deallocate( oqfmsim )
  if( allocated( oqfmsin ) ) deallocate( oqfmsin )

  if( allocated( opp     ) ) deallocate( opp     )
  if( allocated( oqq     ) ) deallocate( oqq     )

  if( allocated( oiota   ) ) deallocate( oiota   )
  if( allocated( oflux   ) ) deallocate( oflux   )
  if( allocated( dflux   ) ) deallocate( dflux   )

  if( allocated( orcmn   ) ) deallocate( orcmn   )
  if( allocated( orsmn   ) ) deallocate( orsmn   )
  if( allocated( otcmn   ) ) deallocate( otcmn   )
  if( allocated( otsmn   ) ) deallocate( otsmn   )

  if( allocated( drcmn   ) ) deallocate( drcmn   )
  if( allocated( drsmn   ) ) deallocate( drsmn   )
  if( allocated( dtcmn   ) ) deallocate( dtcmn   )
  if( allocated( dtsmn   ) ) deallocate( dtsmn   )

  if( allocated( oPxyz   ) ) deallocate( oPxyz   )

  if( allocated( oqfmscplxin  ) ) deallocate( oqfmscplxin  )
  if( allocated( oqfmscplxout ) ) deallocate( oqfmscplxout )

  if( allocated( ijoteta ) ) deallocate( ijoteta )
  if( allocated( ijozeta ) ) deallocate( ijozeta )

  if( allocated( pp      ) ) deallocate( pp      )
  if( allocated( qq      ) ) deallocate( qq      )

  if( allocated( iota    ) ) deallocate( iota    )
  if( allocated( flux    ) ) deallocate( flux    )
  if( allocated( Pxyz    ) ) deallocate( Pxyz    )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  inquire( file = "."//trim(hmnfile)//"."//skperp//".qfms", exist = qfms_exist )
  
  if( qfms_exist ) then
   
   cput = MPI_WTIME()
   
   open( munit, file = "."//trim(hmnfile)//"."//skperp//".qfms", status='unknown', form='unformatted' )
   
   read(munit) oNqfms, oqfmsmn, oqfmsMpol, oqfmsNtor
   
!  write(ounit,1030) cput-tstart, "."//trim(hmnfile)//"."//skperp//".qfms", oNqfms, oqfmsMpol, oqfmsNtor
   
1030 format("readqfms  : ",f10.1,"s : reading ",a," ; oNqfms =",i5," ; ( oM, oN ) = (",i4,",",i4," ) ;")
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmsim ) ) then
    write(0,'("macros    : 0123456789 : oqfmsim already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmsim already allocated ;'
   endif
!#endif

   allocate( oqfmsim(1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmsim ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmsim ;'
   endif
!#endif

   oqfmsim(1:oqfmsmn) = 0 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmsin ) ) then
    write(0,'("macros    : 0123456789 : oqfmsin already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmsin already allocated ;'
   endif
!#endif

   allocate( oqfmsin(1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmsin ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmsin ;'
   endif
!#endif

   oqfmsin(1:oqfmsmn) = 0 


   
   read(munit) oqfmsim, oqfmsin
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( opp     ) ) then
    write(0,'("macros    : 0123456789 : opp     already allocated ;")') 
    stop      'macros    : 0123456789 : opp     already allocated ;'
   endif
!#endif

   allocate( opp    (1:oNqfms ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating opp     ;")') 
    stop      'macros   : 0123456789 : error allocating opp     ;'
   endif
!#endif

   opp    (1:oNqfms ) = 0 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqq     ) ) then
    write(0,'("macros    : 0123456789 : oqq     already allocated ;")') 
    stop      'macros    : 0123456789 : oqq     already allocated ;'
   endif
!#endif

   allocate( oqq    (1:oNqfms ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqq     ;")') 
    stop      'macros   : 0123456789 : error allocating oqq     ;'
   endif
!#endif

   oqq    (1:oNqfms ) = 0 


   
   read(munit) opp    , oqq
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oiota ) ) then
    write(0,'("macros    : 0123456789 : oiota already allocated ;")') 
    stop      'macros    : 0123456789 : oiota already allocated ;'
   endif
!#endif

   allocate( oiota(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oiota ;")') 
    stop      'macros   : 0123456789 : error allocating oiota ;'
   endif
!#endif

   oiota(1:oNqfms) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oflux ) ) then
    write(0,'("macros    : 0123456789 : oflux already allocated ;")') 
    stop      'macros    : 0123456789 : oflux already allocated ;'
   endif
!#endif

   allocate( oflux(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oflux ;")') 
    stop      'macros   : 0123456789 : error allocating oflux ;'
   endif
!#endif

   oflux(1:oNqfms) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dflux ) ) then
    write(0,'("macros    : 0123456789 : dflux already allocated ;")') 
    stop      'macros    : 0123456789 : dflux already allocated ;'
   endif
!#endif

   allocate( dflux(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dflux ;")') 
    stop      'macros   : 0123456789 : error allocating dflux ;'
   endif
!#endif

   dflux(1:oNqfms) = zero 

 ! this is not read in from file, but it is later used for radial interpolation of the flux against iota ; 02/19/2021 ;
   
   read(munit) oiota
   read(munit) oflux
   
   oflux = oflux - oflux(1) ; oflux = oflux / oflux(oNqfms) ! normalizing flux ; 02/19/2021 ;
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( orcmn ) ) then
    write(0,'("macros    : 0123456789 : orcmn already allocated ;")') 
    stop      'macros    : 0123456789 : orcmn already allocated ;'
   endif
!#endif

   allocate( orcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating orcmn ;")') 
    stop      'macros   : 0123456789 : error allocating orcmn ;'
   endif
!#endif

   orcmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( orsmn ) ) then
    write(0,'("macros    : 0123456789 : orsmn already allocated ;")') 
    stop      'macros    : 0123456789 : orsmn already allocated ;'
   endif
!#endif

   allocate( orsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating orsmn ;")') 
    stop      'macros   : 0123456789 : error allocating orsmn ;'
   endif
!#endif

   orsmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( otcmn ) ) then
    write(0,'("macros    : 0123456789 : otcmn already allocated ;")') 
    stop      'macros    : 0123456789 : otcmn already allocated ;'
   endif
!#endif

   allocate( otcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating otcmn ;")') 
    stop      'macros   : 0123456789 : error allocating otcmn ;'
   endif
!#endif

   otcmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( otsmn ) ) then
    write(0,'("macros    : 0123456789 : otsmn already allocated ;")') 
    stop      'macros    : 0123456789 : otsmn already allocated ;'
   endif
!#endif

   allocate( otsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating otsmn ;")') 
    stop      'macros   : 0123456789 : error allocating otsmn ;'
   endif
!#endif

   otsmn(1:oNqfms,1:oqfmsmn) = zero 


   
   read(munit) orcmn
   read(munit) orsmn
   read(munit) otcmn
   read(munit) otsmn
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oPxyz ) ) then
    write(0,'("macros    : 0123456789 : oPxyz already allocated ;")') 
    stop      'macros    : 0123456789 : oPxyz already allocated ;'
   endif
!#endif

   allocate( oPxyz(1:3,0:2,1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oPxyz ;")') 
    stop      'macros   : 0123456789 : error allocating oPxyz ;'
   endif
!#endif

   oPxyz(1:3,0:2,1:oNqfms) = zero 

 ! these are just starting points for Poincare plots; 02/19/2021 ;
   
   read(munit) oPxyz
   
   close(munit)
   
  else ! matches if( qfms_exist ) ; 12/17/20;
   
   oNqfms = 2 ; oqfmsmn = 1 ; oqfmsMpol = 0 ; oqfmsNtor = 0
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmsim ) ) then
    write(0,'("macros    : 0123456789 : oqfmsim already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmsim already allocated ;'
   endif
!#endif

   allocate( oqfmsim(1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmsim ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmsim ;'
   endif
!#endif

   oqfmsim(1:oqfmsmn) = 0 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmsin ) ) then
    write(0,'("macros    : 0123456789 : oqfmsin already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmsin already allocated ;'
   endif
!#endif

   allocate( oqfmsin(1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmsin ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmsin ;'
   endif
!#endif

   oqfmsin(1:oqfmsmn) = 0 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( opp     ) ) then
    write(0,'("macros    : 0123456789 : opp     already allocated ;")') 
    stop      'macros    : 0123456789 : opp     already allocated ;'
   endif
!#endif

   allocate( opp    (1:oNqfms ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating opp     ;")') 
    stop      'macros   : 0123456789 : error allocating opp     ;'
   endif
!#endif

   opp    (1:oNqfms ) = (/ 0, 1 /) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqq     ) ) then
    write(0,'("macros    : 0123456789 : oqq     already allocated ;")') 
    stop      'macros    : 0123456789 : oqq     already allocated ;'
   endif
!#endif

   allocate( oqq    (1:oNqfms ), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqq     ;")') 
    stop      'macros   : 0123456789 : error allocating oqq     ;'
   endif
!#endif

   oqq    (1:oNqfms ) = (/ 1, 1 /) 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oiota ) ) then
    write(0,'("macros    : 0123456789 : oiota already allocated ;")') 
    stop      'macros    : 0123456789 : oiota already allocated ;'
   endif
!#endif

   allocate( oiota(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oiota ;")') 
    stop      'macros   : 0123456789 : error allocating oiota ;'
   endif
!#endif

   oiota(1:oNqfms) = (/ zero, one /) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oflux ) ) then
    write(0,'("macros    : 0123456789 : oflux already allocated ;")') 
    stop      'macros    : 0123456789 : oflux already allocated ;'
   endif
!#endif

   allocate( oflux(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oflux ;")') 
    stop      'macros   : 0123456789 : error allocating oflux ;'
   endif
!#endif

   oflux(1:oNqfms) = (/ zero, one /) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dflux ) ) then
    write(0,'("macros    : 0123456789 : dflux already allocated ;")') 
    stop      'macros    : 0123456789 : dflux already allocated ;'
   endif
!#endif

   allocate( dflux(1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dflux ;")') 
    stop      'macros   : 0123456789 : error allocating dflux ;'
   endif
!#endif

   dflux(1:oNqfms) = (/  one, one /) 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( orcmn ) ) then
    write(0,'("macros    : 0123456789 : orcmn already allocated ;")') 
    stop      'macros    : 0123456789 : orcmn already allocated ;'
   endif
!#endif

   allocate( orcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating orcmn ;")') 
    stop      'macros   : 0123456789 : error allocating orcmn ;'
   endif
!#endif

   orcmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( orsmn ) ) then
    write(0,'("macros    : 0123456789 : orsmn already allocated ;")') 
    stop      'macros    : 0123456789 : orsmn already allocated ;'
   endif
!#endif

   allocate( orsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating orsmn ;")') 
    stop      'macros   : 0123456789 : error allocating orsmn ;'
   endif
!#endif

   orsmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( otcmn ) ) then
    write(0,'("macros    : 0123456789 : otcmn already allocated ;")') 
    stop      'macros    : 0123456789 : otcmn already allocated ;'
   endif
!#endif

   allocate( otcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating otcmn ;")') 
    stop      'macros   : 0123456789 : error allocating otcmn ;'
   endif
!#endif

   otcmn(1:oNqfms,1:oqfmsmn) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( otsmn ) ) then
    write(0,'("macros    : 0123456789 : otsmn already allocated ;")') 
    stop      'macros    : 0123456789 : otsmn already allocated ;'
   endif
!#endif

   allocate( otsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating otsmn ;")') 
    stop      'macros   : 0123456789 : error allocating otsmn ;'
   endif
!#endif

   otsmn(1:oNqfms,1:oqfmsmn) = zero 


   
   orcmn(1,     1) = zero
   orcmn(oNqfms,1) =  one
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oPxyz ) ) then
    write(0,'("macros    : 0123456789 : oPxyz already allocated ;")') 
    stop      'macros    : 0123456789 : oPxyz already allocated ;'
   endif
!#endif

   allocate( oPxyz(1:3,0:2,1:oNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oPxyz ;")') 
    stop      'macros   : 0123456789 : error allocating oPxyz ;'
   endif
!#endif

   oPxyz(1:3,0:2,1:oNqfms) = zero 


   
  endif ! end of if( qfms_exist ) ; 12/17/20;
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( drcmn ) ) then
    write(0,'("macros    : 0123456789 : drcmn already allocated ;")') 
    stop      'macros    : 0123456789 : drcmn already allocated ;'
   endif
!#endif

   allocate( drcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating drcmn ;")') 
    stop      'macros   : 0123456789 : error allocating drcmn ;'
   endif
!#endif

   drcmn(1:oNqfms,1:oqfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( drsmn ) ) then
    write(0,'("macros    : 0123456789 : drsmn already allocated ;")') 
    stop      'macros    : 0123456789 : drsmn already allocated ;'
   endif
!#endif

   allocate( drsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating drsmn ;")') 
    stop      'macros   : 0123456789 : error allocating drsmn ;'
   endif
!#endif

   drsmn(1:oNqfms,1:oqfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dtcmn ) ) then
    write(0,'("macros    : 0123456789 : dtcmn already allocated ;")') 
    stop      'macros    : 0123456789 : dtcmn already allocated ;'
   endif
!#endif

   allocate( dtcmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dtcmn ;")') 
    stop      'macros   : 0123456789 : error allocating dtcmn ;'
   endif
!#endif

   dtcmn(1:oNqfms,1:oqfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( dtsmn ) ) then
    write(0,'("macros    : 0123456789 : dtsmn already allocated ;")') 
    stop      'macros    : 0123456789 : dtsmn already allocated ;'
   endif
!#endif

   allocate( dtsmn(1:oNqfms,1:oqfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating dtsmn ;")') 
    stop      'macros   : 0123456789 : error allocating dtsmn ;'
   endif
!#endif

   dtsmn(1:oNqfms,1:oqfmsmn) = zero 


  
  do ii = 1, oqfmsmn
   ie01bef = 0 ; call E01BEF( oNqfms, oflux(1:oNqfms), orcmn(1:oNqfms,ii), drcmn(1:oNqfms,ii), ie01bef ) ! construct radial interpolation ; 12/17/20;
   ie01bef = 0 ; call E01BEF( oNqfms, oflux(1:oNqfms), orsmn(1:oNqfms,ii), drsmn(1:oNqfms,ii), ie01bef )
   ie01bef = 0 ; call E01BEF( oNqfms, oflux(1:oNqfms), otcmn(1:oNqfms,ii), dtcmn(1:oNqfms,ii), ie01bef )
   ie01bef = 0 ; call E01BEF( oNqfms, oflux(1:oNqfms), otsmn(1:oNqfms,ii), dtsmn(1:oNqfms,ii), ie01bef )
  enddo
  
  ie01bef = 0 ; call E01BEF( oNqfms, oiota(1:oNqfms), oflux(1:oNqfms), dflux(1:oNqfms), ie01bef ) ! this is used for initial guesses ; 02/19/2021 ;
  
  oqfmsNteta = 8 * max(oqfmsMpol,1) ; oqfmsNzeta = 8 * max(oqfmsNtor,1) ; oqfmsNtz = oqfmsNteta * oqfmsNzeta
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmscplxin  ) ) then
    write(0,'("macros    : 0123456789 : oqfmscplxin  already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmscplxin  already allocated ;'
   endif
!#endif

   allocate( oqfmscplxin (1:oqfmsNteta,1:oqfmsNzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmscplxin  ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmscplxin  ;'
   endif
!#endif

   oqfmscplxin (1:oqfmsNteta,1:oqfmsNzeta) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( oqfmscplxout ) ) then
    write(0,'("macros    : 0123456789 : oqfmscplxout already allocated ;")') 
    stop      'macros    : 0123456789 : oqfmscplxout already allocated ;'
   endif
!#endif

   allocate( oqfmscplxout(1:oqfmsNteta,1:oqfmsNzeta), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating oqfmscplxout ;")') 
    stop      'macros   : 0123456789 : error allocating oqfmscplxout ;'
   endif
!#endif

   oqfmscplxout(1:oqfmsNteta,1:oqfmsNzeta) = zero 


  
  oqfmsplanb = fftw_plan_dft_2d( oqfmsNzeta, oqfmsNteta, oqfmscplxin, oqfmscplxout, fftw_backward, fftw_measure + fftw_destroy_input )
  oqfmsplanf = fftw_plan_dft_2d( oqfmsNzeta, oqfmsNteta, oqfmscplxin, oqfmscplxout, fftw_forward , fftw_measure + fftw_destroy_input )
  
   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ijoteta ) ) then
    write(0,'("macros    : 0123456789 : ijoteta already allocated ;")') 
    stop      'macros    : 0123456789 : ijoteta already allocated ;'
   endif
!#endif

   allocate( ijoteta(1:oqfmsNtz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ijoteta ;")') 
    stop      'macros   : 0123456789 : error allocating ijoteta ;'
   endif
!#endif

   ijoteta(1:oqfmsNtz) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ijozeta ) ) then
    write(0,'("macros    : 0123456789 : ijozeta already allocated ;")') 
    stop      'macros    : 0123456789 : ijozeta already allocated ;'
   endif
!#endif

   allocate( ijozeta(1:oqfmsNtz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ijozeta ;")') 
    stop      'macros   : 0123456789 : error allocating ijozeta ;'
   endif
!#endif

   ijozeta(1:oqfmsNtz) = zero 


  
  do kk = 0, oqfmsNzeta-1 ; zeta = kk * pi2 / oqfmsNzeta
   do jj = 0, oqfmsNteta-1 ; teta = jj * pi2 / oqfmsNteta ; jk = 1 + jj + kk * oqfmsNteta ; ijoteta(jk) = teta
    ;                      ;                              ;                               ; ijozeta(jk) = zeta
   enddo
  enddo
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  inquire( file = "."//trim(hmnfile)//"."//skperp//".pq", exist = pq_exist )  

  if( .not.pq_exist ) NFarey = 0
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  select case( NFarey )
    
  case(  : -1 ) ! read selection of (p,q)'s from file ; 12/17/20;
   
   open( munit, file = "."//trim(hmnfile)//"."//skperp//".pq", status='unknown', form='formatted' )   
   Nqfms = 0
   do
    read( munit, *, iostat=iostat) lpp, lqq, lorder
    if( iostat.eq.0 ) Nqfms = Nqfms + 1
    if( iostat.ne.0 ) exit
   enddo
   close( munit )
   
   if( ipqadd.gt.0 ) Nqfms = Nqfms + 1
   
   open( munit, file = "."//trim(hmnfile)//"."//skperp//".pq", status='unknown', form='formatted' )
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( pp ) ) then
    write(0,'("macros    : 0123456789 : pp already allocated ;")') 
    stop      'macros    : 0123456789 : pp already allocated ;'
   endif
!#endif

   allocate( pp(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating pp ;")') 
    stop      'macros   : 0123456789 : error allocating pp ;'
   endif
!#endif

   pp(1:Nqfms) = 0 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qq ) ) then
    write(0,'("macros    : 0123456789 : qq already allocated ;")') 
    stop      'macros    : 0123456789 : qq already allocated ;'
   endif
!#endif

   allocate( qq(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qq ;")') 
    stop      'macros   : 0123456789 : error allocating qq ;'
   endif
!#endif

   qq(1:Nqfms) = 0 


   
   if( ipqadd.eq.0 ) lNqfms = Nqfms
   if( ipqadd.gt.0 ) lNqfms = Nqfms-1
   
   do ipq = 1, lNqfms ; read(munit,*) lpp, lqq ; pp(ipq) = lpp ; qq(ipq) = lqq
   enddo
   
   close( munit )
   
   cput = MPI_WTIME()
   
   if( ipqadd.gt.0 ) then
    pp(ipqadd+1:Nqfms) = pp(ipqadd:Nqfms-1)
    qq(ipqadd+1:Nqfms) = qq(ipqadd:Nqfms-1)
    pp(ipqadd+1) = pp(ipqadd) + pp(ipqadd+2)
    qq(ipqadd+1) = qq(ipqadd) + qq(ipqadd+2)
    if( ireadqfms.eq.1 ) write(ounit,1000) cput-tstart, pp(ipqadd+0), qq(ipqadd+0), pp(ipqadd+1), qq(ipqadd+1), pp(ipqadd+2), qq(ipqadd+2)
   endif
   
1000 format("readqfms  : "f10.1"s : adding "3("("i3","i4" ), ")" ; ")

!   write(ounit,'("readqfms  : "10x"s : ("i3","i4" ) ; ")') ( pp(ipq), qq(ipq), ipq = 1, Nqfms )
   
   do ipq = 2, Nqfms
    if( pp(ipq)*one/qq(ipq).le.pp(ipq-1)*one/qq(ipq-1) ) write(ounit,'("readqfms  : " 10x"s : "4i6)') pp(ipq-1), qq(ipq-1), pp(ipq), qq(ipq)
    if(pp(ipq)*one/qq(ipq).le.pp(ipq-1)*one/qq(ipq-1) ) then
          write(ounit,'("readqfms  :    fatal    : pp(ipq)*one/qq(ipq).le.pp(ipq-1)*one/qq(ipq-1)  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
   enddo
   
  case( 0:    ) ! construct Farey tree ; 12/17/20;
   
   Nqfms = 2
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( pp ) ) then
    write(0,'("macros    : 0123456789 : pp already allocated ;")') 
    stop      'macros    : 0123456789 : pp already allocated ;'
   endif
!#endif

   allocate( pp(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating pp ;")') 
    stop      'macros   : 0123456789 : error allocating pp ;'
   endif
!#endif

   pp(1:Nqfms) = (/ 0, 1 /) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qq ) ) then
    write(0,'("macros    : 0123456789 : qq already allocated ;")') 
    stop      'macros    : 0123456789 : qq already allocated ;'
   endif
!#endif

   allocate( qq(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qq ;")') 
    stop      'macros   : 0123456789 : error allocating qq ;'
   endif
!#endif

   qq(1:Nqfms) = (/ 1, 1 /) 


   
   do iFarey = 1, NFarey
    
    lNqfms = 2 * Nqfms-1
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( lp ) ) then
    write(0,'("macros    : 0123456789 : lp already allocated ;")') 
    stop      'macros    : 0123456789 : lp already allocated ;'
   endif
!#endif

   allocate( lp(1:lNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating lp ;")') 
    stop      'macros   : 0123456789 : error allocating lp ;'
   endif
!#endif

   lp(1:lNqfms) = 0 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( lq ) ) then
    write(0,'("macros    : 0123456789 : lq already allocated ;")') 
    stop      'macros    : 0123456789 : lq already allocated ;'
   endif
!#endif

   allocate( lq(1:lNqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating lq ;")') 
    stop      'macros   : 0123456789 : error allocating lq ;'
   endif
!#endif

   lq(1:lNqfms) = 0 


    
    lp(1:lNqfms:2) = pp(1:Nqfms) ; lp(2:lNqfms-1:2) = pp(1:Nqfms-1) + pp(2:Nqfms)
    lq(1:lNqfms:2) = qq(1:Nqfms) ; lq(2:lNqfms-1:2) = qq(1:Nqfms-1) + qq(2:Nqfms)
    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( pp  ) ) then
    write(0,'("macros   : 0123456789 : pp  not already allocated ;")') 
    stop      'macros   : 0123456789 : pp  not already allocated ;'
   endif
!#endif

   deallocate( pp , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating pp  ;")') 
    stop      'macros   : 0123456789 : error de-allocating pp  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( qq  ) ) then
    write(0,'("macros   : 0123456789 : qq  not already allocated ;")') 
    stop      'macros   : 0123456789 : qq  not already allocated ;'
   endif
!#endif

   deallocate( qq , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating qq  ;")') 
    stop      'macros   : 0123456789 : error de-allocating qq  ;'
   endif
!#endif


    
    Nqfms = lNqfms
    
     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( pp ) ) then
    write(0,'("macros    : 0123456789 : pp already allocated ;")') 
    stop      'macros    : 0123456789 : pp already allocated ;'
   endif
!#endif

   allocate( pp(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating pp ;")') 
    stop      'macros   : 0123456789 : error allocating pp ;'
   endif
!#endif

   pp(1:Nqfms) = lp(1:Nqfms) 


     ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qq ) ) then
    write(0,'("macros    : 0123456789 : qq already allocated ;")') 
    stop      'macros    : 0123456789 : qq already allocated ;'
   endif
!#endif

   allocate( qq(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qq ;")') 
    stop      'macros   : 0123456789 : error allocating qq ;'
   endif
!#endif

   qq(1:Nqfms) = lq(1:Nqfms) 


    
     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( lp  ) ) then
    write(0,'("macros   : 0123456789 : lp  not already allocated ;")') 
    stop      'macros   : 0123456789 : lp  not already allocated ;'
   endif
!#endif

   deallocate( lp , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating lp  ;")') 
    stop      'macros   : 0123456789 : error de-allocating lp  ;'
   endif
!#endif


     ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( lq  ) ) then
    write(0,'("macros   : 0123456789 : lq  not already allocated ;")') 
    stop      'macros   : 0123456789 : lq  not already allocated ;'
   endif
!#endif

   deallocate( lq , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating lq  ;")') 
    stop      'macros   : 0123456789 : error de-allocating lq  ;'
   endif
!#endif


    
   enddo ! end of do iFarey ; 12/17/20;
   
  end select
  
  cput = MPI_WTIME()
  
!  write(ounit,'("readqfms  : ",f10.1,"s : Nqfms =",i4," ;")') cput-tstart, Nqfms
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( Laction )
   
  case( 0 ) ! Laction = 0 ; nothing will be done regarding constructing QFM surfaces ; 02/19/2021 ;
   
    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( pp  ) ) then
    write(0,'("macros   : 0123456789 : pp  not already allocated ;")') 
    stop      'macros   : 0123456789 : pp  not already allocated ;'
   endif
!#endif

   deallocate( pp , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating pp  ;")') 
    stop      'macros   : 0123456789 : error de-allocating pp  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( qq  ) ) then
    write(0,'("macros   : 0123456789 : qq  not already allocated ;")') 
    stop      'macros   : 0123456789 : qq  not already allocated ;'
   endif
!#endif

   deallocate( qq , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating qq  ;")') 
    stop      'macros   : 0123456789 : error de-allocating qq  ;'
   endif
!#endif


   
   Nqfms = oNqfms
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Pxyz ) ) then
    write(0,'("macros    : 0123456789 : Pxyz already allocated ;")') 
    stop      'macros    : 0123456789 : Pxyz already allocated ;'
   endif
!#endif

   allocate( Pxyz(1:3,0:2,1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Pxyz ;")') 
    stop      'macros   : 0123456789 : error allocating Pxyz ;'
   endif
!#endif

   Pxyz(1:3,0:2,1:Nqfms) = oPxyz 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( pp ) ) then
    write(0,'("macros    : 0123456789 : pp already allocated ;")') 
    stop      'macros    : 0123456789 : pp already allocated ;'
   endif
!#endif

   allocate( pp(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating pp ;")') 
    stop      'macros   : 0123456789 : error allocating pp ;'
   endif
!#endif

   pp(1:Nqfms) = opp(1:oNqfms) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( qq ) ) then
    write(0,'("macros    : 0123456789 : qq already allocated ;")') 
    stop      'macros    : 0123456789 : qq already allocated ;'
   endif
!#endif

   allocate( qq(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating qq ;")') 
    stop      'macros   : 0123456789 : error allocating qq ;'
   endif
!#endif

   qq(1:Nqfms) = oqq(1:oNqfms) 


   
  case default ! Laction.ne.0 ; 
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( iota ) ) then
    write(0,'("macros    : 0123456789 : iota already allocated ;")') 
    stop      'macros    : 0123456789 : iota already allocated ;'
   endif
!#endif

   allocate( iota(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating iota ;")') 
    stop      'macros   : 0123456789 : error allocating iota ;'
   endif
!#endif

   iota(1:Nqfms) = pp(1:Nqfms) * one / qq(1:Nqfms) 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( flux ) ) then
    write(0,'("macros    : 0123456789 : flux already allocated ;")') 
    stop      'macros    : 0123456789 : flux already allocated ;'
   endif
!#endif

   allocate( flux(1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating flux ;")') 
    stop      'macros   : 0123456789 : error allocating flux ;'
   endif
!#endif

   flux(1:Nqfms) = pp(1:Nqfms) * one / qq(1:Nqfms) 


   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Pxyz ) ) then
    write(0,'("macros    : 0123456789 : Pxyz already allocated ;")') 
    stop      'macros    : 0123456789 : Pxyz already allocated ;'
   endif
!#endif

   allocate( Pxyz(1:3,0:2,1:Nqfms), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Pxyz ;")') 
    stop      'macros   : 0123456789 : error allocating Pxyz ;'
   endif
!#endif

   Pxyz(1:3,0:2,1:Nqfms) = zero 

 ! Nqfms was set in global:readinput ; this will be assigned in action ; 20/10/21 ;
   
  end select
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( oPxyz  ) ) then
    write(0,'("macros   : 0123456789 : oPxyz  not already allocated ;")') 
    stop      'macros   : 0123456789 : oPxyz  not already allocated ;'
   endif
!#endif

   deallocate( oPxyz , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating oPxyz  ;")') 
    stop      'macros   : 0123456789 : error de-allocating oPxyz  ;'
   endif
!#endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine readqfms

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
