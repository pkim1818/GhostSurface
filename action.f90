














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine action( ipqadd, iaction )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, one, two, pi2, small, ounit, munit, myid, tstart, &
                      hmnfile, &
                      Nqfms, pp, qq, iota, flux, &
                      oNqfms, oiota, oflux, dflux, &
                      qfmsNteta, qfmsNzeta, qfmsNtz, &
                      qfmsMpol, qfmsNtor, qfmsmn, qfmsim, qfmsin, &
                      oqfmsMpol, oqfmsNtor, orcmn, orsmn, otcmn, otsmn, opp, oqq, &
                      pqMpol, pqNtor, NDpq, &
                      Pxyz, &
                     !iin, jin, npint, opint, Ni, Nj, Tijk, & ! 02/19/2021 ;
                     !chi, PX, PY, lamda, mu, cc, e02defrwrk, e02defiwrk, & ! 02/19/2021 ;
                      Linterpol, skperp

  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  LOGICAL              :: qfms_exist, Lsetold, Lsetpq
  integer              :: astat, ipq, lpp, lqq, jpq, qN, Ndof, iuser(1:3), ic05qbf, ii, jj, kk, nn, mm, jk, fM, qfM, ic06faf, npmodq, opq
  integer              :: ij, ib, jb, kb, ie02def, ie01bgf, iaction, ipqadd
  real                 :: cput, cpul, xtol, liota, ruser(1:1), dforce, teta, alpha, zeta, nzq, loflux, ldflux
  real                 :: rjk(1:qfmsNtz), tjk(1:qfmsNtz), stz(1:3), xyz(1:3,0:3), lx, ly, win(1:16)
  real                 :: efmn(1:qfmsmn), ofmn(1:qfmsmn), cfmn(1:qfmsmn), sfmn(1:qfmsmn), ijreal(1:qfmsNtz), ijimag(1:qfmsNtz)
  real                 :: temp, lxx(1:1), lyy(1:1), lf(1:1), lflux
  real   , allocatable :: rcn(:,:), rsn(:,:), tcn(:,:), tsn(:,:), nuu(:), xx(:), ff(:)
  real   , allocatable :: rcmn(:,:), rsmn(:,:), tcmn(:,:), tsmn(:,:), rajk(:,:), tajk(:,:), rwork(:)

  EXTERNAL             :: actiongradient

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cput = MPI_WTIME()
 
  cpul = cput
 
  inquire( file = "."//trim(hmnfile)//".qfms", exist = qfms_exist )

  Lsetold = .false. ! not presently used; 05/21/21;

  if( qfms_exist .and. oqfmsMpol.eq.qfmsMpol .and. oqfmsNtor.eq.qfmsNtor ) Lsetold = .true. ! not presently used; 05/21/21;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  xtol = small ; fM = NDpq * pqMpol
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
 
! write(ounit,'("action  : "10x"s : Nqfms ="i3", oNqfms ="i3", qfsmn ="i3", oqfsmn ="i3" ;")') Nqfms, oNqfms, qfsmn, oqfsmn

   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rcmn ) ) then
    write(0,'("macros    : 0123456789 : rcmn already allocated ;")') 
    stop      'macros    : 0123456789 : rcmn already allocated ;'
   endif
!#endif

   allocate( rcmn(1:Nqfms,1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rcmn ;")') 
    stop      'macros   : 0123456789 : error allocating rcmn ;'
   endif
!#endif

   rcmn(1:Nqfms,1:qfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rsmn ) ) then
    write(0,'("macros    : 0123456789 : rsmn already allocated ;")') 
    stop      'macros    : 0123456789 : rsmn already allocated ;'
   endif
!#endif

   allocate( rsmn(1:Nqfms,1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rsmn ;")') 
    stop      'macros   : 0123456789 : error allocating rsmn ;'
   endif
!#endif

   rsmn(1:Nqfms,1:qfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( tcmn ) ) then
    write(0,'("macros    : 0123456789 : tcmn already allocated ;")') 
    stop      'macros    : 0123456789 : tcmn already allocated ;'
   endif
!#endif

   allocate( tcmn(1:Nqfms,1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating tcmn ;")') 
    stop      'macros   : 0123456789 : error allocating tcmn ;'
   endif
!#endif

   tcmn(1:Nqfms,1:qfmsmn) = zero 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( tsmn ) ) then
    write(0,'("macros    : 0123456789 : tsmn already allocated ;")') 
    stop      'macros    : 0123456789 : tsmn already allocated ;'
   endif
!#endif

   allocate( tsmn(1:Nqfms,1:qfmsmn), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating tsmn ;")') 
    stop      'macros   : 0123456789 : error allocating tsmn ;'
   endif
!#endif

   tsmn(1:Nqfms,1:qfmsmn) = zero 


  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  do ipq = 1, Nqfms ; lpp = pp(ipq) ; lqq = qq(ipq) ; liota = iota(ipq) ; qN = lqq * pqNtor ; qfM = lqq * fM
   
!   if( Lsetold ) then
!    Lsetpq = .false.
!    do opq = 1, oNqfms
!     if( opp(opq).eq.lpp .and. oqq(opq).eq.lqq ) then
!      write(ounit,'("action   : "10x"s : setting to old ;")') 
!      rcmn(ipq,1:qfmsmn) = orcmn(opq,1:qfmsmn)
!      rsmn(ipq,1:qfmsmn) = orsmn(opq,1:qfmsmn)
!      tcmn(ipq,1:qfmsmn) = otcmn(opq,1:qfmsmn)
!      tsmn(ipq,1:qfmsmn) = otsmn(opq,1:qfmsmn)
!      Lsetpq = .true.
!     endif
!    enddo
!    if( Lsetpq ) cycle
!   endif

    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rcn ) ) then
    write(0,'("macros    : 0123456789 : rcn already allocated ;")') 
    stop      'macros    : 0123456789 : rcn already allocated ;'
   endif
!#endif

   allocate( rcn(0:qN,0:fM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rcn ;")') 
    stop      'macros   : 0123456789 : error allocating rcn ;'
   endif
!#endif

   rcn(0:qN,0:fM-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rsn ) ) then
    write(0,'("macros    : 0123456789 : rsn already allocated ;")') 
    stop      'macros    : 0123456789 : rsn already allocated ;'
   endif
!#endif

   allocate( rsn(0:qN,0:fM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rsn ;")') 
    stop      'macros   : 0123456789 : error allocating rsn ;'
   endif
!#endif

   rsn(0:qN,0:fM-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( tcn ) ) then
    write(0,'("macros    : 0123456789 : tcn already allocated ;")') 
    stop      'macros    : 0123456789 : tcn already allocated ;'
   endif
!#endif

   allocate( tcn(0:qN,0:fM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating tcn ;")') 
    stop      'macros   : 0123456789 : error allocating tcn ;'
   endif
!#endif

   tcn(0:qN,0:fM-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( tsn ) ) then
    write(0,'("macros    : 0123456789 : tsn already allocated ;")') 
    stop      'macros    : 0123456789 : tsn already allocated ;'
   endif
!#endif

   allocate( tsn(0:qN,0:fM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating tsn ;")') 
    stop      'macros   : 0123456789 : error allocating tsn ;'
   endif
!#endif

   tsn(0:qN,0:fM-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( nuu ) ) then
    write(0,'("macros    : 0123456789 : nuu already allocated ;")') 
    stop      'macros    : 0123456789 : nuu already allocated ;'
   endif
!#endif

   allocate( nuu(     0:fM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating nuu ;")') 
    stop      'macros   : 0123456789 : error allocating nuu ;'
   endif
!#endif

   nuu(     0:fM-1) = zero 


   
   mm = 1 ; ie01bgf = 1 ; call E01BGF( oNqfms, oiota(1:oNqfms), oflux(1:oNqfms), dflux(1:oNqfms), mm, liota, loflux, ldflux, ie01bgf )

   rcn(0,0) = loflux ! initial guess for radial coordinate (toroidal flux) from flux as a function of rotational-transform ; 05/21/21;
   tcn(0,0) =  zero - (pi2/lqq) / fM
   
   Ndof = 1+qN + qN + 1+qN + qN + 1
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( xx ) ) then
    write(0,'("macros    : 0123456789 : xx already allocated ;")') 
    stop      'macros    : 0123456789 : xx already allocated ;'
   endif
!#endif

   allocate( xx(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating xx ;")') 
    stop      'macros   : 0123456789 : error allocating xx ;'
   endif
!#endif

   xx(1:Ndof) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( ff ) ) then
    write(0,'("macros    : 0123456789 : ff already allocated ;")') 
    stop      'macros    : 0123456789 : ff already allocated ;'
   endif
!#endif

   allocate( ff(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating ff ;")') 
    stop      'macros   : 0123456789 : error allocating ff ;'
   endif
!#endif

   ff(1:Ndof) = zero 


   
   do jpq = 0, fM-1 ; iuser(1:3) = (/ ipq, jpq, qN /)

    ruser(1) = jpq * (pi2/lqq) / fM ! ruser is target area ;
    
    if( jpq.gt.0 ) then ; rcn(0:qN,jpq) = rcn(0:qN,jpq-1) ! copy across initial guess; 12/17/20;
     ;                  ; rsn(0:qN,jpq) = rsn(0:qN,jpq-1)
     ;                  ; tcn(0:qN,jpq) = tcn(0:qN,jpq-1)
     ;                  ; tsn(0:qN,jpq) = tsn(0:qN,jpq-1)
     ;                  ; nuu(     jpq) = nuu(     jpq-1)
    endif

    tcn(0,jpq) = tcn(0,jpq) + (pi2/lqq) / fM ! estimated shift; 12/17/20;
    
    xx(1                :1+qN           ) = rcn(0:qN,jpq) + one
    xx(1+qN+1           :1+qN+qN        ) = rsn(1:qN,jpq) + one
    xx(1+qN+qN+1        :1+qN+qN+1+qN   ) = tcn(0:qN,jpq) + one
    xx(1+qN+qN+1+qN+1   :1+qN+qN+1+qN+qN) = tsn(1:qN,jpq) + one
    xx(1+qN+qN+1+qN+qN+1                ) = nuu(     jpq) + one
    
    ic05qbf = 1 ; call C05QBF( actiongradient, Ndof, xx(1:Ndof), ff(1:Ndof), xtol, iuser(1:3), ruser(1:1), ic05qbf )

    dforce = sqrt(sum(ff*ff)/Ndof)

    rcn(0:qN,jpq) = xx(1                :1+qN           ) - one
    rsn(1:qN,jpq) = xx(1+qN+1           :1+qN+qN        ) - one
    tcn(0:qN,jpq) = xx(1+qN+qN+1        :1+qN+qN+1+qN   ) - one
    tsn(1:qN,jpq) = xx(1+qN+qN+1+qN+1   :1+qN+qN+1+qN+qN) - one
    nuu(     jpq) = xx(1+qN+qN+1+qN+qN+1                ) - one
    
    cput = MPI_WTIME()

    select case( ic05qbf )
    case( 0 )    ;!write(ounit,1000) cput-tstart, lpp, lqq, jpq, ic05qbf, dforce, nuu(jpq)
    case default ;!write(ounit,1000) cput-tstart, lpp, lqq, jpq, ic05qbf, dforce, nuu(jpq)
    end select    

1000 format("action    : ",f10.1,"s : (",i5,",",i5," ) ; ",i3," : ic05qbf =",i3," ; |df| =",es8.1," ; nuu =",es13.5," ;")

   enddo ! end of do jpq ; 12/17/20;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rajk ) ) then
    write(0,'("macros    : 0123456789 : rajk already allocated ;")') 
    stop      'macros    : 0123456789 : rajk already allocated ;'
   endif
!#endif

   allocate( rajk(0:qfM-1,0:qfmsNzeta-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rajk ;")') 
    stop      'macros   : 0123456789 : error allocating rajk ;'
   endif
!#endif

   rajk(0:qfM-1,0:qfmsNzeta-1) = zero 


    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( tajk ) ) then
    write(0,'("macros    : 0123456789 : tajk already allocated ;")') 
    stop      'macros    : 0123456789 : tajk already allocated ;'
   endif
!#endif

   allocate( tajk(0:qfM-1,0:qfmsNzeta-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating tajk ;")') 
    stop      'macros   : 0123456789 : error allocating tajk ;'
   endif
!#endif

   tajk(0:qfM-1,0:qfmsNzeta-1) = zero 


   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   do jpq = 0, fM-1
    
    do mm = 0, lqq-1
     do kk = 0, qfmsNzeta-1 ; zeta = mm * pi2 + kk * (pi2/qfmsNzeta)
      
      ii = mod( jpq + mm * lpp * fM, qfM ) ; npmodq = ( mm * lpp ) / lqq
      
      ;  nn = 0     ;                       ; rajk(ii,kk) =               rcn(nn,jpq)
      ;                                     ; tajk(ii,kk) =               tcn(nn,jpq) + liota * zeta - npmodq * pi2
      do nn = 1, qN ; nzq = nn * zeta / lqq ; rajk(ii,kk) = rajk(ii,kk) + rcn(nn,jpq) * cos(nzq) + rsn(nn,jpq) * sin(nzq)
       ;            ;                       ; tajk(ii,kk) = tajk(ii,kk) + tcn(nn,jpq) * cos(nzq) + tsn(nn,jpq) * sin(nzq)
      enddo
      
     enddo ! end of do kk ; 12/17/20;
    enddo ! end of do mm ; 12/17/20;
    
   enddo ! end of do jpq ; 12/17/20;
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   do jpq = 0, 2 ; Pxyz(1:3,jpq,ipq) = (/ rajk(jpq*pqMpol,0), tajk(jpq*pqMpol,0), zero /) ! this is used in poincare.h to initialize trajectories ; 12/17/20;
   enddo
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   do ii = 0, qfM-1
    do kk = 0, qfmsNzeta-1 ; tajk(ii,kk) = tajk(ii,kk) - ( ii*pi2/qfM + liota * kk*pi2/qfmsNzeta ) ! hereafter tajk is really lambda ; 12/17/20;
    enddo
   enddo

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rwork ) ) then
    write(0,'("macros    : 0123456789 : rwork already allocated ;")') 
    stop      'macros    : 0123456789 : rwork already allocated ;'
   endif
!#endif

   allocate( rwork(0:qfM-1), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rwork ;")') 
    stop      'macros   : 0123456789 : error allocating rwork ;'
   endif
!#endif

   rwork(0:qfM-1) = zero 

 ! one-dimensional FFT workspace; interpolation in alpha; 12/17/20;
   
   do kk = 0, qfmsNzeta-1 ; zeta = kk * pi2 / qfmsNzeta
    
    call C06FAF( rajk(0:qfM-1,kk), qfM, rwork, ic06faf ) ; rajk(0:qfM-1,kk) = rajk(0:qfM-1,kk) / sqrt(one*qfM) ; rajk(1:qfM-1,kk) = rajk(1:qfM-1,kk) * two
    call C06FAF( tajk(0:qfM-1,kk), qfM, rwork, ic06faf ) ; tajk(0:qfM-1,kk) = tajk(0:qfM-1,kk) / sqrt(one*qfM) ; tajk(1:qfM-1,kk) = tajk(1:qfM-1,kk) * two
    
    do jj = 0, qfmsNteta-1 ; teta = jj * pi2 / qfmsNteta ; alpha = teta - liota * zeta ; jk = 1 + jj + kk*qfmsNteta
     
     ;  mm = 0             ; rjk(jk) =           rajk(mm,kk)
     ;                     ; tjk(jk) =           tajk(mm,kk)
     do mm = 1, lqq*pqMpol ; rjk(jk) = rjk(jk) + rajk(mm,kk) * cos(mm*alpha) - rajk(qfM-mm,kk) * sin(mm*alpha)
      ;                    ; tjk(jk) = tjk(jk) + tajk(mm,kk) * cos(mm*alpha) - tajk(qfM-mm,kk) * sin(mm*alpha)
     enddo
     
     stz(1:3) = (/ rjk(jk), teta + tjk(jk), zeta /) ; call stzxyz( stz(1:3), xyz(1:3,0:3), Linterpol ) ; rjk(jk) = xyz(1,0) ; tjk(jk) = xyz(2,0) - teta
     
    enddo ! end of do jj ; 12/17/20;
    
   enddo ! end of do kk ; 12/17/20;
   
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



   call qfmstfft( qfmsNteta, qfmsNzeta, rjk(1:qfmsNtz), tjk(1:qfmsNtz), &
                  qfmsmn, qfmsim, qfmsin, efmn(1:qfmsmn), ofmn(1:qfmsmn), cfmn(1:qfmsmn), sfmn(1:qfmsmn) )

   do ii = 1, qfmsmn ; rcmn(ipq,ii) = efmn(ii)
    ;                ; rsmn(ipq,ii) = ofmn(ii)
    ;                ; tcmn(ipq,ii) = cfmn(ii)
    ;                ; tsmn(ipq,ii) = sfmn(ii)
   enddo

    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rajk  ) ) then
    write(0,'("macros   : 0123456789 : rajk  not already allocated ;")') 
    stop      'macros   : 0123456789 : rajk  not already allocated ;'
   endif
!#endif

   deallocate( rajk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rajk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rajk  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( tajk  ) ) then
    write(0,'("macros   : 0123456789 : tajk  not already allocated ;")') 
    stop      'macros   : 0123456789 : tajk  not already allocated ;'
   endif
!#endif

   deallocate( tajk , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating tajk  ;")') 
    stop      'macros   : 0123456789 : error de-allocating tajk  ;'
   endif
!#endif


   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rcn  ) ) then
    write(0,'("macros   : 0123456789 : rcn  not already allocated ;")') 
    stop      'macros   : 0123456789 : rcn  not already allocated ;'
   endif
!#endif

   deallocate( rcn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rcn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rcn  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rsn  ) ) then
    write(0,'("macros   : 0123456789 : rsn  not already allocated ;")') 
    stop      'macros   : 0123456789 : rsn  not already allocated ;'
   endif
!#endif

   deallocate( rsn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rsn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rsn  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( tcn  ) ) then
    write(0,'("macros   : 0123456789 : tcn  not already allocated ;")') 
    stop      'macros   : 0123456789 : tcn  not already allocated ;'
   endif
!#endif

   deallocate( tcn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating tcn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating tcn  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( tsn  ) ) then
    write(0,'("macros   : 0123456789 : tsn  not already allocated ;")') 
    stop      'macros   : 0123456789 : tsn  not already allocated ;'
   endif
!#endif

   deallocate( tsn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating tsn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating tsn  ;'
   endif
!#endif


    ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( nuu  ) ) then
    write(0,'("macros   : 0123456789 : nuu  not already allocated ;")') 
    stop      'macros   : 0123456789 : nuu  not already allocated ;'
   endif
!#endif

   deallocate( nuu , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating nuu  ;")') 
    stop      'macros   : 0123456789 : error de-allocating nuu  ;'
   endif
!#endif


   
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
   if( .not.allocated( ff  ) ) then
    write(0,'("macros   : 0123456789 : ff  not already allocated ;")') 
    stop      'macros   : 0123456789 : ff  not already allocated ;'
   endif
!#endif

   deallocate( ff , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating ff  ;")') 
    stop      'macros   : 0123456789 : error de-allocating ff  ;'
   endif
!#endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   cfmn(1:qfmsmn) = + qfmsim(1:qfmsmn) * sfmn(1:qfmsmn)
   sfmn(1:qfmsmn) = - qfmsim(1:qfmsmn) * cfmn(1:qfmsmn)

   call qfmsinvfft( qfmsmn, qfmsim(1:qfmsmn), qfmsin(1:qfmsmn), efmn(1:qfmsmn), ofmn(1:qfmsmn), cfmn(1:qfmsmn), sfmn(1:qfmsmn), &
                    qfmsNteta, qfmsNzeta, ijreal(1:qfmsNtz), ijimag(1:qfmsNtz) )
   
   flux(ipq) = sum( ijreal * ( one + ijimag ) ) / qfmsNtz
   
   cput = MPI_WTIME()

   if( iaction.eq.1 ) write(ounit,1020) cput-tstart, lpp, lqq, Ndof, cput-cpul

1020 format("action    : ",f10.1,"s : constructed (",i5,",",i5," ).qfms ; Ndof =",i6," ;",f10.2,"s ;")

   cpul = cput

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  enddo ! end of do ipq; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  cput = MPI_WTIME()

  if( iaction.eq.1 ) write(ounit,1010) cput-tstart, "."//trim(hmnfile)//"."//skperp//".qfms", Nqfms, qfmsMpol, qfmsNtor

1010 format("action    : ",f10.1,"s : writing ",a," ;  Nqfms =",i5," ; (  M,  N ) = (",i4,",",i4," ) ;")

  open( munit, file = "."//trim(hmnfile)//"."//skperp//".qfms", status='unknown', form='unformatted' )
  write(munit) Nqfms, qfmsmn, qfmsMpol, qfmsNtor
  write(munit) qfmsim, qfmsin
  write(munit) pp, qq
  write(munit) iota
  write(munit) flux
  write(munit) rcmn
  write(munit) rsmn
  write(munit) tcmn
  write(munit) tsmn
  write(munit) Pxyz
  close(munit)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rcmn  ) ) then
    write(0,'("macros   : 0123456789 : rcmn  not already allocated ;")') 
    stop      'macros   : 0123456789 : rcmn  not already allocated ;'
   endif
!#endif

   deallocate( rcmn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rcmn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rcmn  ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( rsmn  ) ) then
    write(0,'("macros   : 0123456789 : rsmn  not already allocated ;")') 
    stop      'macros   : 0123456789 : rsmn  not already allocated ;'
   endif
!#endif

   deallocate( rsmn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating rsmn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating rsmn  ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( tcmn  ) ) then
    write(0,'("macros   : 0123456789 : tcmn  not already allocated ;")') 
    stop      'macros   : 0123456789 : tcmn  not already allocated ;'
   endif
!#endif

   deallocate( tcmn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating tcmn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating tcmn  ;'
   endif
!#endif


   ! macro expansion of dallocate = de- allocate;

!#ifdef DEBUG
   if( .not.allocated( tsmn  ) ) then
    write(0,'("macros   : 0123456789 : tsmn  not already allocated ;")') 
    stop      'macros   : 0123456789 : tsmn  not already allocated ;'
   endif
!#endif

   deallocate( tsmn , stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error de-allocating tsmn  ;")') 
    stop      'macros   : 0123456789 : error de-allocating tsmn  ;'
   endif
!#endif



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine action

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine actiongradient( Ndof, xx, ff, iuser, ruser, iflag )

  use globals, only : zero, one, two, half, pi2, pi, ounit, myid, tstart, &
                      pp, qq, iota, &
                      pqMpol, pqNtor, &
                      Bzeta
                                
  implicit none
  
  include "mpif.h"

  integer            :: Ndof, iuser(1:3), iflag
  real               :: xx(1:Ndof), ff(1:Ndof), ruser(1:1)

  integer, parameter :: MM = 8 

  integer            :: ipq, jpq, lpp, lqq, qN, Nfft, ii, nn, ic06faf
  real               :: rcn(0:iuser(3)), rsn(1:iuser(3)), tcn(0:iuser(3)), tsn(1:iuser(3)), nu, cnzq, snzq
  real               :: liota, targetarea, area, zz, dz, nzq, xy(1:6), Bxy(1:6)
  real               :: rho(0:MM*iuser(3)-1), the(0:MM*iuser(3)-1), rd(0:MM*iuser(3)-1), td(0:MM*iuser(3)-1), rwork(0:MM*iuser(3)-1)

  ipq = iuser(1) ; jpq = iuser(2)

  lpp = pp(ipq) ; lqq = qq(ipq) ; liota = iota(ipq) ; qN = lqq * pqNtor ; Nfft = MM * qN ; dz = pi2 / ( MM * pqNtor )

  targetarea = ruser(1)

  rcn(0:qN) = xx(1                :1+qN           ) - one
  rsn(1:qN) = xx(1+qN+1           :1+qN+qN        ) - one
  tcn(0:qN) = xx(1+qN+qN+1        :1+qN+qN+1+qN   ) - one
  tsn(1:qN) = xx(1+qN+qN+1+qN+1   :1+qN+qN+1+qN+qN) - one
  nu        = xx(1+qN+qN+1+qN+qN+1                ) - one
  
  do ii = 0, Nfft - 1
   
   zz = ii * dz

   rho(ii) = rcn(0)
   the(ii) = tcn(0) + liota * zz 
    
   do nn = 1, qN
    
    nzq = nn * zz / lqq ; cnzq = cos(nzq) ; snzq = sin(nzq)
    
    rho(ii) = rho(ii) + rcn(nn) * cnzq + rsn(nn) * snzq
    the(ii) = the(ii) + tcn(nn) * cnzq + tsn(nn) * snzq

   enddo

   xy(1:6) = (/ the(ii), rho(ii), zero, zero, zero, zero /)

!  call bfield( zz, xy(1:6), Bxy(1:6) )
   call cfield( zz, xy(1:6), Bxy(1:6) ) ! the calculation takes place in the pre-existing chaotic coordinates ; 12/17/20;

   td(ii) = Bxy(1)
   rd(ii) = Bxy(2) - nu / Bzeta

  enddo ! end of do ii ; 12/17/20;
  
  area = ( sum( the(0:Nfft-1) ) + lpp*pi ) * dz / (lqq*pi2) - lpp * pi

  call C06FAF( td(0:Nfft-1), Nfft, rwork(0:Nfft-1), ic06faf ) ; td = td / sqrt(one*Nfft) ; td(1:Nfft-1) = td(1:Nfft-1) * two
  call C06FAF( rd(0:Nfft-1), Nfft, rwork(0:Nfft-1), ic06faf ) ; rd = rd / sqrt(one*Nfft) ; rd(1:Nfft-1) = rd(1:Nfft-1) * two

  ;  nn = 0     ; ff(1                ) =                            - rd(      0)
  do nn = 1, qN ; ff(1+nn             ) = rsn(nn) * (   nn*one/lqq ) - rd(     nn)
   ;            ; ff(1+qN+nn          ) = rcn(nn) * ( - nn*one/lqq ) + rd(Nfft-nn)
  enddo
  ;  nn = 0     ; ff(1+qN+qN+1        ) = liota                      - td(     nn)
  do nn = 1, qN ; ff(1+qN+qN+1+nn     ) = tsn(nn) * (   nn*one/lqq ) - td(     nn)
   ;            ; ff(1+qN+qN+1+qN+nn  ) = tcn(nn) * ( - nn*one/lqq ) + td(Nfft-nn)
  enddo
  ;             ; ff(1+qN+qN+1+qN+qN+1) = area - targetarea

! if( something illegal ) iflag = -1

  return
  
end subroutine actiongradient

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

subroutine cfield( zeta, xys, Bxys )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  use globals, only : zero, one, tstart, ounit, &
                      Bzeta, &
                      Linterpol
  
  implicit none  
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  integer, parameter   :: Node = 6
  real   , intent(in)  :: zeta, xys(1:Node)
  real   , intent(out) :: Bxys(1:Node)
  
  integer              :: lLinterpol
  real                 :: xy(1:Node), Bxy(1:Node), Delta, AA(1:3,1:3), Brtp(1:3), stz(1:3), xyz(1:3,0:3)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  stz(1:3) = (/ xys(2), xys(1), zeta /)

! lLinterpol = Linterpol
  lLinterpol = 2         ! using C1 coordinate transformation ; 02/19/2021 ;

  call stzxyz( stz(1:3), xyz(1:3,0:3), lLinterpol ) ! coordinate transformation ;

  Delta = xyz(2,2) * xyz(1,1) - xyz(2,1) * xyz(1,2)
  
  AA(1,1) = + xyz(2,2) / Delta  ; AA(1,2) = - xyz(1,2) / Delta ; AA(1,3) = (   xyz(2,3) * xyz(1,2) - xyz(2,2) * xyz(1,3) ) / Delta ! vector transformation ; 
  AA(2,1) = - xyz(2,1) / Delta  ; AA(2,2) = + xyz(1,1) / Delta ; AA(2,3) = ( - xyz(2,3) * xyz(1,1) + xyz(2,1) * xyz(1,3) ) / Delta
  AA(3,1) =   zero              ; AA(3,2) =   zero             ; AA(3,3) =     one

  xy(1:Node) = (/ xyz(2,0), xyz(1,0), one, zero, zero, one /)

  call bfield( zeta, xy(1:Node), Bxy(1:Node) )

  Brtp(1:3) = matmul( AA(1:3,1:3), (/ Bxy(2), Bxy(1), Bzeta /) )
    
  Bxys(1:Node) = (/ Brtp(2), Brtp(1), zero, zero, zero, zero /)

  Bzeta = one

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
  
end subroutine cfield

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
