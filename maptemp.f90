














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine maptemp
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, pi2, ounit, tunit, tstart, &
                      chi, hmnfile, metrixsuff, ksuff, &
                      Ndof, Ni, Nj, Nk, dx, dy, dz, Nteta, Nzeta, &
                      imethod, &
                      PX, PY, lamda, mu, cc, e02defrwrk, e02defiwrk, &
                      mn, im, in, &
                      Linterpol
                                
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, parameter   :: MM = 1
  integer              :: ii, jj, kk, ie02def, astat
  real                 :: srho, teta, stz(1:3), xyz(1:3,0:3), lxx(1:MM), lyy(1:MM), lff(1:MM)
  real   , allocatable :: Mijk(:,:)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  select case( imethod )
   
  case( 0 )
   
  case( 1 )
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Mijk ) ) then
    write(0,'("macros    : 0123456789 : Mijk already allocated ;")') 
    stop      'macros    : 0123456789 : Mijk already allocated ;'
   endif
!#endif

   allocate( Mijk(0:Ni,0:Nj), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Mijk ;")') 
    stop      'macros   : 0123456789 : error allocating Mijk ;'
   endif
!#endif

   Mijk(0:Ni,0:Nj) = zero 


   
   do ii = 0, Ni ; teta = ii * dx
    
    do jj = 0, Nj ; srho = chi%low + jj * dy
     
     stz(1:3) = (/ srho, teta, zero /)
     
     call stzxyz( stz, xyz, Linterpol )
     
     lxx(1) = min(pi2,max(zero,xyz(2,0))) ; lyy(1) = min(chi%upp,max(chi%low,xyz(1,0)))
     
     ie02def = 0
     
     call E02DEF( MM, PX, PY, lxx(1:MM), lyy(1:MM), lamda(1:PX), mu(1:PY), cc(1:(PX-4)*(PY-4)), lff(1:MM), e02defrwrk(1:PY-4), e02defiwrk(1:PY-4), ie02def )
     
     Mijk(ii,jj) = lff(1)
     
    enddo
    
   enddo
   
   open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Me", status='replace', form='unformatted' )
   write(tunit) Ni, Nj
   write(tunit) Mijk(0:Ni,0:Nj)
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



  case( 2 ) 
   
    ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Mijk ) ) then
    write(0,'("macros    : 0123456789 : Mijk already allocated ;")') 
    stop      'macros    : 0123456789 : Mijk already allocated ;'
   endif
!#endif

   allocate( Mijk(0:Ni,0:Nk), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Mijk ;")') 
    stop      'macros   : 0123456789 : error allocating Mijk ;'
   endif
!#endif

   Mijk(0:Ni,0:Nk) = zero 


   
   do ii = 0, Ni ; teta = ii * dx
    
    do kk = 0, Nk ; srho = chi%low + kk * dz
     
     stz(1:3) = (/ srho, teta, zero /)
     
     call stzxyz( stz(1:3), xyz(1:3,0:3), Linterpol )
     
     lxx(1) = min(pi2,max(zero,xyz(2,0))) ; lyy(1) = min(chi%upp,max(chi%low,xyz(1,0)))
     
     ie02def = 0
     
     call E02DEF( MM, PX, PY, lxx(1:MM), lyy(1:MM), lamda(1:PX), mu(1:PY), cc(1:(PX-4)*(PY-4)), lff(1:MM), e02defrwrk(1:PY-4), e02defiwrk(1:PY-4), ie02def )
     
     Mijk(ii,kk) = lff(1)
     
    enddo
    
   enddo
   
   open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".M", status='replace', form='unformatted' )
   write(tunit) Ni, Nk
   write(tunit) Mijk(0:Ni,0:Nk)
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


   
  case( 3 ) ! imethod = 3; 20/10/21 ;
   
!   SALLOCATE( Mijk, (0:Nteta,0:Ni), zero )
!   
!   do ii = 0, Nteta ; teta = ii * pi2 / Nteta
!    
!    do jj = 0, Ni    ; srho = chi%low + ii * ( chi%upp - chi%low) / Ni
!     
!     stz(1:3) = (/ srho, teta, zero /)
!     
!     call stzxyz( stz(1:3), xyz(1:3,0:3), Linterpol )
!     
!     Mijk(ii,jj) = zero
!     do kk = 1, mn
!      Mijk(ii,jj) = Mijk(ii,jj) + fTT(kk) * cos( im(kk) * xyz(2,0) + in(kk) * xyz(3,0) )
!     enddo
!     
!    enddo ! end of do jj ; 20/10/21 ;
!    
!   enddo ! end of do ii ; 20/10/21 ;
!   
!   open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".M", status='replace', form='unformatted' )
!   write(tunit) Nteta, Ni
!   write(tunit) Mijk(0:Nteta,0:Ni)
!   close(tunit)
!   
!   DALLOCATE( Mijk )
   
  case default
   
  end select
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine maptemp

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !


