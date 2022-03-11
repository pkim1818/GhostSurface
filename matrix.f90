














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine matrix
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals, only : zero, ounit, munit, ten, four, twelve, myid, tstart, &
                            kpara, tkperp, Tlow, Tupp, &
                            Ni, Nj, Nk, dx, dy, dz, metrixsuff, ksuff, hmnfile, &
                            idiff, kdiff, opint, &
                            dxydabf, Bfh, &
                            Ndof, Mnz, nz, irow, jcol, Acs, rhs, matrix_exist, &
                            chi, nperp, npint, ifd, jfd, wfd, iin, jin, wpd
  
  implicit none
  
  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  real    :: cpui, cpuo, cput
  integer :: ierr, astat, ii, jj, kk, nn, mm, ia, ja, ka, ib, jb, kb, li, lj
  real    :: lx, ly
! INTEGER :: ij, nperp(0:2), lperp, ifd(0:2,1:9), jfd(0:2,1:9), npint(1:2), lpint, iin(1:2,1:16), jin(1:2,1:16), lpara
  integer :: ij,             lperp,                                         lpint,                               lpara
! REAL    :: wfd(0:2,1:9), win(1:16), wpd(1:2,-2:2,-2:2), fpara
  real    ::               win(1:16)
  real    :: fpara

! REAL    :: xx, yy, zz, SS ! include source; only if lsource.eq.1;

  integer :: ios, inz

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cpui = MPI_WTIME() ; cpuo = cpui
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! nperp(0:2) = (/ 5, 9, 9 /) ! #points used in perpendicular finite differencing for idiff = 0, 1, 2; ! globals; 12/17/20;

! npint(1:2) = (/ 4, 16 /) ! #points used in perpendicular interpolation for opint = 1, 2; ! globals ; 12/17/20;

  Mnz = ( nperp(idiff) + 2 * kdiff * npint(opint) ) * Ndof ! upper bound on #non-zero elements in matrix;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( rhs ) ) then
    write(0,'("macros    : 0123456789 : rhs already allocated ;")') 
    stop      'macros    : 0123456789 : rhs already allocated ;'
   endif
!#endif

   allocate( rhs(1:Ndof), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating rhs ;")') 
    stop      'macros   : 0123456789 : error allocating rhs ;'
   endif
!#endif

   rhs(1:Ndof) = zero 

 ! boundary condition; Ndof is set in global; Ndof counts degrees of freedom;

   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( Acs ) ) then
    write(0,'("macros    : 0123456789 : Acs already allocated ;")') 
    stop      'macros    : 0123456789 : Acs already allocated ;'
   endif
!#endif

   allocate( Acs(1:Mnz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating Acs ;")') 
    stop      'macros   : 0123456789 : error allocating Acs ;'
   endif
!#endif

   Acs(1:Mnz) = zero 



   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( irow ) ) then
    write(0,'("macros    : 0123456789 : irow already allocated ;")') 
    stop      'macros    : 0123456789 : irow already allocated ;'
   endif
!#endif

   allocate( irow(1:Mnz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating irow ;")') 
    stop      'macros   : 0123456789 : error allocating irow ;'
   endif
!#endif

   irow(1:Mnz) = 0 


   ! macro expansion of sallocate = set allocate;

!#ifdef DEBUG
   if( allocated( jcol ) ) then
    write(0,'("macros    : 0123456789 : jcol already allocated ;")') 
    stop      'macros    : 0123456789 : jcol already allocated ;'
   endif
!#endif

   allocate( jcol(1:Mnz), stat=astat )

!#ifdef DEBUG
   if( astat.ne.0 ) then
    write(0,'("macros   : 0123456789 : error allocating jcol ;")') 
    stop      'macros   : 0123456789 : error allocating jcol ;'
   endif
!#endif

   jcol(1:Mnz) = 0 



! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  inquire( file = "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", exist=matrix_exist )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if( matrix_exist ) then
   
   if( myid.eq.0 ) then
    
    open( munit, file = "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", status='old', form='unformatted' )
    
    read(munit) nz ! probably should put some error flags on the read; 12/17/20;
    read(munit) irow
    read(munit) jcol
    read(munit) Acs
    read(munit) rhs
    
    close(munit)
    
   endif
   
! should broadcast
   
   cput = MPI_WTIME()
   
   write(ounit,1101) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", cput-cpuo ; cpuo = cput

1101 format("matrix    : "f10.2"s : sparse coupling read from ",a," ; time =",f10.2," ;")
   
   return

  endif ! end of if( matrix_exist ) ; 12/17/20;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  cput = MPI_WTIME()

  write(ounit,'("matrix    : "f10.2"s : computing "a" ;")') cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix"

! call flush(ounit)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ifd(0,1:5) = (/   0, - 1,   0, + 1,   0 /) ! global ; 12/17/20;
! jfd(0,1:5) = (/ - 1,   0,   0,   0,   1 /) ! global ; 12/17/20;
  
! wfd(0,1:5) = (/   0,   1, - 2,   1,   0 /) / dx / dx & ! global ; 12/17/20;
!            + (/   1,   0, - 2,   0,   1 /) / dy / dy  
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ifd(1,1:9) = (/ - 1,   0,   1, - 1,   0,   1, - 1,   0,   1 /) ! global ; 12/17/20;

  Mnz = ( nperp(idiff) + 2 * kdiff * npint(opint) ) * Ndof ! upper bound on #non-zero elements in matrix;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! SALLOCATE( rhs, (1:Ndof), zero ) ! boundary condition; Ndof is set in global; Ndof counts degrees of freedom;

! SALLOCATE( Acs, (1:Mnz), zero )

! SALLOCATE( irow, (1:Mnz), 0 )
! SALLOCATE( jcol, (1:Mnz), 0 )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! inquire( file = "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", exist=matrix_exist )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
! if( matrix_exist ) then
   
!  if( myid.eq.0 ) then
    
!   open( munit, file = "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", status='old', form='unformatted' )
    
!   read(munit) nz ! probably should put some error flags on the read; 12/17/20;
!   read(munit) irow
!   read(munit) jcol
!   read(munit) Acs
!   read(munit) rhs
    
!   close(munit)
    
!  endif
   
! should broadcast
   
!  cput = MPI_WTIME()
   
!  write(ounit,1101) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", cput-cpuo ; cpuo = cput

!101 format("matrix    : "f10.2"s : sparse coupling read from ",a," ; time =",f10.2," ;")
   
!  return

! endif ! end of if( matrix_exist ) ; 12/17/20;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! cput = MPI_WTIME()

! write(ounit,'("matrix    : "f10.2"s : computing "a" ;")') cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix"

! call flush(ounit)

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ifd(0,1:5) = (/   0, - 1,   0, + 1,   0 /) ! global ; 12/17/20;
! jfd(0,1:5) = (/ - 1,   0,   0,   0,   1 /) ! global ; 12/17/20;
  
! wfd(0,1:5) = (/   0,   1, - 2,   1,   0 /) / dx / dx & ! global ; 12/17/20;
!            + (/   1,   0, - 2,   0,   1 /) / dy / dy  
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ifd(1,1:9) = (/ - 1,   0,   1, - 1,   0,   1, - 1,   0,   1 /) ! global ; 12/17/20;
! jfd(1,1:9) = (/ - 1, - 1, - 1,   0,   0,   0,   1,   1,   1 /) ! global ; 12/17/20;
  
! wfd(1,1:9) = (/   1, - 2,   1,   2, - 4,   2,   1, - 2,   1 /) / dx / dx / four & ! global ; 12/17/20;
!            + (/   1,   2,   1, - 2, - 4, - 2,   1,   2,   1 /) / dy / dy / four

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! ifd(2,1:9) = (/   0,   0, - 2, - 1,   0,   1,   2,   0,   0 /) ! global ; 12/17/20;
! jfd(2,1:9) = (/ - 2, - 1,   0,   0,   0,   0,   0,   1,   2 /) ! global ; 12/17/20;
  
! wfd(2,1:9) = (/   0,   0, - 1,  16, -30,  16, - 1,   0,   0 /) / dx / dx / twelve & ! global ; 12/17/20;
!            + (/ - 1,  16,   0,   0, -30,   0,   0,  16, - 1 /) / dy / dy / twelve ! global ; 12/17/20;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! iin(1,1: 4) = (/  0, 1, 0, 1 /) ! linear interpolation? ; 12/17/20; ! global ; 12/17/20;
! jin(1,1: 4) = (/  0, 0, 1, 1 /) ! global ; 12/17/20;

! iin(2,1:16) = (/ -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2 /) ! cubic interpolation? ; 12/17/20; ! global ; 12/17/20;
! jin(2,1:16) = (/ -1,-1,-1,-1,  0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2 /) ! global ; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! wpd(1,-2,-2:2) = (/   0,   0,   0,   0,   0 /) * kpara / dz**2 ! global ; 12/17/20;
! wpd(1,-1,-2:2) = (/   0,   1,   0,   0,   0 /) * kpara / dz**2 ! global ; 12/17/20;
! wpd(1, 0,-2:2) = (/   0, - 1,   0, - 1,   0 /) * kpara / dz**2 ! global ; 12/17/20;
! wpd(1, 1,-2:2) = (/   0,   0,   0,   1,   0 /) * kpara / dz**2 ! global ; 12/17/20;
! wpd(1, 2,-2:2) = (/   0,   0,   0,   0,   0 /) * kpara / dz**2 ! global ; 12/17/20;
  
! wpd(2,-2,-2:2) = (/ - 1,   0,   0,   0,   0 /) * kpara / dz**2 / 12 ! global ; 12/17/20;
! wpd(2,-1,-2:2) = (/   0,  16,   0,   0,   0 /) * kpara / dz**2 / 12 ! global ; 12/17/20;
! wpd(2, 0,-2:2) = (/   1, -16,   0, -16,   1 /) * kpara / dz**2 / 12 ! global ; 12/17/20;
! wpd(2, 1,-2:2) = (/   0,   0,   0,  16,   0 /) * kpara / dz**2 / 12 ! global ; 12/17/20;
! wpd(2, 2,-2:2) = (/   0,   0,   0,   0, - 1 /) * kpara / dz**2 / 12 ! global ; 12/17/20;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! now build the sparse matrix

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  nz = 0 ! initialize counter ; 12/17/20;

  do kk = 0, Nk-1 ! loop over toroidal planes;
   do jj = 1, Nj-1 ! loop over grid points; radial ; 12/17/20;
    do ii = 0, Ni-1 ! loop over grid points; poloidal ; 12/17/20;
     
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
     
     nn = 1 + ii + (jj-1)*Ni + kk*(Nj-1)*Ni

     if(nn.lt.1 .or. nn.gt.Ndof ) then
          write(ounit,'("matrix :    fatal    : nn.lt.1 .or. nn.gt.Ndof  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

     lpara = kdiff

     if( jj.le.1 .or. jj.ge.Nj-1 ) then ; lperp =   1   ! must use first order differencing near boundary ; 12/17/20;
     else                               ; lperp = idiff
     endif
     
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

     do ij = 1, nperp(lperp) ; jb = jj + jfd(lperp,ij) ; ib = ii + ifd(lperp,ij) ! construct perpendicular diffusion;
      
      if( ib.lt.0  ) ib = ib + Ni ! periodicity in x;
      if( ib.ge.Ni ) ib = ib - Ni ! periodicity in x;
      
      if    ( jb.le. 0 ) then ;             ; rhs(nn) = tkperp * wfd(lperp,ij) * Tlow + rhs(nn)                      ! lower boundary;
      elseif( jb.ge.Nj ) then ;             ; rhs(nn) = tkperp * wfd(lperp,ij) * Tupp + rhs(nn)                      ! upper boundary;
      else                    ; nz = nz + 1 ; Acs(nz) = tkperp * wfd(lperp,ij)                                       ! degree of freedom;
       ;                      ;             ; mm = 1 + ib + (jb-1)*Ni + kk*(Nj-1)*Ni ; irow(nz) = nn ; jcol(nz) = mm ! element identification;
       if(mm.lt.1 .or. mm.gt.Ndof ) then
          write(ounit,'("matrix :    fatal    : mm.lt.1 .or. mm.gt.Ndof  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
      endif
      
      if( ib.eq.ii .and. jb.eq.jj ) then

       fpara = sum( wpd(lpara,0,-lpara:lpara) * Bfh(-lpara:lpara,ii,jj,kk) )
       Acs(nz) = fpara + Acs(nz) ! parallel;

      endif
      
     enddo ! end of ij loop;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
     
     do ka = -lpara, lpara ; kb = kk + ka ! loop over backward/forward mapped points;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

      if( ka.eq.0 ) cycle ! parallel self-coupling included above;
      
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      
      do
       if( kb.lt.  0 ) kb = kb + Nk ! enforce periodicity in z;
       if( kb.ge. Nk ) kb = kb - Nk ! enforce periodicity in z;
       if( kb.ge.0 .and. kb.lt.Nk ) exit
      enddo
      
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

      call nearest( dxydabf(1,ka*2,ii,jj,kk), dxydabf(2,ka*2,ii,jj,kk), li, lj, lx, ly )

      if( lj.le.1 .or. lj.ge.Nj-1 ) then ; lpint = 1     ! must use first-order differencing near boundary ; 12/17/20;
      else                               ; lpint = opint
      endif

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
      
      if    ( lpint.eq.1 ) then ; call biliner( li, lj, lx, ly, win(1: 4) ) ! use bi-linear weights;
      elseif( lpint.eq.2 ) then ; call bicubic( li, lj, lx, ly, win(1:16) ) ! use bi-cubic weights;
      endif
      
      do ij = 1, npint(lpint) ; jb = lj + jin(lpint,ij) ; ib = li + iin(lpint,ij)
       
       if( ib.lt.0  ) ib = ib + Ni ! periodicity in x;
       if( ib.ge.Ni ) ib = ib - Ni ! periodicity in x;
       
       fpara = sum( wpd(lpara,ka,-lpara:lpara) * Bfh(-lpara:lpara,ii,jj,kk) )

       if    ( jb.le. 0 ) then ;             ; rhs(nn) = fpara * win(ij) * Tlow + rhs(nn)
       elseif( jb.ge.Nj ) then ;             ; rhs(nn) = fpara * win(ij) * Tupp + rhs(nn)
       else                    ; nz = nz + 1 ; Acs(nz) = fpara * win(ij) ; mm = 1 + ib + (jb-1)*Ni + kb*(Nj-1)*Ni ; irow(nz) = nn ; jcol(nz) = mm
       if(mm.lt.1 .or. mm.gt.Ndof ) then
          write(ounit,'("matrix :    fatal    : mm.lt.1 .or. mm.gt.Ndof  ;")')
          call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
          stop
         endif
       endif
       
      enddo ! end of ij loop;
       
     enddo ! end of ka loop;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
     
    enddo ! end of ii loop;
   enddo ! end of jj loop;
  enddo ! end of kk loop;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  rhs = - rhs
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  cput = MPI_WTIME()
  
  write( ounit, 1001 ) cput-tstart, idiff, kdiff, opint, cput-cpuo

1001 format("matrix    : "f10.2"s : constructed sparse coupling ; idiff ="i2" ; kdiff ="i2" ; opint ="i2" ; time ="f10.2" ;")

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return
  
end subroutine matrix

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

