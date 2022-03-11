














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine readdata( iok )
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  use globals,only : ounit, tunit, tstart, myid, infile, hmnfile, &
                             kperp, iupdate, lintol, &
                             iNi, iNj, iNk, Ni, Nj, Nk, Tijk, Tlow, Tupp, metrixsuff, ksuff
  
  implicit none

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  integer, intent(out) :: iok
  
  LOGICAL              :: exist
  integer              :: oiNi, oiNj, oiNk, oNi, oNj, oNk, iostat, ii, jj, kk, nn, okperp
  real                 :: cpuo, cput, oerr, ornorm, oTlow, oTupp, oyylow, oyyupp
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 

  cput = MPI_WTIME() ; cpuo = cput

  iok = 0

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  inquire( file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Te", exist=exist )
  
  if( .not.exist ) then ; iok = 2 ; goto 1000
  endif
   
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

! cput = MPI_WTIME()

! write(ounit,'("readdata  : "f10.2"s : "a" ; reading ;")') cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Te"

  open( tunit, file="."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Te", status='old', form='unformatted', iostat=iostat )

  if( iostat.ne.0 ) then

   cput = MPI_WTIME()

   write( ounit, '("readdata  : "f10.2"s : error opening a file that exists ? ;")') cput-tstart

   iok = 2

   close(tunit)

   goto 1000

  endif

  read( tunit, iostat=iostat ) oiNi, oiNj, oiNk, okperp, oerr, ornorm, oTlow, oTupp, oyylow, oyyupp

  if( iostat.ne.0 ) then

   cput = MPI_WTIME()

   write( ounit,'("readdata  : "f10.2"s : error reading oiNi, oiNj, oiNk, . . . ;")') cput-tstart

   iok = 2

   close(tunit)

   goto 1000

  endif
  
  if( oiNi.ne.iNi .or. oiNj.ne.iNj .or. oiNk.ne.iNk ) then

   cput = MPI_WTIME()

   write( ounit, 2001 ) cput-tstart, myid

   iok = 2

   close(tunit)

   goto 1000

  endif
  
  if( okperp  .ne. kperp  ) iok = 1 ! the .Te file has probably been copied to provide a good initial guess;
  if( oerr    .gt. lintol ) iok = 1 ! unacceptable convergence;
 !if( ornorm  .gt. lintol ) iok = 1 ! unacceptable convergence;
  if( iupdate .eq.   1    ) iok = 1 ! force additional iterations;
  if( iupdate .eq.   2    ) iok = 0 ! unconditionally accept previous solution;

  read( tunit, iostat=iostat ) Tijk ! ( 0:Ni-1, 1:Nj-1, 0:Nk-1 ) ! read in previous data;

  if( iostat.ne.0 ) then

   cput = MPI_WTIME()

   write(ounit,'("readdata  : "f10.2"s : error reading Tijk ;")') cput-tstart

   iok = 2

  endif

  close(tunit)
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

1000 continue

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  if( iok.eq.0 .or. iok.eq.1 ) then ! old data correctly read from file;
   
   cput = MPI_WTIME()

   write( ounit, 1001 ) cput-tstart, okperp, oerr, ornorm, lintol, iupdate, iok, cput-cpuo

   cpuo = cput

1001 format("readdata  : ",f10.2,"s : reading Tijk ; old kperp = ",i3.2," ; oerr =",es9.2, &
  " ; ornorm =",es9.2," ; lintol =",es9.2," ; iupdate =",i2" ; iok =",i2," ; time =",f10.2," ;")
   
  else ! old data could not be read from file; initialize with linear profile;

   do kk = 0, Nk-1
    do ii = 0, Ni-1
     do jj = 1, Nj-1 ; Tijk( 1 + ii + (jj-1)*Ni + kk*Ni*(Nj-1) ) = Tlow + jj * (Tupp-Tlow) / Nj
     enddo
    enddo
   enddo
   
   cput = MPI_WTIME()

   write(ounit,1002) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".Te", iok, cput-cpuo ; cpuo = cput
   
1002 format("readdata  : "f10.2"s : "a" ; not exist (or corrupted) ; initialize linear ; iok ="i2" ; time ="f10.2" ;")
   
  endif

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  return

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

2001  format("readdata  : "f10.2"s : resolution error ;")
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine readdata

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

