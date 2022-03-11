














! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

subroutine relax
  
  use globals, only : ounit, munit, myid, tstart, &
                              Ni, Nj, Nk, Ndof, Tijk, &
                              metrixsuff, ksuff, hmnfile, &
                              irow, jcol, Acs, nz, rhs, matrix_exist, &
                              lintol, Mits, &
                              lbasis, method, precon, omega, &
                              chi
  implicit none

  include "mpif.h"
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  integer              :: astat, ifail, if11def, Mbasis, extra, lwork, its, ii, jj, kk, nn, ios, inz
  integer, allocatable :: iwork(:), istr(:)
  real                 :: cpui, cput, cpuo, erri, errf, rnorm
  real   , allocatable :: work(:)  
  character            :: cero, dup, trans, check
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if(allocated(iwork)) stop 'iwork allocated'
	allocate(iwork(2*Ndof+1) ,stat=astat)
	if(astat.ne.0) stop 'error allocating iwork'

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cpuo = MPI_WTIME()
  
  if( .not.matrix_exist ) then ! the re-ordered representation has been saved;
   
   if(allocated(istr)) stop 'istr allocated'
	allocate(istr(Ndof+1) ,stat=astat)
	if(astat.ne.0) stop 'error allocating istr'
   
   dup = 'S' ; cero = 'R'
   ifail = 0 ; call F11ZAF( Ndof, nz, Acs(1:nz), irow(1:nz), jcol(1:nz), dup, cero, istr(1:Ndof+1), iwork(1:Ndof), ifail )
   
   cput = MPI_WTIME()
   
   write(ounit,1001) cput-tstart, ifail, cput-cpuo
   
1001 format("relax     : "f10.2"s : called F11ZAF (re-ordering) ; ifail ="i3" ; time ="f10.2" ;")
   
   if(.not.allocated(istr)) stop 'istr not allocated'
	deallocate(istr,stat=astat)
	if(astat.ne.0) stop 'error deallocating istr'
   
   open( munit, file = "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", status='unknown', form='unformatted' )   
   write(munit) nz
   write(munit) irow
   write(munit) jcol
   write(munit) Acs  
   write(munit) rhs  
   close(munit)
   
   cput = MPI_WTIME()
   
   write( ounit, 1011 ) cput-tstart, "."//trim(hmnfile)//"."//metrixsuff//"."//ksuff//".matrix", cput-cpuo
   
1011 format("relax     : "f10.2"s : "a" ; written ; time ="f10.2" ;") 
   
  endif

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  Mbasis = 0

  if( method.eq.'RGMRES'   ) Mbasis = min( Ndof, 50, lbasis )
  if( method.eq.'BICGSTAB' ) Mbasis = min( Ndof, 10, lbasis )
  
  extra = 0
  
  if( precon.ne.'J' .and. precon.ne.'S' ) precon = 'N'
  
  if( precon.eq.'J' .or. precon.eq.'S') extra = Ndof
  
  if( method.eq.'RGMRES' )   lwork =  4 * Ndof + Mbasis * ( Mbasis + Ndof + 5 )           + extra + 101
  if( method.eq.'CGS' )      lwork =  8 * Ndof                                            + extra + 100
  if( method.eq.'BICGSTAB' ) lwork =  2 * Ndof * ( Mbasis + 3 ) + Mbasis * ( Mbasis + 2 ) + extra + 100
  if( method.eq.'TFQMR'    ) lwork = 11 * Ndof                                            + extra + 100

  if(allocated(work)) stop 'work allocated'
	allocate(work(lwork) ,stat=astat)
	if(astat.ne.0) stop 'error allocating work'

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  trans = 'N' ; check = 'C'

  ifail = 0

  call F11XAF( trans, Ndof, nz, Acs(1:nz), irow(1:nz), jcol(1:nz), check, Tijk(1:Ndof), work(1:Ndof), ifail )

  erri = sqrt( sum( (rhs(1:Ndof)-work(1:Ndof))**2 ) / Ndof ) ! error initial;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  cput = MPI_WTIME()

  write(ounit,1012) cput-tstart, method, precon, Mbasis, Ndof, Mits

1012 format("relax     : "f10.2"s : calling F11DEF ; method = "a" ; precon = "a1" ; Mbasis ="i3" ; Ndof = "i9" ; Mits ="i8" ;")

  if11def = 1

  call F11DEF( method, precon, Ndof, nz, Acs(1:nz), irow(1:nz), jcol(1:nz), omega, rhs(1:Ndof), &
               Mbasis, lintol, Mits, Tijk(1:Ndof), rnorm, its, work(1:lwork), lwork, iwork(1:2*Ndof+1), if11def )

  cput = MPI_WTIME()
  
  select case( if11def )
  case( 0 )    ; write(ounit,'("relax     : "f10.2"s : F11DEF success ;")') cput-tstart
  case( 1 )    ; write(ounit,'("relax     : "f10.2"s : F11DEF input error ;")') cput-tstart
  case( 5 )    ; write(ounit,'("relax     : "f10.2"s : F11DEF required accuracy not obtained in allowed iterations ;")') cput-tstart
  case default ; write(ounit,'("relax     : "f10.2"s : F11DEF ifail not recognized ; if11def =",i," ;")') cput-tstart, if11def
  end select

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  trans = 'N' ; check = 'C'

  ifail = 0

  call F11XAF( trans, Ndof, nz, Acs(1:nz), irow(1:nz), jcol(1:nz), check, Tijk(1:Ndof), work(1:Ndof), ifail )

  errf = sqrt( sum( (rhs(1:Ndof)-work(1:Ndof))**2 ) / Ndof ) ! error final;

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  if(.not.allocated(work )) stop 'work  not allocated'
	deallocate(work ,stat=astat)
	if(astat.ne.0) stop 'error deallocating work '
  if(.not.allocated(iwork )) stop 'iwork  not allocated'
	deallocate(iwork ,stat=astat)
	if(astat.ne.0) stop 'error deallocating iwork '

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  cput = MPI_WTIME()

  write(ounit,1000) cput-tstart, method, precon, Mbasis, omega, lintol, if11def, erri, errf, rnorm, its, cput-cpuo

  cpuo = cput
  
1000 format("relax     : "f10.2"s : "a8" ; precon = "a1" ; Mbasis ="i3" ; omega ="f5.2" ; lintol ="es9.2" ; ":&
   "if11def ="i2" ; err ="es9.2" -->"es9.2" ; rnorm ="es9.2" ; its ="i12" ; time ="f10.2"s ;")
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  call savedata( Ni, Nj, Nk, Tijk, errf, rnorm )

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
  if(.not.allocated(rhs  )) stop 'rhs   not allocated'
	deallocate(rhs  ,stat=astat)
	if(astat.ne.0) stop 'error deallocating rhs  ' ! allocated in matrix;
  if(.not.allocated(irow )) stop 'irow  not allocated'
	deallocate(irow ,stat=astat)
	if(astat.ne.0) stop 'error deallocating irow ' ! allocated in matrix;
  if(.not.allocated(jcol )) stop 'jcol  not allocated'
	deallocate(jcol ,stat=astat)
	if(astat.ne.0) stop 'error deallocating jcol ' ! allocated in matrix;
  if(.not.allocated(Acs  )) stop 'Acs   not allocated'
	deallocate(Acs  ,stat=astat)
	if(astat.ne.0) stop 'error deallocating Acs  ' ! allocated in matrix;
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  return
  
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  
end subroutine relax

! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
