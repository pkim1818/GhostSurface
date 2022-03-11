#!/bin/sh

###############################################################################################################################################################

#files=bfield        sumpolynomial metrics matrix nearest biliner bicubic poincare pqfms coords relax savedata readdata freemem cputime
#afiles=bfield stzxyz cfield sumpolynomial metrics volume matrix nearest biliner bicubic poincare pqfms viscous action coords relax savedata readdata
 afiles=bfield stzxyz        sumpolynomial metrics volume matrix nearest biliner bicubic poincare       viscous action        relax savedata readdata
 files=$(afiles) readqfms preset maptemp tprofile diffusion fdiffusion tsolve

 allfiles=main global $(files)

 extras=macros
 
###############################################################################################################################################################

 MACROS=macros
 F90=mpif90

 RFLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback
 DFLAGS=#-check all -check noarg_temp_created -debug full -D DEBUG

 NAGFLAGS=-L${NAG_ROOT}/lib -lnag_nag
#FFTWFLAGS=-I$(FFTWHOME)/include -L$(FFTWHOME)/lib -lfftw3
 FFTWFLAGS=-I$(FFTW_HOME)/include -L$(FFTW_HOME)/lib -lfftw3

# NETCDFFLAGS=-I$(NETCDF_C_HOME)/include -L$(NETCDF_C_HOME)/lib -lnetcdf
# NETCDF_FORT_FLAGS=-I$(NETCDF_FORTRAN_HOME)/include -L$(NETCDF_FORTRAN_HOME)/lib -lnetcdff
# EZSPLINEFLAGS=-I$(PSPLINE_HOME)/include -L$(PSPLINE_HOME)/lib -lpspline -lpthread -lz -lm
# EZCDFFLAGS=                      -L$(NTCC_HOME)/lib -lezcdf

 LFLAGS=$(NAGFLAGS) $(FFTWFLAGS) $(EZCDFFLAGS) $(NETCDF_FORT_FLAGS) $(NETCDFFLAGS) $(EZSPLINEFLAGS)

###############################################################################################################################################################

 date:=$(shell date)
 suff:=$(shell date | awk '{print $$6$$2$$3}')

###############################################################################################################################################################

xheat: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(F90) $(RFLAGS) $(DFLAGS) -o xheat $(addsuffix .o,$(allfiles)) $(LFLAGS)

dheat: $(addsuffix .o,$(allfiles)) $(MACROS) Makefile
	$(F90) $(RFLAGS) $(DFLAGS) -o dheat $(addsuffix .o,$(allfiles)) $(LFLAGS)

###############################################################################################################################################################

global.o: global.h $(MACROS) Makefile
	m4 -P $(MACROS) global.h > global.f90
	$(F90) $(RFLAGS) $(DFLAGS) -o global.o -c global.f90 $(LFLAGS)

###############################################################################################################################################################

%.o: %.h global.o $(MACROS) Makefile
	m4 -P $(MACROS) $*.h > $*.f90
	$(F90) $(RFLAGS) $(DFLAGS) -o $*.o -c $*.f90 $(LFLAGS)

###############################################################################################################################################################

main.o: main.h global.o $(addsuffix .o,$(files)) $(MACROS) Makefile
	@awk -v date='$(date)' -v pwd='$(PWD)' -v macros='$(MACROS)' -v f90='$(F90)' -v FLAGS='$(RFLAGS) $(DFLAGS)' -v nag='$(NAGFLAGS)' \
	'{if($$2=="compilation") {print "  if( myid.eq.0 ) then " ; \
	                          print "    write(ounit,*)\"         : "$$2" : date    = "date" ; \"" ; \
	                          print "    write(ounit,*)\"         : "$$2" : dir     = "pwd" ; \"" ; \
	                          print "    write(ounit,*)\"         : "$$2" : macros  = "macros" ; \"" ; \
	                          print "    write(ounit,*)\"         : "$$2" : f90     = "f90" ; \"" ; \
	                          print "    write(ounit,*)\"         : "$$2" : flags   = "flags" ; \"" ; \
	                          print "    write(ounit,*)\"         : "$$2" : nag     = "nag" ; \"" ; \
	                          print "  endif " } \
	  else {print}}' main.h > mmain.h
	m4 -P $(MACROS) mmain.h > main.f90
	$(F90) $(RFLAGS) $(DFLAGS) -o main.o -c main.f90 $(LFLAGS)
	@rm -f mmain.h

###############################################################################################################################################################

clean:
	rm -f *.o ; rm -f *.mod ; rm -f *.f90 ; rm -f *.pdf ; rm -f *.dvi

###############################################################################################################################################################

