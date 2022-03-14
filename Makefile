# Compiler
F90=mpiifort


# Compile flags
RFLAGS=-r8 -mcmodel=large -O3 -m64 -unroll0 -fno-alias -ip -traceback
DFLAGS=#-check all -check noarg_temp_created -debug full -D DEBUG

NAGFLAGS=-L${NAG_ROOT}/lib -lnag_nag
FFTWFLAGS=-I$(FFTW_HOME)/include -L$(FFTW_HOME)/lib -lfftw3
LFLAGS=$(NAGFLAGS) $(FFTWFLAGS)


# Objects
OBJS := main.o global.o \
	bfield.o stzxyz.o sumpolynomial.o metrics.o volume.o matrix.o \
	nearest.o biliner.o bicubic.o poincare.o viscous.o action.o \
	relax.o savedata.o readdata.o \
	readqfms.o preset.o maptemp.o tprofile.o \
	diffusion.o fdiffusion.o tsolve.o


# Programs
xheat : $(OBJS)
	$(F90) $(RFLAGS) $(DFLAGS) $(OBJS) $(LFLAGS) -o $@


dheat : $(OBJS)
	$(F90) $(RFLAGS) $(DFLAGS) $(OBJS) $(LFLAGS) -o $@



global.o: global.f90
	$(F90) $(RFLAGS) $(DFLAGS) -o global.o -c global.f90 $(LFLAGS)

%.o: %.f90 global.o 
	$(F90) $(RFLAGS) $(DFLAGS) -o $*.o -c $*.f90 $(LFLAGS)

main.o: main.f90 $(OBJS)
	$(F90) $(RFLAGS) $(DFLAGS) -o $@ -c main.f90 $(LFLAGS)

.PHONY: clean
clean:
	-rm -f *.o *.mod

