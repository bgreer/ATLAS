### NEED CFITSIO, FFTW, MKL

### LIBRARIES
FITSLOC=/shared/cfitsio
#FITSLOC=/home8/begr7169/SOFTWARE/cfitsio
FITSLIB=$(FITSLOC)/lib
FITSINC=$(FITSLOC)/include

FFTWLOC=/shared/fftw-3.2.1
FFTWLIB=$(FFTWLOC)/lib
FFTWINC=$(FFTWLOC)/include

MKL = -L/central/intel/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
FFTW = -I/shared/fftw/include -L/shared/fftw/lib -lfftw3
CFITSIO = -L/shared/cfitsio-3.24/lib -lcfitsio 

CLIBS = $(MKL) $(CFITSIO) $(FFTW) -lm
FLIBS = -I$(FFTWINC) -L$(FFTWLIB) -lfftw3 -L$(FITSLIB) -lcfitsio -I$(FITSINC) $(FLAGS) -I/home8/begr7169/SOFTWARE/openmpi-1.6.5/include -L/home8/begr7169/SOFTWARE/openmpi-1.6.5/lib

# FORTRAN COMPILE FLAGS
FFLAGS = -O2 -ip -ipo -g -CB -free -openmp
# C COMPILE FLAGS
CFLAGS = -O2 -ip -ipo -g
# CROSS-COMPILE FLAGS
ALLFLAGS = -O2 -ip -ipo -openmp

# FORTRAN COMPILER
#F90 = /home8/begr7169/SOFTWARE/openmpi-1.6.5/bin/mpif90
F90 = mpif90
# C COMPILER
ICC = icc


EXECUTABLE = atlas

FOBJECTS = ParseInput.o \
		  Communication_Library.o \
		  File_Management.o \
		  Interpolation.o \
		  Timing.o \
		  Pspec.o \
		  Track.o \
		  Grid.o \
		  Fit.o \
		  Main.o

COBJECTS = fit.o \
		  io.o \
		  function.o \
		  errbars.o \
		  ridge.o \
		  mpfit.o \
		  mrf_fit.o \
		  cart_to_polar.o \
		  parse_input.o \
		  fit_wrapper.o

# compile fortran and C parts together
$(EXECUTABLE) : $(FOBJECTS) $(COBJECTS)
	$(info ************ CAUTION: TRY COMPILING TWICE OR MORE BEFORE BELIEVING ANY ERROR MESSAGES ************)
	$(info ************ COMPILING ALL TOGETHER ************)
	$(F90) $(ALLFLAGS) -o $(EXECUTABLE) $(FOBJECTS) $(COBJECTS) $(FLIBS) $(CLIBS)


# compile C parts
$(COBJECTS) : *.c
	$(info ************ CAUTION: TRY COMPILING TWICE OR MORE BEFORE BELIEVING ANY ERROR MESSAGES ************)
	$(info ************ COMPILING C CODE ************)
	$(ICC) -c $(CFLAGS) *.c $(CLIBS)

# compile FORTRAN parts
$(FOBJECTS) : *.F90
	$(info ************ CAUTION: TRY COMPILING TWICE OR MORE BEFORE BELIEVING ANY ERROR MESSAGES ************)
	$(info ************ COMPILING FORTRAN CODE ************)
	$(F90) -c $(FFLAGS) *.F90 $(FLIBS)


clean : 
	rm -f $(COBJECTS) $(FOBJECTS) *.mod $(EXECUTABLE)

