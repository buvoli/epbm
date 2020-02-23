# BSD Liscence
# ============================================================================================================
# Makefile for Fortran Implementation of Numerical Code in Buvoli T., "ETD Spectral Deferred Correction", 2014
# Complied using gfortran 4.9.0 and gnu make 3.81
#
# External Dependencies:
# 	BLAS 	(http://www.netlib.org/blas/)
# 	LAPACK 	(http://www.netlib.org/lapack/)
# 	FFTW 	(http://www.fftw.org)
# ============================================================================================================

# Optional Variables
ifndef $(DEBUG)
DEBUG=FALSE
endif

ifndef $(OMP)
OMP=FALSE
endif

# Include Appropriate Directories
VPATH = tools:equations

# Compiler
FC = gfortran

# Flags
LDFLAGS = -L /usr/local/lib/ -lfftw3 -llapack -lblas
#FCFLAGS = -o3 -fexternal-blas

ifeq ($(DEBUG), TRUE)
FCFLAGS += -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
endif

ifeq ($(OMP), TRUE)
FCFLAGS += -fopenmp
endif

# Objects
EQNS = quasigeostrophic_mod.o nikolaevskiy_mod.o kuramoto_mod.o kdv_mod.o nls_mod.o
METHODS = etdrk4_mod.o etdsdc_mod.o imexsdc_mod.o etdab_mod.o etdpbm_mod.o etdpbm_F1_mod.o etdpbm_omp_mod.o etdpbm_F1_omp_mod.o
PHI = phi_mod.o
TOOLS = tools_mod.o

PROGRAMS = sisc-experiments.exe

all: $(PROGRAMS)

$(EQNS): $(TOOLS)
$(METHODS): $(TOOLS) $(PHI)

fft_omp.o: $(TOOLS)
fft_omp: $(TOOLS)

matprod_omp.o: $(TOOLS)
matprod_omp: $(TOOLS)

singlerun_omp.o: $(TOOLS) $(METHODS) $(EQNS)
singlerun_omp.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

singlerun.o: $(TOOLS) $(METHODS) $(EQNS)
singlerun.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

experimentrun.o: $(TOOLS) $(METHODS) $(EQNS)
experimentrun.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

experimentrun_mp.o: $(TOOLS) $(METHODS) $(EQNS)
experimentrun_mp.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

sisc-experiments.o: $(TOOLS) $(METHODS) $(EQNS)
sisc-experiments.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

sisc-experiments-epbm-a.o: $(TOOLS) $(METHODS) $(EQNS)
sisc-experiments-epbm-a.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

sisc-experiments-epbm-m.o: $(TOOLS) $(METHODS) $(EQNS)
sisc-experiments-epbm-m.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

sisc-experiments-epbm-eab.o: $(TOOLS) $(METHODS) $(EQNS)
sisc-experiments-epbm-eab.exe: $(TOOLS) $(PHI) $(METHODS) $(EQNS)

# ======================================================================
# General rules. Taken From
# (http://www.webalice.it/o.drofa/davide/makefile-fortran/makefile-fortran.html)
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%.exe: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD
