# Compiler
FC = ifort

# Optimization and debug flags
FFLAGS = -O3
# for gfortran
#FFLAGS = -O3 -ffree-line-length-none

# Uncomment for debugging:
# FFLAGS = -O0 -g -traceback -fbounds-check -DDEBUG -pg

# Library paths (EDIT THESE!)
#LIBPATH = -L/path/to/blas \
          -L/path/to/lapack \
          -L/path/to/fftw3

# Libraries (EDIT IF NEEDED!)
LIBS = -lblas -llapack -lfftw3

BIN_DIR = ../bin
UTL_DIR = ../utility

# Source files
SRC = global_module.f90 \
      modelvar_module.f90 \
      parse_module.f90 \
      sysdef_module.f90 \
      initial_module.f90 \
      rpmd_module.f90 \
      models_module.f90 \
      print_module.f90 \
      dboc_module.f90 \
      propagation_module.f90 \
      nonadiabatic_module.f90 \
      runtraj_module.f90 

# Main program
MAIN = rpmdsh.f90
# Executable
XOBJ = $(BIN_DIR)/sharp.x
# Object files
OBJS = ${SRC:.f90=.o}

#AverageCode
CAVG1 = $(UTL_DIR)/average.f90
XAVG1 = $(BIN_DIR)/average.x

CAVG2 = $(UTL_DIR)/average2db.f90
XAVG2 = $(BIN_DIR)/average2db.x

all: $(XOBJ) $(XAVG1) $(XAVG2)
	@echo "==> Build complete. Executables are in $(BIN_DIR)"

# Compilation rules
%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

%.mod : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# Linking rule
$(XOBJ): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(MAIN) -o $(XOBJ) $(LIBPATH) $(LIBS)

$(XAVG1): $(CAVG1)
	$(FC) $(FFLAGS) $< -o $@

$(XAVG2): $(CAVG2)
	$(FC) $(FFLAGS) $< -o $@

# Cleaning
clean:
	@echo "==> Cleaning ALL intermediate and executable files"
	rm -f *.mod *.o $(XOBJ) $(XAVG1) $(XAVG2)

