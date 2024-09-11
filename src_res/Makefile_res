################################################
#	COMPILER
################

#GCC
H5FC = h5pfc
FC = gfortran
H5PFC = h5pfc
MPIFC = mpif90

#inputs para hacer pruebas
HDF=1
OPENMP=1
MPI=1
DBG=1
MN=0

# Preprocessing options/regions

ifeq ($(HDF),1)
	UHDF = HDF
        ifeq ($(MPI),1)
                FC = $(H5PFC)
		UPARALELO=PARALELO
	else
                FC = $(H5FC)
		UPARALELO=NOPARALELO	
        endif
else
	UHDF = NOHDF
        ifeq ($(MPI),1)
                FC = $(MPIFC)
		UPARALELO=PARALELO
        else
		ifeq ($(IBM),1)
			FC = $(IFC)
		endif
                UPARALELO=NOPARALELO
        endif
endif



ifeq ($(OPENMP),1)
        UOPENMP = OPENMP
else
        UOPENMP = NOOPENMP
endif

	ifeq ($(OPENMP),1)
           OMPF = -fopenmp
	endif

	ifeq ($(DBG),1)
            FCOPT = -ffpe-trap=invalid,zero,overflow -ffree-line-length-170 -m64 -march=native -fimplicit-none -O3 -fcheck=all -g -fimplicit-none  
            #-Wall #-Wextra -Wuninitialized
        else
            FCOPT = -ffree-line-length-170 -m64 -march=native -fimplicit-none -O3 #-xHost #-r8
        endif

# Defining preprocessing and flags 

        DEFINES =  -D$(UOPENMP) -D$(UHDF) -D$(UPARALELO) # -DNONSTCPP
        FLAGS = -ffree-form #-x f95-cpp-input
        FFLAGS = $(FCOPT) $(OMPF) $(FLAGS) $(DEFINES)

# DIRECTORIES
SRCDIR     := src
BINDIR     := bin
VPATH := $(SRCDIR) $(dir $(wildcard $(SRCDIR)/*/.))

# SOURCES
MODULE_SRC 	:=  modulos.F main3bu.F bndrybu.F coef.F eos.F filnam.F geom.F        \
      getprfq3d.F gridbu.F hydrow3.F initbu.F init_ext.F input.F pltout.F psgrav.F \
      riem3dp.F tdelay.F tstep_nocooling.F paralelo.F stmass0.F sweeps_nocooling.F misc.F \
      restrt_hdf5.F rec_hdf5.F rst_hdf5.F  

MAIN_SRC 	:= main3bu.F

# EXECUTABLES
MAIN_EXE := RATPENAT
FFLAGS+= -J $(BINDIR)

# CREATE OBJECTS
MODULE_SRC := $(foreach file,$(MODULE_SRC),$(SRCDIR)/$(file))
MAIN_SRC   := $(foreach file,$(MAIN_SRC),$(SRCDIR)/$(file))

MODULE_OBJ := $(patsubst %,$(BINDIR)/%,$(notdir $(MODULE_SRC:.F=.o)))
MAIN_OBJ   := $(patsubst %,$(BINDIR)/%,$(notdir $(MAIN_SRC:.F=.o)))

ifneq ($(BINDIR),)
$(shell test -d $(BINDIR) || mkdir -p $(BINDIR))
endif

# PROGRAM
default: $(MAIN_EXE)

$(MAIN_EXE): $(MODULE_OBJ) $(MAIN_OBJ)
	$(FC) $(FFLAGS) $^ -o $@

$(BINDIR)/%.o : %.F
	$(FC) $(FFLAGS) -c -o $@ $<

# CLEAN
clean:
	@rm -rf $(BINDIR)/*.o $(BINDIR)/*.mod
	@rm -rf $(MAIN_EXE)
	@echo "bin executables deleted"
