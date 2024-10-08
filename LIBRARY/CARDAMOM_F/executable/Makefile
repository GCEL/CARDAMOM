#### Settings..

compiler       = ifort    # ifort or gfortran
compiler_flags = normal   # normal, debug, or timing
linker_style   = dynamic  # dynamic or static

###################################################

#### Declaration of source files..

general_SOURCES = ../misc/math_functions.f90  \
                  ../misc/oksofar.f90  \
                  ../misc/brent_zero/brent_zero.f90  \
                  ../general/cardamom_structures.f90 \
                  ../method/MHMCMC/MCMC_FUN/MHMCMC_STRUCTURES.f90 \
                  ../model/DALEC_CDEA_FR/src/DALEC_CDEA_FR_PARS.f90 \
                  ../general/cardamom_io.f90   \
                  ../method/MHMCMC/MCMC_FUN/MHMCMC.f90  \
                  ../model/DALEC_CDEA_FR/likelihood/MODEL_LIKELIHOOD.f90       \
                  ../model/DALEC_CDEA_FR/src/DALEC_CDEA_FR.f90 \
                  ../general/cardamom_main.f90 

###################################################

#### Declaration of general compiler/linker flags..
# ( which one actually gets used depends upon what )
# (  the user specifies at the top of this file.   )

timing = -pg

# For the Intel fortran compiler..(v11 or newer!)
ifort_normal  = -fpp -vec-report=0 -O2
ifort_debug   = $(ifort_normal) -check all -debug-parameters all -ftrapuv -g -fpe0 \
                -implicitnone -p -stand f03 -traceback 
ifort_timing  = $(ifort_normal) $(timing) 
ifort_dynamic = -i-dynamic                   # use the dynamic intel libs
ifort_static  = -static -i-static -I/usr/lib # only use the static intel libraries

# For the GNU fortran compiler..(v4.4.5 or newer!)
gfortran_normal  = -cpp -ffree-line-length-none -O2
gfortran_debug   = $(gfortran_normal) -fbounds-check -fimplicit-none -frange-check \
                   -ftree-vectorizer-verbose=0 -ggdb -pedantic -std=f2003 -Wall
gfortran_timing  = $(gfortran_normal) $(timing)
gfortran_dynamic =                           #
gfortran_static  = -static # -static-libgcc  # only use the static gfortran libraries

#############################################
# PLEASE DO NOT TOUCH BELOW HERE UNLESS YOU #
# ARE HAPPY YOU KNOW HOW MAKEFILES WORK.    #
#############################################

#### Trim any white-space off the user's entries..

FC=$(strip $(compiler))
COMP=$(strip $(compiler_flags))
LINK=$(strip $(linker_style))

#### Select the compiler/linker flags the user asked for...

ifeq ($(FC),gfortran)
	NC_INC_FLAGS = -I$(ncpath)/include -I/usr/lib64/gfortran/modules/
else
	NC_INC_FLAGS = -I$(ncpath)/include
endif

FC_FLAGS   = $($(FC)_$(COMP)) $(NC_PREPROCESSOR_FLAG)
INC_FLAGS  = $(NC_INC_FLAGS)
LINK_FLAGS = $($(FC)_$(LINK)) $(NC_LINK_FLAGS)

#### Select the sources files the user asked for..

	cardamom_SOURCES = $(general_SOURCES)

#### Instructions for building CARDAMOM...

cardamom.exe:	$(cardamom_SOURCES:%.f90=%.o)
	$(FC) -o $@ $^ $($(COMP)) $(LINK_FLAGS)

clean:
	rm -f *.mod ../misc/*.o ../model/AT_DALEC/src/*.o ../model/AT_DALEC/likelihood/*.o ../general/*.o ../method/MHMCMC/MCMC_FUN/*.o make_dep cardamom.exe gmon.out

%.o:	%.f90
	$(FC) $(FC_FLAGS) $(INC_FLAGS) -c $< -o $@

Makefile:	make_dep

make_dep:	$(cardamom_SOURCES)
		./f90_dependencies.py $^ > $@

-include make_dep
