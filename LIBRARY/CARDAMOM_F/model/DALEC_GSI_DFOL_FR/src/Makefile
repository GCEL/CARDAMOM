#### Settings..

compiler       = ifort    # ifort or gfortran
compiler_flags = normal   # normal, debug, or timing
linker_style   = dynamic  # dynamic or static

###################################################

#### Declaration of source files..

general_SOURCES = DALEC_GSI_DFOL_FR.f90 \
	          DALEC_GSI_DFOL_FR_CROP.f90 \
                  DALEC_CROP_GLOBAL_VARIABLES.f90 \
		  DALEC_CROP_IO.f90 \

general_SOURCES_path = ../src


MAIN_SOURCE = MAIN_DALEC_CROP.f90


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
#ifort_dynamic = -i-dynamic                   # use the dynamic intel libs
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

SOURCES = $(general_SOURCES) $(MAIN_SOURCE) 

#### Instructions for building DALEC_crop...

DALEC_CROP.exe:	$(SOURCES:%.f90=%.o)
		$(FC) -o $@ $^ $($(COMP)) $(LINK_FLAGS)
		mv DALEC_CROP.exe ../

clean:
	rm -f *.mod *.o make_dep DALEC_CROP.exe gmon.out

%.o:	%.f90
	$(FC) $(FC_FLAGS) $(INC_FLAGS) -I$(general_SOURCES_path) -c $< -o $@

Makefile:	make_dep

make_dep:	$(general_SOURCES)
		./f90_dependencies.py $^ > $@


-include make_dep


