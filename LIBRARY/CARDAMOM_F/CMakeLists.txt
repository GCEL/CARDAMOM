set(FILES1 misc/math_functions.f90
                  misc/brent_zero/brent_zero.f90
                  general/cardamom_structures.f90
                )
# Model to compile. Overridable from the command line with -DMODEL=<modelname>;
# the three model source files are derived from it by the fixed directory layout.
set(MODEL "DALEC.A1.C1.D2.F2.H2.P1.004" CACHE STRING "DALEC model to compile")
set(MODELPARS
                  model/${MODEL}/src/${MODEL}_PARS.f90
                )
set(FILES2
                  model/model_shared.f90
                  general/cardamom_io.f90
                  general/model_likelihood_wrapper.f90
                  general/cardamom_main_utils.f90
                  general/MHMCMC_StressTests.f90
                )
if(${CMAKE_BUILD_TYPE} STREQUAL "RELEASE")
  list(APPEND FILES2 general/wrappers_for_R.f90)
endif()

set(MODELLIKELIHOODFILE
                  model/${MODEL}/likelihood/MODEL_LIKELIHOOD.f90
                )
set(MODELFILE
                  model/${MODEL}/src/${MODEL}.f90
                )

set(ALLMODELNAMES
#DALEC_1005 
#DALEC_1005a
DALEC.A1.C1.D2.F2.H1.P1.003
DALEC.A1.C1.D2.F2.H2.P1.004
DALEC.A1.C1.D2.F2.H2.P1.R1.005
DALEC.A1.C1.D2.F2.H2.P2.018
DALEC.A1.C1.D2.F2.H2.P5.021
DALEC.A1.C1.D2.F2.H2.P6.022
#DALEC.A1.C1.D2.F2.H3.P1.029
DALEC.A1.C2.D2.F2.H1.P3.R1.008
DALEC.A1.C2.D2.F2.H1.P4.R2.010
DALEC.A4.C6.D2.F2.H3.P10.026
DALEC.A1.C2.D2.F2.H2.P1.R1.006
DALEC.A1.C2.D2.F2.H2.P2.R1.007
DALEC.A1.C2.D2.F2.H2.P2.R3.019
DALEC.A1.C2.D2.F2.H2.P3.R1.009
DALEC.A1.C2.D2.F2.H2.P4.R2.011
DALEC.A1.C2.D2.F2.H2.P7.R2.023
DALEC.A2.C1.D2.F2.H2.P1.020
DALEC.A3.C1.D2.F2.H2.P1.030
DALEC.A3.C3.H2.M1.015
DALEC.A3.H2.M2.017
DALEC.A4.C6.D2.F2.H2.P11.031
DALEC.C1.D1.F2.P1.002
#DALEC.C3.M1.014
DALEC.C4.D1.F2.012
DALEC.C5.D1.F2.P1.013
DALEC.D1.F2.001
DALEC.M2.016
)

add_library(cardamom-samplers-lib
              method/math.f90
              method/random_uniform.f90
              method/samplers_io/samplers_io.f90
              method/samplers_shared.f90
              method/APMCMC/APMCMC.f90
              method/APMCMC/MHMCMC.f90
              method/DEMCMC/DEMCz.f90
)

set_target_properties(
  "cardamom-samplers-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)
target_include_directories(
  "cardamom-samplers-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  make_directory("${PROJECT_BINARY_DIR}/include")
endif()

if(USE_OMP)
  find_package(OpenMP)
  if(OpenMP_Fortran_FOUND)
    target_link_libraries(cardamom-samplers-lib PUBLIC OpenMP::OpenMP_Fortran)
  endif()
endif()

if(COMPILE_ALL_MODELS)
  #  loop through compiling with each model in list above, to check that they compile
  foreach(n IN LISTS ALLMODELNAMES)
    set(FILES ${FILES1})
    list(APPEND FILES "model/${n}/src/${n}_PARS.f90")
    list(APPEND FILES ${FILES2})
    list(APPEND FILES "model/${n}/src/${n}.f90")
    list(APPEND FILES "model/${n}/likelihood/MODEL_LIKELIHOOD.f90")
    add_library("cardamom-lib-${n}" STATIC ${FILES})
    add_dependencies("cardamom-lib-${n}" cardamom-samplers-lib) 
    target_link_libraries("cardamom-lib-${n}" PUBLIC cardamom-samplers-lib)
    add_executable("cardamom-${n}" general/cardamom_main.f90)
    target_link_libraries("cardamom-${n}" PUBLIC "cardamom-lib-${n}")
  endforeach()
endif()

# now with the model we chose to run with
add_library(cardamom-lib ${FILES1} ${MODELPARS} ${FILES2} ${MODELLIKELIHOODFILE} ${MODELFILE})
add_dependencies(cardamom-lib cardamom-samplers-lib) 
target_link_libraries(cardamom-lib PUBLIC cardamom-samplers-lib)
add_executable(cardamom general/cardamom_main.f90)
target_link_libraries(cardamom PUBLIC cardamom-lib)
#add_executable(cardamom-diffev general/cardamom_main_diffev.f90)
#target_link_libraries(cardamom-diffev PUBLIC cardamom-lib)

# Place the runnable executables in the legacy LIBRARY/CARDAMOM_F/executable
# directory and give them the legacy .exe file names so downstream tooling
# (cardamom_project_setup.r) can copy them from a fixed, known location.
#set_target_properties(cardamom cardamom-diffev
#  PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/executable")
set_target_properties(cardamom        PROPERTIES OUTPUT_NAME "cardamom.exe")
#set_target_properties(cardamom-diffev PROPERTIES OUTPUT_NAME "cardamom-diffev.exe")

# Shared library used by the R post-processing for forward DALEC runs. Emitted
# as dalec.so (no lib prefix) into the executable directory, with the same MODEL
# selection. Source set mirrors the historical standalone build (brent_zero +
# model + R interface); <MODEL>_R_interface.f90 only `use`s CARBON_MODEL_MOD, so
# no further modules are required here.
#
# dalec.so MUST ALWAYS be built with GNU gfortran, independent of the compiler
# used for the executables / cardamom-lib (which may be ifx). CMake uses one
# Fortran compiler per build tree, so we cannot make an add_library() target use
# a different compiler; instead we invoke gfortran directly via a custom command.
# The sources are self-contained and compiled into their own module directory
# (include_dalec), so there is no dependence on the project's (possibly
# Intel-format) .mod files.
find_program(GFORTRAN_EXECUTABLE gfortran REQUIRED)

set(DALEC_MODDIR "${PROJECT_BINARY_DIR}/include_dalec")
set(DALEC_SO     "${CMAKE_CURRENT_SOURCE_DIR}/executable/dalec.so")
set(DALEC_SRCS
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/brent_zero/brent_zero.f90"
    "${CMAKE_CURRENT_SOURCE_DIR}/${MODELFILE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/model/${MODEL}/src/${MODEL}_R_interface.f90")

# gfortran flags for dalec.so
if(CMAKE_BUILD_TYPE STREQUAL "DEGUG")
  set(DALEC_FFLAGS -O0 -g -fcheck=all -ffree-line-length-0 -fPIC -shared)
else ()
  set(DALEC_FFLAGS -O2 -ffree-line-length-0 -fPIC -shared)
endif()

file(MAKE_DIRECTORY "${DALEC_MODDIR}")
add_custom_command(
  OUTPUT  "${DALEC_SO}"
  COMMAND "${GFORTRAN_EXECUTABLE}"
          ${DALEC_FFLAGS}
          ${DALEC_SRCS}
          -o "${DALEC_SO}"
  WORKING_DIRECTORY "${DALEC_MODDIR}"
  DEPENDS ${DALEC_SRCS}
  COMMENT "Building dalec.so with gfortran (forced GNU, independent of project compiler)"
  VERBATIM)

add_custom_target(dalec DEPENDS "${DALEC_SO}")


