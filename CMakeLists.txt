cmake_minimum_required(VERSION 3.22)

project(CARDAMOM VERSION 2.0 
	DESCRIPTION "University of Edinburgh / NCEO (UK) CARbon DAta MOdel FraMework" 
  LANGUAGES Fortran)
enable_language(Fortran)

include(CTest)

option(USE_OMP "Use openmp" ON)
option(COMPILE_ALL_MODELS "Compile all models, resulting in an executable for each.  For checking that all compile." OFF)
option(ENABLE_TESTING "Compile tests" OFF)

# Build-type flags, per compiler. Selected via -DCMAKE_Fortran_COMPILER
# (gfortran -> GNU, ifx -> IntelLLVM). The GNU flags are not valid for ifx, so
# each toolchain gets its own equivalent set.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG "-ffree-line-length-0 -O2 -g -pg -Wall")
  set(CMAKE_Fortran_FLAGS_RELEASE "-ffree-line-length-0 -O3")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -warn all")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
else()
  message(WARNING "Unrecognised Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'; using CMake default flags.")
endif()

# Explicitly mark the stack non-executable at link time so the intended,
# secure behaviour is asserted rather than inferred.
add_link_options(-Wl,-z,noexecstack)

# Make Release the default CMAKE_BUILD_TYPE
IF(NOT CMAKE_BUILD_TYPE)
    get_property(
        CMAKE_BUILD_TYPE_HELPSTRING
        CACHE CMAKE_BUILD_TYPE
        PROPERTY HELPSTRING
    )
    set(
      CMAKE_BUILD_TYPE "RELEASE" CACHE STRING
        "${CMAKE_BUILD_TYPE_HELPSTRING}" FORCE
    )
ENDIF()

IF(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -g")
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF) #for gprof
ELSE()
  option(BUILD_SHARED_LIBS "Build using shared libraries" ON) #For R wrappers
ENDIF()
add_subdirectory(LIBRARY/CARDAMOM_F)

set_target_properties(
  "cardamom-lib"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include"
)

target_include_directories(
  "cardamom-lib"
  PUBLIC
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)
if(NOT EXISTS "${PROJECT_BINARY_DIR}/include")
  make_directory("${PROJECT_BINARY_DIR}/include")
endif()

#install(TARGETS cardamom DESTINATION "bin")

if(ENABLE_TESTING)
  # add the testsuite
  # Code coverage analysis # TODO improve
  target_compile_options(cardamom-samplers-lib PRIVATE --coverage -fprofile-arcs -ftest-coverage) 
  target_link_options(cardamom-samplers-lib PRIVATE -lgcov --coverage) 
  target_compile_options(cardamom-lib PRIVATE --coverage -fprofile-arcs -ftest-coverage) 
  target_link_options(cardamom-lib PRIVATE -lgcov --coverage) 
  enable_testing()
  add_subdirectory("test")
endif()
