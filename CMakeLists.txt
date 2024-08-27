# Environment
#
#     FC            path to MPI/Fortran compiler
#     CC            path to MPI/C compiler
#     CXX           path to MPI/C++ compiler
#
#     FFTW_ROOT     root path of FFTW3 (if not in system paths)
#
# Command-line options for cmake
#
#     -D CMAKE_BUILD_TYPE=⟨Release|Check|Debug⟩
#         Release:  default build type
#         Check:    similar to Release, with additional runtime checks
#         Debug:    enables debugging, may fail with some case
#
#     -D OPT=⟨n⟩
#         overrides preset optimization level to -On
#
#     -D REPORT=⟨0|1⟩
#         0:  do not generate optimization report (default)
#         1:  generate optimization report (only with Intel, so far)
#
#     -D SKIP_LIBXSMM=⟨0|1⟩
#         0:  use LIBXSMM if available (default)
#         1:  do not use LIBXSMM
#
#     -D OpenMP=⟨0|1⟩
#         0:  compile without OpenMP (default)
#         1:  compile with OpenMP

cmake_minimum_required(VERSION 3.12)

project(HiSPEET)

#===============================================================================
# CMake

enable_language(Fortran)
enable_testing()

include(FindPackageHandleStandardArgs)
include(ExternalProject)

# path to CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# directories for generated libraries and archives
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#===============================================================================
# Prerequisites

#-------------------------------------------------------------------------------
# OpenMP

if (OpenMP OR NOT SKIP_LIBXSMM)
#   find_package(OpenMP REQUIRED)
   find_package(OpenMP)
endif ()

#-------------------------------------------------------------------------------
# MPI

set(MPI_SKIP_COMPILER_WRAPPER TRUE)
set(MPI_SKIP_GUESSING         TRUE)
find_package(MPI)
if (NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "MPI does not have MPI_F08 module")
endif ()

#===============================================================================
# Flags

# Set default build type to "Release"
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# Flags for checking
if (CMAKE_BUILD_TYPE MATCHES "Check")
    add_compile_definitions(CHECK)
endif ()

# Preprocessor settings for debugging
if (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_compile_definitions(DEBUG)
endif ()

# Preprocessor setting for using the INS grad-div penalty instead bulk viscosity
if (INS_PENALTY)
    add_compile_definitions(__INS_PENALTY__)
endif ()

#-------------------------------------------------------------------------------
# GNU

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")

    set(CMAKE_Fortran_FLAGS_CHECK "-Wall -Wno-c-binding-type -fcheck=all -fbacktrace")
    set(CMAKE_Fortran_FLAGS_DEBUG "-Wall -Wno-c-binding-type -fcheck=all -fbacktrace -g")

    if (OPT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O${OPT}")
    else (OPT)
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    endif (OPT)

#-------------------------------------------------------------------------------
# Intel

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -heap-arrays")
    set(CMAKE_Fortran_FLAGS_CHECK   "-warn all -CA -traceback -fpe0")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-warn all -CA -traceback -fpe0 -g -gsplit-dwarf")
    set(CMAKE_Fortran_FLAGS_RELEASE "-ip" )

    if (OPT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O${OPT}")
    else (OPT)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3" )
    endif (OPT)

    if (REPORT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qopt-report=4")
    endif (REPORT)

    if (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable=cpu-dispatch")

    elseif (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} -fno-common")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size,0x40000000")

    endif ()

endif ()

#-------------------------------------------------------------------------------
# OpenMP

if (OpenMP)
    set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
elseif (NOT SKIP_LIBXSMM)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OpenMP)

#===============================================================================
# External libraries and utilities

#-------------------------------------------------------------------------------
# BLAS and LAPACK

if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(LAPACK_LIBRARIES "-qmkl=sequential")
    set(LAPACK_FOUND TRUE)
    message(STATUS "Found LAPACK: using MKL")
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    message(STATUS "Found LAPACK libraries: ${LAPACK_LIBRARIES}") 
endif ()

#-------------------------------------------------------------------------------
# VTK XML export

find_package(VTK REQUIRED COMPONENTS IOXML)

# provides: 
#   - VTK_INCLUDE_DIRS 
#   - VTK_LIBRARIES

#-------------------------------------------------------------------------------
# HDF5

find_package(HDF5 REQUIRED COMPONENTS Fortran HL)

# provides: 
#   - HDF5_Fortran_INCLUDE_DIRS 
#   - HDF5_Fortran_LIBRARIES

#-------------------------------------------------------------------------------
# KD-Tree

include(external/cmake/kdtree.cmake)

#-------------------------------------------------------------------------------
# LIBXSMM

include(external/cmake/libxsmm.cmake)

# provides: 
#   - LIBXSMM_INCLUDE_DIRS 
#   - LIBXSMM_LIBRARIES

#-------------------------------------------------------------------------------
# FFTW, ParMeTiS, ppFUnit, SLURM

find_package(FFTW)         # FFTW_INCLUDE_DIRS    , FFTW_LIBRARIES
find_package(ParMETIS)     # ParMETIS_INCLUDE_DIRS, ParMETIS_LIBRARIES
find_package(pFUnit QUIET)
find_package(Slurm QUIET)

#===============================================================================
# Subprojects

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(app)
add_subdirectory(develop)

# experimental programs excluded from make all
add_subdirectory(dirty EXCLUDE_FROM_ALL)
