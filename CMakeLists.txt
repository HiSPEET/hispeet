# Recommended environment
#
#     FC      set to MPI/Fortran compiler
#     CC      set to MPI/C compiler
#     CXX     set to MPI/C++ compiler
#     CMAKE_Fortran_COMPILER

cmake_minimum_required(VERSION 3.5)

project(HiSPEET)

#===============================================================================
# CMake

enable_language(Fortran)
enable_testing()

include(FindPackageHandleStandardArgs)

# path to CMake modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# path to Python packages
set(PYTHONPATH ${PROJECT_SOURCE_DIR}/bin)

# directories for generated libraries and archives
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

#===============================================================================
# Prerequisites

#-------------------------------------------------------------------------------
# MPI

set(MPI_SKIP_COMPILER_WRAPPER TRUE)
set(MPI_SKIP_GUESSING         TRUE)
find_package(MPI REQUIRED)
if (NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "MPI does not have MPI_F08 module")
endif ()

#-------------------------------------------------------------------------------
# External libraries and utilities

find_package(LAPACK REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(Slurm)
find_package(ParMETIS)
find_package(Exodus)
find_package(pFUnit PATHS $ENV{PFUNIT})

#===============================================================================
# Flags

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

#-------------------------------------------------------------------------------
# GNU

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -Wall")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fbacktrace -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native" )

    # compiler environment module needed for generating slurm scripts
    set(ENV_MODULE_COMPILER "compiler/gnu")

#-------------------------------------------------------------------------------
# Intel

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -heap-arrays -warn all")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-traceback -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -xHost" )

    if (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable=cpu-dispatch")

    elseif (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} -fno-common")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size,0x40000000")

    endif ()

    # do not warn about using Intel compiler directives
    set(Fortran_FLAGS_TPO "-diag-disable=6477")

    # compiler environment module needed for generating slurm scripts
    set(ENV_MODULE_COMPILER "compiler/intel")

#-------------------------------------------------------------------------------
# PGI

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")

    string(JOIN \  CMAKE_Fortran_FLAGS
           ${CMAKE_Fortran_FLAGS}
           -Mfree -Minform=inform
           -Mbackslash
           -Minfo=inline
           -Minline=maxsize:250
           -Mvect=simd
           -Munroll
          )

    set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast -O3" )

    set(Fortran_FLAGS_TPO "-Minline=reshape")

    # compiler environment module needed for generating slurm scripts
    set(ENV_MODULE_COMPILER "compiler/pgi")

endif ()

#-------------------------------------------------------------------------------
# OpenMP

if (OpenMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OpenMP)

#-------------------------------------------------------------------------------
# OpenACC

if (OpenACC)
    find_package(OpenACC REQUIRED)
    if (NOT DEFINED OpenACC_ACCEL_TARGET)
        set(OpenACC_ACCEL_TARGET tesla:kepler+)
    endif()
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
    set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} ${OpenACC_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")

    # linker flags
    if (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Mcuda -lcublas")
    endif()

endif ()

#===============================================================================
# Subprojects

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
add_subdirectory(develop)

add_subdirectory(dirty)
