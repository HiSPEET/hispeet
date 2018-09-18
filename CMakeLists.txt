# Required environment
#
#     CC     C compiler
#     CXX    C++ compiler
#     FC     Fortran compiler
#
# must support MPI

cmake_minimum_required(VERSION 3.11)

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

# directories for generated libraries, archives and module files
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

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
# OpenMP

if (OpenMP)
    find_package(OpenMP REQUIRED)
endif ()

#-------------------------------------------------------------------------------
# External libraries and utilities

find_package(LAPACK REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)
find_package(PythonInterp 3 REQUIRED)
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

#-------------------------------------------------------------------------------
# PGI  *** incomplete ***

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")

    set(Fortran_FLAGS_TPO "-Minline=reshape")

endif ()

#-------------------------------------------------------------------------------
# OpenMP

if (OpenMP)
    set(CMAKE_C_FLAGS          "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OpenMP)

#===============================================================================
# Subprojects

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
add_subdirectory(develop)

#add_subdirectory(dirty)
