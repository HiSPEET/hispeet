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
#         Debug:    enables debugging, may fail with some compilers (Intel)
#
#     -D OPT=⟨n⟩
#         overrides preset optimization level to -On
#
#     -D SKIP_LIBXSMM=⟨0|1⟩
#         0:  use LIBXSMM if available (default)
#         1:  do not use LIBXSMM
#
#     -D OpenMP=⟨0|1⟩
#         0:  compile without OpenMP (default)
#         1:  compile with OpenMP
#
#     -D OpenACC=⟨0|1⟩
#         0:  compile without OpenACC (default)
#         1:  compile with OpenACC

cmake_minimum_required(VERSION 3.10)

project(HiSPEET)

#===============================================================================
# CMake

enable_language(Fortran)
enable_testing()

include(FindPackageHandleStandardArgs)
include(ExternalProject)

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
#find_package(MPI REQUIRED)
find_package(MPI)
if (NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "MPI does not have MPI_F08 module")
endif ()

#-------------------------------------------------------------------------------
# External libraries and utilities

find_package(LAPACK REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(Slurm)
find_package(FFTW)
find_package(ParMETIS)
find_package(Exodus)
find_package(pFUnit QUIET)

#===============================================================================
# Flags

# Set default build type to "Release"
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# Flags for checkingh
if (CMAKE_BUILD_TYPE MATCHES "Check")
    add_definitions(-DCHECK)
    add_definitions(-DOMP_SERIAL_SCALARPRODUCT)
endif ()

# Preprocessor settings for debugging
if (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DDEBUG)
endif ()


#-------------------------------------------------------------------------------
# GNU

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form -Wall")
    set(CMAKE_Fortran_FLAGS_CHECK   "-fcheck=all -fbacktrace -fpe0")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-fcheck=all -fbacktrace -g")

    if (OPT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O${OPT}")
    else (OPT)
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3" )
    endif (OPT)

#-------------------------------------------------------------------------------
# Intel

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free -heap-arrays -warn all")
    #set(CMAKE_Fortran_FLAGS_CHECK   "-CA -traceback -fpe0 -ftrapuv")
    set(CMAKE_Fortran_FLAGS_CHECK   "-CA -traceback -fpe0")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-CA -traceback -fpe0 -ftrapuv -g -gsplit-dwarf")
    set(CMAKE_Fortran_FLAGS_RELEASE "-ip" )

    if (OPT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O${OPT}")
    else (OPT)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3" )
    endif (OPT)

    if (CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -diag-disable=cpu-dispatch")

    elseif (CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
        set(CMAKE_Fortran_FLAGS    "${CMAKE_Fortran_FLAGS} -fno-common")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-stack_size,0x40000000")

    endif ()

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
    set(CMAKE_Fortran_FLAGS_CHECK   "-Mbounds")
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast" )

    if (OPT)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O${OPT}")
    else (OPT)
        set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3" )
    endif (OPT)

    set(Fortran_FLAGS_TPO "-Minline=reshape")

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

    # LIBXSMM is not used
    set(SKIP_LIBXSMM "1")

    ## using CMake module -- not working :(
    ## find_package(OpenACC REQUIRED)

    if (NOT DEFINED OpenACC_ACCEL_TARGET)
        set(OpenACC_ACCEL_TARGET "tesla:kepler+")
    endif ()

    # Fortran and C compiler IDs must match to avoid problems with linker
    if ( CMAKE_Fortran_COMPILER_ID MATCHES "PGI" AND
         CMAKE_C_COMPILER_ID       MATCHES "PGI" )

        set( OpenACC_C_FLAGS
             "-acc -Minfo=accel -Minline=3 -ta=${OpenACC_ACCEL_TARGET}" )

        set( OpenACC_Fortran_FLAGS
             "-acc -Minfo=accel -Minline=3 -ta=${OpenACC_ACCEL_TARGET}" )

        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Mcuda -lcublas")

        set(OpenACC_FOUND TRUE)
        message (STATUS "Found OpenACC (PGI)")

    else ()
        message(FATAL_ERROR "Could NOT find OpenACC")
    endif ()

    # Add OpenACC compiler flags ()
    set(CMAKE_C_FLAGS       "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
    set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS} ${OpenACC_CXX_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenACC_Fortran_FLAGS}")

endif (OpenACC)

#===============================================================================
# External libraries

include(external/cmake/libxsmm.cmake)

#===============================================================================
# Subprojects

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(app)
add_subdirectory(develop)

add_subdirectory(dirty)
