# What is needed to be set
#
# CC  C compiler/wrapper supporting MPI
# FC  Fortran compiler/wrapper supporting MPI

cmake_minimum_required(VERSION 3.11)

project(HiSPEET)

enable_language(Fortran)
enable_testing()
include(FindPackageHandleStandardArgs)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# MPI
set(MPI_SKIP_COMPILER_WRAPPER TRUE)
set(MPI_SKIP_GUESSING         TRUE)
find_package(MPI REQUIRED)
if (NOT MPI_Fortran_HAVE_F08_MODULE)
    message(FATAL_ERROR "MPI does not have MPI_F08 module")
endif ()

# Compilers
if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    list(APPEND CMAKE_Fortran_FLAGS "-ffree-form -O3")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    list(APPEND CMAKE_Fortran_FLAGS "-free -stand:f08 -O3")
    set(Fortran_FLAGS_TPO "-diag-disable=6477")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(Fortran_FLAGS_TPO -Minline=reshape)
endif ()

# path to own Python packages
set(PYTHONPATH ${PROJECT_SOURCE_DIR}/bin)

# directories for generated libraries, archives and module files
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# External libraries and utilities .............................................

find_package(LAPACK REQUIRED)
find_package(VTK COMPONENTS vtkIOXML REQUIRED)
find_package(PythonInterp 3 REQUIRED)
find_package(ParMETIS)
find_package(Exodus)
find_package(pFUnit PATHS $ENV{PFUNIT})

# subprojects ................................................................

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
add_subdirectory(develop)

#add_subdirectory(dirty)
