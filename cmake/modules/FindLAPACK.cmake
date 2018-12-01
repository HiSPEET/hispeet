# BLAS and LAPACK
#
# This module sets
#
#     LAPACK_LIBRARIES
#
# or, in case of Intel, the compiler flags for using MKL

# With Intel use MKL
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential")
    set(LAPACK_FOUND TRUE)

# With PGI use libraries shipped with the compiler
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lblas -llapack")
    set(LAPACK_FOUND TRUE)

else ()

    # Try OpenBLAS first
    find_package(OpenBLAS QUIET)
    if (OpenBLAS_FOUND)
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
        set(LAPACK_FOUND TRUE)

    else ()

        # Try to found standalone BLAS and LAPACK libraries
        find_library(BLAS_LIBRARY blas)
        find_library(LAPACK_LIBRARY lapack)
        if (BLAS_LIBRARY AND LAPACK_LIBRARY)
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
            set(LAPACK_FOUND TRUE)
        endif ()

    endif ()

endif ()

if (LAPACK_FOUND)
    if (NOT LAPACK_FIND_QUIETLY)
        message(STATUS "Found LAPACK")
    endif()
else ()
    set(LAPACK_FOUND FALSE)
    if (LAPACK_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find LAPACK")
    endif()
endif ()

