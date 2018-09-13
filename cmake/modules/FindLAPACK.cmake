# BLAS and LAPACK
#
# This module sets
#
#     LAPACK_LIBRARIES
#
# or, in case of Intel, the compiler flags for using MKL

# With Intel use MKL
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    list(APPEND CMAKE_Fortran_FLAGS "-mkl=sequential")
    set(LAPACK_IS_AVAILABLE TRUE)

# With PGI use libraries shipped with the compiler
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    # TBD

else ()

    # Try OpenBLAS first
    find_package(OpenBLAS QUIET)
    if (OpenBLAS_FOUND)
        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})
        set(LAPACK_IS_AVAILABLE TRUE)

    else ()

        # Try to found standalone BLAS and LAPACK libraries
        find_library(BLAS_LIBRARY blas)
        find_library(LAPACK_LIBRARY lapack)
        if (BLAS_LIBRARY AND LAPACK_LIBRARY)
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
            set(LAPACK_IS_AVAILABLE TRUE)
        endif ()

    endif ()


endif ()

find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_IS_AVAILABLE)
