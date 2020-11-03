# BLAS and LAPACK

# With Intel use MKL
if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential")
    set(LAPACK_FOUND TRUE)
    message(STATUS "Found LAPACK: using MKL")

# With PGI use libraries shipped with the compiler
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -lblas -llapack")
    set(LAPACK_FOUND TRUE)
    message(STATUS "Found LAPACK: using libraries shipped with PGI compiler")

else ()

    # Try OpenBLAS first
    find_package(OpenBLAS QUIET)

    if (OpenBLAS_FOUND)

        set(LAPACK_LIBRARIES ${OpenBLAS_LIBRARIES})

    else ()

        # Try to found standalone BLAS and LAPACK libraries
        find_library(BLAS_LIBRARY blas)
        find_library(LAPACK_LIBRARY lapack)
        if (BLAS_LIBRARY AND LAPACK_LIBRARY)
            set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
        endif ()

    endif ()

    find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARIES)

endif ()
