# FFTW

find_library(FFTW3 fftw3)

set(FFTW_LIBRARIES "${FFTW3}")

find_path(FFTW_INCLUDE_DIRS fftw3.f03)

find_package_handle_standard_args( FFTW DEFAULT_MSG 
                                   FFTW_LIBRARIES 
                                   FFTW_INCLUDE_DIRS )

mark_as_advanced(FFTW_LIBRARIES FFTW_INCLUDE_DIRS)
