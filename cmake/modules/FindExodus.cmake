# Exodus

find_library(Exodus exodus)
find_library(NetCDF netcdf)

set(Exodus_LIBRARIES "${Exodus};${NetCDF}")

find_path(Exodus_INCLUDE_DIRS exodusII.h)

find_package_handle_standard_args( Exodus DEFAULT_MSG 
                                   Exodus_LIBRARIES 
                                   Exodus_INCLUDE_DIRS )

mark_as_advanced(Exodus_LIBRARIES Exodus_INCLUDE_DIRS)
