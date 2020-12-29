# Exodus

find_library(Exodus exodus)
find_library(NetCDF netcdf)

set(Exodus_LIBRARIES "${Exodus};${NetCDF}")

find_path(Exodus_INCLUDE_DIRS exodusII.h)

find_package_handle_standard_args( Exodus DEFAULT_MSG 
                                   Exodus_LIBRARIES 
                                   Exodus_INCLUDE_DIRS )

if (NOT Exodus_FOUND)
    set(Exodus_LIBRARIES "")
    set(Exodus_INCLUDE_DIRS "")
else ()
    add_compile_definitions(__EXODUS__)
endif ()

mark_as_advanced(Exodus_LIBRARIES Exodus_INCLUDE_DIRS)
