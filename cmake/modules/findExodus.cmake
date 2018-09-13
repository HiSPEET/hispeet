# Exodus

if (DEFINED ENV{EXODUS})

    find_library(EXODUS exodus PATHS $ENV{EXODUS}/lib)
    find_library(NETCDF netcdf PATHS $ENV{EXODUS}/lib)

    set(Exodus_LIBRARIES "${EXODUS};${NETCDF}")
    set(Exodus_INCLUDE_DIRS "$ENV{EXODUS}/include")

endif ()

find_package_handle_standard_args(Exodus DEFAULT_MSG Exodus_LIBRARIES)
