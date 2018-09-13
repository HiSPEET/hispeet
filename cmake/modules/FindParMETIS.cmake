# ParMETIS

if (DEFINED ENV{PARMETIS})

    find_library(METIS    metis    PATHS $ENV{PARMETIS}/lib)
    find_library(PARMETIS parmetis PATHS $ENV{PARMETIS}/lib)

    set(ParMETIS_LIBRARIES "${PARMETIS};${METIS}")
    set(ParMETIS_INCLUDE_DIRS "$ENV{PARMETIS}/include")

endif ()

find_package_handle_standard_args(ParMETIS DEFAULT_MSG ParMETIS_LIBRARIES)
