# ParMETIS

find_library(METIS metis)
find_library(ParMETIS parmetis)

set(ParMETIS_LIBRARIES "${ParMETIS};${METIS}")

find_path(ParMETIS_INCLUDE_DIRS parmetis.h)

find_package_handle_standard_args( ParMETIS DEFAULT_MSG 
                                   ParMETIS_LIBRARIES
                                   ParMETIS_INCLUDE_DIRS )

if (ParMETIS_FOUND)
    add_compile_definitions(__ParMETIS__)
else ()
    set(ParMETIS_LIBRARIES "")
    set(ParMETIS_INCLUDE_DIRS "")
endif ()

mark_as_advanced(ParMETIS_LIBRARIES ParMETIS_INCLUDE_DIRS)

