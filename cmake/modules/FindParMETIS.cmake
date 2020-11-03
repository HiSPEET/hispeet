# ParMETIS

find_library(METIS metis)
find_library(ParMETIS parmetis)

set(ParMETIS_LIBRARIES "${ParMETIS};${METIS}")

find_path(ParMETIS_INCLUDE_DIRS parmetis.h)

find_package_handle_standard_args( ParMETIS DEFAULT_MSG 
                                   ParMETIS_LIBRARIES
                                   ParMETIS_INCLUDE_DIRS )

mark_as_advanced(ParMETIS_LIBRARIES ParMETIS_INCLUDE_DIRS)
