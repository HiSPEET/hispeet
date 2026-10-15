#-------------------------------------------------------------------------------
# LIBXSMM integration
#
# Sources:  git submodule external/libxsmm
#           -> the commit recorded in the superproject is the ONLY version pin.
#              Update:  git submodule update --remote external/libxsmm
#                       git add external/libxsmm && git commit
#
# Build:    LIBXSMM's own Makefile, run out-of-tree in the CMake build tree
#           (source tree stays clean, each build directory has its own copy)
#           using the same C/C++/Fortran compilers as HiSPEET.
#
# Provides: imported target LIBXSMM::LIBXSMM
#             - Fortran module libxsmm.mod (include directory)
#             - static libraries libxsmmf.a, [libxsmmext.a,] libxsmm.a
#               + system libraries
#               (libxsmmext is linked only if the checked-out version builds
#                it by default: older versions, e.g. 50c6702, need it because
#                the Fortran module references libxsmm_gemm_batch_omp;
#                newer ones, e.g. c14cbc6, no longer build it)
#             - compile definition __LIBXSMM__

option(SKIP_LIBXSMM "Build HiSPEET without LIBXSMM" OFF)

if (SKIP_LIBXSMM)
  message(STATUS "LIBXSMM: disabled (SKIP_LIBXSMM=ON)")
  return()
endif ()

# directories ..................................................................

set(LIBXSMM_SOURCE_DIR  ${PROJECT_SOURCE_DIR}/external/libxsmm)
set(LIBXSMM_BINARY_DIR  ${PROJECT_BINARY_DIR}/external/libxsmm)
set(LIBXSMM_INCLUDE_DIR ${LIBXSMM_BINARY_DIR}/include)
set(LIBXSMM_LIB_DIR     ${LIBXSMM_BINARY_DIR}/lib)

if (NOT EXISTS ${LIBXSMM_SOURCE_DIR}/Makefile)
  message(FATAL_ERROR
    "LIBXSMM submodule not found in ${LIBXSMM_SOURCE_DIR}.\n"
    "Run  'git submodule update --init external/libxsmm'  "
    "or configure with  -D SKIP_LIBXSMM=ON")
endif ()

find_package(Git QUIET)
if (GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
                  WORKING_DIRECTORY ${LIBXSMM_SOURCE_DIR}
                  OUTPUT_VARIABLE   LIBXSMM_GIT_VERSION
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_QUIET)
endif ()
message(STATUS "LIBXSMM: using submodule ${LIBXSMM_GIT_VERSION}")

# build options ................................................................

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(LIBXSMM_DBG 1)
else ()
  set(LIBXSMM_DBG 0)
endif ()

include(ProcessorCount)
ProcessorCount(LIBXSMM_NPROC)
if (LIBXSMM_NPROC EQUAL 0)
  set(LIBXSMM_NPROC 1)
endif ()

find_program(GNU_MAKE NAMES gmake make)
if (NOT GNU_MAKE)
  message(FATAL_ERROR "LIBXSMM: GNU make is required to build LIBXSMM")
endif ()

# re-configure whenever the submodule version (i.e. its Makefile) changes
set_property(DIRECTORY APPEND PROPERTY
             CMAKE_CONFIGURE_DEPENDS ${LIBXSMM_SOURCE_DIR}/Makefile)

# does this LIBXSMM version build libxsmmext by default ("libs: ... elib")?
file(STRINGS ${LIBXSMM_SOURCE_DIR}/Makefile LIBXSMM_LIBS_RULE
     REGEX "^libs:.*[ \t]elib([ \t]|$)")

set(LIBXSMM_LIBRARIES ${LIBXSMM_LIB_DIR}/libxsmmf.a)
if (LIBXSMM_LIBS_RULE)
  list(APPEND LIBXSMM_LIBRARIES ${LIBXSMM_LIB_DIR}/libxsmmext.a)
endif ()
list(APPEND LIBXSMM_LIBRARIES ${LIBXSMM_LIB_DIR}/libxsmm.a)

# external project .............................................................

ExternalProject_Add(LIBXSMM
  PREFIX            ${PROJECT_BINARY_DIR}/external
  SOURCE_DIR        ${LIBXSMM_SOURCE_DIR}
  BINARY_DIR        ${LIBXSMM_BINARY_DIR}
  DOWNLOAD_COMMAND  ""
  UPDATE_COMMAND    ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND     ${GNU_MAKE} -f ${LIBXSMM_SOURCE_DIR}/Makefile
                    -j${LIBXSMM_NPROC}
                    CC=${CMAKE_C_COMPILER}
                    CXX=${CMAKE_CXX_COMPILER}
                    FC=${CMAKE_Fortran_COMPILER}
                    STATIC=1 OMP=0 BLAS=0 DBG=${LIBXSMM_DBG}
  BUILD_ALWAYS      TRUE      # make is incremental; picks up submodule updates
  BUILD_BYPRODUCTS  ${LIBXSMM_LIBRARIES}
                    ${LIBXSMM_INCLUDE_DIR}/libxsmm.mod
  INSTALL_COMMAND   ""
  )

# imported target ..............................................................

# INTERFACE_INCLUDE_DIRECTORIES must exist at configure time
file(MAKE_DIRECTORY ${LIBXSMM_INCLUDE_DIR})

find_package(Threads REQUIRED)

add_library(LIBXSMM::LIBXSMM INTERFACE IMPORTED GLOBAL)
set_target_properties(LIBXSMM::LIBXSMM PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${LIBXSMM_INCLUDE_DIR}"
  INTERFACE_COMPILE_DEFINITIONS "__LIBXSMM__"
  INTERFACE_LINK_LIBRARIES      "${LIBXSMM_LIBRARIES};Threads::Threads;${CMAKE_DL_LIBS};m"
  )
add_dependencies(LIBXSMM::LIBXSMM LIBXSMM)
