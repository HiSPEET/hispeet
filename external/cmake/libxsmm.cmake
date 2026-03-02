#-------------------------------------------------------------------------------
# LIBXSMM integration
#
# NOTE:
# GIT_TAG is fixed for stability, but should be updated from time to time

if (NOT SKIP_LIBXSMM)

  # directories, modules and libraries .........................................

  set( LIBXSMM_SRC_DIR       ${CMAKE_CURRENT_BINARY_DIR}/external/src/LIBXSMM )
  set( LIBXSMM_LIB_DIR       "${LIBXSMM_SRC_DIR}/lib"                 )
  set( LIBXSMM_MODULE        "${LIBXSMM_SRC_DIR}/include/libxsmm.mod" )
  set( LIBXSMM_LIBRARIES     "${LIBXSMM_SRC_DIR}/lib/libxsmmf.a"
                             "${LIBXSMM_SRC_DIR}/lib/libxsmm.a"       ) 
###                             "${LIBXSMM_SRC_DIR}/lib/libxsmmext.a"
  set( LIBXSMM_INCLUDE_DIRS  "${LIBXSMM_SRC_DIR}/include"             )    

  # compiler and flags .........................................................

###  if ( CMAKE_SYSTEM_NAME    MATCHES "Linux" AND
###       CMAKE_C_COMPILER_ID  MATCHES "GNU"     )
###
###      set( LIBXSMM_LIBRARIES ${LIBXSMM_LIBRARIES} "-lpthread -lrt -ldl -lm -lc" )
###
###  endif ()

  if (CMAKE_BUILD_TYPE MATCHES "Debug")
      set(DBG "1")
  else ()
      set(DBG "0")
  endif ()

  ExternalProject_Add(LIBXSMM 
      PREFIX               ${CMAKE_CURRENT_BINARY_DIR}/external
      GIT_REPOSITORY       ${CMAKE_CURRENT_SOURCE_DIR}/external/libxsmm
#     GIT_TAG              "main"
###     GIT_TAG              "50c6702"
      GIT_TAG              "c14cbc6"
#     UPDATE_DISCONNECTED  TRUE
      CONFIGURE_COMMAND    ""
      BUILD_COMMAND        test -d ${LIBXSMM_LIB_DIR} || make OMP=1 BLAS=0 DBG=${DBG}
      SOURCE_DIR           ${LIBXSMM_SRC_DIR}
      BUILD_IN_SOURCE      TRUE
      BUILD_BYPRODUCTS     ${LIBXSMM_MODULE} ${LIBXSMM_LIBRARIES}
      INSTALL_COMMAND      ""
      )

  add_compile_definitions(__LIBXSMM__)

endif ()
