# libxsmm integration ..........................................................

if ( NOT DEFINED OpenACC )

  set( LibXSMM_SOURCE_DIR    ${CMAKE_CURRENT_BINARY_DIR}/external/src/LibXSMM )
  set( LibXSMM_INCLUDE_DIR  "${LibXSMM_SOURCE_DIR}/include"      )    
  set( LibXSMM_LIB_DIR      "${LibXSMM_SOURCE_DIR}/lib"          )
  set( LIBXSMM_MODULE       "${LibXSMM_INCLUDE_DIR}/libxsmm.mod" )
  set( LibXSMM_LIBRARIES    "${LibXSMM_LIB_DIR}/libxsmmf.a"
                            "${LibXSMM_LIB_DIR}/libxsmm.a"       ) 

  if ( CMAKE_SYSTEM_NAME    MATCHES "Linux" AND
       CMAKE_C_COMPILER_ID  MATCHES "GNU"     )

      set( LibXSMM_LIBRARIES ${LibXSMM_LIBRARIES} "-lpthread -lrt -ldl -lm -lc" )

  endif ()

  ExternalProject_Add(LibXSMM 
      PREFIX               ${CMAKE_CURRENT_BINARY_DIR}/external
      GIT_REPOSITORY       ${CMAKE_CURRENT_SOURCE_DIR}/external/libxsmm
      UPDATE_DISCONNECTED  TRUE
      CONFIGURE_COMMAND    ""
      BUILD_COMMAND        test -d ${LibXSMM_LIB_DIR} || make
      SOURCE_DIR           ${LibXSMM_SOURCE_DIR}
      BUILD_IN_SOURCE      TRUE
      BUILD_BYPRODUCTS     ${LIBXSMM_MODULE} ${LibXSMM_LIBRARIES}
      INSTALL_COMMAND      ""
      )

  add_definitions(-D__LIBXSMM__)

endif  ( NOT DEFINED OpenACC )
