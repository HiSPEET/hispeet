# libxsmm integration ..........................................................

    ExternalProject_Add(LibXSMM 
        PREFIX            ${CMAKE_CURRENT_BINARY_DIR}/external
        GIT_REPOSITORY    ${CMAKE_CURRENT_SOURCE_DIR}/external/libxsmm
        CONFIGURE_COMMAND ""
        BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM}
        BUILD_IN_SOURCE   TRUE
        )

#        BUILD_COMMAND     ${CMAKE_MAKE_PROGRAM}

    set(LibXSMM_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/src/LibXSMM)

    set(LibXSMM_INCLUDE_DIRS "${LibXSMM_INSTALL_DIR}/include" )    
    set(LibXSMM_LIBRARIES    "${LibXSMM_INSTALL_DIR}/lib/libxsmmf.a"
                             "${LibXSMM_INSTALL_DIR}/lib/libxsmm.a" ) 
