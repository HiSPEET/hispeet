#===============================================================================
# essentials

file(GLOB ESSENTIALS essentials/*.f)

#===============================================================================
# bindings

set(BINDINGS bindings/c_binding.f
             bindings/mpi_binding.F
             bindings/vtk_wrapper.c
             bindings/vtk_binding.f
             )

if (ParMETIS_FOUND)
    list(APPEND BINDINGS bindings/parmetis_wrapper.c
                         bindings/parmetis_binding.f )
endif ()

if (Exodus_FOUND)
    list(APPEND BINDINGS bindings/exodus_wrapper.c
                         bindings/exodus_binding.f )
endif ()

#===============================================================================
# xmpi

file(GLOB XMPI xmpi/*.f)

#===============================================================================
# linear_algebra

set(LINEAR_ALGEBRA linear_algebra/array_assignments.F
                   linear_algebra/array_reductions.F
                   linear_algebra/eigenproblems.f
                   linear_algebra/linear_equations.f
                   linear_algebra/matrix_operators.f
                   linear_algebra/singular_value_decomposition.f
                   )

#===============================================================================
# tensor_product

# pre-existing sources
set(TENSOR_PRODUCT tensor_product/tpo_sddd.f
                   tensor_product/tpo_sddd__gen.f )

# configuration script
set(TPO_CONFIG ${PROJECT_SOURCE_DIR}/bin/tpo_configure.py)

# target directory of generated sources
set(DEST ${CMAKE_CURRENT_BINARY_DIR})

# parametrized operators
include(cmake/tpo_aaa.cmake)

# activate extra compile flags
if (DEFINED Fortran_FLAGS_TPO)
    set_source_files_properties(${TENSOR_PRODUCT}
        PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_TPO})
endif ()

#===============================================================================
# utilities

file(GLOB UTILITIES utilities/*.f)

#===============================================================================
# spectral_element

file(GLOB SPECTRAL_ELEMENT spectral_element/*.f)

#===============================================================================
# time_integration

file(GLOB TIME_INTEGRATION time_integration/*.f)

#===============================================================================
# mesh

file(GLOB MESH mesh/*.f)

#===============================================================================
# base library

set(BASE_SRC ${ESSENTIALS}
             ${BINDINGS}
             ${XMPI}
             ${LINEAR_ALGEBRA}
             ${TENSOR_PRODUCT}
             ${SPECTRAL_ELEMENT}
             ${TIME_INTEGRATION}
             ${UTILITIES}
             ${MESH}
             )

add_library(HiBase ${BASE_SRC})

set_property(TARGET HiBase
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_include_directories(HiBase
    PUBLIC  ${PROJECT_BINARY_DIR}/modules
    PRIVATE ${VTK_INCLUDE_DIRS}
            ${ParMETIS_INCLUDE_DIRS}
            ${Exodus_INCLUDE_DIRS}
            )

target_link_libraries(HiBase
    PUBLIC  ${LAPACK_LIBRARIES}
            ${VTK_LIBRARIES}
            ${ParMETIS_LIBRARIES}
            ${Exodus_LIBRARIES}
            )
