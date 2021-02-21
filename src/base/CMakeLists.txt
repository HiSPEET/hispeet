#===============================================================================
# essentials

file(GLOB ESSENTIALS essentials/*.f)

#===============================================================================
# bindings

set(BINDINGS bindings/c_binding.f
             bindings/openmp_binding.F
             bindings/mpi_binding.F
             bindings/vtk_wrapper.c
             bindings/vtk_binding.f
             bindings/fftw_binding.f
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

# add primitive directories to include path
include_directories( tensor_product/primitives/IxIxQt
                     tensor_product/primitives/IxQtxI
                     tensor_product/primitives/QtxIxI
                     tensor_product/primitives/IxIxAt 
                     tensor_product/primitives/IxBtxI
                     tensor_product/primitives/CtxIxI )

# parametrized operators
include(cmake/tpo__aaa_3d.cmake)
include(cmake/tpo__diagonal_3d.cmake)
include(cmake/tpo__elliptic_3d.cmake)
include(cmake/tpo__div_3d.cmake)
include(cmake/tpo__rot_3d.cmake)
include(cmake/tpo__grad_3d.cmake)
include(cmake/tpo__rotrot_3d.cmake)
include(cmake/tpo__spectral_3d.cmake)
include(cmake/tpo__schwarz_3d.cmake)

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
# mesh_3d

set ( MESH_3D mesh_3d/generic_mesh_3d.f
              mesh_3d/generic_mesh_3d__mp_createcylinder.f
              mesh_3d/generic_mesh_3d__mp_createannulargap.f
              mesh_3d/mesh_3d__boundary.f
              mesh_3d/mesh_3d__element_indexing.f
              mesh_3d/mesh_3d__element.f
              mesh_3d/mesh_3d__face.f
              mesh_3d/mesh_3d__link.f
              mesh_3d/mesh_3d__partition.f
              mesh_3d/mesh_3d__partition__mp_buildfaces.f
              mesh_3d/mesh_3d__partition__mp_buildlinks.f
              mesh_3d/mesh_3d__partition__mp_buildghosts.f
              mesh_3d/mesh_3d__partition__mp_getpoints.f
              mesh_3d/mesh_3d__partition__mp_getpointvalency.f
              mesh_3d/mesh_3d__partition__mp_identifyedges.f
              mesh_3d/mesh_3d__partition__mp_identifyprimaries.f
              mesh_3d/mesh_3d__partition__mp_importgenericmesh.f
              mesh_3d/mesh_3d__structured_indexing.f
              mesh_3d/mesh_3d__generate_regular_mesh.f
              )

#===============================================================================
# util_3d

set( UTIL_3D util_3d/export/export_vtk_3d__volume_data.f
             util_3d/transfer/element_transfer_buffer_3d.f
             util_3d/transfer/element_transfer_buffer_3d__sm_integer_idk.F
             util_3d/transfer/element_transfer_buffer_3d__sm_integer_ixs.F
             util_3d/transfer/element_transfer_buffer_3d__sm_integer_ixl.F
             util_3d/transfer/element_transfer_buffer_3d__sm_real_rnp.F
             util_3d/transfer/trace_transfer_buffer_3d.f
             util_3d/verification/verify_mesh_3d.f
             )

if (Exodus_FOUND)
    set (UTIL_3D ${UTIL_3D} util_3d/import/import_exodus_3d.f)
endif ()

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
             ${MESH_3D}
             ${UTIL_3D}
             )

add_library(HiBase ${BASE_SRC})

set_property(TARGET HiBase
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

if (NOT SKIP_LIBXSMM)
    add_dependencies(HiBase LIBXSMM)
endif ()

target_include_directories(HiBase
    PUBLIC  ${PROJECT_BINARY_DIR}/modules
            ${LIBXSMM_INCLUDE_DIRS}
    PRIVATE ${VTK_INCLUDE_DIRS}
            ${ParMETIS_INCLUDE_DIRS}
            ${FFTW_INCLUDE_DIRS}
            ${Exodus_INCLUDE_DIRS}
            )

target_link_libraries(HiBase
    PUBLIC  ${LAPACK_LIBRARIES}
            ${VTK_LIBRARIES}
            ${FFTW_LIBRARIES}
            ${ParMETIS_LIBRARIES}
            ${Exodus_LIBRARIES}
            ${LIBXSMM_LIBRARIES}
            )
