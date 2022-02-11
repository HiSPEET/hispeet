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
include(cmake/tpo__aaa__3d.cmake)
include(cmake/tpo__diagonal__3d.cmake)
include(cmake/tpo__div__3d.cmake)
include(cmake/tpo__diffusion__3d.cmake)
include(cmake/tpo__grad__3d.cmake)
include(cmake/tpo__ins_convection__3d.cmake)
include(cmake/tpo__rot__3d.cmake)
include(cmake/tpo__rotrot__3d.cmake)
include(cmake/tpo__schwarz__3d.cmake)
include(cmake/tpo__spectral__3d.cmake)
include(cmake/tpo__stress__3d.cmake)

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

set ( MESH_3D mesh_3d/generic_mesh__3d.f
              mesh_3d/generic_mesh__3d__mp_createannulus.f
              mesh_3d/generic_mesh__3d__mp_createcylinder.f
              mesh_3d/generic_mesh__3d__mp_creatediamonds.f
              mesh_3d/generic_mesh__3d__mp_createonerotated.f
              mesh_3d/mesh_boundary__3d.f
              mesh_3d/mesh_element__3d.f
              mesh_3d/mesh_element_indexing__3d.f
              mesh_3d/mesh_face__3d.f
              mesh_3d/mesh_link__3d.f
              mesh_3d/mesh__3d.f
              mesh_3d/mesh__3d__mp_buildfaces.f
              mesh_3d/mesh__3d__mp_buildghosts.f
              mesh_3d/mesh__3d__mp_buildlinks.f
              mesh_3d/mesh__3d__mp_getcuboids.f
              mesh_3d/mesh__3d__mp_getpoints.f
              mesh_3d/mesh__3d__mp_getpointvalency.f
              mesh_3d/mesh__3d__mp_identifyedges.f
              mesh_3d/mesh__3d__mp_identifyranks.f
              mesh_3d/mesh__3d__mp_identifyelementneighbors.f
              mesh_3d/mesh__3d__mp_importgenericmesh.f
              mesh_3d/structured_mesh_indexing__3d.f
              mesh_3d/generate_regular_mesh__3d.f
              mesh_3d/mesh_metrics__3d.f
              mesh_3d/spectral_element_mesh__3d.f
              )

#===============================================================================
# util_3d

set( UTIL_3D util_3d/export/export_vtk_volume_data__3d.f
             util_3d/transfer/element_transfer_buffer__3d.f
             util_3d/transfer/element_transfer_buffer__3d__sm_integer_idk.F
             util_3d/transfer/element_transfer_buffer__3d__sm_integer_ixs.F
             util_3d/transfer/element_transfer_buffer__3d__sm_integer_ixl.F
             util_3d/transfer/element_transfer_buffer__3d__sm_real_rdp.F
             util_3d/transfer/element_transfer_buffer__3d__sm_real_rsp.F
             util_3d/transfer/element_face_transfer_buffer__3d.f
             util_3d/transfer/element_face_transfer_buffer__3d__sm_real_rdp.F
             util_3d/transfer/trace_transfer_buffer__3d.f
             util_3d/ops/assembly__3d.f
             util_3d/ops/assembly__3d__sm_0u.f
             util_3d/ops/convert_traces__3d.f
             util_3d/verification/verify_mesh__3d.f
             )

if (Exodus_FOUND)
    set (UTIL_3D ${UTIL_3D} util_3d/import/import_exodus__3d.f)
endif ()

#===============================================================================
# variables_3d

set( VARIABLES_3D variables_3d/spectral_element_variable__3d.f
                  variables_3d/spectral_element_variable__3d__mp_gettraces.f
                  variables_3d/spectral_element_scalar__3d.f
                  variables_3d/spectral_element_vector__3d.f
                  variables_3d/spectral_element_boundary_variable__3d.f
                  )

#===============================================================================
# variables_3d

set( DIFFUSION_3D diffusion_3d/dg__diffusion_operator__3d.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_setdiffusivity.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_apply.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_apply_rc.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_apply_rv.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_apply_dc.F
                  diffusion_3d/dg__diffusion_operator__3d__mp_addfluxes_dc_gen.f
                  diffusion_3d/dg__diffusion_operator__3d__mp_addfluxes_dc_xsmm.F
                  )

#                  diffusion_3d/dg__diffusion_operator__3d__mp_addbc_rc.f
#                  diffusion_3d/dg__diffusion_operator__3d__mp_addbc_rv.f
#                  diffusion_3d/dg__diffusion_operator__3d__mp_addbc_dc.f

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
             ${VARIABLES_3D}
             ${DIFFUSION_3D}
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
