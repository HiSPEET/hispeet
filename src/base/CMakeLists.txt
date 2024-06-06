#===============================================================================
# essentials

file(GLOB ESSENTIALS essentials/*.f)

#===============================================================================
# bindings

set(BINDINGS bindings/c_binding.f
             bindings/hdf5_binding.f
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
include(cmake/tpo__average__3d.cmake)
include(cmake/tpo__diagonal__3d.cmake)
include(cmake/tpo__div__3d.cmake)
include(cmake/tpo__elliptic__3d.cmake)
include(cmake/tpo__grad__3d.cmake)
include(cmake/tpo__ins_convection__3d.cmake)
include(cmake/tpo__rot__3d.cmake)
include(cmake/tpo__rotrot__3d.cmake)
include(cmake/tpo__schwarz__3d.cmake)
include(cmake/tpo__spectral__3d.cmake)
include(cmake/tpo__stress__3d.cmake)

# legacy Schwarz TPO for CART modules
include(cmake/tpo__schwarz__3d_cart.cmake)

# activate extra compile flags
if (DEFINED Fortran_FLAGS_TPO)
    set_source_files_properties(${TENSOR_PRODUCT}
        PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_TPO})
endif ()

#===============================================================================
# utilities

set ( UTILITIES 
      utilities/cache_blocking.f
      utilities/execution_control.f
      utilities/logging_levels.f
      utilities/quick_sort.f
      utilities/schwarz_weighting.f
    )

#===============================================================================
# spectral_element

set ( SPECTRAL_ELEMENT 
      spectral_element/cg__element_operators__1d.f
      spectral_element/dg__element_operators__1d.f
      spectral_element/embedded_interpolation__1d.f
      spectral_element/embedded_interpolation__3d.f
      spectral_element/coarse_to_fine_interpolation__1d.f
      spectral_element/fine_to_coarse_projection__1d.f
      spectral_element/projection_operator__1d.f
      spectral_element/projection_operator__3d.f
      spectral_element/standard_operators__1d.f

      spectral_element/parent_to_child_interpolation__1d.f # obsolete

    )

#===============================================================================
# time_integration

set ( TIME_INTEGRATION 
      time_integration/imex_runge_kutta_method.f
      time_integration/spectral_deferred_correction.f
    )

#===============================================================================
# mesh

file(GLOB MESH mesh/*.f)

#===============================================================================
# mesh_3d

set ( MESH_3D 
      mesh_3d/generic_mesh__3d.f
      mesh_3d/generic_mesh__3d__mp_createannulus.f
      mesh_3d/generic_mesh__3d__mp_createcylinder.f
      mesh_3d/generic_mesh__3d__mp_creatediamonds.f
      mesh_3d/generic_mesh__3d__mp_createonerotated.f
      mesh_3d/mesh_boundary__3d.f
      mesh_3d/mesh_boundary__3d__sm_hdf5.f
      mesh_3d/mesh_element__3d.F
      mesh_3d/mesh_element__3d__sm_hdf5.f
      mesh_3d/mesh_element__3d__sm_mpi.f
      mesh_3d/mesh_element_indexing__3d.f
      mesh_3d/mesh_face__3d.f
      mesh_3d/mesh_link__3d.f
      mesh_3d/mesh_map_to_child__3d.f
      mesh_3d/mesh_map_to_parent__3d.f
      mesh_3d/mesh__3d.f
      mesh_3d/mesh__3d__mp_buildboundaryfaces.f
      mesh_3d/mesh__3d__mp_buildcommunicator.f
      mesh_3d/mesh__3d__mp_buildcuboids.f
      mesh_3d/mesh__3d__mp_buildfaces.f
      mesh_3d/mesh__3d__mp_buildghosts.f
      mesh_3d/mesh__3d__mp_buildlinks.f
      mesh_3d/mesh__3d__mp_buildmaptochild.f
      mesh_3d/mesh__3d__mp_buildmaptoparent.f
      mesh_3d/mesh__3d__mp_getcuboids.f
      mesh_3d/mesh__3d__mp_getpoints.f
      mesh_3d/mesh__3d__mp_getpointvalency.f
      mesh_3d/mesh__3d__mp_identifyclusters.f
      mesh_3d/mesh__3d__mp_identifyedges.f
      mesh_3d/mesh__3d__mp_identifyelementneighbors.f
      mesh_3d/mesh__3d__mp_identifyranks.f
      mesh_3d/mesh__3d__mp_identifyvertices.f
      mesh_3d/mesh__3d__mp_importgenericmesh.f
      mesh_3d/mesh__3d__mp_writehdf5
      mesh_3d/mesh__3d__sm_hdf5.f
      mesh_3d/structured_mesh_indexing__3d.f
      mesh_3d/generate_regular_mesh__3d.f
      mesh_3d/mesh_metrics__3d.f
      mesh_3d/spectral_element_mesh__3d.f
      mesh_3d/boundary_variable__3d.f
    )

#===============================================================================
# util_1d

set( UTIL_1D 
     util_1d/cg__condensed_solver__1d.f
     util_1d/cg__utilities__1d.f
     util_1d/dg__utilities__1d.f
     util_1d/export/export_vtk_spacetime_data__1d.f
    )

#===============================================================================
# util_3d

set( UTIL_3D util_3d/verification/verify_mesh__3d.f )

set( UTIL_3D ${UTIL_3D}
     util_3d/export/export_vtk_volume_data__3d.f
     util_3d/export/export_vtk_schwarz_domain__3d.f
    )

set( UTIL_3D ${UTIL_3D}
     util_3d/meshgen/create_annulus.f
     util_3d/meshgen/create_cuboid_cartesian.f
     util_3d/meshgen/create_cuboid_diamonds.f
     util_3d/meshgen/create_cuboid_onerotated.f
     util_3d/meshgen/create_cylinder.f
    )

set( UTIL_3D ${UTIL_3D}
     util_3d/import/import_gmsh__3d.f
     util_3d/import/import_gmsh__3d__mp_mapshellnodes.f
   )

set( UTIL_3D ${UTIL_3D}
     util_3d/transfer/data_exchange__3d.f
     util_3d/transfer/element_transfer_buffer__3d.f
     util_3d/transfer/element_transfer_buffer__3d__sm_integer_idk.F
     util_3d/transfer/element_transfer_buffer__3d__sm_integer_ixs.F
     util_3d/transfer/element_transfer_buffer__3d__sm_integer_ixl.F
     util_3d/transfer/element_transfer_buffer__3d__sm_real_rdp.F
     util_3d/transfer/element_transfer_buffer__3d__sm_real_rsp.F
     util_3d/transfer/element_face_transfer_buffer__3d.f
     util_3d/transfer/element_face_transfer_buffer__3d__sm_real_rdp.F
    )

set( UTIL_3D ${UTIL_3D}
     util_3d/ops/assembly__3d.f
     util_3d/ops/assembly__3d__sm_0u.f
     util_3d/ops/trace_operators__3d.f
     util_3d/ops/volume_integrals__3d.f
     util_3d/ops/surface_integrals__3d.f
    )

set( UTIL_3D ${UTIL_3D}
     util_3d/adaptation/child_distribution_map__3d.f
     util_3d/adaptation/child_mesh_adaptation__3d.F
     util_3d/adaptation/child_mesh_adaptation__3d__mp_buildchilddata.f
     util_3d/adaptation/child_mesh_adaptation__3d__mp_buildchildmap.f
     util_3d/adaptation/child_mesh_adaptation__3d__mp_buildconnections.f
     util_3d/adaptation/element_distribution_map__3d.f
     util_3d/adaptation/element_distribution_data__3d.f
     util_3d/adaptation/globalize_adaptation_pattern__3d.f
     util_3d/adaptation/partitioner_interface__3d.f
     util_3d/adaptation/process_adaptation_pattern__3d.f
     util_3d/adaptation/restrict_adaptation_pattern__3d.f
     util_3d/adaptation/root_mesh_partitioning__3d.F
    )

#===============================================================================
# elliptic_1d

set( ELLIPTIC_1D 
     elliptic_1d/dg__elliptic_operator__1d.f
     elliptic_1d/dg__elliptic_operator__1d__mp_cg_method_rx.f
     elliptic_1d/dg__elliptic_operator__1d__mp_eval_rc.f
     elliptic_1d/dg__elliptic_operator__1d__mp_eval_rv.f
     elliptic_1d/dg__elliptic_operator__1d__mp_schwarz_method_rx.f
     elliptic_1d/dg__elliptic_operator__1d__mp_schwarzpcg_method_rx.f
     elliptic_1d/dg__elliptic_operator__1d__mp_hybrid_solver_rc.f
     elliptic_1d/dg__schwarz_operator__1d.f
   )

#===============================================================================
# elliptic_3d

set( ELLIPTIC_3D 
     elliptic_3d/dg__elliptic_operator__3d.f
     elliptic_3d/dg__elliptic_operator__3d__mp_eval_rc.f
     elliptic_3d/dg__elliptic_operator__3d__mp_eval_rv.f
     elliptic_3d/dg__elliptic_operator__3d__mp_eval_dc.F
     elliptic_3d/dg__elliptic_operator__3d__mp_addfluxes_dc_gen.f
     elliptic_3d/dg__elliptic_operator__3d__mp_addfluxes_dc_xsmm.F
     elliptic_3d/dg__elliptic_operator__3d__mp_cg_method_x.f
     elliptic_3d/dg__elliptic_operator__3d__mp_schwarz_method_x.f
     elliptic_3d/dg__elliptic_operator__3d__mp_schwarzpcg_method_x.f
     elliptic_3d/dg__schwarz_operator__3d.f
     elliptic_3d/dg__schwarz_operator__3d__rdp__mp_merge_core.F
     elliptic_3d/dg__schwarz_operator__3d__rdp__mp_merge_structured.F
     elliptic_3d/dg__schwarz_operator__3d__rdp__mp_merge_unstructured.F
     elliptic_3d/dg__schwarz_operator__3d__rdp__mp_restrict_structured.F
     elliptic_3d/dg__schwarz_operator__3d__rdp__mp_restrict_unstructured.F
     elliptic_3d/dg__schwarz_operator__3d__rsp__mp_merge_core.F
     elliptic_3d/dg__schwarz_operator__3d__rsp__mp_merge_structured.F
     elliptic_3d/dg__schwarz_operator__3d__rsp__mp_merge_unstructured.F
     elliptic_3d/dg__schwarz_operator__3d__rsp__mp_restrict_structured.F
     elliptic_3d/dg__schwarz_operator__3d__rsp__mp_restrict_unstructured.F
   )

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
             ${UTIL_1D}
             ${UTIL_3D}
             ${VARIABLES_3D}
             ${ELLIPTIC_1D}
             ${ELLIPTIC_3D}
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
            ${HDF5_Fortran_INCLUDE_DIRS}
            ${ParMETIS_INCLUDE_DIRS}
            ${FFTW_INCLUDE_DIRS}
            )

target_link_libraries(HiBase
    PUBLIC  ${LAPACK_LIBRARIES}
            ${VTK_LIBRARIES}
            ${HDF5_Fortran_LIBRARIES}
            ${FFTW_LIBRARIES}
            ${ParMETIS_LIBRARIES}
            ${LIBXSMM_LIBRARIES}
            )
