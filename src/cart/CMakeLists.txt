#===============================================================================
# tensor_product

# configuration script
set(TPO_CONFIG ${PROJECT_SOURCE_DIR}/bin/tpo_configure.py)

# shared directory for non/uniform operators
set(UNIFORM    ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/shared/uniform)
set(NONUNIFORM ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/shared/nonuniform)

# target directory of generated sources
set(DEST ${CMAKE_CURRENT_BINARY_DIR})

# operators, adding themselves to TENSOR_PRODUCT
include(cmake/cart__tpo_grad.cmake)
include(cmake/cart__tpo_div.cmake)
include(cmake/cart__tpo_rot.cmake)
include(cmake/cart__tpo_rotrot.cmake)
include(cmake/cart__tpo_elliptic_ci.cmake)
include(cmake/cart__tpo_elliptic_vi.cmake)
include(cmake/cart__tpo_spectral.cmake)
include(cmake/cart__tpo_spectral_iso.cmake)
include(cmake/cart__tpo_schwarz.cmake)
include(cmake/cart__tpo_schwarz_iso.cmake)
include(cmake/cart__tpo_schwarz_old.cmake)

if (DEFINED Fortran_FLAGS_TPO)
    set_source_files_properties(${TENSOR_PRODUCT}
        PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_TPO})
endif ()

#===============================================================================
# mesh

set(MESH mesh/cart__mesh_face.f
         mesh/cart__mesh_element.f
         mesh/cart__mesh_boundary.f
         mesh/cart__mesh_link.f
         mesh/cart__mesh_partition.f
         mesh/cart__mesh_partition__mp_getpoints.f
         mesh/cart__mesh_partition__mp_getpointvalency.f
         mesh/cart__element_connectivity.f
         mesh/cart__generate_mesh_links.f
         mesh/cart__generate_structured_mesh.f
         )

#===============================================================================
# variables

set(VARIABLES variables/cart__boundary_variable.f)

#===============================================================================
# intra_level

set(INTRA_LEVEL intra_level/cart__face_transfer_buffer.f
                intra_level/cart__trace_transfer_buffer.f
                intra_level/cart__element_transfer_buffer.f
                )

#===============================================================================
# inter_level

set(INTER_LEVEL inter_level/cart__dg_pmg_transfer_operators.f)

#===============================================================================
# spectral_element

file(GLOB SPECTRAL_ELEMENT spectral_element/*.f)

#===============================================================================
# operators

set(OPERATORS operators/cart__assembly_operator.f
              operators/cart__trace_operator.f
              operators/cart__normal_trace_operator.f
              operators/cart__dg_weak_gradient.f
              operators/cart__dg_weak_divergence.f
              )

#===============================================================================
# elliptic

set(ELLIPTIC elliptic/cart__schwarz_operator.f
             elliptic/cart__schwarz_operator__mp_buildeigensystems_ip.f
             elliptic/cart__schwarz_operator__mp_buildsubdomains_ci.f
             elliptic/cart__schwarz_operator__mp_buildsubdomains_vi.f
             elliptic/cart__elliptic_operator.f
             elliptic/cart__elliptic_operator_ip.f
             elliptic/cart__elliptic_operator_ip__mp_setproblem.f
             elliptic/cart__elliptic_operator_ip__mp_apply.f
             elliptic/cart__elliptic_operator_ip__mp_apply_ci.f
             elliptic/cart__elliptic_operator_ip__mp_apply_vi.f
             elliptic/cart__elliptic_operator_ip__mp_bctorhs.f
             elliptic/cart__elliptic_operator_ip__mp_residual.f
             elliptic/cart__elliptic_operator_ip__mp_schwarzmethod.f
             elliptic/cart__elliptic_pmg_level.f
             elliptic/cart__elliptic_pmg.f
             elliptic/cart__elliptic_pmg__mp_init.f
             elliptic/cart__elliptic_pmg__mp_setproblem.f
             elliptic/cart__cg_elliptic_ci_bc.f
             elliptic/cart__dg_elliptic_cis_bc.f
             elliptic/cart__dg_elliptic_ciu_bc.f
             elliptic/cart__dg_elliptic_cis_operator.f
             elliptic/cart__dg_elliptic_ciu_operator.f
             elliptic/cart__dg_elliptic_ci_operator.f
             elliptic/cart__dg_elliptic_ci_residual.f
             elliptic/cart__dg_elliptic_ci_conj_grad.f
             elliptic/cart__dg_elliptic_ci_schwarz.f
             elliptic/cart__dg_elliptic_ci_pmg_level.f
             elliptic/cart__dg_elliptic_ci_pmg.f
             )

#            elliptic/cart__elliptic_operator_cg.f


#===============================================================================
# isp_flow

set(ISP_FLOW isp_flow/cart__isp_flow__boundary_values.f
             isp_flow/cart__isp_flow__convection__eq.f
             isp_flow/cart__isp_flow__convection__iq.f
             isp_flow/cart__isp_flow__convection.f
             isp_flow/cart__isp_flow__diffusion.f
             isp_flow/cart__isp_flow__operators.f
             isp_flow/cart__isp_flow__pressure__eo.f
             isp_flow/cart__isp_flow__pressure__mo.f
             isp_flow/cart__isp_flow__pressure.f
             isp_flow/cart__isp_flow__projection.f
             isp_flow/cart__isp_flow__time_derivative.f
             isp_flow/cart__isp_flow__time_integrator.f
             isp_flow/cart__isp_flow__time_integrator__euler.f
             isp_flow/cart__isp_flow__time_integrator__runge_kutta.f
             isp_flow/cart__isp_flow__euler.f
             isp_flow/cart__isp_flow__sdc.f
             isp_flow/cart__isp_flow__sdc_method.f
             isp_flow/cart__isp_flow__sdc_corrector.f
             isp_flow/cart__isp_flow__sdc_corrector__euler.f
             )

#===============================================================================
# cart library

list(APPEND CART_SRC ${TENSOR_PRODUCT}
                     ${SPECTRAL_ELEMENT}
                     ${MESH}
                     ${VARIABLES}
                     ${INTRA_LEVEL}
                     ${INTER_LEVEL}
                     ${OPERATORS}
                     ${ELLIPTIC}
                     ${ISP_FLOW}
                     )

add_library(HiCart ${CART_SRC})

set_property(TARGET HiCart
             PROPERTY Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

target_link_libraries(HiCart PUBLIC HiBase HiPhysics)
