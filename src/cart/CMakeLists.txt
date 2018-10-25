#===============================================================================
# tensor_product

# configuration script
set(TPO_CONFIG ${PROJECT_SOURCE_DIR}/bin/tpo_configure.py)

# shared directory for uniform operators
set(UNIFORM ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/shared/uniform)

# target directory of generated sources
set(DEST ${CMAKE_CURRENT_BINARY_DIR})

# operators
include(cmake/cart__tpo_grad.cmake)
include(cmake/cart__tpo_div.cmake)
include(cmake/cart__tpo_rot.cmake)
include(cmake/cart__tpo_rotrot.cmake)
include(cmake/cart__tpo_elliptic_ci.cmake)
include(cmake/cart__tpo_spectral.cmake)
include(cmake/cart__tpo_spectral_iso.cmake)
include(cmake/cart__tpo_schwarz_iso.cmake)

list(APPEND TENSOR_PRODUCT ${CART__TPO_DIV_SRC}
                           ${CART__TPO_GRAD_SRC}
                           ${CART__TPO_ROT_SRC}
                           ${CART__TPO_ROTROT_SRC}
                           ${CART__TPO_ELLIPTIC_CI_SRC}
                           ${CART__TPO_SPECTRAL_SRC}
                           ${CART__TPO_SPECTRAL_ISO_SRC}
                           ${CART__TPO_SCHWARZ_ISO_SRC}
                           )

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
         mesh/cart__mesh_partition__get_points.f
         mesh/cart__mesh_partition__get_valency.f
         mesh/cart__element_connectivity.f
         mesh/cart__generate_mesh_links.f
         mesh/cart__generate_structured_mesh.f
         )

#===============================================================================
# spectral_element

file(GLOB SPECTRAL_ELEMENT spectral_element/*.f)

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
# operators

set(OPERATORS operators/cart__trace_operator.f
              operators/cart__normal_trace_operator.f
              operators/cart__dg_weak_gradient.f
              operators/cart__dg_weak_divergence.f
              )

#===============================================================================
# elliptic

set(ELLIPTIC elliptic/cart__cg_elliptic_ci_bc.f
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


#===============================================================================
# isp_flow

set(ISP_FLOW isp_flow/cart__isp_flow__boundary_values.f
             isp_flow/cart__isp_flow__operators.f
             isp_flow/cart__isp_flow__convection__eq.f
             isp_flow/cart__isp_flow__convection__iq.f
             isp_flow/cart__isp_flow__convection.f
             isp_flow/cart__isp_flow__diffusion.f
             isp_flow/cart__isp_flow__time_derivative.f
             isp_flow/cart__isp_flow__pressure__dg_eo.f
             isp_flow/cart__isp_flow__pressure.f
             isp_flow/cart__isp_flow__euler.f
             isp_flow/cart__isp_flow__sdc.f
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
