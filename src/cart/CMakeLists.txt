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
include(cmake/cart__tpo_diffusion.cmake)
include(cmake/cart__tpo_spectral.cmake)
include(cmake/cart__tpo_spectral_iso.cmake)
include(cmake/cart__tpo_schwarz_iso.cmake)

list(APPEND TENSOR_PRODUCT ${CART__TPO_DIV_SRC}
                           ${CART__TPO_GRAD_SRC}
                           ${CART__TPO_ROT_SRC}
                           ${CART__TPO_ROTROT_SRC}
                           ${CART__TPO_DIFFUSION_SRC}
                           ${CART__TPO_SPECTRAL_SRC}
                           ${CART__TPO_SPECTRAL_ISO_SRC}
                           ${CART__TPO_SCHWARZ_ISO_SRC}
                           )

if (DEFINED Fortran_FLAGS_TPO)
    set_source_files_properties(${TENSOR_PRODUCT}
        PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_TPO})
endif ()

#===============================================================================
# spectral_element

file(GLOB SPECTRAL_ELEMENT spectral_element/*.f)

#===============================================================================
# mesh

set(MESH mesh/cart__mesh_face.f
         mesh/cart__mesh_element.f
         mesh/cart__mesh_boundary.f
         mesh/cart__mesh_link.f
         mesh/cart__mesh_partition.f
         mesh/cart__mesh_partition__get_points.f
         mesh/cart__element_connectivity.f
         mesh/cart__generate_mesh_links.f
         mesh/cart__generate_structured_mesh.f
         )

#===============================================================================
# variables

set(VARIABLES variables/cart__boundary_variable.f)

#===============================================================================
# intra_level_transfer

set(INTRA_LEVEL intra_level_transfer/cart__face_transfer_buffer.f
                intra_level_transfer/cart__trace_transfer_buffer.f
                intra_level_transfer/cart__element_transfer_buffer.f
                )

#===============================================================================
# inter_level_transfer

set(INTER_LEVEL inter_level_transfer/cart__dg_pmg_transfer_operators.f)

#===============================================================================
# operators

set(OPERATORS operators/cart__trace_operator.f
              operators/cart__dg_weak_gradient.f
              operators/cart__dg_weak_divergence.f
              )

#===============================================================================
# diffusion

set(DIFFUSION diffusion/cart__dg_diffusion_cis_bc.f
              diffusion/cart__dg_diffusion_cis_operator.f
              diffusion/cart__dg_diffusion_ci_operator.f
              diffusion/cart__dg_diffusion_ci_residual.f
              diffusion/cart__dg_diffusion_ci_conj_grad.f
              diffusion/cart__dg_diffusion_ci_schwarz.f
              diffusion/cart__dg_diffusion_ci_pmg_level.f
              diffusion/cart__dg_diffusion_ci_pmg.f
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
                     ${DIFFUSION}
                     )

add_library(HiCart ${CART_SRC})

target_link_libraries(HiCart PUBLIC HiBase HiPhysics)
