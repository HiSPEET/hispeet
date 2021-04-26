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
             )

#            elliptic/cart__elliptic_operator_cg.f

#===============================================================================
# isp_flow

set(ISP_FLOW isp_flow/cart__isp_flow__boundary_values.F
             isp_flow/cart__isp_flow__convection__eq.F
             isp_flow/cart__isp_flow__convection__iq.F
             isp_flow/cart__isp_flow__convection.F
             isp_flow/cart__isp_flow__diffusion.F
             isp_flow/cart__isp_flow__operators.F
             isp_flow/cart__isp_flow__pressure__eo.F
             isp_flow/cart__isp_flow__pressure__mo.F
             isp_flow/cart__isp_flow__pressure.F
             isp_flow/cart__isp_flow__projection.F
             isp_flow/cart__isp_flow__time_derivative.F
             isp_flow/cart__isp_flow__time_integrator.F
             isp_flow/cart__isp_flow__time_integrator__euler.F
             isp_flow/cart__isp_flow__time_integrator__bdf2.F
             isp_flow/cart__isp_flow__time_integrator__tr.F
             isp_flow/cart__isp_flow__time_integrator__runge_kutta.F
             isp_flow/cart__isp_flow__euler.F
             isp_flow/cart__isp_flow__sdc.F
             isp_flow/cart__isp_flow__sdc_method.F
             isp_flow/cart__isp_flow__sdc_method__euler.F
             )

#===============================================================================
# cart library

list(APPEND CART_SRC ${SPECTRAL_ELEMENT}
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
