# cart__tpo_diffusion ..........................................................

set(CART__TPO_DIFFUSION_SRC cart__tpo_diffusion.f
                            cart__tpo_diffusion__gen.f
                            cart__tpo_diffusion__gen_2.f
                            cart__tpo_diffusion__gen_3.f
                            cart__tpo_diffusion__gen_4.f )

if(DEFINED ENV{TUNING_FLAG})
  list(APPEND CART__TPO_DIFFUSION_SRC cart__tpo_diffusion__par.F)
endif()

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_diffusion)

# target for updating the cart__tpo_diffusion configuration
add_custom_target(config__cart__tpo_diffusion ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Diffusion"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_diffusion.var
    COMMENT "Updating cart__tpo_diffusion.var"
    )

# command for generating the cart__tpo_diffusion sources
add_custom_command(OUTPUT  ${CART__TPO_DIFFUSION_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_diffusion.py
                                     --tmpl ${TMPL}
                                     --dest ${DEST}
                                     --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_diffusion.var
    DEPENDS config__cart__tpo_diffusion
    )

