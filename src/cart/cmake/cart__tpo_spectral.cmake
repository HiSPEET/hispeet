# cart__tpo_spectral ...........................................................

set(CART__TPO_SPECTRAL_SRC cart__tpo_spectral.f
                           cart__tpo_spectral__gen.f)

if(DEFINED ENV{TUNING_FLAG})
  list(APPEND CART__TPO_SPECTRAL_SRC cart__tpo_spectral__par.F)
endif()

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_spectral)

# target for updating the cart__tpo_spectral configuration
add_custom_target(config__cart__tpo_spectral ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Spectral"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_spectral.var
    COMMENT "Updating cart__tpo_spectral.var"
    )

# command for generating the cart__tpo_spectral sources
add_custom_command(OUTPUT  ${CART__TPO_SPECTRAL_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_spectral.py
                                          --tmpl ${TMPL}
                                          --dest ${DEST}
                                          --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_spectral.var
    DEPENDS config__cart__tpo_spectral
    )

