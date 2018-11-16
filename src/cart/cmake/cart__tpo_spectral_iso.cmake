# cart__tpo_spectral_iso .......................................................

set(CART__TPO_SPECTRAL_ISO_SRC cart__tpo_spectral_iso.f
                               cart__tpo_spectral_iso__var.F)

list(APPEND TENSOR_PRODUCT ${CART__TPO_SPECTRAL_ISO_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_spectral_iso)

# target for updating the cart__tpo_spectral_iso configuration
add_custom_target(config__cart__tpo_spectral_iso ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Spectral_Iso"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_spectral_iso.var
    COMMENT "Updating cart__tpo_spectral_iso.var"
    )

# command for generating the cart__tpo_spectral_iso sources
add_custom_command(OUTPUT  ${CART__TPO_SPECTRAL_ISO_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_spectral_iso.py
                                              --tmpl ${TMPL}
                                              --dest ${DEST}
                                              --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_spectral_iso.var
    DEPENDS config__cart__tpo_spectral_iso
    )

