# cart__tpo_schwarz_cubic ......................................................

set(CART__TPO_SCHWARZ_CUBIC_SRC cart__tpo_schwarz_cubic.f
                                 cart__tpo_schwarz_cubic__var.F)

list(APPEND TENSOR_PRODUCT ${CART__TPO_SCHWARZ_CUBIC_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_schwarz_cubic)

# target for updating the cart__tpo_schwarz_cubic configuration
add_custom_target(config__cart__tpo_schwarz_cubic ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Schwarz_Cubic"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_schwarz_cubic.var
    COMMENT "Updating cart__tpo_schwarz_cubic.var"
    )

# command for generating the cart__tpo_schwarz_cubic sources
add_custom_command(OUTPUT  ${CART__TPO_SCHWARZ_CUBIC_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_schwarz_cubic.py
                                              --tmpl ${TMPL}
                                              --dest ${DEST}
                                              --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_schwarz_cubic.var
    DEPENDS config__cart__tpo_schwarz_cubic
    )

