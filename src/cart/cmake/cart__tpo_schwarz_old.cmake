# cart__tpo_schwarz_old ........................................................

set(CART__TPO_SCHWARZ_OLD_SRC cart__tpo_schwarz_old.f
                              cart__tpo_schwarz_old__var.F)

list(APPEND TENSOR_PRODUCT ${CART__TPO_SCHWARZ_OLD_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_schwarz_old)

# target for updating the cart__tpo_schwarz_old configuration
add_custom_target(config__cart__tpo_schwarz_old ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Schwarz_Old"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_schwarz_old.var
    COMMENT "Updating cart__tpo_schwarz_old.var"
    )

# command for generating the cart__tpo_schwarz_old sources
add_custom_command(OUTPUT  ${CART__TPO_SCHWARZ_OLD_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_schwarz_old.py
                                              --tmpl ${TMPL}
                                              --dest ${DEST}
                                              --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_schwarz_old.var
    DEPENDS config__cart__tpo_schwarz_old
    )

