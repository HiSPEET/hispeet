# cart__tpo_schwarz ...........................................................

set(CART__TPO_SCHWARZ_SRC cart__tpo_schwarz.f cart__tpo_schwarz__var.F)

list(APPEND TENSOR_PRODUCT ${CART__TPO_SCHWARZ_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_schwarz)

# target for updating the cart__tpo_schwarz configuration
add_custom_target(config__cart__tpo_schwarz ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Schwarz"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_schwarz.var
    COMMENT "Updating cart__tpo_schwarz.var"
    )

# command for generating the cart__tpo_schwarz sources
add_custom_command(OUTPUT  ${CART__TPO_SCHWARZ_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_schwarz.py
                                          --tmpl ${TMPL}
                                          --dest ${DEST}
                                          --shared ${NONUNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_schwarz.var
    DEPENDS config__cart__tpo_schwarz
    )

