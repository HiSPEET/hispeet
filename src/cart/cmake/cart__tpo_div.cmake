# cart__tpo_div ................................................................

set(CART__TPO_DIV_SRC cart__tpo_div.f cart__tpo_div__var.F)

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_div)

# target for updating the cart__tpo_div configuration
add_custom_target(config__cart__tpo_div ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Div"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_div.var
    COMMENT "Updating cart__tpo_div.var"
    )

# command for generating the cart__tpo_div sources
add_custom_command(OUTPUT  ${CART__TPO_DIV_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_div.py
                                     --tmpl ${TMPL}
                                     --dest ${DEST}
                                     --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_div.var
    DEPENDS config__cart__tpo_div
    )

