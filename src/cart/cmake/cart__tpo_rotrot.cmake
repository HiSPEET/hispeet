# cart__tpo_rotrot .............................................................

set(CART__TPO_ROTROT_SRC cart__tpo_rotrot.f cart__tpo_rotrot__var.F)

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_rotrot)

# target for updating the cart__tpo_rotrot configuration
add_custom_target(config__cart__tpo_rotrot ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_RotRot"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_rotrot.var
    COMMENT "Updating cart__tpo_rotrot.var"
    )

# command for generating the cart__tpo_rotrot sources
add_custom_command(OUTPUT  ${CART__TPO_ROTROT_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_rotrot.py
                                        --tmpl ${TMPL}
                                        --dest ${DEST}
                                        --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_rotrot.var
    DEPENDS config__cart__tpo_rotrot
    )

