# cart__tpo_grad ...............................................................

set(CART__TPO_GRAD_SRC cart__tpo_grad.f cart__tpo_grad__var.F)

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_grad)

# target for updating the cart__tpo_grad configuration
add_custom_target(config__cart__tpo_grad ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Grad"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_grad.var
    COMMENT "Updating cart__tpo_grad.var"
    )

# command for generating the cart__tpo_grad sources
add_custom_command(OUTPUT  ${CART__TPO_GRAD_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_grad.py
                                      --tmpl ${TMPL}
                                      --dest ${DEST}
                                      --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_grad.var
    DEPENDS config__cart__tpo_grad
    )

