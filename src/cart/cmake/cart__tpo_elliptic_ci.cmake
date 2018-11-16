# cart__tpo_elliptic_ci ........................................................

set(CART__TPO_ELLIPTIC_CI_SRC
    cart__tpo_elliptic_ci.f
    cart__tpo_elliptic_ci__var.F
    )

list(APPEND TENSOR_PRODUCT ${CART__TPO_ELLIPTIC_CI_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_elliptic_ci)

# target for updating the cart__tpo_elliptic_ci configuration
add_custom_target(config__cart__tpo_elliptic_ci ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Elliptic_CI"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_elliptic_ci.var
    COMMENT "Updating cart__tpo_elliptic_ci.var"
    )

# command for generating the cart__tpo_elliptic_ci sources
add_custom_command(OUTPUT  ${CART__TPO_ELLIPTIC_CI_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_elliptic_ci.py
                                     --tmpl ${TMPL}
                                     --dest ${DEST}
                                     --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_elliptic_ci.var
    DEPENDS config__cart__tpo_elliptic_ci
    )

