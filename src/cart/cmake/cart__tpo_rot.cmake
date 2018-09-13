# cart__tpo_rot ................................................................

set(CART__TPO_ROT_SRC cart__tpo_rot.f cart__tpo_rot__gen.f)

if(DEFINED ENV{TUNING_FLAG})
  list(APPEND CART__TPO_ROT_SRC cart__tpo_rot__par.F)
endif()

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/cart__tpo_rot)

# target for updating the cart__tpo_rot configuration
add_custom_target(config__cart__tpo_rot ALL
    COMMAND ${TPO_CONFIG} --op "CART__TPO_Rot"
                          --tmpl ${TMPL}
                          --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/cart__tpo_rot.var
    COMMENT "Updating cart__tpo_rot.var"
    )

# command for generating the cart__tpo_rot sources
add_custom_command(OUTPUT  ${CART__TPO_ROT_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/cart__tpo_rot.py
                                     --tmpl ${TMPL}
                                     --dest ${DEST}
                                     --shared ${UNIFORM}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/cart__tpo_rot.var
    DEPENDS config__cart__tpo_rot
    )

