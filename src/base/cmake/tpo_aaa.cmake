# tpo_aaa ......................................................................

set(TPO_AAA_SRC tpo_aaa.f tpo_aaa__var.F)
list(APPEND TENSOR_PRODUCT ${TPO_AAA_SRC})

# template directory
set(TMPL ${CMAKE_CURRENT_SOURCE_DIR}/tensor_product/tpo_aaa)

# target for updating the tpo_aaa configuration
add_custom_target(config__tpo_aaa ALL
    COMMAND ${TPO_CONFIG} --op "TPO_AAA" --tmpl ${TMPL} --dest ${DEST}
    WORKING_DIRECTORY ${TMPL}
    BYPRODUCTS ${DEST}/tpo_aaa.var
    COMMENT "Updating tpo_aaa.var"
    )

# command for generating the tpo_aaa sources
add_custom_command(OUTPUT  ${TPO_AAA_SRC}
    COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${PYTHONPATH}
            ${PYTHON_EXECUTABLE} ${TMPL}/tpo_aaa.py --tmpl ${TMPL} --dest ${DEST}
    WORKING_DIRECTORY ${TMPL_}
    MAIN_DEPENDENCY ${DEST}/tpo_aaa.var
    DEPENDS config__tpo_aaa
    )
