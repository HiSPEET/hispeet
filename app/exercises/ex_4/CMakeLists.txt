add_executable(exact_integration EXCLUDE_FROM_ALL exact_integration.f)
target_link_libraries(exact_integration HiBase)

add_executable(tensor_products EXCLUDE_FROM_ALL tensor_products.f)
target_link_libraries(tensor_products HiBase)

configure_file(tensor_products.prm tensor_products.prm COPYONLY)

configure_file(exact_integration.py exact_integration.py COPYONLY)

configure_file(tensor_products.py tensor_products.py COPYONLY)

configure_file(README.md README.md COPYONLY)
