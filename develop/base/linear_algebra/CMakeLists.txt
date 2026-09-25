#===============================================================================
# test_scalar_product

add_executable( test_scalar_product
                test_scalar_product.f )

target_link_libraries(test_scalar_product HiBase)

if (Slurm_FOUND)
    configure_file( test_scalar_product.slurm
                    test_scalar_product.slurm COPYONLY )
endif (Slurm_FOUND)
