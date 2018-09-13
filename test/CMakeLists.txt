if (PFUNIT_FOUND)

    # base .....................................................................

    set(TEST_BASE_SOURCES gauss_jacobi__test.pf)
    add_pfunit_test(unit_test_base "${TEST_BASE_SOURCES}" "" "")

    target_link_libraries(unit_test_base hibase)

    # parmetis .................................................................

    if (ParMETIS_FOUND)

        set(TEST_ParMETIS_SOURCES parmetis_binding__test.pf)
        add_pfunit_test(unit_test_parmetis "${TEST_ParMETIS_SOURCES}" "" "" NPES 3)

        target_compile_options(unit_test_parmetis PRIVATE "-DUSE_MPI")
        target_link_libraries (unit_test_parmetis hibase)

    endif (ParMETIS_FOUND)

endif (PFUNIT_FOUND)
