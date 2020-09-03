set( TPO__SCHWARZ_3D_DIR tensor_product/tpo__schwarz_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d.f
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ci.F
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ci__gen.f
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ci__hand.F
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ci__xsmm.F
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ca.F
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ca__gen.f
                   ${TPO__SCHWARZ_3D_DIR}/tpo__schwarz_3d_ca__xsmm.F
                   )
