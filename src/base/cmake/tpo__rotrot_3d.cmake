set( TPO__ROTROT_3D_DIR tensor_product/tpo__rotrot_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__ROTROT_3D_DIR}/tpo__rotrot_3d_r.F
                   ${TPO__ROTROT_3D_DIR}/tpo__rotrot_3d_r__gen.f
                   ${TPO__ROTROT_3D_DIR}/tpo__rotrot_3d_r__hand.F
                   ${TPO__ROTROT_3D_DIR}/tpo__rotrot_3d_r__xsmm.F
                   )
