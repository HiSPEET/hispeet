set( TPO__DIV_3D_DIR tensor_product/tpo__div_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__DIV_3D_DIR}/tpo__div_3d_r.F
                   ${TPO__DIV_3D_DIR}/tpo__div_3d_r__gen.f
                   ${TPO__DIV_3D_DIR}/tpo__div_3d_r__hand.F
                   ${TPO__DIV_3D_DIR}/tpo__div_3d_r__xsmm.F
                   )
