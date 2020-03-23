set( TPO__ROT_3D_DIR tensor_product/tpo__rot_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__ROT_3D_DIR}/tpo__rot_3d_r.F
                   ${TPO__ROT_3D_DIR}/tpo__rot_3d_r__gen.f
                   ${TPO__ROT_3D_DIR}/tpo__rot_3d_r__hand.F
                   ${TPO__ROT_3D_DIR}/tpo__rot_3d_r__xsmm.F
                   )
