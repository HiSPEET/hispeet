set( TPO__GRAD_3D_DIR tensor_product/tpo__grad_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__GRAD_3D_DIR}/tpo__grad_3d_r.F
                   ${TPO__GRAD_3D_DIR}/tpo__grad_3d_r__gen.f
                   ${TPO__GRAD_3D_DIR}/tpo__grad_3d_r__hand.F
                   ${TPO__GRAD_3D_DIR}/tpo__grad_3d_r__xsmm.F
                   )
