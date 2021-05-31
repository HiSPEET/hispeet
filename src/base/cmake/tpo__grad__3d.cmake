set( TPO__GRAD__3D_DIR tensor_product/tpo__grad__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__gen.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__hand.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__xsmm.F
                   )
