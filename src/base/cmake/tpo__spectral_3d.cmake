set( TPO__SPECTRAL_3D_DIR tensor_product/tpo__spectral_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__SPECTRAL_3D_DIR}/tpo__spectral_3d.f
                   ${TPO__SPECTRAL_3D_DIR}/tpo__spectral_3d_i.F
                   ${TPO__SPECTRAL_3D_DIR}/tpo__spectral_3d_i__gen.f
                   ${TPO__SPECTRAL_3D_DIR}/tpo__spectral_3d_i__hand.F
                   ${TPO__SPECTRAL_3D_DIR}/tpo__spectral_3d_i__xsmm.F
                   )
