set( TPO__SPECTRAL__3D_DIR tensor_product/tpo__spectral__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d.f
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ci.F
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ci__gen.f
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ci__hand.F
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ci__xsmm.F
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ca.F
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ca__gen.f
                   ${TPO__SPECTRAL__3D_DIR}/tpo__spectral__3d_ca__xsmm.F
                   )
