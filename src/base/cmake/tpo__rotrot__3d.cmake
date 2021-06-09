set( TPO__ROTROT__3D_DIR tensor_product/tpo__rotrot__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r.F
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__gen.f
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__gen_rdp.F
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__hand.f
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__hand_rdp.F
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__xsmm.f
                   ${TPO__ROTROT__3D_DIR}/tpo__rotrot__3d_r__xsmm_rdp.F
                   )
