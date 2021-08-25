set( TPO__DIV__3D_DIR tensor_product/tpo__div__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r.F
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__gen.f
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__gen_rdp.F
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_d__gen.f
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_d__gen_rdp.F
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__hand.f
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__hand_rdp.F
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__xsmm.f
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_r__xsmm_rdp.F
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_d__xsmm.f
                   ${TPO__DIV__3D_DIR}/tpo__div__3d_d__xsmm_rdp.F
                   )
