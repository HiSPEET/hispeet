set( TPO__ROT__3D_DIR tensor_product/tpo__rot__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r.F
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__gen.f
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__gen_rdp.F
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__hand.f
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__hand_rdp.F
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__xsmm.f
                   ${TPO__ROT__3D_DIR}/tpo__rot__3d_r__xsmm_rdp.F
                   )
