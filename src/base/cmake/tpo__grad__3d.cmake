set( TPO__GRAD__3D_DIR tensor_product/tpo__grad__3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__gen.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__gen_rdp.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__gen.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__gen_rdp.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__hand.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__hand_rdp.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__hand.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__hand_rdp.F
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__xsmm.f
                   ${TPO__GRAD__3D_DIR}/tpo__grad__3d_r__xsmm_rdp.F
		               ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__xsmm.f
		               ${TPO__GRAD__3D_DIR}/tpo__grad__3d_d__xsmm_rdp.F
		               )
