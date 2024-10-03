set(TPO__STRESS__3D_DIR tensor_product/tpo__stress__3d)

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d.f

                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlci.F
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlci__gen.f
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlci__gen_rdp.F
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlci__xsmm.f
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlci__xsmm_rdp.F

                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlvi.F
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlvi__gen.f
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlvi__gen_rdp.F
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlvi__xsmm.f
                   ${TPO__STRESS__3D_DIR}/tpo__stress__3d_dlvi__xsmm_rdp.F

                   )
