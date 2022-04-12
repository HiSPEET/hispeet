set(TPO__INS_CONVECTION__3D_DIR tensor_product/tpo__ins_convection__3d)

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d.f
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d_d.F
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d_d__gen.f
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d_d__gen_rdp.F
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d_d__xsmm.f
                   ${TPO__INS_CONVECTION__3D_DIR}/tpo__ins_convection__3d_d__xsmm_rdp.F
                   )
