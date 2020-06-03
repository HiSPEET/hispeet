set( TPO__ELLIPTIC_3D_DIR tensor_product/tpo__elliptic_3d )

set(TENSOR_PRODUCT ${TENSOR_PRODUCT}
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d.f
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlci.F
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlci__gen.f
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlci__hand.F
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlci__xsmm.F
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlvi.F
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlvi__gen.f
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlvi__hand.F
                   ${TPO__ELLIPTIC_3D_DIR}/tpo__elliptic_3d_rlvi__xsmm.F
                   )
