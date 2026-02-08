!> summary:  3D curvilinear convection operator for INS
!> author:   Jerome Michel
!> date:     2022/01/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module TPO__INS_Convection__3D_D
  use Kind_Parameters, only: RDP
  use TPO__INS_Convection__3D_D__Gen
  use TPO__INS_Convection__3D_D__XSMM
  implicit none
  private

  public :: TPO_INS_Convection

  interface TPO_INS_Convection
    module procedure TPO_INS_Convection_D_RDP
  end interface

contains

  subroutine TPO_INS_Convection_D_RDP &
                 (nv, nq, ne, D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c)

    integer, intent(in) :: nv !< number of velocity points per element direction
    integer, intent(in) :: nq !< number of quadrature points per direction
    integer, intent(in) :: ne !< number of elements

    real(RDP), intent(in) :: D_v (nv,nv) !< 1D diff operator on v-points
    real(RDP), intent(in) :: I_vq(nq,nv) !< 1D interpolation from to q-points
    real(RDP), intent(in) :: w_q (nq)    !< 1D quadrature weights

    real(RDP), intent(in)  :: Jd_q(nq,nq,nq,ne)     !< Jacobian determinant
    real(RDP), intent(in)  :: Ji_q(nq,nq,nq,ne,3,3) !< inverse Jacobian matrix
    real(RDP), intent(in)  :: a_q(nq,nq,6,ne)       !< face area coefficients
    real(RDP), intent(in)  :: n_q(nq,nq,6,ne,3)     !< face unit normal vectors
    real(RDP), intent(in)  :: v(nv,nv,nv,ne,3)      !< velocity
    real(RDP), intent(in)  :: vp(nv,nv,6,ne,3)      !< velocity v⁺ at faces
    real(RDP), intent(out) :: F_c(nv,nv,nv,ne,3)    !< convective flux integrals

    real(RDP) :: alpha = 0.0

#ifdef __LIBXSMM__
      call TPO_INS_Convection_D_XSMM(nv, nq, ne, alpha, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c)
#else
      call TPO_INS_Convection_D_Gen(nv, nq, ne, alpha, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c)
#endif

  end subroutine TPO_INS_Convection_D_RDP

  !=============================================================================

end module TPO__INS_Convection__3D_D
