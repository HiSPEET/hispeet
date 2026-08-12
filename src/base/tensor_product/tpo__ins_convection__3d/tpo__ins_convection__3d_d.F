!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  3D curvilinear convection operator for INS
!> author:   Jerome Michel, Jörg Stiller
!> date:     2022/01/21 - 2026/03/08
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

  subroutine TPO_INS_Convection_D_RDP( nv, nq, ne, bc, D_v, I_vq, w_q   &
                                     , Jd_q, Ji_q, a_q, n_q, v, vp, F_c &
                                     , form, e_max )

    integer, intent(in) :: nv !< number of velocity points per element direction
    integer, intent(in) :: nq !< number of quadrature points per direction
    integer, intent(in) :: ne !< number of elements

    character, intent(in) :: bc(6,ne) !< element face BC

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

    character(4), optional, intent(in) :: form  !< {flux,skew,conv} [flux]
    integer,      optional, intent(in) :: e_max !< last element to evaluate [ne]

    character(4) :: form_
    real(RDP)    :: alpha

    if (present(form)) then
      form_ = form
    else
      form_ = 'flux'
    end if

    select case(form_)
    case('skew')
      alpha = 0.5
    case('conv')
      alpha = 1.0
    case default
      alpha = 0.0
    end select

#ifdef __LIBXSMM__
    call TPO_INS_Convection_D_XSMM( nv, nq, ne, bc, alpha, D_v, I_vq, w_q   &
                                  , Jd_q, Ji_q, a_q, n_q, v, vp, F_c, e_max )
#else
    call TPO_INS_Convection_D_Gen( nv, nq, ne, bc, alpha, D_v, I_vq, w_q   &
                                 , Jd_q, Ji_q, a_q, n_q, v, vp, F_c, e_max )
#endif

  end subroutine TPO_INS_Convection_D_RDP

  !=============================================================================

end module TPO__INS_Convection__3D_D
