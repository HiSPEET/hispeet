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

!> summary:  Sigma subgrid-scale model for incompressible flow
!> author:   Yang Liu, Martina Grotteschi, Joerg Stiller
!> date:     2020/11/06
!>
!> This module computes the turbulent viscosity of a single point according to
!> the sigma model from Nicoud et al. (2011). Input parameters are the model
!> constant `c_sigma`, subgrid characteristic length `delta` and the velocity
!> gradient matrix `grad_v`. The algorithm derives from Hasan et al. (2001).
!===============================================================================

module ISP_Flow_SGS__Sigma_Model

  use Kind_Parameters, only: RNP
  use Constants      , only: ZERO, ONE, TWO, HALF, THIRD, PI
  implicit none
  private

  public :: SigmaModel
  public :: DSigma

contains

  !-----------------------------------------------------------------------------
  !> Evaluation on the sigma model for given velocity gradient

  pure subroutine SigmaModel(c_sigma, delta, grad_v, nu_t)
    real(RNP), intent(in)  :: c_sigma     !< model constant
    real(RNP), intent(in)  :: delta       !< subgrid characteristic length
    real(RNP), intent(in)  :: grad_v(3,3) !< velocity gradient,
                                          !! grad_v(i,j) = ∂v_j /∂ x_i
    real(RNP), intent(out) :: nu_t        !< SGS turbulence viscosity

    real(RNP) :: d_sigma

    call DSigma(grad_v,d_sigma)

    nu_t = sign(ONE, c_sigma) * (c_sigma * delta)**2 * d_sigma ! keep sign

  end subroutine SigmaModel

  !-----------------------------------------------------------------------------
  !> Evaluation on the differential operator d_sigma

  pure subroutine DSigma(grad_v, d_sigma)
    real(RNP), intent(in)  :: grad_v(3,3) !< velocity gradient,
                                          !! grad_v(i,j) = ∂v_j /∂ x_i
    real(RNP), intent(out) :: d_sigma     !< differential operator

    real(RNP) :: gg(3,3)        ! gg = grad_vᵀ·grad_v
    real(RNP) :: gg_square(3,3) ! gg_square = gg·gg

    ! the following are intermediate variables to calculate sigular values
    real(RNP) :: alpha_1, alpha_2, alpha_3
    real(RNP) :: sigma_1, sigma_2, sigma_3
    real(RNP) :: i_1, i_2, i_3
    real(RNP) :: a
    integer   :: i, j, k

    ! calculate gg and gg_square matrix
    gg = 0
    do j = 1,3
    do i = 1,3
    do k = 1,3
      gg(i,j) = gg(i,j) + grad_v(k,i)*grad_v(k,j)
    end do
    end do
    end do

    ! calculate gg and gg_square matrix
    gg_square = 0
    do j = 1,3
    do i = 1,3
    do k = 1,3
      gg_square(i,j) = gg_square(i,j) + gg(i,k)*gg(k,j)
    end do
    end do
    end do

    ! calculate the singular values sigma
    i_1 = Tr(gg)
    i_2 = HALF * (Tr(gg)**2 - Tr(gg_square))
    i_3 = Det(gg)

    alpha_1 = i_1**2 / 9 - i_2*THIRD

    if (alpha_1 <= 0) then

      ! By alpha_1 <= 0 there exists no meaningful sigular values
      d_sigma = 0

    else

      alpha_2 = i_1**3 / 27 - i_1*i_2/6 + i_3*HALF

      a = alpha_2 / sqrt(alpha_1**3)      ! compute the argument of alpha_3
      a = sign(ONE, a) * min(abs(a), ONE) ! avoide invalid input for acos function
      alpha_3 = acos(a) * THIRD

      sigma_1 = sqrt( max(i_1*THIRD + TWO*sqrt(alpha_1)*cos(alpha_3),ZERO) )
      sigma_2 = sqrt( max(i_1*THIRD - TWO*sqrt(alpha_1)*cos(PI*THIRD+alpha_3),ZERO) )
      sigma_3 = sqrt( max(i_1*THIRD - TWO*sqrt(alpha_1)*cos(PI*THIRD-alpha_3),ZERO) )
      ! Here use the max function to avoid invalid input for sqrt function

      d_sigma = sigma_3 * (sigma_1 - sigma_2) * (sigma_2 - sigma_3) / sigma_1**2

    end if

  end subroutine DSigma

  !-----------------------------------------------------------------------------
  !> Trace of a 3*3 a

  pure real(RNP) function Tr(a) result(trace)
    real(RNP),intent(in) :: a(3,3)

    trace = a(1,1) + a(2,2) + a(3,3)

  end function Tr

  !-----------------------------------------------------------------------------
  !> Determinant of a 3*3 a

  pure real(RNP) function Det(a) result(d)
    real(RNP),intent(in) :: a(3,3)

    d = a(1,1) * a(2,2) * a(3,3)  &
      + a(1,2) * a(2,3) * a(3,1)  &
      + a(1,3) * a(2,1) * a(3,2)  &
      - a(3,1) * a(2,2) * a(1,3)  &
      - a(1,1) * a(2,3) * a(3,2)  &
      - a(1,2) * a(2,1) * a(3,3)

  end function Det

  !=============================================================================

end module ISP_Flow_SGS__Sigma_Model
