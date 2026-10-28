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

!> summary:  Lambda2 vortex sensor
!> author:   Joerg Stiller
!> date:     2026/10/28
!>
!> Computation of the vortex sensor λ₂ according to Jeong & Hussain (1995).
!> Given the 3D vector field v, λ₂ is the second eigenvalue of the tensor S²+Ω²,
!> where S represents the symmetric part of ∇v, and Ω the antisymmetric part.
!>
!> J. Jeong and F. Hussain. On the Identification of a Vortex.
!> J. Fluid Mechanics, 285:69-94, 1995.
!===============================================================================

module Vortex_Sensors__3D
  use Kind_Parameters!, only: RNP
  use Constants!, only: HALF
  use Vector_Algebra__3D
  implicit none
  private

  public :: Lambda2

contains

  !-----------------------------------------------------------------------------
  !> Evaluates the λ₂ vortex sensor based on the velocity gradient tensor

  pure real(RNP) function Lambda2(grad_v)
    real(RNP), intent(in) :: grad_v(3,3) !< vector gradient ∇v = [∂v_j/∂x_i]

    real(RNP) :: S(3,3)    ! symmetric part of grad_v
    real(RNP) :: O(3,3)    ! antisymmetric part of grad_v
    real(RNP) :: A(3,3)    ! A = S² + Ω²
    real(RNP) :: lambda(3) ! eigenvalues of A

    S = HALF * (grad_v + transpose(grad_v))
    O = HALF * (grad_v - transpose(grad_v))
    A = matmul(S,S) + matmul(O,O)

    lambda = SymEigenvalues(A)

    Lambda2 = lambda(2)

  end function Lambda2

  !=============================================================================

end module Vortex_Sensors__3D
