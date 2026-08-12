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

!> summary:   Singular value decomposition
!> author:    Joerg Stiller
!> date:      2014/12/05
!>
!>### Singular value decomposition
!==============================================================================

module Singular_Value_Decomposition
  use Kind_Parameters
  implicit none
  private

  public :: SingularValueDecomposition

contains

  !---------------------------------------------------------------------------
  !> Singular value decomposition of a matrix
  !>
  !> Computes the singular value decomposition of the real matrix A(m,n)
  !> such that
  !>
  !>     A = U S V^T
  !>
  !> where
  !>
  !>   *  U is an unitary matrix of dimension (m,m), i.e. U Uᵀ = I.
  !>   *  S a non-negative rectangular diagonal matrix of dimension (m,n).
  !>      Its diagonal entries are the singular values of A.
  !>   *  V is an unitary matrix of dimension (n,n)
  !>
  !> The current implementation employs DGESVD from LAPACK.
  !>
  !> @note
  !> Thread-safe, can be used with OpenMP.
  !> @endnote

  subroutine SingularValueDecomposition(A, s, U, Vt)
    real(RDP), intent(in)  :: A(:,:)  !< matrix of shape [m,n]
    real(RDP), intent(out) :: s(:)    !< singular values of A, size(s) = min(m,n)
    real(RDP), intent(out) :: U(:,:)  !< left unitary matrix, U
    real(RDP), intent(out) :: Vt(:,:) !< transposed right unitary matrix, Vt = Vᵀ

    real(RDP) :: C(size(A,1), size(A,2))
    real(RDP) :: work(5*max(size(A,1),size(A,2)))
    integer   :: m, n, info

    m = size(A,1)
    n = size(A,2)
    C = A
    call DGESVD('A','A', m, n, C, m, s, U, m, Vt, n, work, size(work), info)

    if (info /= 0) then
      error stop 'SingularValueDecomposition: DGESVD failed'
    end if

  end subroutine SingularValueDecomposition

  !============================================================================

end module Singular_Value_Decomposition
