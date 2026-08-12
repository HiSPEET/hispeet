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

!> summary:   Matrix operations
!> author:    Joerg Stiller
!> date:      2014/12/05
!>
!>### Matrix operations
!==============================================================================

module Matrix_Operators
  use Kind_Parameters
  use Linear_Equations,             only: LU_Decomposition, LU_Solver
  use Singular_Value_Decomposition, only: SingularValueDecomposition
  implicit none
  private

  public :: Inverse
  public :: PseudoInverse

  interface PseudoInverse
    module procedure PseudoInverse_from_SVD
    module procedure PseudoInverse_from_Matrix
  end interface PseudoInverse

  !-----------------------------------------------------------------------------
  ! private parameter

  real(RNP), parameter ::  eps0 = 1E6 * tiny(1.0_RNP)

contains

  !-----------------------------------------------------------------------------
  !> Returns the inverse of a regular matrix, or 0 in case of failure.

  pure function Inverse(A) result(AI)
    real(RNP), intent(in) :: A(:,:)  !< regular matrix
    real(RNP) :: AI(size(A,1),size(A,1))

    real(RHP) :: C(size(A,1),size(A,1)) ! work array
    integer   :: p(size(A,1))           ! permutation vector
    integer   :: i

    AI = 0

    call LU_Decomposition(A, C, p)

    forall(i=1:size(p)) AI(i,i) = 1

    do i = 1, size(p)
      call LU_Solver(C, p, AI(:,i))
    end do

  end function Inverse

  !===========================================================================
  ! Moore-Penrose pseudoinverse

  !---------------------------------------------------------------------------
  !> Computes the Moore-Penrose pseudoinverse from given SVD

  function PseudoInverse_from_SVD(s, U, Vt) result(Ai)
    real(RDP), intent(in) :: s(:)    !< singular values, size(s) = min(m,n)
    real(RDP), intent(in) :: U(:,:)  !< left unitary matrix of shape [m,m]
    real(RDP), intent(in) :: Vt(:,:) !< right unitary matrix of shape [n,n]

    real(RDP) :: Ai(size(Vt,2),size(U,1))

    real(RDP) :: V(size(Vt,2),size(Vt,1))
    integer   :: j, l

    l = min(size(U,2), size(Vt,1), count(s > eps0))
    do j = 1, l
       V(:,j) = 1/s(j) * Vt(j,:)
    end do
    Ai = matmul(V(:,1:l), transpose(U(:,1:l)))

  end function PseudoInverse_from_SVD

  !---------------------------------------------------------------------------
  !> Computes the Moore-Penrose pseudoinverse for given matrix

  function PseudoInverse_from_Matrix(A) result(Ai)
    real(RDP), intent(in)  ::  A(:,:) !< matrix of shape [m,n]

    real(RDP) :: Ai(size(A,2),size(A,1))

    real(RDP) :: s(max(size(A,1),size(A,2)))
    real(RDP) :: U(size(A,1),size(A,1))
    real(RDP) :: V(size(A,2),size(A,2))
    integer   :: j, l, m, n

    ! intialization
    m = size(A,1)
    n = size(A,2)

    ! SVD of A
    call SingularValueDecomposition(A, s, U, V)

    ! compose pseudeinverse
    l = min(m, n, count(s > eps0))
    U = transpose(U)
    V = transpose(V)
    do j = 1, l
       V(:,j) = 1/s(j) * V(:,j)
    end do
    Ai = matmul(V(:,1:l), U(1:l,:))

  end function PseudoInverse_from_Matrix

  !============================================================================

end module Matrix_Operators
