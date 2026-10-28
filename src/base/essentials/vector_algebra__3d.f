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

!> summary:  Elementary linear operations on 3-vectors and 3x3 matrices
!> author:   Joerg Stiller
!> date:     2016/08/15
!===============================================================================

module Vector_Algebra__3D
  use Kind_Parameters, only: RNP
  use Constants
  implicit none
  private

  public ::  operator(.o.)
  public ::  operator(.x.)

  public ::  Length
  public ::  UnitVector
  public ::  Projection
  public ::  SymEigenvalues

  !-----------------------------------------------------------------------------
  !> Scalar/inner product between vectors and/or matrices

  interface operator(.o.)
     module procedure  InnerProductVV
     module procedure  InnerProductVM
     module procedure  InnerProductMV
     module procedure  InnerProductMM
  end interface

  !-----------------------------------------------------------------------------
  !> Cross product between vectors

  interface operator(.x.)
     module procedure  Cross
  end interface

contains

  !=============================================================================
  ! Inner product

  !-----------------------------------------------------------------------------
  !> Scalar product of two 3-vectors

  pure function InnerProductVV(a, b) result(c)
    real(RNP), intent(in) :: a(3) !< first operand
    real(RNP), intent(in) :: b(3) !< second operand
    real(RNP) ::  c

    c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

  end function InnerProductVV

  !-----------------------------------------------------------------------------
  !> Dot product of 3-vector and 3x3 matrix

  pure function InnerProductVM(l, m) result(c)
    real(RNP), intent(in), dimension(3)   :: l !< first operand
    real(RNP), intent(in), dimension(3,3) :: m !< second operand
    real(RNP), dimension(3)   :: c

    c(1) = l(1)*m(1,1) + l(2)*m(2,1) + l(3)*m(3,1)
    c(2) = l(1)*m(1,2) + l(2)*m(2,2) + l(3)*m(3,2)
    c(3) = l(1)*m(1,3) + l(2)*m(2,3) + l(3)*m(3,3)

  end function InnerProductVM

  !-----------------------------------------------------------------------------
  !> Dot product of 3x3 matrix and 3-vector

  pure function InnerProductMV(m, r) result(c)
    real(RNP), intent(in), dimension(3,3) :: m !< first operand
    real(RNP), intent(in), dimension(3)   :: r !< second operand
    real(RNP), dimension(3)   :: c

    c(1) = m(1,1)*r(1) + m(1,2)*r(2) + m(1,3)*r(3)
    c(2) = m(2,1)*r(1) + m(2,2)*r(2) + m(2,3)*r(3)
    c(3) = m(3,1)*r(1) + m(3,2)*r(2) + m(3,3)*r(3)

  end function InnerProductMV

  !-----------------------------------------------------------------------------
  !> Inner product of two 3x3 matrices

  pure function InnerProductMM(l, r) result(c)
    real(RNP), intent(in), dimension(3,3) :: l !< first operand
    real(RNP), intent(in), dimension(3,3) :: r !< second operand
    real(RNP), dimension(3,3) :: c

    c(1,1) = l(1,1) * r(1,1) &
           + l(1,2) * r(2,1) &
           + l(1,3) * r(3,1)
    c(2,1) = l(2,1) * r(1,1) &
           + l(2,2) * r(2,1) &
           + l(2,3) * r(3,1)
    c(3,1) = l(3,1) * r(1,1) &
           + l(3,2) * r(2,1) &
           + l(3,3) * r(3,1)
    c(1,2) = l(1,1) * r(1,2) &
           + l(1,2) * r(2,2) &
           + l(1,3) * r(3,2)
    c(2,2) = l(2,1) * r(1,2) &
           + l(2,2) * r(2,2) &
           + l(2,3) * r(3,2)
    c(3,2) = l(3,1) * r(1,2) &
           + l(3,2) * r(2,2) &
           + l(3,3) * r(3,2)
    c(1,3) = l(1,1) * r(1,3) &
           + l(1,2) * r(2,3) &
           + l(1,3) * r(3,3)
    c(2,3) = l(2,1) * r(1,3) &
           + l(2,2) * r(2,3) &
           + l(2,3) * r(3,3)
    c(3,3) = l(3,1) * r(1,3) &
           + l(3,2) * r(2,3) &
           + l(3,3) * r(3,3)

  end function InnerProductMM

  !=============================================================================

  !-----------------------------------------------------------------------------
  !> Cross product of two vectors

  pure function Cross(a, b) result(c)
    real(RNP), intent(in) :: a(3) !< first operand
    real(RNP), intent(in) :: b(3) !< second operand
    real(RNP) :: c(3)

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)

  end function Cross

  !-----------------------------------------------------------------------------
  !> Length of a vector

  pure function Length(a) result(l)
    real(RNP), intent(in) ::  a(3) !< given vector
    real(RNP) ::  l

    l = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

  end function Length

  !-----------------------------------------------------------------------------
  !> Returns normalized (unit) vector to a, or zero if a is too short

  pure function UnitVector(a) result(b)
    real(RNP), intent(in) :: a(3) !< given vector
    real(RNP) :: b(3)

    real(RNP) :: l

    l = Length(a)
    if (l > epsilon(l)) then
       b = a / l
    else
       b = 0
    end if

  end function UnitVector

  !-----------------------------------------------------------------------------
  !> Projection of vector a to vector b: c = (a·b)/(b·b) b

  pure function Projection(a, b) result(c)
    real(RNP), intent(in) :: a(3) !< vector which is to be projected
    real(RNP), intent(in) :: b(3) !< vector to which a is to be projected
    real(RNP) ::  c(3)

    c = UnitVector(b)
    c = (a .o. c) * c

  end function Projection

  !-----------------------------------------------------------------------------
  !> Eigenvalues of a symmetric matrix in ascending order, λ₁ ≤ λ₂ ≤ λ₃

  pure function SymEigenvalues(A) result(lmb)
    real(RNP), intent(in) :: A(3,3)
    real(RNP)             :: lmb(3)

    real(RNP) :: p1, p2, p, q, r, phi
    real(RNP) :: b11, b22, b33, b12, b13, b23, detB

    ! sum of squares of the off-diagonal elements
    p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

    ! special case: diagonal matrix
    if (p1 == ZERO) then
       lmb = [A(1,1), A(2,2), A(3,3)]
       call SortAscending(lmb)
       return
    end if

    q  = (A(1,1) + A(2,2) + A(3,3)) * THIRD
    p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + TWO*p1
    p  = sqrt(p2 / 6)

    ! B = (A - q*I) / p  (symmetric, only 6 elements needed)
    b11 = (A(1,1) - q) / p
    b22 = (A(2,2) - q) / p
    b33 = (A(3,3) - q) / p
    b12 = A(1,2) / p
    b13 = A(1,3) / p
    b23 = A(2,3) / p

    detB = b11*(b22*b33 - b23*b23) &
         - b12*(b12*b33 - b23*b13) &
         + b13*(b12*b23 - b22*b13)

    ! r lies in [-1,1] in theory; guard against rounding errors
    r   = max(-ONE, min(ONE, HALF*detB))
    phi = acos(r) * THIRD

    lmb(1) = q + TWO*p*cos(phi + TWO*THIRD*real(PI, RNP))
    lmb(3) = q + TWO*p*cos(phi)
    lmb(2) = THREE*q - lmb(1) - lmb(3)

  end function SymEigenvalues

  !-----------------------------------------------------------------------------
   !> Sorts three values in ascending order
   pure subroutine SortAscending(x)
      real(RNP), intent(inout) :: x(3)
      if (x(1) > x(2)) call Swap(x(1), x(2))
      if (x(2) > x(3)) call Swap(x(2), x(3))
      if (x(1) > x(2)) call Swap(x(1), x(2))
   end subroutine SortAscending

   pure subroutine Swap(a, b)
      real(RNP), intent(inout) :: a, b
      real(RNP) :: t
      t = a;  a = b;  b = t
   end subroutine Swap

  !=============================================================================

end module Vector_Algebra__3D
