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
  implicit none
  private

  public ::  operator(.o.)
  public ::  operator(.x.)

  public ::  Length
  public ::  UnitVector
  public ::  Projection

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

!===============================================================================
! Inner product

!-------------------------------------------------------------------------------
!> Scalar product of two 3-vectors

pure function InnerProductVV(a, b) result(c)
  real(RNP), intent(in) :: a(3) !< first operand
  real(RNP), intent(in) :: b(3) !< second operand
  real(RNP) ::  c

  c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

end function InnerProductVV

!-------------------------------------------------------------------------------
!> Dot product of 3-vector and 3x3 matrix

pure function InnerProductVM(l, m) result(c)
  real(RNP), intent(in), dimension(3)   :: l !< first operand
  real(RNP), intent(in), dimension(3,3) :: m !< second operand
  real(RNP), dimension(3)   :: c

  c(1) = l(1)*m(1,1) + l(2)*m(2,1) + l(3)*m(3,1)
  c(2) = l(1)*m(1,2) + l(2)*m(2,2) + l(3)*m(3,2)
  c(3) = l(1)*m(1,3) + l(2)*m(2,3) + l(3)*m(3,3)

end function InnerProductVM

!-------------------------------------------------------------------------------
!> Dot product of 3x3 matrix and 3-vector

pure function InnerProductMV(m, r) result(c)
  real(RNP), intent(in), dimension(3,3) :: m !< first operand
  real(RNP), intent(in), dimension(3)   :: r !< second operand
  real(RNP), dimension(3)   :: c

  c(1) = m(1,1)*r(1) + m(1,2)*r(2) + m(1,3)*r(3)
  c(2) = m(2,1)*r(1) + m(2,2)*r(2) + m(2,3)*r(3)
  c(3) = m(3,1)*r(1) + m(3,2)*r(2) + m(3,3)*r(3)

end function InnerProductMV

!-------------------------------------------------------------------------------
!> Inner product of two 3x3 matrices

pure function InnerProductMM(l, r) result(c)
  real(RNP), intent(in), dimension(3,3) :: l !< first operand
  real(RNP), intent(in), dimension(3,3) :: r !< second operand
  real(RNP), dimension(3,3) :: c

  c(1,1) = l(1,1) *r (1,1) &
         + l(1,2) *r (2,1) &
         + l(1,3) *r (3,1)
  c(2,1) = l(2,1) *r (1,1) &
         + l(2,2) *r (2,1) &
         + l(2,3) *r (3,1)
  c(3,1) = l(3,1) *r (1,1) &
         + l(3,2) *r (2,1) &
         + l(3,3) *r (3,1)
  c(1,2) = l(1,1) *r (1,2) &
         + l(1,2) *r (2,2) &
         + l(1,3) *r (3,2)
  c(2,2) = l(2,1) *r (1,2) &
         + l(2,2) *r (2,2) &
         + l(2,3) *r (3,2)
  c(3,2) = l(3,1) *r (1,2) &
         + l(3,2) *r (2,2) &
         + l(3,3) *r (3,2)
  c(1,3) = l(1,1) *r (1,3) &
         + l(1,2) *r (2,3) &
         + l(1,3) *r (3,3)
  c(2,3) = l(2,1) *r (1,3) &
         + l(2,2) *r (2,3) &
         + l(2,3) *r (3,3)
  c(3,3) = l(3,1) *r (1,3) &
         + l(3,2) *r (2,3) &
         + l(3,3) *r (3,3)

end function InnerProductMM

!===============================================================================

!-------------------------------------------------------------------------------
!> Cross product of two vectors

pure function Cross(a, b) result(c)
  real(RNP), intent(in) :: a(3) !< first operand
  real(RNP), intent(in) :: b(3) !< second operand
  real(RNP) :: c(3)

  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)

end function Cross

!-------------------------------------------------------------------------------
!> Length of a vector

pure function Length(a) result(l)
  real(RNP), intent(in) ::  a(3) !< given vector
  real(RNP) ::  l

  l = sqrt(a(1)**2 + a(2)**2 + a(3)**2)

end function Length

!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
!> Projection of vector a to vector b: c = (a·b)/(b·b) b

pure function Projection(a, b) result(c)
  real(RNP), intent(in) :: a(3) !< vector which is to be projected
  real(RNP), intent(in) :: b(3) !< vector to which a is to be projected
  real(RNP) ::  c(3)

  c = UnitVector(b)
  c = (a .o. c) * c

end function Projection

!===============================================================================

end module Vector_Algebra__3D
