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

!> summary:  Support for affine transformations in 3D
!> author:   Joerg Stiller
!> date:     2025/06/18
!===============================================================================

module Affine_Transformation__3D
  use Kind_Parameters, only: RNP
  use Execution_Control
  implicit none
  private

  public :: InverseAffineMap_3D

  real(RNP), parameter, public :: &
      AFFINE_IDENTITY_MAP_3D(4,4) = reshape( [ 1, 0, 0, 0 &
                                             , 0, 1, 0, 0 &
                                             , 0, 0, 1, 0 &
                                             , 0, 0, 0, 1 ], [4,4] )

contains

  !-----------------------------------------------------------------------------
  !> Returns the inverse affine map

  function InverseAffineMap_3D(map) result(inv_map)
    real(RNP), intent(in) :: map(4,4)
    real(RNP) :: inv_map(4,4)

    real(RNP) :: A(3,3), A_inv(3,3), b(3), b_inv(3), det_A

    A = map(1:3,1:3)  ! rotation ...
    b = map(1:3,4)    ! translation

    det_A = A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2)) &
          + A(1,2) * (A(2,3) * A(3,1) - A(2,1) * A(3,3)) &
          + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))

    if (abs(det_A) < epsilon(1.0)) then
      call Error('InverseAffinityMatrix','singular matrix A','Import_GMSH__3D')
    end if

    A_inv(1,1) =  A(2,2)*A(3,3) - A(3,2)*A(2,3)
    A_inv(1,2) =  A(3,2)*A(1,3) - A(1,2)*A(3,3)
    A_inv(1,3) =  A(1,2)*A(2,3) - A(2,2)*A(1,3)

    A_inv(2,1) =  A(2,3)*A(3,1) - A(3,3)*A(2,1)
    A_inv(2,2) =  A(3,3)*A(1,1) - A(1,3)*A(3,1)
    A_inv(2,3) =  A(1,3)*A(2,1) - A(2,3)*A(1,1)

    A_inv(3,1) =  A(2,1)*A(3,2) - A(3,1)*A(2,2)
    A_inv(3,2) =  A(3,1)*A(1,2) - A(1,1)*A(3,2)
    A_inv(3,3) =  A(1,1)*A(2,2) - A(2,1)*A(1,2)

    A_inv = 1/det_A * A_inv

    b_inv(1) = -(A_inv(1,1)*b(1) + A_inv(1,2)*b(2) + A_inv(1,3)*b(3))
    b_inv(2) = -(A_inv(2,1)*b(1) + A_inv(2,2)*b(2) + A_inv(2,3)*b(3))
    b_inv(3) = -(A_inv(3,1)*b(1) + A_inv(3,2)*b(2) + A_inv(3,3)*b(3))

    inv_map(1:3,1:3) = A_inv
    inv_map(1:3, 4 ) = b_inv
    inv_map( 4 ,1:3) = 0
    inv_map( 4 , 4 ) = 1

  end function InverseAffineMap_3D

!=============================================================================

end module Affine_Transformation__3D
