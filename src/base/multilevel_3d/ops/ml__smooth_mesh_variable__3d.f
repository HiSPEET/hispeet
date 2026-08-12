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

!> summary:  Utility for removing jumps and smoothing multilevel mesh variables
!> author:   Joerg Stiller
!> date:     2026/08/10
!===============================================================================

module ML__Smooth_Mesh_Variable__3D
  use Smooth_Mesh_Data__3D
  use Child_To_Parent_Projection__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  implicit none
  private

  public :: ML_SmoothMeshVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Remove jumps and smooth multilevel mesh data

  subroutine ML_SmoothMeshVariable_3D(ml_op, u, filter, order)
    class(ML_MeshOperators_3D), intent(in)    :: ml_op !< ML mesh operators
    class(ML_MeshVariable_3D),  intent(inout) :: u     !< ML mesh variable
    integer, intent(in) :: filter !< 0/1/2/3: none/cut-off/erfc-log/exponential
    integer, intent(in) :: order  !< filter order (0: auto)

    integer :: l

    do l = 1, size(u%level)

      ! smooth current mesh level
      call SmoothMeshData_3D( mesh   = ml_op % sem(l) % mesh   &
                            , eop    = ml_op % sem(l) % std_op &
                            , u      = u % level(l) % val      &
                            , filter = filter                  &
                            , order  = order                   )

      if (l == 1) exit

      ! project to parent level (currently using embedded interpolation)
      call ChildToParentProjection_3D( child  = ml_op % sem(l  ) % mesh &
                                     , parent = ml_op % sem(l-1) % mesh &
                                     , pop    = ml_op % iop_fc_x(l)     &
                                     , v_c    = u % level(l  ) % val    &
                                     , v_p    = u % level(l-1) % val    )

    end do

  end subroutine ML_SmoothMeshVariable_3D

  !=============================================================================

end module ML__Smooth_Mesh_Variable__3D
