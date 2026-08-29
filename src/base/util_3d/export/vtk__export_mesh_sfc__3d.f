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

!> summary:  Export SFC through 3d mesh element centers to legacy VTK file
!> author:   Joerg Stiller
!> date:     2026/01/22
!===============================================================================

module VTK__Export_Mesh_SFC__3D
  use Execution_Control
  use VTK_Binding
  use Mesh__3D
  implicit none
  private

  public :: VTK_ExportMeshSFC_3D

contains

  !-----------------------------------------------------------------------------
  !> Export space filling curve through element centers to legacy VTK file.

  subroutine VTK_ExportMeshSFC_3D(mesh, file)
    class(Mesh_3D),    intent(in) :: mesh    !< local mesh partition
    character(len=*),  intent(in) :: file    !< base name of VTK file

    ! VTK data .................................................................

    integer(VTK_INT32) :: cell_type = VTK_POLY_LINE
    integer(VTK_INT32), allocatable :: cells(:,:)
    real(VTK_FLOAT64),  allocatable :: points(:,:)
    integer(VTK_INT32), allocatable :: attrib(:,:)

    ! internal variables .......................................................

    character(len=8), parameter :: attrib_names(2) = [ 'sfc_rank', 'elem_id ' ]

    integer, allocatable :: map(:)
    integer :: e, k, ne, ns, rk_min, rk_max, sfc_rank_0

    ! initialization :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    if (mesh % part < 0 .or. .not. mesh % has_sfc) then
      return
    end if

    ne = mesh % n_elem
    ns = size(attrib_names)

    rk_min = minval(mesh % element % sfc_rank)
    rk_max = maxval(mesh % element % sfc_rank)
    if (rk_max - rk_min + 1 /= ne) then
      block
        character(len = 80) :: msg
        write(msg,'(3(A,X,I0))') &
          'SFC fragmented: ne =', ne, ', rk_min/max =', rk_min, ',', rk_max
        call Error('VTK_ExportMeshSFC_3D',trim(msg),'VTK__Export_Mesh_SFC__3D')
      end block
    end if
    sfc_rank_0 = 1 - rk_min

    ! map from elements to SFC vertices
    allocate(map(ne))
    do e = 1, ne
      map(mesh%element(e)%sfc_rank + sfc_rank_0) = e
    end do

    ! data :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! points ...................................................................

    allocate(points(3,ne))

    ! element midpoints points ordered according to their SFC rank
    do k = 1, ne
      points(1:3,k) = mesh % element(map(k)) % geometry % x_c(0,1:3)
    end do

    ! cells ....................................................................

    allocate(cells(ne,1))

    do k = 1, ne
      cells(k,1) = k - 1
    end do

    ! scalars ..................................................................

    allocate(attrib(ne,ns))

    ! point data
    do k = 1, ne
      attrib(k,1) = k - sfc_rank_0  ! SFC rank
      attrib(k,2) = map(k)          ! element ID
    end do

    ! Export :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    call VTK_WriteXML_Unstructured( points, cells, cell_type  &
                                  , pa       = attrib         &
                                  , pa_names = attrib_names   &
                                  , file     = file           &
                                  , piece    = mesh % part    &
                                  , n_pieces = mesh % n_parts )

  end subroutine VTK_ExportMeshSFC_3D

  !=============================================================================

end module VTK__Export_Mesh_SFC__3D
