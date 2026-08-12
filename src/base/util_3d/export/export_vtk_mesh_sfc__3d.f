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

module Export_VTK_Mesh_SFC__3D
  use Execution_Control
  use C_Binding
  use VTK_Binding
  use Mesh__3D
  implicit none
  private

  public :: ExportVTK_MeshSFC

contains

  !-----------------------------------------------------------------------------
  !> Export space filling curve through element centers to legacy VTK file.

  subroutine ExportVTK_MeshSFC(mesh, file)
    class(Mesh_3D),   intent(in) :: mesh !< local mesh partition
    character(len=*), intent(in) :: file !< base name of VTK file

    ! VTK data .................................................................

    integer(C_INT) :: vtk     ! writer handle
    integer(C_INT) :: success ! success flag

    integer(C_VTK_ID), allocatable :: cell(:,:)
    real(C_DOUBLE),    allocatable :: x(:,:)
    integer(C_INT),    allocatable :: s(:,:)

    ! internal variables .......................................................

    character(len=80) :: tag
    integer, allocatable :: map(:)
    integer :: e, k, ne, rk_min, rk_max, sfc_rank_0
    integer :: pvtu

    ! initialization ...........................................................

    if (mesh % part < 0 .or. .not. mesh % has_sfc) then
      return
    else
      write(tag, fmt='(A2,I0)') '_p', mesh % part
    end if

    ne = mesh % n_elem

    rk_min = minval(mesh % element % sfc_rank)
    rk_max = maxval(mesh % element % sfc_rank)
    if (rk_max - rk_min + 1 /= ne) then
      block
        character(len = 80) :: msg
        write(msg,'(3(A,X,I0))') &
          'SFC fragmented: ne =', ne, ', rk_min/max =', rk_min, ',', rk_max
        call Error('ExportVTK_MeshSFC', trim(msg), 'Export_VTK_Mesh_SFC__3D')
      end block
    end if
    sfc_rank_0 = 1 - rk_min

    ! map from elements to SFC vertices
    allocate(map(ne))
    do e = 1, ne
      map(mesh%element(e)%sfc_rank + sfc_rank_0) = e
    end do

    ! set up VTK file ..........................................................

    call VTK_XMLWriter_New(vtk)
    call VTK_XMLWriter_SetDataObjectType(vtk, VTK_UNSTRUCTURED_GRID)
    call VTK_XMLWriter_SetDataModeType(vtk, VTK_APPENDED)
    call VTK_XMLWriter_SetFileName(vtk, trim(file)//trim(tag)//'.vtu')

    ! points ...................................................................

    allocate(x(3,ne))

    ! element midpoints points ordered according to their SFC rank
    do k = 1, ne
      x(1:3,k) = mesh % element(map(k)) % geometry % x_c(0,1:3)
    end do

    call VTK_XMLWriter_SetPoints(vtk, x)

    ! cells ....................................................................

    allocate(cell(0:ne,1))

    ! number of polyline points
    cell(0,1) = ne

    ! point indices
    do k = 1, ne
      cell(k,1) = k - 1
    end do

    call VTK_XMLWriter_SetCellsWithType(vtk, VTK_POLY_LINE, cell)

    ! scalars ..................................................................

    allocate(s(ne,2))

    ! point data
    do k = 1, ne
      s(k,1) = k - sfc_rank_0  ! SFC rank
      s(k,2) = map(k)          ! element ID
    end do

    call VTK_XMLWriter_SetPointData(vtk, 'sfc_rank', s(:,1), '')
    call VTK_XMLWriter_SetPointData(vtk, 'elem_id' , s(:,2), '')

    ! write VTK file ...........................................................

    call VTK_XMLWriter_Write(vtk, success)
    call VTK_XMLWriter_Delete(vtk)

    ! write PVTU file ..........................................................

    if (mesh % part == 0) then

      open(newunit=pvtu, file=trim(file)//'.pvtu')

      ! header
      write(pvtu,'(1A)') '<?xml version="1.0"?>'
      write(pvtu,'(3A)') '<VTKFile type="PUnstructuredGrid" version="0.1" ', &
                         'byte_order="LittleEndian" ',                       &
                         'compressor="vtkZLibDataCompressor">'

      write(pvtu,'(2X,A)') '<PUnstructuredGrid GhostLevel="0">'

      ! point data section
      write(pvtu,'(4X,A)')  '<PPointData>'
      write(pvtu,'(6X,4A)') '<PDataArray type="Int32" ','Name="sfc_rank"/>'
      write(pvtu,'(6X,4A)') '<PDataArray type="Int32" ','Name="elem_id"/>'
      write(pvtu,'(4X,A)')  '</PPointData>'

      ! points section
      write(pvtu,'(4X,A)') '<PPoints>'
      write(pvtu,'(6X,A)') '<PDataArray type="Float64" NumberOfComponents="3"/>'
      write(pvtu,'(4X,A)') '</PPoints>'

      ! piece sources
      do k = 0, mesh%n_parts-1
        write(pvtu,'(4X,3A,I0,A)') '<Piece Source="',trim(file),'_p',k,'.vtu"/>'
      end do

      ! trailer
      write(pvtu,'(2X,1A)') '</PUnstructuredGrid>'
      write(pvtu,'(A)') '</VTKFile>'

      ! close file
      close(pvtu)

    end if

  end subroutine ExportVTK_MeshSFC

  !=============================================================================

end module Export_VTK_Mesh_SFC__3D
