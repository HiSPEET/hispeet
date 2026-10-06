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

!> summary:  Read a multilevel mesh variable from HDF5
!> author:   Joerg Stiller
!> date:     2025/07/29
!===============================================================================

submodule(ML__Mesh_Variable__3D) MP_ReadHDF5
  use, intrinsic :: ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Read multilevel mesh variable from HDF5 file
  !>
  !> The variable must be initialized with a sufficient number of components.
  !> The components stored in `file` are updated, including their names.

  module subroutine ReadHDF5(this, file)
    class(ML_MeshVariable_3D), intent(inout) :: this !< multilevel mesh variable
    character(len=*),          intent(in)    :: file !< name of HDF5 file

    character(len=:), allocatable :: file_pr
    character(len=:), allocatable, target :: name(:), name_old(:)
    character(len=80) :: tag
    integer(HID_T)    :: data_id, file_id, group_id, type_id
    type(C_Ptr)       :: buf
    integer, target   :: l_top, ln, nc
    integer           :: err, l
    integer           :: ln_old, nc_old
    logical           :: exists

    associate( comm => this % level(1) % mesh % comm_world &
             , proc => this % level(1) % mesh % proc       )

      ! preliminaries ..........................................................

      ! safeguard
      call Init_HDF5_Binding()

      ! append process rank to file name and check if it exists
      write(tag,'(I0)') proc
      file_pr = trim(file)//'_'//trim(tag)//'.h5'
      inquire(file=file_pr, exist=exists)

      if (exists) then
        call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
      end if

      ! attributes .............................................................

      ln_old = len (this%name)
      nc_old = size(this%name)

      if (proc == 0) then

        call H5Gopen_f(file_id, '/ml_variable/attrib', group_id, err)

        ! l_top
        buf = C_Loc(l_top)
        call H5Dopen_f(group_id, 'l_top', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        ! ln
        buf = C_Loc(ln)
        call H5Dopen_f(group_id, 'ln', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        ! nc
        buf = C_Loc(nc)
        call H5Dopen_f(group_id, 'nc', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        ! name
        allocate(character(ln) :: name(nc))
        buf = C_Loc(name(1)(1:1))
        call H5Dopen_f(group_id, 'name', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        call H5Gclose_f(group_id, err)

      end if

      ! globalize read attributes
      call XMPI_Bcast(l_top, 0, comm)
      call XMPI_Bcast(ln   , 0, comm)
      call XMPI_Bcast(nc   , 0, comm)
      if (proc > 0) then
        allocate(character(ln) :: name(nc))
      end if
      call XMPI_Bcast(name, 0, comm)

      ! check dimensions
      if (l_top > size(this % level) .or. nc > nc_old) then
        call Error( 'ReadHDF5' &
                  , 'mismatch between given and stored variable dimensions' &
                  , 'MP_WriteHDF5 @ ML__Mesh_Variable__3D')
      end if

      ! update names
      if (nc == nc_old) then
        this%name = name
      else
        if (ln > ln_old) then
          call move_alloc(this%name, name_old)
          allocate(character(ln) :: this%name(nc_old))
          this%name(nc+1:nc_old) = name_old(nc+1:nc_old)
          deallocate(name_old)
        end if
        this%name(1:nc) = name
      end if

    end associate

    ! levels ...................................................................

    if (exists) then
      do l = 1, l_top
        if (size(this%level(l)%val) == 0) cycle
        buf = C_Loc(this % level(l) % val(0,0,0,1,1))
        write(tag,'(I0)') l
        call H5Gopen_f(file_id, '/ml_variable/level_'//trim(tag), group_id, err)
        call H5Dopen_f(group_id, 'val', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)
        call H5Gclose_f(group_id, err)
      end do
    end if

    ! finalize .................................................................

    call H5Fclose_f(file_id, err)

  end subroutine ReadHDF5

  !=============================================================================

end submodule MP_ReadHDF5
