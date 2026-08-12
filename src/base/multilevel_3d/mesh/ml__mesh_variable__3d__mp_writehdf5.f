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

!> summary:  Writing a multilevel mesh variable to HDF5
!> author:   Joerg Stiller
!> date:     2025/07/29
!===============================================================================

submodule(ML__Mesh_Variable__3D) MP_WriteHDF5
  use, intrinsic :: ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Write multilevel mesh partition into HDF5 file

  module subroutine WriteHDF5(this, file)
    class(ML_MeshVariable_3D), target, intent(in) :: this !< ML mesh variable
    character(len=*),                  intent(in) :: file !< name of HDF5 file

    character(len=:), allocatable :: file_pr
    character(len=80) :: tag
    integer(HID_T)    :: file_id, group_id, ml_variable_id
    integer(HID_T)    :: data_id, space_id, type_id
    integer(HSIZE_T)  :: dims(5)
    integer, target   :: l_top, ln, nc
    integer           :: err, l, n_dims

    ! preliminaries ............................................................

    ! safeguard
    call Init_HDF5_Binding()

    ! append process rank to file name
    write(tag,'(I0)') this % level(1) % mesh % proc
    file_pr = trim(file)//'_'//trim(tag)//'.h5'

    l_top = size(this % level)
    ln    = len (this % name) ! length of component names
    nc    = size(this % name) ! number of components

    ! HDF5 file and base group .................................................

    call H5Fcreate_f(file_pr, H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, '/ml_variable', ml_variable_id, err)

    ! attributes ...............................................................

    n_dims = 1

    call H5Gcreate_f(file_id, '/ml_variable/attrib', group_id, err)

    ! l_top
    dims = 1
    call H5Screate_simple_f(n_dims, dims, space_id, err)
    call H5Dcreate_f(group_id, 'l_top', H5T_INTEGER, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_INTEGER, C_Loc(l_top), err)
    call H5Dclose_f(data_id, err)
    call H5Sclose_f(space_id, err)

    ! ln
    dims = 1
    call H5Screate_simple_f(n_dims, dims, space_id, err)
    call H5Dcreate_f(group_id, 'ln', H5T_INTEGER, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_INTEGER, C_Loc(ln), err)
    call H5Dclose_f(data_id, err)
    call H5Sclose_f(space_id, err)

    ! nc
    dims = 1
    call H5Screate_simple_f(n_dims, dims, space_id, err)
    call H5Dcreate_f(group_id, 'nc', H5T_INTEGER, space_id, data_id, err)
    call H5Dwrite_f(data_id, H5T_INTEGER, C_Loc(nc), err)
    call H5Dclose_f(data_id, err)
    call H5Sclose_f(space_id, err)

    ! name
    dims = nc
    call H5Tcopy_f(H5T_CHARACTER, type_id, err)
    call H5Tset_size_f(type_id, int(ln, SIZE_T), err)
    call H5Screate_simple_f(n_dims, dims, space_id, err)
    call H5Dcreate_f(group_id, 'name', type_id, space_id, data_id, err)
    call H5Dwrite_f(data_id, type_id, C_Loc(this%name(1)(1:1)), err)
    call H5Dclose_f(data_id, err)
    call H5Sclose_f(space_id, err)
    call H5Tclose_f(type_id, err)

    call H5Gclose_f(group_id, err)

    ! levels ...................................................................

    n_dims = 5

    do l = 1, l_top
      dims = shape(this % level(l) % val)
      write(tag,'(I0)') l
      call H5Gcreate_f(file_id, '/ml_variable/level_'//trim(tag), group_id, err)
      call H5Screate_simple_f(n_dims, dims, space_id, err)
      call H5Dcreate_f(group_id, 'val', H5T_REAL_RNP, space_id, data_id, err)
      call H5Dwrite_f(data_id, H5T_REAL_RNP, &
                      C_Loc(this%level(l)%val(0,0,0,1,1)), err)
      call H5Dclose_f(data_id, err)
      call H5Sclose_f(space_id, err)
      call H5Gclose_f(group_id, err)
    end do

    ! finalize .................................................................

    call H5Gclose_f(ml_variable_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine WriteHDF5

  !=============================================================================

end submodule MP_WriteHDF5
