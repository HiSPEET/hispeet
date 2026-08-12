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

!> summary:  Reading a multilevel mesh partition from HDF5
!> author:   Joerg Stiller
!> date:     2024/10/29
!===============================================================================

submodule(ML__Mesh__3D) MP_ReadHDF5
  use, intrinsic :: ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Read multilevel mesh partition from HDF5 file

  module subroutine ReadHDF5(this, file, comm)
    class(ML_Mesh_3D), intent(inout) :: this !< multilevel mesh partition
    character(len=*),  intent(in)    :: file !< name of HDF5 file
    type(MPI_Comm),    intent(in)    :: comm !< MPI "mesh_world" communicator

    character(len=:), allocatable :: file_pr
    character(len=80) :: tag
    integer(HID_T)    :: file_id, group_id
    integer, target   :: l_top
    integer           :: err, rank, l
    logical           :: exists

    ! preliminaries ............................................................

    ! safeguard
    call Init_HDF5_Binding()

    ! append process rank to file name and check if it exists
    call MPI_Comm_rank(comm, rank)
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'
    inquire(file=file_pr, exist=exists)

    if (exists) then
      call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
    end if

    ! attributes ...............................................................

    block
      integer(HID_T) :: data_id, type_id
      type(C_Ptr)    :: buf

      if (rank == 0) then

        call H5Gopen_f(file_id, '/ml_mesh/attrib', group_id, err)

        ! l_top
        buf = C_Loc(l_top)
        call H5Dopen_f(group_id, 'l_top', data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        call H5Gclose_f(group_id, err)

      end if

      call XMPI_Bcast(l_top, 0, comm)

    end block

    ! mesh levels ..............................................................

    allocate(this % mesh(l_top))

    do l = 1, l_top

      if (exists) then
        write(tag,'(I0)') l
        call H5Gopen_f(file_id, '/ml_mesh/level_'//trim(tag), group_id, err)
      else
        group_id = H5I_INVALID_HID_F
      end if

      call this % mesh(l) % ReadHDF5(group_id, comm)

      if (exists) then
        call H5Gclose_f(group_id, err)
      end if

    end do

    ! finalize .................................................................

    call H5Fclose_f(file_id, err)

  end subroutine ReadHDF5

  !=============================================================================

end submodule MP_ReadHDF5
