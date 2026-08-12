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

!> summary:  Reading a mesh partition from HDF5
!> author:   Joerg Stiller, Erik Pfister
!> date:     2024/06/07
!>
!> @todo
!> Remove deallocate statements once HDF5 is fixed
!===============================================================================

submodule(Mesh__3D) MP_ReadHDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Read mesh partition from given HDF5 file

  module subroutine ReadHDF5_F(mesh, file, comm)
    class(Mesh_3D),   intent(inout) :: mesh !< mesh partition
    character(len=*), intent(in)    :: file !< name of HDF5 file
    type(MPI_Comm),   intent(in)    :: comm !< MPI "world" communicator

    integer(HID_T)    :: file_id, group_id
    integer           :: err, rank
    logical           :: exists
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    ! append process rank to file name and check if it exists
    call MPI_Comm_rank(comm, rank)
    write(tag,'(I0)') rank
    file_pr = trim(file)//'_'//trim(tag)//'.h5'
    inquire(file=file_pr, exist=exists)

    if (exists) then
      ! open HDF5 file and group for reading
      call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
      call H5Gopen_f(file_id, '/mesh', group_id, err)
    else
      group_id = H5I_INVALID_HID_F
    end if

    ! read mesh partition
    call ReadHDF5_G(mesh, group_id, comm)

    if (exists) then
      ! close HDF5 group and file
      call H5Gclose_f(group_id, err)
      call H5Fclose_f(file_id, err)
    end if

  end subroutine ReadHDF5_F

  !-----------------------------------------------------------------------------
  !> Read mesh partition from given HDF5 group
  !>
  !> Pass `group_id = H5I_INVALID_HID_F` if no input data is available for the
  !> present MPI process.

  module subroutine ReadHDF5_G(mesh, group_id, comm)
    class(Mesh_3D), target, intent(inout) :: mesh  !< mesh partition
    integer(HID_T), intent(in) :: group_id !< ID of related HDF5 group
    type(MPI_Comm), intent(in) :: comm     !< MPI "world" communicator

    ! internal data ............................................................

    ! dynamical data
    type(MeshAttributes_3D),                      target :: attrib
    type(MeshElementNeighbor_3D),    allocatable, target :: neighbor(:)
    real(RNP),                       allocatable, target :: xc(:)
    integer,                         allocatable, target :: mesh_dim(:)

    ! HDF5 datatype, dataspace, dimensions and buffer adress pointer
    integer(HID_T)   :: type_id
    integer(HID_T)   :: space_id
    integer(HID_T)   :: data_id
    integer(HSIZE_T) :: dims(1), maxdims(1)
    type(C_Ptr)      :: buf

    ! HDF5 dataset names
    character(len=*), parameter :: name_mp  = 'mesh_part'
    character(len=*), parameter :: name_ma  = 'mesh_attributes'
    character(len=*), parameter :: name_mba = 'mesh_boundary_attributes'
    character(len=*), parameter :: name_md  = 'mesh_dimensions'
    character(len=*), parameter :: name_me  = 'mesh_elements'
    character(len=*), parameter :: name_men = 'mesh_element_neighbors'
    character(len=*), parameter :: name_mec = 'mesh_element_coordinates'

    integer :: e, i, nc, nd, nn, np, po, rank
    integer :: err

    ! safeguard ................................................................

    call Init_HDF5_Binding()

    call mesh % Delete_Mesh_3D()

    ! preliminaries ............................................................

    call MPI_Comm_rank(comm, rank)

    if (group_id /= H5I_INVALID_HID_F) then

      ! get mesh part ..........................................................

      buf = C_Loc(mesh%part)
      call H5Dopen_f(group_id, name_mp, data_id, err)
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Dclose_f(data_id, err)
      call H5Tclose_f(type_id, err)

      ! get mesh attributes ....................................................

      buf = C_Loc(attrib)
      call H5Dopen_f(group_id, name_ma, data_id, err)
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Dclose_f(data_id, err)
      call H5Tclose_f(type_id, err)

      ! get mesh boundary attributes ...........................................

      allocate(attrib % boundary(attrib%n_bound))

      if (attrib%n_bound > 0) then
        buf = C_Loc(attrib%boundary)
        call H5Dopen_f(group_id, name_mba, data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)
      end if

      ! initialize mesh ........................................................

      call mesh % Init_Mesh_3D(attrib, comm)

      ! get mesh dimensions ....................................................

      call H5Dopen_f(group_id, name_md, data_id, err)
      call H5Dget_space_f(data_id, space_id, err)
      call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)
      call H5Sclose_f(space_id, err)

      nd = dims(1)
      allocate(mesh_dim(nd))

      buf = C_Loc(mesh_dim)
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Dclose_f(data_id, err)
      call H5Tclose_f(type_id, err)

      mesh % n_vert        = mesh_dim( 1)
      mesh % n_edge        = mesh_dim( 2)
      mesh % n_face        = mesh_dim( 3)
      mesh % n_elem        = mesh_dim( 4)
      mesh % n_elem_active = mesh_dim( 5)
      mesh % n_elem_frozen = mesh_dim( 6)
      mesh % n_cluster     = mesh_dim( 7)
      mesh % n_ghost       = mesh_dim( 8)
      mesh % n_link        = mesh_dim( 9)
      mesh % n_child       = mesh_dim(10)
      mesh % n_parent      = mesh_dim(11)
      mesh % n_elem_1      = mesh_dim(12)
      mesh % n_elem_2      = mesh_dim(13)
      mesh % n_elem_3      = mesh_dim(14)

      ! get elements ...........................................................

      call H5Dopen_f(group_id, name_me, data_id, err)
      call H5Dget_space_f(data_id, space_id, err)
    ! call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)
      call H5Sclose_f(space_id, err)

    ! mesh % n_elem = dims(1)               ! retain read-in value
      allocate(mesh % element(mesh%n_elem))

      if (mesh%n_elem > 0) then
        buf = C_Loc(mesh%element)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Tclose_f(type_id, err)
      end if

      call H5Dclose_f(data_id, err)

      ! get element neighbors ..................................................

      if (mesh%n_elem > 0) then

        call H5Dopen_f(group_id, name_men, data_id, err)
        call H5Dget_space_f(data_id, space_id, err)
        call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)
        call H5Sclose_f(space_id, err)

        nn = dims(1)
        allocate(neighbor(nn))

        buf = C_Loc(neighbor)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        i = 1
        do e = 1, mesh % n_elem
          associate(element => mesh % element(e))
            ! number of element neighbors
            nn = sum(element % vertex % n_neighbor) &
               + sum(element % edge   % n_neighbor) &
               + sum(element % face   % n_neighbor)
            allocate(element%neighbor(nn), source = neighbor(i:i+nn-1))
            i = i + nn
          end associate
        end do

      end if

      ! get element coordinates ................................................

      if (mesh%n_elem > 0) then

        call H5Dopen_f(group_id, name_mec, data_id, err)
        call H5Dget_space_f(data_id, space_id, err)
        call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)
        call H5Sclose_f(space_id, err)

        nc = dims(1)
        allocate(xc(nc))

        buf = C_Loc(xc)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Dclose_f(data_id, err)
        call H5Tclose_f(type_id, err)

        i = 1
        do e = 1, mesh % n_elem
          associate(geometry => mesh % element(e) % geometry)
            po = geometry % po
            np = po + 1
            nc = np * np * np * 3
            allocate(geometry % x_e(0:po,0:po,0:po,3), &
                     source = reshape(xc(i:i+nc-1), [np,np,np,3]))
            i = i + nc
          end associate
        end do

      end if

      ! complete mesh ..........................................................

      if (mesh%n_elem > 0) then

        call mesh % IdentifyEdges()
        call mesh % BuildFaces()
        call mesh % BuildBoundaryFaces()
        call mesh % BuildLinks()
        call mesh % BuildGhosts()
      ! call mesh % BuildCuboids()            ! retain read-in values
      ! call mesh % IdentifyRanks()           ! retain read-in values
        call mesh % BuildMapToParent()
        call mesh % BuildMapToChild()

      end if

      ! release resources ......................................................

      if (allocated(neighbor)) deallocate(neighbor)
      if (allocated(xc))       deallocate(xc)

    end if

    ! initialize new empty partitions ..........................................

    call attrib % Bcast(0, comm)

    if (group_id == H5I_INVALID_HID_F) then
      call mesh % Init_Mesh_3D(attrib, comm)
    end if

    ! build intracommunicator .................................................

    call mesh % BuildCommunicator()

    ! clean-up ................................................................

    if (allocated(attrib%boundary))  deallocate(attrib%boundary)
    if (allocated(mesh_dim))         deallocate(mesh_dim)

  end subroutine ReadHDF5_G

  !=============================================================================

end submodule MP_ReadHDF5
