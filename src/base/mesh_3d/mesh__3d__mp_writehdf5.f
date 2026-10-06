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

!> summary:  Writing a mesh partition to HDF5
!> author:   Joerg Stiller, Erik Pfister
!> date:     2024/06/05
!===============================================================================

submodule(Mesh__3D) MP_WriteHDF5
  use, intrinsic :: ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Write mesh partition into HDF5 file

  module subroutine WriteHDF5_F(mesh, file)
    class(Mesh_3D),   intent(in) :: mesh  !< mesh partition
    character(len=*), intent(in) :: file  !< name of HDF5 file

    integer(HID_T)    :: file_id, group_id
    integer           :: err
    character(len=80) :: tag
    character(len=:), allocatable :: file_pr

    ! append process rank to file name
    write(tag,'(I0)') mesh % proc
    file_pr = trim(file)//'_'//trim(tag)//'.h5'

    ! create HDF5 file and group
    call H5Fcreate_f(file_pr, H5F_ACC_TRUNC_F, file_id, err)
    call H5Gcreate_f(file_id, '/mesh', group_id, err)

    ! write mesh partition
    call WriteHDF5_G(mesh, group_id)

    ! close HDF5 file and group
    call H5Gclose_f(group_id, err)
    call H5Fclose_f(file_id, err)

  end subroutine WriteHDF5_F

  !-----------------------------------------------------------------------------
  !> Write mesh partition into given HDF5 group

  module subroutine WriteHDF5_G(mesh, group_id)
    class(Mesh_3D), intent(in) :: mesh     !< mesh partition
    integer(HID_T), intent(in) :: group_id !< ID of related HDF5 group

    ! internal data ............................................................

    ! HDF5 datatypes
    integer(HID_T) :: H5T_MeshAttributes_3D
    integer(HID_T) :: H5T_MeshBoundaryAttributes_3D
    integer(HID_T) :: H5T_MeshElement_3D
    integer(HID_T) :: H5T_MeshElementNeighbor_3D

    ! HDF5 dataspace IDs
    integer(HID_T) :: space_mp    ! space ID of mesh part
    integer(HID_T) :: space_ma    ! space ID of mesh attributes
    integer(HID_T) :: space_mba   ! space ID of mesh boundary attributes
    integer(HID_T) :: space_md    ! space ID of mesh dimensions
    integer(HID_T) :: space_me    ! space ID of mesh elements
    integer(HID_T) :: space_men   ! space ID of mesh element neighbor data
    integer(HID_T) :: space_mec   ! space ID of mesh element coordinates

    ! HDF5 dataspace dimensions
    integer(HSIZE_T) :: dim_mp (1)  ! dimension of mesh part
    integer(HSIZE_T) :: dim_ma (1)  ! dimension of mesh attributes
    integer(HSIZE_T) :: dim_mba(1)  ! dimension of mesh boundary attributes
    integer(HSIZE_T) :: dim_md (1)  ! dimension of mesh dimensions
    integer(HSIZE_T) :: dim_me (1)  ! dimension of mesh elements
    integer(HSIZE_T) :: dim_men(1)  ! dimension of mesh element neighbor data
    integer(HSIZE_T) :: dim_mec(1)  ! dimension of mesh element coordinates

    ! HDF5 dataset names
    character(len=*), parameter :: name_mp  = 'mesh_part'
    character(len=*), parameter :: name_ma  = 'mesh_attributes'
    character(len=*), parameter :: name_mba = 'mesh_boundary_attributes'
    character(len=*), parameter :: name_md  = 'mesh_dimensions'
    character(len=*), parameter :: name_me  = 'mesh_elements'
    character(len=*), parameter :: name_men = 'mesh_element_neighbors'
    character(len=*), parameter :: name_mec = 'mesh_element_coordinates'

    ! HDF5 dataset IDs
    integer(HID_T) :: data_mp    ! dataset ID of mesh part
    integer(HID_T) :: data_ma    ! dataset ID of mesh attributes
    integer(HID_T) :: data_mba   ! dataset ID of mesh boundary attributes
    integer(HID_T) :: data_md    ! dataset ID of mesh dimensions
    integer(HID_T) :: data_me    ! dataset ID of mesh elements
    integer(HID_T) :: data_men   ! dataset ID of mesh element neighbor data
    integer(HID_T) :: data_mec   ! dataset ID of mesh element coordinates

    class(Mesh_3D), allocatable, target :: copy_mesh

    type(MeshAttributes_3D),                      target :: attrib
    type(MeshBoundaryAttributes_3D), allocatable, target :: attrib_bound(:)
    type(MeshElementNeighbor_3D),    allocatable, target :: neighbor(:)
    real(RNP),                       allocatable, target :: xc(:)
    integer,                         allocatable, target :: mesh_dim(:)

    integer, allocatable :: nn_elem(:)
    integer, allocatable :: nc_elem(:)

    integer :: e, i, j, nc, nn
    integer :: err

    ! preliminaries ............................................................

    ! safeguard
    call Init_HDF5_Binding()

    ! auxiliary data ...........................................................

    copy_mesh = mesh

    ! mesh and mesh boundary attributes
    attrib = MeshAttributes_3D(copy_mesh)
    call move_alloc(attrib%boundary, attrib_bound)

    ! mesh dimensions and max valencies
    mesh_dim = [ copy_mesh % n_vert        &
               , copy_mesh % n_edge        &
               , copy_mesh % n_face        &
               , copy_mesh % n_elem        &
               , copy_mesh % n_elem_active &
               , copy_mesh % n_elem_frozen &
               , copy_mesh % n_cluster     &
               , copy_mesh % n_ghost       &
               , copy_mesh % n_link        &
               , copy_mesh % n_child       &
               , copy_mesh % n_parent      &
               , copy_mesh % n_elem_1      &
               , copy_mesh % n_elem_2      &
               , copy_mesh % n_elem_3      &
               , copy_mesh % max_vert_val  &
               , copy_mesh % max_edge_val  ]

    ! count number of neighbors and coordinates
    allocate(nn_elem(copy_mesh%n_elem), nc_elem(copy_mesh%n_elem))
    nn = 0
    nc = 0
    do e = 1, copy_mesh % n_elem
      nn_elem(e) = size(copy_mesh % element(e) % neighbor)
      nc_elem(e) = size(copy_mesh % element(e) % geometry % x_e)
      nn = nn + nn_elem(e)
      nc = nc + nc_elem(e)
    end do

    ! collect neighbor data
    allocate(neighbor(nn))
    i = 1
    do e = 1, copy_mesh % n_elem
      j = i + nn_elem(e) - 1
      neighbor(i:j) = copy_mesh % element(e) % neighbor
      deallocate(copy_mesh % element(e) % neighbor)
      i = i + nn_elem(e)
    end do

    ! collect element coordinates
    allocate(xc(nc))
    i = 1
    do e = 1, copy_mesh % n_elem
      j = i + nc_elem(e) - 1
      xc(i:j) = reshape(copy_mesh % element(e) % geometry % x_e, [nc_elem(e)])
      deallocate(copy_mesh % element(e) % geometry % x_e)
      i = i + nc_elem(e)
    end do

    ! datatypes ................................................................

    ! get datatypes
    call Get_H5T_MeshElement_3D(H5T_MeshElement_3D)
    call Get_H5T_MeshElementNeighbor_3D(H5T_MeshElementNeighbor_3D)
    call Get_H5T_MeshAttributes_3D(H5T_MeshAttributes_3D)
    call Get_H5T_MeshBoundaryAttributes_3D(H5T_MeshBoundaryAttributes_3D)

    ! dataspaces ...............................................................

    ! dimensions
    dim_mp  = 1
    dim_ma  = 1
    dim_mba = copy_mesh % n_bound
    dim_md  = size(mesh_dim)
    dim_me  = copy_mesh % n_elem
    dim_men = nn
    dim_mec = nc

    ! spaces
    call H5Screate_simple_f(1, dim_mp , space_mp , err)
    call H5Screate_simple_f(1, dim_ma , space_ma , err)
    call H5Screate_simple_f(1, dim_mba, space_mba, err)
    call H5Screate_simple_f(1, dim_md , space_md , err)
    call H5Screate_simple_f(1, dim_me , space_me , err)
    call H5Screate_simple_f(1, dim_men, space_men, err)
    call H5Screate_simple_f(1, dim_mec, space_mec, err)

    ! datasets .................................................................

    ! create mesh part dataset
    call H5Dcreate_f(group_id, name_mp, H5T_INTEGER, space_mp, data_mp, err)

    ! write mesh part dataset
    call H5Dwrite_f(data_mp, H5T_INTEGER, C_Loc(copy_mesh%part), err)

    ! create mesh attribute dataset
    call H5Dcreate_f(group_id, name_ma, H5T_MeshAttributes_3D, &
                     space_ma, data_ma, err)

    ! write mesh attribute dataset
    call H5Dwrite_f(data_ma, H5T_MeshAttributes_3D, C_Loc(attrib), err)

    ! create mesh boundary attribute dataset
    call H5Dcreate_f(group_id, name_mba, H5T_MeshBoundaryAttributes_3D, &
                     space_mba, data_mba, err)

    ! write mesh boundary attribute dataset (may be empty!)
    call H5Dwrite_f(data_mba, H5T_MeshBoundaryAttributes_3D, &
                    C_Loc(attrib_bound), err)

    ! create mesh dimensions dataset
    call H5Dcreate_f(group_id, name_md, H5T_INTEGER, space_md, data_md, err)

    ! write mesh dimensions dataset
    call H5Dwrite_f(data_md, H5T_INTEGER, C_Loc(mesh_dim), err)

    ! create mesh element dataset
    call H5Dcreate_f(group_id, name_me, H5T_MeshElement_3D, &
                     space_me, data_me, err)

    ! write mesh element dataset (may be empty!)
    call H5Dwrite_f(data_me, H5T_MeshElement_3D, C_Loc(copy_mesh%element), err)

    ! create mesh element neighbor dataset
    call H5Dcreate_f(group_id, name_men, H5T_MeshElementNeighbor_3D, &
                     space_men, data_men, err)

    ! write mesh element neighbor dataset (may be empty!)
    call H5Dwrite_f(data_men, H5T_MeshElementNeighbor_3D, C_Loc(neighbor), err)

    ! create mesh element coordinates dataset
    call H5Dcreate_f(group_id, name_mec, H5T_REAL_RNP, space_mec, data_mec, err)

    ! write mesh element coordinates dataset (may be empty!)
    call H5Dwrite_f(data_mec, H5T_REAL_RNP, C_Loc(xc), err)

    ! release resources ........................................................

    deallocate(nn_elem, nc_elem, neighbor, xc)
    deallocate(attrib_bound, mesh_dim)
    deallocate(copy_mesh)

    call H5Sclose_f(space_mp , err)
    call H5Sclose_f(space_ma , err)
    call H5Sclose_f(space_mba, err)
    call H5Sclose_f(space_md , err)
    call H5Sclose_f(space_me , err)
    call H5Sclose_f(space_men, err)
    call H5Sclose_f(space_mec, err)

    call H5Dclose_f(data_mp , err)
    call H5Dclose_f(data_ma , err)
    call H5Dclose_f(data_mba, err)
    call H5Dclose_f(data_md , err)
    call H5Dclose_f(data_me , err)
    call H5Dclose_f(data_men, err)
    call H5Dclose_f(data_mec, err)

  end subroutine WriteHDF5_G

  !=============================================================================

end submodule MP_WriteHDF5
