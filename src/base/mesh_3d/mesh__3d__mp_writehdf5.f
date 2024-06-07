!> summary:  Writing a mesh partition into HDF5 group
!> author:   Joerg Stiller, Erik Pfister, Moritz Kreuseler
!> date:     2024/06/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_WriteHDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Write mesh partition into given HDF5 group

  module subroutine WriteHDF5(mesh, group)
    class(Mesh_3D), target, intent(in) :: mesh  !< mesh partition
    integer(hid_t),         intent(in) :: group !< ID of related HDF5 group

    ! internal data ............................................................

    ! dynamical data
    integer, allocatable :: nn_elem(:)  ! number of element neighbors
    integer, allocatable :: nc_elem(:)  ! number of element coordinates
    type(MeshAttributes_3D), target :: attrib ! mesh attributes
    real(RNP), allocatable,  target :: xc(:)       ! mesh element coordinates
    type(MeshElementNeighbor_3D), allocatable, target :: neighbor(:)

    ! HDF5 datatypes
    integer(hid_t) :: H5T_MeshAttributes_3D
    integer(hid_t) :: H5T_MeshBoundaryAttributes_3D
    integer(hid_t) :: H5T_MeshElement_3D
    integer(hid_t) :: H5T_MeshElementNeighbor_3D

    ! HDF5 dataspace IDs
    integer(hid_t) :: space_ma    ! space ID of mesh attributes
    integer(hid_t) :: space_mba   ! space ID of mesh boundary attributes
    integer(hid_t) :: space_me    ! space ID of mesh elements
    integer(hid_t) :: space_men   ! space ID of mesh element neighbor data
    integer(hid_t) :: space_mec   ! space ID of mesh element coordinates

    ! HDF5 dataspace dimensions
    integer(hsize_t) :: dim_ma (1)  ! dimension of mesh attributes
    integer(hsize_t) :: dim_mba(1)  ! dimension of mesh boundary attributes
    integer(hsize_t) :: dim_me (1)  ! dimension of mesh elements
    integer(hsize_t) :: dim_men(1)  ! dimension of mesh element neighbor data
    integer(hsize_t) :: dim_mec(1)  ! dimension of mesh element coordinates

    ! HDF5 dataset names
    character(len=*), parameter :: name_ma  = 'mesh_attributes'
    character(len=*), parameter :: name_mba = 'mesh_boundary_attributes'
    character(len=*), parameter :: name_me  = 'mesh_elements'
    character(len=*), parameter :: name_men = 'mesh_element_neighbors'
    character(len=*), parameter :: name_mec = 'mesh_element_coordinates'

    ! HDF5 dataset IDs
    integer(hid_t) :: data_ma    ! dataset ID of mesh attributes
    integer(hid_t) :: data_mba   ! dataset ID of mesh boundary attributes
    integer(hid_t) :: data_me    ! dataset ID of mesh elements
    integer(hid_t) :: data_men   ! dataset ID of mesh element neighbor data
    integer(hid_t) :: data_mec   ! dataset ID of mesh element coordinates

    integer :: e, i, j, nc, nn
    integer :: err

    ! auxiliary data ...........................................................

    ! mesh attributes
    attrib = MeshAttributes_3D(mesh)

    ! count number of neighbors and coordinates
    allocate(nn_elem(mesh%n_elem), nc_elem(mesh%n_elem))
    nn = 0
    nc = 0
    do e = 1, mesh % n_elem
      nn_elem(e) = size(mesh % element(e) % neighbor)
      nc_elem(e) = size(mesh % element(e) % geometry % x_e)
      nn = nn + nn_elem(e)
      nc = nc + nc_elem(e)
    end do

    ! collect neighbor data
    allocate(neighbor(nn))
    i = 1
    do e = 1, mesh % n_elem
      j = i + nn_elem(e) - 1
      neighbor(i:j) = mesh % element(e) % neighbor
      i = i + nn_elem(e)
    end do

    ! collect element coordinates
    allocate(xc(nc))
    i = 1
    do e = 1, mesh % n_elem
      j = i + nc_elem(e) - 1
      xc(i:j) = reshape(mesh % element(e) % geometry % x_e, [nc_elem(e)])
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
    dim_ma  = 1
    dim_mba = mesh % n_bound
    dim_me  = mesh % n_elem
    dim_men = nn
    dim_mec = nc

    ! spaces
    call H5Screate_simple_f(1, dim_ma , space_ma , err)
    call H5Screate_simple_f(1, dim_mba, space_mba, err)
    call H5Screate_simple_f(1, dim_me , space_me , err)
    call H5Screate_simple_f(1, dim_men, space_men, err)
    call H5Screate_simple_f(1, dim_mec, space_mec, err)

    ! datasets .................................................................

    ! create mesh attribute dataset
    call H5Dcreate_f(group, name_ma, H5T_MeshAttributes_3D, &
                     space_ma, data_ma, err)

    ! write mesh attribute dataset
    call H5Dwrite_f(data_ma, H5T_MeshAttributes_3D, C_Loc(attrib), err)

    ! create mesh boundary attribute dataset
    call H5Dcreate_f(group, name_mba, H5T_MeshBoundaryAttributes_3D, &
                     space_mba, data_mba, err)

    ! write mesh boundary attribute dataset (may be empty!)
    call H5Dwrite_f(data_mba, H5T_MeshBoundaryAttributes_3D, &
                    C_Loc(attrib%boundary), err)

    ! create mesh element dataset
    call H5Dcreate_f(group, name_me, H5T_MeshElement_3D, space_me, data_me, err)

    ! write mesh element dataset (may be empty!)
    call H5Dwrite_f(data_me, H5T_MeshElement_3D, C_Loc(mesh%element), err)

    ! create mesh element neighbor dataset
    call H5Dcreate_f(group, name_men, H5T_MeshElementNeighbor_3D, &
                     space_men, data_men, err)

    ! write mesh element neighbor dataset (may be empty!)
    call H5Dwrite_f(data_men, H5T_MeshElementNeighbor_3D, C_Loc(neighbor), err)

    ! create mesh element neighbor dataset
    call H5Dcreate_f(group, name_mec, H5T_REAL_RNP, space_mec, data_mec, err)

    ! write mesh element neighbor dataset (may be empty!)
    call H5Dwrite_f(data_mec, H5T_REAL_RNP, C_Loc(xc), err)

    ! release resources ........................................................

    deallocate(nn_elem, nc_elem, neighbor, xc)

    if (allocated(attrib%boundary)) then
      deallocate(attrib%boundary)
    end if

    call H5Sclose_f(space_ma , err)
    call H5Sclose_f(space_mba, err)
    call H5Sclose_f(space_me , err)
    call H5Sclose_f(space_men, err)
    call H5Sclose_f(space_mec, err)

    call H5Dclose_f(data_ma , err)
    call H5Dclose_f(data_mba, err)
    call H5Dclose_f(data_me , err)
    call H5Dclose_f(data_men, err)
    call H5Dclose_f(data_mec, err)

  end subroutine WriteHDF5

  !=============================================================================

end submodule MP_WriteHDF5
