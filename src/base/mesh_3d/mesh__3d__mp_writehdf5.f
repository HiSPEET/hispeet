!> summary:  ...
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
  !> ...

  module subroutine WriteHDF5(this, group_id)
    class(Mesh_3D), intent(in) :: mesh
    integer(hid_t), intent(in) :: group_id ! ID of HDF5 group

    ! internal data ............................................................

    ! dynamical data buffers
    type(MeshAttributes_3D) :: attrib     ! mesh attributes
    integer, allocatable    :: nn_elem(:) ! number of element neighbors
    integer, allocatable    :: nc_elem(:) ! number of element coordinates
    real(RNP), allocatable  :: xc(:)      ! coordinates of element points
    type(MeshElementNeighbor_3D), allocatable :: neighbor(:) ! neighbor data

    ! HDF5 datatypes
    integer(hid_t) :: H5T_MeshElement_3D
    integer(hid_t) :: H5T_MeshElementNeighbor_3D
    integer(hid_t) :: H5T_MeshAttributes_3D
    integer(hid_t) :: H5T_MeshBoundaryAttributes_3D

    ! HDF5 datatype sizes
    integer(size_t) :: size_rnp
    integer(size_t) :: size_int
    integer(size_t) :: size_me
    integer(size_t) :: size_ma
    integer(size_t) :: size_mba

    integer :: e, i, j, nn, nc
    integer :: err

    ! auxiliary data ...........................................................

    ! mesh attributes
    attrib = MeshAttributes_3D(mesh)

    if (mesh % n_ele)

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

    ! get datatype sizes
    call H5Tget_size_f(H5T_INTEGER, size_int, err)
    call H5Tget_size_f(H5T_REAL_RNP, size_rnp, err)
    call H5Tget_size_f(H5T_MeshElement_3D, size_me, err)
    call H5Tget_size_f(H5T_MeshElementNeighbor_3D, size_men, err)
    call H5Tget_size_f(H5T_MeshAttributes_3D, size_ma, err)
    call H5Tget_size_f(H5T_MeshBoundaryAttributes_3D, size_mba, err)

! TBD: spaces
! TBD: datasets

  end subroutine WriteHDF5

  !=============================================================================

end submodule MP_WriteHDF5
