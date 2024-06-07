!> summary:  Reading a mesh partition from HDF5 group
!> author:   Joerg Stiller, Erik Pfister, Moritz Kreuseler
!> date:     2024/06/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh__3D) MP_ReadHDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Read mesh partition from given HDF5 group

  module subroutine ReadHDF5(mesh, group, comm)
    class(Mesh_3D), target, intent(inout) :: mesh  !< mesh partition
    integer(hid_t),         intent(in)    :: group !< ID of related HDF5 group
    type(MPI_Comm),         intent(in)    :: comm  !< "world" communicator

    ! internal data ............................................................

    ! dynamical data
    type(MeshAttributes_3D), target :: attrib ! mesh attributes
    real(RNP), allocatable,  target :: xc(:)       ! mesh element coordinates
    type(MeshElementNeighbor_3D), allocatable, target :: neighbor(:)

    ! HDF5 datatype, dataspace and dimensions
    integer(hid_t)   :: type_id
    integer(hid_t)   :: space_id
    integer(hid_t)   :: data_id
    integer(hsize_t) :: dims(1), maxdims(1)

    ! HDF5 dataset names
    character(len=*), parameter :: name_ma  = 'mesh_attributes'
    character(len=*), parameter :: name_mba = 'mesh_boundary_attributes'
    character(len=*), parameter :: name_me  = 'mesh_elements'
    character(len=*), parameter :: name_men = 'mesh_element_neighbors'
    character(len=*), parameter :: name_mec = 'mesh_element_coordinates'

    integer :: e, i, j, nc, nn, np, po
    integer :: err

    ! preliminaries ............................................................

    ! ???

    ! get mesh attributes ......................................................

    call H5Dopen_f(group, name_ma, data_id, err)
    call H5Dget_type_f(data_id, type_id, err)
    call H5Dread_f(data_id, type_id, C_Loc(attrib), err)
    call H5Tclose_f(type_id, err)
    call H5Dclose_f(data_id, err)

    ! get mesh boundary attributes .............................................

    allocate(attrib % boundary(attrib%n_bound))

    if (attrib%n_bound > 0) then
      call H5Dopen_f(group, name_mba, data_id, err)
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, C_Loc(attrib%boundary), err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)
    end if

    ! initialize mesh ..........................................................

    call mesh % Init_Mesh_3D(attrib, comm)

    ! get elements .............................................................

    call H5Dopen_f(group, name_me, data_id, err)
    call H5Sget_space_f(data_id, space_id, err)
    call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)

    mesh % n_elem = dims(1)
    allocate(mesh % element(mesh%n_elem))

    if (mesh%n_elem > 0) then
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, C_Loc(mesh%element), err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)
    end if

    ! get element neighbors ....................................................

    if (mesh%n_elem > 0) then

      call H5Dopen_f(group, name_men, data_id, err)
      call H5Sget_space_f(data_id, space_id, err)
      call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)

      nn = dims(1)
      allocate(neighbor(nn))

      if (nn > 0) then
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, C_Loc(neighbor), err)
        call H5Tclose_f(type_id, err)
        call H5Dclose_f(data_id, err)
      end if

      i = 1
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))
          ! number of element neighbors
          nn = sum(element % vertex % n_neighbor) &
             + sum(element % edge   % n_neighbor) &
             + sum(element % face   % n_neighbor)
          if (nn > 0) then
            allocate(element%neighbor(nn), source = neighbor(i:i+nn-1))
          else
            allocate(element%neighbor(0))
          end if
          i = i + nn
        end associate
      end do

    end if

    ! get element coordinates ..................................................

    if (mesh%n_elem > 0) then

      call H5Dopen_f(group, name_mec, data_id, err)
      call H5Sget_space_f(data_id, space_id, err)
      call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)

      nc = dims(1)
      allocate(xc(nc))

      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, C_Loc(xc), err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)

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

    ! complete mesh ............................................................

    call mesh % IdentifyEdges()
    call mesh % BuildFaces()
    call mesh % BuildBoundaryFaces()
    call mesh % BuildLinks()
    call mesh % BuildGhosts()
  ! call mesh % BuildCuboids()        ! retain read-in values
  ! call mesh % IdentifyRanks()       ! retain read-in values
    call mesh % BuildMapToParent()
    call mesh % BuildMapToChild()

    ! intracommunicator between active partitions
    call mesh % BuildCommunicator()

    ! release resources ........................................................

    deallocate(neighbor, xc)

    if (allocated(attrib%boundary)) then
      deallocate(attrib%boundary)
    end if

  end subroutine ReadHDF5

  !=============================================================================

end submodule MP_ReadHDF5
