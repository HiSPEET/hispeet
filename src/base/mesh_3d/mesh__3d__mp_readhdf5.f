!> summary:  Reading a mesh partition from HDF5
!> author:   Joerg Stiller, Erik Pfister
!> date:     2024/06/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
!### CHECK
print '(99G0)', '$p',rank,', ReadHDF5_F: file_pr = "',file_pr,'"'
!### CHECK END

    if (exists) then
      ! open HDF5 file and group for reading
      call H5Fopen_f(file_pr, H5F_ACC_RDWR_F, file_id, err)
      call H5Gopen_f(file_id, 'mesh', group_id, err)
    else
      group_id = H5I_INVALID_HID_F
    end if
!### CHECK
print '(99G0)', '$p',rank,', ReadHDF5_F: group_id = ',group_id
!### CHECK END

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
    integer(hid_t), intent(in) :: group_id !< ID of related HDF5 group
    type(MPI_Comm), intent(in) :: comm     !< MPI "world" communicator

    ! internal data ............................................................

    ! dynamical data
    type(MeshAttributes_3D), target :: attrib
    type(MeshBoundaryAttributes_3D), allocatable, target :: attrib_boundary(:)
    type(MeshElementNeighbor_3D), allocatable, target :: neighbor(:)
    real(RNP), allocatable, target :: xc(:)

    ! HDF5 datatype, dataspace, dimensions and buffer adress pointer
    integer(hid_t)   :: type_id
    integer(hid_t)   :: space_id
    integer(hid_t)   :: data_id
    integer(hsize_t) :: dims(1), maxdims(1)
    type(C_Ptr)      :: buf

    ! HDF5 dataset names
    character(len=*), parameter :: name_mp  = 'mesh_part'
    character(len=*), parameter :: name_ma  = 'mesh_attributes'
    character(len=*), parameter :: name_mba = 'mesh_boundary_attributes'
    character(len=*), parameter :: name_me  = 'mesh_elements'
    character(len=*), parameter :: name_men = 'mesh_element_neighbors'
    character(len=*), parameter :: name_mec = 'mesh_element_coordinates'

    integer :: e, i, nc, nn, np, po, rank
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
!### CHECK
print '(99G0)', '$p',rank,', ReadHDF5_G: data_id = ',data_id,', err = ', err
!### CHECK END
      call H5Dget_type_f(data_id, type_id, err)
!### CHECK
print '(99G0)', '$p',rank,', ReadHDF5_G: data_id = ',type_id,', err = ', err
!### CHECK END
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)
!### CHECK
print '(99G0)', '$p',rank,', ReadHDF5_G: part = ',mesh%part
!### CHECK END

      ! get mesh attributes ....................................................

      buf = C_Loc(attrib)
      call H5Dopen_f(group_id, name_ma, data_id, err)
      call H5Dget_type_f(data_id, type_id, err)
      call H5Dread_f(data_id, type_id, buf, err)
      call H5Tclose_f(type_id, err)
      call H5Dclose_f(data_id, err)

      ! get mesh boundary attributes ...........................................

      allocate(attrib % boundary(attrib%n_bound))

      if (attrib%n_bound > 0) then
        buf = C_Loc(attrib%boundary)
        call H5Dopen_f(group_id, name_mba, data_id, err)
        call H5Dget_type_f(data_id, type_id, err)
        call H5Dread_f(data_id, type_id, buf, err)
        call H5Tclose_f(type_id, err)
        call H5Dclose_f(data_id, err)
      end if

      ! initialize mesh ........................................................

      call mesh % Init_Mesh_3D(attrib, comm)

      ! get elements ...........................................................

      call H5Dopen_f(group_id, name_me, data_id, err)
      call H5Dget_space_f(data_id, space_id, err)
      call H5Sget_simple_extent_dims_f(space_id, dims, maxdims, err)
      call H5Sclose_f(space_id, err)

      mesh % n_elem = dims(1)
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
        call H5Tclose_f(type_id, err)
        call H5Dclose_f(data_id, err)

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

      ! complete mesh ..........................................................

      call mesh % IdentifyEdges()
      call mesh % BuildFaces()
      call mesh % BuildBoundaryFaces()
      call mesh % BuildLinks()
      call mesh % BuildGhosts()
    ! call mesh % BuildCuboids()        ! retain read-in values
    ! call mesh % IdentifyRanks()       ! retain read-in values
      call mesh % BuildMapToParent()
      call mesh % BuildMapToChild()

      ! release resources ......................................................

      deallocate(neighbor, xc)

    end if

    ! initialize new empty partitions ..........................................

    call attrib % Bcast(0, comm)

    if (group_id == H5I_INVALID_HID_F) then
      call mesh % Init_Mesh_3D(attrib, comm)
    end if

    ! build intracommunicator .................................................

    call mesh % BuildCommunicator()

    ! clean-up ................................................................

    if (allocated(neighbor)) deallocate(neighbor)
    if (allocated(xc))       deallocate(xc)

    if (allocated(attrib%boundary)) then
      deallocate(attrib%boundary)
    end if

  end subroutine ReadHDF5_G

  !=============================================================================

end submodule MP_ReadHDF5
