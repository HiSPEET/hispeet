program HDF5_Test
  use Kind_Parameters
  USE ISO_C_BINDING   ! required for C_Loc
  use HDF5_Binding
  implicit none

  ! simplified element type ....................................................

  type ElementVertex
    integer      :: id         = -1 !< local mesh vertex ID
    integer(IXS) :: n_neighbor =  0 !< number of neighbor elements
  end type ElementVertex

  type ElementNeighbor
    integer      :: id        = -1 !< local ID of neighbor element
    integer(IXS) :: component = -1 !< coupled component
  end type ElementNeighbor

  type ElementGeometry
    integer                :: po           !< polynomial order
    real(RNP), allocatable :: x_e(:,:,:,:) !< element Lobatto points
    real(RNP)              :: x_c(0:3,1:3) !< corresponding cuboid
  end type ElementGeometry

  type MeshElement
    integer :: id     = -1                            !< local element ID
    logical :: frozen = .false.                       !< switch for elements
    type(ElementVertex)  :: vertex(8)                 !< vertex data
    type(ElementNeighbor), allocatable :: neighbor(:) !< neighbor data
    type(ElementGeometry) :: geometry                 !< geometry data
  end type MeshElement

  type(MeshElement) :: element_orig(2)
  type(MeshElement) :: element_read(2)

  ! HDF5 types .................................................................

  integer(hid_t) :: H5T_ElementVertex   = -1
  integer(hid_t) :: H5T_ElementGeometry = -1
  integer(hid_t) :: H5T_MeshElement     = -1

  ! auxiliary variables ........................................................

  integer :: err

  ! HDF5 initialization ........................................................

  call H5open_f(err)

  call Init_HDF5_Binding()
  print '(A,I0)', 'H5T_INTEGER         =  ', H5T_INTEGER
  print '(A,I0)', 'H5T_INTEGER_IXS     =  ', H5T_INTEGER_IXS
  print '(A,I0)', 'H5T_INTEGER_IXL     =  ', H5T_INTEGER_IXL
  print '(A,I0)', 'H5T_REAL_RSP        =  ', H5T_REAL_RSP
  print '(A,I0)', 'H5T_REAL_RDP        =  ', H5T_REAL_RDP
  print '(A,I0)', 'H5T_REAL_RHP        =  ', H5T_REAL_RHP
  print '(A,I0)', 'H5T_REAL_RNP        =  ', H5T_REAL_RNP
  print '(A,I0)', 'H5T_CHARACTER       =  ', H5T_CHARACTER
  print '(A,I0)', 'H5T_LOGICAL         =  ', H5T_LOGICAL

  call Init_H5T_ElementVertex()
  print '(A,I0)', 'H5T_ElementVertex   =  ', H5T_ElementVertex

  call Init_H5T_ElementGeometry()
  print '(A,I0)', 'H5T_ElementGeometry =  ', H5T_ElementGeometry

  call Init_H5T_Element()
  print '(A,I0)', 'H5T_MeshElement     =  ', H5T_MeshElement

  ! create dummy data ..........................................................

  ! TBD

  ! save data ..................................................................

  ! TBD

  ! read data ..................................................................

  ! TBD

  ! check data .................................................................

  ! TBD

  call H5close_f(err)

contains

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementVertex

  subroutine Init_H5T_ElementVertex()

    type(ElementVertex), target :: vertex(2)
    integer(size_t) :: offset
    integer :: err

    ! create HDF5 type
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementVertex, err)

    ! insert id
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%id))
    call H5Tinsert_f(H5T_ElementVertex, 'id', offset, H5T_INTEGER, err)

    ! insert n_neighbor
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%n_neighbor))
    call H5Tinsert_f(H5T_ElementVertex, 'n_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

  end subroutine Init_H5T_ElementVertex

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementGeometry

  subroutine Init_H5T_ElementGeometry()

    type(ElementGeometry), target :: geo(2)
    integer(size_t)  :: offset
    integer          :: rank    = 2
    integer(hsize_t) :: dims(2) = [ 4, 3 ]
    integer(hid_t)   :: tid
    integer :: err

    ! create HDF5 type
    offset = H5offsetof(C_Loc(geo(1)), C_Loc(geo(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementGeometry, err)

    ! insert po
    offset = H5offsetof(C_Loc(geo(1)), C_Loc(geo(1)%po))
    call H5Tinsert_f(H5T_ElementGeometry, 'po', offset, H5T_INTEGER, err)

    ! insert x_c
    call H5Tarray_create_f(H5T_REAL_RNP, rank, dims, tid, err)
    offset = H5offsetof(C_Loc(geo(1)), C_Loc(geo(1)%x_c(0,1)))
    call H5Tinsert_f(H5T_ElementGeometry, 'x_c', offset, tid, err)

  end subroutine Init_H5T_ElementGeometry

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_MeshElement

  subroutine Init_H5T_Element()

    type(MeshElement), target :: element(2)
    integer(size_t)  :: offset
    integer          :: rank = 1
    integer(hsize_t) :: dims(1)
    integer(hid_t)   :: tid
    integer :: err

    ! create HDF5 type
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_MeshElement, err)

    ! insert id
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%id))
    call H5Tinsert_f(H5T_MeshElement, 'id', offset, H5T_INTEGER, err)

    ! insert frozen
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%frozen))
    call H5Tinsert_f(H5T_MeshElement, 'frozen', offset, H5T_LOGICAL, err)

    ! insert vertex
    dims(1) = size(element(1)%vertex)
    call H5Tarray_create_f(H5T_ElementVertex, rank, dims, tid, err)
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%vertex(1)))
    call H5Tinsert_f(H5T_MeshElement, 'vertex', offset, tid, err)

    ! insert geometry
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%geometry))
    call H5Tinsert_f(H5T_MeshElement, 'geometry', offset, &
                     H5T_ElementGeometry, err)

  end subroutine Init_H5T_Element

  !=============================================================================

end program HDF5_Test
