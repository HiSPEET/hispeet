!> summary:  Generation of MPI datatype for mesh elements
!> author:   Joerg Stiller, Erik Pfister, Moritz Kreuseler
!> date:     2024/06/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_Element__3D) SM_HDF5
  use, intrinsic ::  ISO_C_Binding
  use HDF5_Binding
  implicit none

  !-----------------------------------------------------------------------------
  !> HDF5 datatype for mesh elements

  integer(HID_T) :: H5T_Element = -1

  !-----------------------------------------------------------------------------
  !> HDF5 datatype for mesh element neighbor data

  integer(HID_T) :: H5T_ElementNeighbor = -1

  !-----------------------------------------------------------------------------
  ! Auxiliary HDF5 datatypes for element components

  integer(HID_T) :: H5T_ElementVertex     = -1
  integer(HID_T) :: H5T_ElementEdge       = -1
  integer(HID_T) :: H5T_ElementFace       = -1
  integer(HID_T) :: H5T_ElementAdaptation = -1
  integer(HID_T) :: H5T_ElementGeometry   = -1

contains

  !-----------------------------------------------------------------------------
  !> Get HDF5 datatype for essential static components of MeshElement_3D

  module subroutine Get_H5T_MeshElement_3D( H5T_MeshElement_3D )
    integer(HID_T), intent(out) :: H5T_MeshElement_3D

    !$omp master

    if (H5T_Element < 0) then
      call Init_H5T_Element()
    end if

    H5T_MeshElement_3D = H5T_Element

    !$omp end master

  end subroutine Get_H5T_MeshElement_3D

  !-----------------------------------------------------------------------------
  !> Initialization of HDF5_MeshElement

  subroutine Init_H5T_Element()

    type(MeshElement_3D), target :: element(2)
    integer(SIZE_T)  :: offset
    integer(HSIZE_T) :: dims(1)
    integer(HID_T)   :: tid
    integer :: err

    ! preliminaries ............................................................

    ! initialize element datatype
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_Element, err)

    ! initialize component datatypes
    call Init_H5T_ElementVertex()
    call Init_H5T_ElementEdge()
    call Init_H5T_ElementFace()
    call Init_H5T_ElementAdaptation()
    call Init_H5T_ElementGeometry()

    ! insert components ........................................................

    ! insert vertex
    dims   = size(element(1)%vertex)
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%vertex(1)))
    call H5Tarray_create_f(H5T_ElementVertex, 1, dims, tid, err)
    call H5Tinsert_f(H5T_Element, 'vertex', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert edge
    dims   = size(element(1)%edge)
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%edge(1)))
    call H5Tarray_create_f(H5T_ElementEdge, 1, dims, tid, err)
    call H5Tinsert_f(H5T_Element, 'edge', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert face
    dims   = size(element(1)%face)
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%face(1)))
    call H5Tarray_create_f(H5T_ElementFace, 1, dims, tid, err)
    call H5Tinsert_f(H5T_Element, 'face', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert adapdation
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%adaptation))
    call H5Tinsert_f(H5T_Element, 'adaptation', offset, H5T_ElementAdaptation, &
                     err)

    ! insert geometry
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%geometry))
    call H5Tinsert_f(H5T_Element, 'geometry', offset, H5T_ElementGeometry, err)

    ! release datatypes ........................................................

    call H5Tclose_f(H5T_ElementVertex, err)
    call H5Tclose_f(H5T_ElementEdge, err)
    call H5Tclose_f(H5T_ElementFace, err)
    call H5Tclose_f(H5T_ElementAdaptation, err)
    call H5Tclose_f(H5T_ElementGeometry, err)

  end subroutine Init_H5T_Element

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementVertex

  subroutine Init_H5T_ElementVertex()

    type(MeshElementVertex_3D), target :: vertex(2)
    integer(SIZE_T) :: offset
    integer :: err

     ! create datatype
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementVertex, err)

    ! insert id
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%id))
    call H5Tinsert_f(H5T_ElementVertex, 'id', offset, H5T_INTEGER, err)

    ! insert n_neighbor
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%n_neighbor))
    call H5Tinsert_f(H5T_ElementVertex, 'n_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert i_neighbor
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%i_neighbor))
    call H5Tinsert_f(H5T_ElementVertex, 'i_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert rank
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%rank))
    call H5Tinsert_f(H5T_ElementVertex, 'rank', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert val
    offset = H5offsetof(C_Loc(vertex(1)), C_Loc(vertex(1)%val))
    call H5Tinsert_f(H5T_ElementVertex, 'val', offset, &
                     H5T_INTEGER_IXS, err)

  end subroutine Init_H5T_ElementVertex

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementEdge

  subroutine Init_H5T_ElementEdge()

    type(MeshElementEdge_3D), target :: edge(2)
    integer(SIZE_T) :: offset
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementEdge, err)

    ! insert id
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%id))
    call H5Tinsert_f(H5T_ElementEdge, 'id', offset, H5T_INTEGER, err)

    ! insert orientation
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%orientation))
    call H5Tinsert_f(H5T_ElementEdge, 'orientation', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert n_neighbor
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%n_neighbor))
    call H5Tinsert_f(H5T_ElementEdge, 'n_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert i_neighbor
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%i_neighbor))
    call H5Tinsert_f(H5T_ElementEdge, 'i_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert rank
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%rank))
    call H5Tinsert_f(H5T_ElementEdge, 'rank', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert val
    offset = H5offsetof(C_Loc(edge(1)), C_Loc(edge(1)%val))
    call H5Tinsert_f(H5T_ElementEdge, 'val', offset, &
                     H5T_INTEGER_IXS, err)

  end subroutine Init_H5T_ElementEdge

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementFace

  subroutine Init_H5T_ElementFace()

    type(MeshElementFace_3D), target :: face(2)
    integer(SIZE_T) :: offset
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementFace, err)

    ! insert id
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%id))
    call H5Tinsert_f(H5T_ElementFace, 'id', offset, H5T_INTEGER, err)

    ! insert boundary
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%boundary))
    call H5Tinsert_f(H5T_ElementFace, 'boundary', offset, H5T_INTEGER, err)

    ! insert normal
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%normal))
    call H5Tinsert_f(H5T_ElementFace, 'normal', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert rotation
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%rotation))
    call H5Tinsert_f(H5T_ElementFace, 'rotation', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert n_neighbor
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%n_neighbor))
    call H5Tinsert_f(H5T_ElementFace, 'n_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert i_neighbor
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%i_neighbor))
    call H5Tinsert_f(H5T_ElementFace, 'i_neighbor', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert rank
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%rank))
    call H5Tinsert_f(H5T_ElementFace, 'rank', offset, &
                     H5T_INTEGER_IXS, err)

    ! insert val
    offset = H5offsetof(C_Loc(face(1)), C_Loc(face(1)%val))
    call H5Tinsert_f(H5T_ElementFace, 'val', offset, &
                     H5T_INTEGER_IXS, err)


  end subroutine Init_H5T_ElementFace

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementAdaptation

  subroutine Init_H5T_ElementAdaptation()

    type(MeshElementAdaptation_3D), target :: adaptation(2)
    integer(SIZE_T) :: offset
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementAdaptation, err)

    ! insert parent_proc
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%parent_proc))
    call H5Tinsert_f(H5T_ElementAdaptation, 'parent_proc', offset, &
                     H5T_INTEGER, err)

    ! insert parent_id
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%parent_id))
    call H5Tinsert_f(H5T_ElementAdaptation, 'parent_id', offset, &
                     H5T_INTEGER, err)

    ! insert refinement
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%refinement))
    call H5Tinsert_f(H5T_ElementAdaptation, 'refinement', offset, &
                     H5T_INTEGER, err)

    ! create child_proc
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%child_proc))
    call H5Tinsert_f(H5T_ElementAdaptation, 'child_proc', offset, &
                     H5T_INTEGER, err)

    ! insert sublevels
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%sublevels))
    call H5Tinsert_f(H5T_ElementAdaptation, 'sublevels', offset, &
                     H5T_INTEGER, err)

    ! insert mark
    offset = H5offsetof(C_Loc(adaptation(1)), C_Loc(adaptation(1)%mark))
    call H5Tinsert_f(H5T_ElementAdaptation, 'mark', offset, H5T_INTEGER, err)

  end subroutine Init_H5T_ElementAdaptation

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementGeometry

  subroutine Init_H5T_ElementGeometry()

    type(MeshElementGeometry_3D), target :: geometry(2)
    integer(SIZE_T) :: offset
    integer(HID_T)  :: tid
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementGeometry, err)

    ! insert po
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%po))
    call H5Tinsert_f(H5T_ElementGeometry, 'po', offset, H5T_INTEGER, err)

    ! insert x_c
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%x_c(0,1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 2, int([4,3], HSIZE_T), tid, err)
    call H5Tinsert_f(H5T_ElementGeometry, 'x_c', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert dx_m
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%dx_m(1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 1, int([6], HSIZE_T), tid, err)
    call H5Tinsert_f(H5T_ElementGeometry, 'dx_m', offset, tid, err)
    call H5Tclose_f(tid, err)

  end subroutine Init_H5T_ElementGeometry

  !-----------------------------------------------------------------------------
  !> Get HDF5 datatype of MeshElementNeighbor_3D

  module subroutine Get_H5T_MeshElementNeighbor_3D( H5T_MeshElementNeighbor_3D )
    integer(HID_T), intent(out) :: H5T_MeshElementNeighbor_3D

    !$omp master

    if (H5T_ElementNeighbor < 0) then
      call Init_H5T_ElementNeighbor()
    end if

    H5T_MeshElementNeighbor_3D = H5T_ElementNeighbor

    !$omp end master

  end subroutine Get_H5T_MeshElementNeighbor_3D

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementNeighbor

  subroutine Init_H5T_ElementNeighbor()

    type(MeshElementNeighbor_3D), target :: neighbor(2)
    integer(SIZE_T) :: offset
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(neighbor(1)), C_Loc(neighbor(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementNeighbor, err)

    ! insert id
    offset = H5offsetof(C_Loc(neighbor(1)), C_Loc(neighbor(1)%id))
    call H5Tinsert_f(H5T_ElementNeighbor, 'id', offset, H5T_INTEGER, err)

    ! insert part
    offset = H5offsetof(C_Loc(neighbor(1)), C_Loc(neighbor(1)%part))
    call H5Tinsert_f(H5T_ElementNeighbor, 'part', offset, H5T_INTEGER, err)

    ! insert component
    offset = H5offsetof(C_Loc(neighbor(1)), C_Loc(neighbor(1)%component))
    call H5Tinsert_f(H5T_ElementNeighbor, 'component', offset &
                   , H5T_INTEGER_IXS, err)

    ! insert orientation
    offset = H5offsetof(C_Loc(neighbor(1)), C_Loc(neighbor(1)%orientation))
    call H5Tinsert_f(H5T_ElementNeighbor, 'orientation', offset &
                   , H5T_INTEGER_IXS, err)

  end subroutine Init_H5T_ElementNeighbor

  !=============================================================================

end submodule SM_HDF5
