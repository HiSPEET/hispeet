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

  integer(hid_t) :: H5T_Element = -1

  !-----------------------------------------------------------------------------
  ! HDF5 datatypes element components and auxiliary data

  integer(hid_t) :: H5T_ElementVertex     = -1
  integer(hid_t) :: H5T_ElementEdge       = -1
  integer(hid_t) :: H5T_ElementFace       = -1
  integer(hid_t) :: H5T_ElementAdaptation = -1
  integer(hid_t) :: H5T_ElementGeometry   = -1

  logical :: initialized = .false.

contains

  !-----------------------------------------------------------------------------
  !> Get HDF5 datatype for essential static components of MeshElement_3D

  module subroutine Get_H5T_MeshElement_3D( H5T_MeshElement_3D )
    integer(hid_t), intent(out) :: H5T_MeshElement_3D

    !$omp master

    if (.not. initialized) then
      call Init_H5T_Element()
    end if

    H5T_MeshElement_3D = H5T_Element

    !$omp end master

  end subroutine Get_H5T_MeshElement_3D

  !-----------------------------------------------------------------------------
  !> Initialization of HDF5_MeshElement

  subroutine Init_H5T_Element()

    type(MeshElement_3D), target :: element(2)
    integer(size_t)  :: offset
    integer(hsize_t) :: dims(1)
    integer(hid_t)   :: tid
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

    ! insert face

    ! insert adapdation

    ! insert geometry
    offset = H5offsetof(C_Loc(element(1)), C_Loc(element(1)%geometry))
    call H5Tinsert_f(H5T_Element, 'geometry', offset, H5T_ElementGeometry, err)

    ! release datatypes ........................................................

    call H5Tclose_f(H5T_ElementVertex, err)
!   call H5Tclose_f(H5T_ElementEdge, err)
!   call H5Tclose_f(H5T_ElementFace, err)
!   call H5Tclose_f(H5T_ElementAdaptation, err)
!   call H5Tclose_f(H5T_ElementGeometry, err)

  end subroutine Init_H5T_Element

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementVertex

  subroutine Init_H5T_ElementVertex()

    type(MeshElementVertex_3D), target :: vertex(2)
    integer(size_t) :: offset
    integer :: err

  end subroutine Init_H5T_ElementVertex

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementEdge

  subroutine Init_H5T_ElementEdge()

    type(MeshElementEdge_3D), target :: edge(2)
    integer(size_t) :: offset
    integer :: err

  end subroutine Init_H5T_ElementEdge

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementFace

  subroutine Init_H5T_ElementFace()

    type(MeshElementFace_3D), target :: face(2)
    integer(size_t) :: offset
    integer :: err

  end subroutine Init_H5T_ElementFace

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementAdaptation

  subroutine Init_H5T_ElementAdaptation()

    type(MeshElementAdaptation_3D), target :: adaptation(2)
    integer(size_t) :: offset
    integer :: err

  end subroutine Init_H5T_ElementAdaptation

  !-----------------------------------------------------------------------------
  !> Initialization of H5T_ElementGeometry

  subroutine Init_H5T_ElementGeometry()

    type(MeshElementGeometry_3D), target :: geometry(2)
    integer(size_t) :: offset
    integer(hid_t)  :: tid
    integer :: err

    ! create datatype
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(2)))
    call H5Tcreate_f(H5T_COMPOUND_F, offset, H5T_ElementGeometry, err)

    ! insert po
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%po))
    call H5Tinsert_f(H5T_ElementGeometry, 'po', offset, H5T_INTEGER, err)

    ! insert x_c
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%x_c(0,1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 2, int([4,3], hsize_t), tid, err)
    call H5Tinsert_f(H5T_ElementGeometry, 'x_c', offset, tid, err)
    call H5Tclose_f(tid, err)

    ! insert dx_m
    offset = H5offsetof(C_Loc(geometry(1)), C_Loc(geometry(1)%dx_m(1)))
    call H5Tarray_create_f(H5T_REAL_RNP, 1, int([6], hsize_t), tid, err)
    call H5Tinsert_f(H5T_ElementGeometry, 'dx_m', offset, tid, err)
    call H5Tclose_f(tid, err)

  end subroutine Init_H5T_ElementGeometry

  !=============================================================================

end submodule SM_HDF5
