!> summary:  Methods for transferring linked element data of type integer
!> author:   Joerg Stiller
!> date:     2021/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Element_Transfer_Buffer__3D) SM_Integer_IDK
  implicit none

  integer, parameter :: IK = kind(1)

contains

  !=============================================================================
  ! ElementTransferBuffer_3D: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single integer variable

  module function New_TransferBuffer_IDK_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    class(Mesh_3D),    intent(in) :: mesh       !< mesh partition
    integer(IK),       intent(in) :: v(:,:,:,:) !< mesh variable
    integer, optional, intent(in) :: nl(3)      !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = 1

    allocate(integer(IK) :: this % master % buf(size(this % master % node)))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_IDK_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of integer variables

  module function New_TransferBuffer_IDK_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    class(Mesh_3D),    intent(in) :: mesh         !< mesh partition
    integer(IK),       intent(in) :: v(:,:,:,:,:) !< mesh variables
    integer, optional, intent(in) :: nl(3)        !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = nc

    allocate(integer(IK) :: this % master % buf(size(this % master % node) * nc))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_IDK_A

  !=============================================================================
  ! ElementTransferBuffer_3D: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer scalar

  module subroutine Transfer_IDK_S(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh       !< mesh partition
    integer(IK),    intent(in) :: v(:,:,:,:) !< mesh variable
    integer,        intent(in) :: tag        !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IDK_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer array

  module subroutine Transfer_IDK_A(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh         !< mesh partition
    integer(IK),    intent(in) :: v(:,:,:,:,:) !< mesh variable
    integer,        intent(in) :: tag          !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IDK_A

  !=============================================================================
  ! ElementTransferBuffer_3D: merge procedures

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer scalar
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng)`, where
  !>
  !>    *  `np` is the number of points per direction, as in this%np
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements

  module subroutine Merge_IDK_S(this, v)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    integer(IK), intent(inout) :: v(:,:,:,:) !< mesh variable

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IDK_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_IX(this, v)

  end subroutine Merge_IDK_S

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer array
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng, nc)`, where
  !>
  !>    *  `np` is the number of points per direction, as in `this%np`
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements
  !>    *  `nc` is the number of components, as in `this%nc`

  module subroutine Merge_IDK_A(this, v)
    class(ElementTransferBuffer_3D)  , intent(inout) :: this
    integer(IK), intent(inout) :: v(:,:,:,:,:) !< mesh variable

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IDK_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_IX(this, v)

  end subroutine Merge_IDK_A

  !=============================================================================
  ! Templated procedures Transfer_IX and Merge_IX

# include "element_transfer_buffer__3d__sm_integer-tmpl.f"

  !=============================================================================

end submodule SM_Integer_IDK
