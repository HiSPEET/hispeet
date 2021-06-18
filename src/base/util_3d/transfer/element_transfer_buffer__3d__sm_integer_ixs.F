!> summary:  Methods for transferring linked element data of type integer(IXS)
!> author:   Joerg Stiller
!> date:     2021/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Element_Transfer_Buffer__3D) SM_Integer_IXS
  use Element_Transfer_Buffer__3D
  implicit none

   integer, parameter :: IK = IXS

contains

  !=============================================================================
  ! ElementTransferBuffer_3D: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single integer variable

  module function New_TransferBuffer_IXS_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    type(MeshPartition_3D),          intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,               optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = 1

    allocate(integer(IK) :: this % master % buf(size(this % master % node)))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_IXS_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of integer variables

  module function New_TransferBuffer_IXS_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    type(MeshPartition_3D),            intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,                 optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = nc

    allocate(integer(IK) :: this % master % buf(size(this % master % node) * nc))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_IXS_A

  !=============================================================================
  ! ElementTransferBuffer_3D: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer scalar

  module subroutine Transfer_IXS_S(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    type(MeshPartition_3D),          intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                         intent(in) :: tag   !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IXS_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer array

  module subroutine Transfer_IXS_A(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    type(MeshPartition_3D),            intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                           intent(in) :: tag   !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IXS_A

  !=============================================================================
  ! ElementTransferBuffer_3D: merge procedures

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer scalar
  !>
  !> Denoting the buffer with `vb`, the following operation will be executed:
  !>
  !>     v  =  alpha * v  +  beta * vb
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng)`, where
  !>
  !>    *  `np` is the number of points per direction, as in this%np
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements

  module subroutine Merge_IXS_S(this, v, alpha, beta)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    integer(IK), dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
    integer,               optional, intent(in)    :: alpha !< coeff of v  [0]
    integer,               optional, intent(in)    :: beta  !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IXS_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IXS_S

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- integer array
  !>
  !> Denoting the buffer with `vb`, the following operation will be executed:
  !>
  !>      v  =  alpha * v  +  beta * vb
  !>
  !> As this operation affects the ghost elements, the mesh variable must
  !> be dimensioned as `v(np(1), np(2), np(3), ne+ng, nc)`, where
  !>
  !>    *  `np` is the number of points per direction, as in `this%np`
  !>    *  `ne` is the number of master elements
  !>    *  `ng` is the number of ghost elements
  !>    *  `nc` is the number of components, as in `this%nc`

  module subroutine Merge_IXS_A(this, v, alpha, beta)
    class(ElementTransferBuffer_3D)  , intent(inout) :: this
    integer(IK), dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
    integer,                 optional, intent(in)    :: alpha !< coeff of v  [0]
    integer,                 optional, intent(in)    :: beta  !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IXS_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IXS_A

  !=============================================================================
  ! Templated procedures Transfer_IX and Merge_IX

# include "element_transfer_buffer__3d__sm_integer-tmpl.f"

  !=============================================================================

end submodule SM_Integer_IXS
