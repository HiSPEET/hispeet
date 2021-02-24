!> summary:  Methods for transferring linked element data of type integer(IXL)
!> author:   Joerg Stiller
!> date:     2021/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Element_Transfer_Buffer_3d) SM_Integer_IXL
  use Element_Transfer_Buffer_3d
  implicit none

   integer, parameter :: IK = IXL

contains

  !=============================================================================
  ! ElementTransferBuffer3d: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single integer variable

  module function New_TransferBuffer_IXL_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer3d) :: this
    type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,               optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer3d_Shared(this, mesh, np, nl)

    this % nc = 1

    allocate(integer(IK) :: this % master % buf(size(this % master % node)))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_IXL_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of integer variables

  module function New_TransferBuffer_IXL_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer3d) :: this
    type(Mesh3d_Partition),            intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,                 optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer3d_Shared(this, mesh, np, nl)

    this % nc = nc

    allocate(integer(IK) :: this % master % buf(size(this % master % node) * nc))
    allocate(integer(IK) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_IXL_A

  !=============================================================================
  ! ElementTransferBuffer3d: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer scalar

  module subroutine Transfer_IXL_S(this, mesh, v, tag)
    class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
    type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                         intent(in) :: tag   !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IXL_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- integer array

  module subroutine Transfer_IXL_A(this, mesh, v, tag)
    class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
    type(Mesh3d_Partition),            intent(in) :: mesh  !< mesh partition
    integer(IK), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                           intent(in) :: tag   !< message tag

    call Transfer_IX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_IXL_A

  !=============================================================================
  ! ElementTransferBuffer3d: merge procedures

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

  module subroutine Merge_IXL_S(this, v, alpha, beta)
    class(ElementTransferBuffer3d) , intent(inout) :: this
    integer(IK), dimension(:,:,:,:), intent(inout) :: v     !< mesh variable
    integer,               optional, intent(in)    :: alpha !< coeff of v  [0]
    integer,               optional, intent(in)    :: beta  !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IXL_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IXL_S

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

  module subroutine Merge_IXL_A(this, v, alpha, beta)
    class(ElementTransferBuffer3d)   , intent(inout) :: this
    integer(IK), dimension(:,:,:,:,:), intent(inout) :: v     !< mesh variable
    integer,                 optional, intent(in)    :: alpha !< coeff of v  [0]
    integer,                 optional, intent(in)    :: beta  !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_IXL_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_IX(this, v, alpha, beta)

  end subroutine Merge_IXL_A

  !=============================================================================
  ! Templated procedures Transfer_IX and Merge_IX

# include "element_transfer_buffer_3d__sm_integer-tmpl.f"

  !=============================================================================

end submodule SM_Integer_IXL
