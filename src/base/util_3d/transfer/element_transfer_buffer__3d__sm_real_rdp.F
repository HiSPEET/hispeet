!> summary:  Methods for transferring linked element data of type real(RDP)
!> author:   Joerg Stiller
!> date:     2021/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Element_Transfer_Buffer__3D) SM_Real_RDP
  use Element_Transfer_Buffer__3D
  implicit none

   integer, parameter :: RK = RDP

contains

  !=============================================================================
  ! ElementTransferBuffer_3D: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single real variable

  module function New_TransferBuffer_RDP_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    type(MeshPartition_3D),        intent(in) :: mesh  !< mesh partition
    real(RDP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,             optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = 1

    allocate(real(RDP) :: this % master % buf(size(this % master % node)))
    allocate(real(RDP) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_RDP_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of real variables

  module function New_TransferBuffer_RDP_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer_3D) :: this
    type(MeshPartition_3D),          intent(in) :: mesh  !< mesh partition
    real(RDP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,               optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer_Shared_3D(this, mesh, np, nl)

    this % nc = nc

    allocate(real(RDP) :: this % master % buf(size(this % master % node) * nc))
    allocate(real(RDP) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_RDP_A

  !=============================================================================
  ! ElementTransferBuffer_3D: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real scalar

  module subroutine Transfer_RDP_S(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    type(MeshPartition_3D),        intent(in) :: mesh  !< mesh partition
    real(RDP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                       intent(in) :: tag   !< message tag

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RDP_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real array

  module subroutine Transfer_RDP_A(this, mesh, v, tag)
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    type(MeshPartition_3D),          intent(in) :: mesh  !< mesh partition
    real(RDP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                         intent(in) :: tag   !< message tag

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RDP_A

  !=============================================================================
  ! ElementTransferBuffer_3D: merge procedures

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real scalar
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

  module subroutine Merge_RDP_S(this, v, alpha, beta)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    real(RDP),   dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RDP),             optional, intent(in)    :: alpha  !< coeff of v  [0]
    real(RDP),             optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RDP_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RDP_S

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real array
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

  module subroutine Merge_RDP_A(this, v, alpha, beta)
    class(ElementTransferBuffer_3D), intent(inout) :: this
    real(RDP), dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RDP),             optional, intent(in)    :: alpha  !< coeff of v  [0]
    real(RDP),             optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RDP_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer__3D'    )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RDP_A

  !=============================================================================
  ! Templated procedures Transfer_RX and Merge_RX

# include "element_transfer_buffer__3d__sm_real-tmpl.f"

  !=============================================================================

end submodule SM_Real_RDP
