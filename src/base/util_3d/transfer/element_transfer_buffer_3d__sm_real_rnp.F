!> summary:  Methods for transferring linked element data of type real(RNP)
!> author:   Joerg Stiller
!> date:     2021/02/20
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Element_Transfer_Buffer_3d) SM_Real_RNP
  use Element_Transfer_Buffer_3d
  implicit none

   integer, parameter :: RK = RNP

contains

  !=============================================================================
  ! ElementTransferBuffer3d: constructor functions

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for a single real variable

  module function New_TransferBuffer_RNP_S(mesh, v, nl) result(this)
    type(ElementTransferBuffer3d) :: this
    type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,             optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3)

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)

    call Init_ElementTransferBuffer3d_Shared(this, mesh, np, nl)

    this % nc = 1

    allocate(real(RNP) :: this % master % buf(size(this % master % node)))
    allocate(real(RNP) :: this % ghost  % buf(size(this % ghost  % node)))

  end function New_TransferBuffer_RNP_S

  !-----------------------------------------------------------------------------
  !> New element transfer buffer for an array of real variables

  module function New_TransferBuffer_RNP_A(mesh, v, nl) result(this)
    type(ElementTransferBuffer3d) :: this
    type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variables
    integer,               optional, intent(in) :: nl(3) !< number of layers

    integer :: np(3), nc

    np(1) = size(v,1)
    np(2) = size(v,2)
    np(3) = size(v,3)
    nc    = size(v,5)

    call Init_ElementTransferBuffer3d_Shared(this, mesh, np, nl)

    this % nc = nc

    allocate(real(RNP) :: this % master % buf(size(this % master % node) * nc))
    allocate(real(RNP) :: this % ghost  % buf(size(this % ghost  % node) * nc))

  end function New_TransferBuffer_RNP_A

  !=============================================================================
  ! ElementTransferBuffer3d: transfer procedures

  !----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real scalar

  module subroutine Transfer_RNP_S(this, mesh, v, tag)
    class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
    type(Mesh3d_Partition),        intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                       intent(in) :: tag   !< message tag

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RNP_S

  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real array

  module subroutine Transfer_RNP_A(this, mesh, v, tag)
    class(ElementTransferBuffer3d), asynchronous, intent(inout) :: this
    type(Mesh3d_Partition),          intent(in) :: mesh  !< mesh partition
    real(RNP), dimension(:,:,:,:,:), intent(in) :: v     !< mesh variable
    integer,                         intent(in) :: tag   !< message tag

    call Transfer_RX(this, mesh, v, size(v,4), tag)

  end subroutine Transfer_RNP_A

  !=============================================================================
  ! ElementTransferBuffer3d: merge procedures

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

  module subroutine Merge_RNP_S(this, v, alpha, beta)
    class(ElementTransferBuffer3d), intent(inout) :: this
    real(RNP),  dimension(:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RNP),            optional, intent(in)    :: alpha  !< coeff of v  [0]
    real(RNP),            optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RNP_S'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RNP_S

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

  module subroutine Merge_RNP_A(this, v, alpha, beta)
    class(ElementTransferBuffer3d) , intent(inout) :: this
    real(RNP), dimension(:,:,:,:,:), intent(inout) :: v      !< mesh variable
    real(RNP),             optional, intent(in)    :: alpha  !< coeff of v  [0]
    real(RNP),             optional, intent(in)    :: beta   !< coeff of vb [1]

    !$omp master
    if (size(v,4) /= this%ne + this%ng) then
      call Error( 'Merge_RNP_A'                    &
                , 'size(v,4) /= this%ne + this%ng' &
                , 'Element_Transfer_Buffer_3d'     )
    end if
    !$omp end master

    call Merge_RX(this, v, alpha, beta)

  end subroutine Merge_RNP_A

  !=============================================================================
  ! Templated procedures Transfer_RX and Merge_RX

# include "element_transfer_buffer_3d__sm_real-tmpl.f"

  !=============================================================================

end submodule SM_Real_RNP
