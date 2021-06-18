!> summary:  Assembly of continuous Lobatto element contributions
!> author:   Joerg Stiller
!> date:     2018/10/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__Assembly_Operator

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use CART__Mesh_Partition
  use CART__Element_Transfer_Buffer

  implicit none
  private

  public :: AssemblyOperator

  !-----------------------------------------------------------------------------
  !> Type for assembling Lobatto element contributions
  !>
  !>
  !>
  !>
  !>
  !>

  type, extends(ElementTransferBuffer) :: AssemblyOperator
  contains
    generic,   public  :: StartAssembly  => StartAssembly_S
    procedure, private :: StartAssembly_S
    generic,   public  :: FinishAssembly => FinishAssembly_S
    procedure, private :: FinishAssembly_S
  end type AssemblyOperator

contains

!-------------------------------------------------------------------------------
!>

subroutine StartAssembly_S(this, mesh, v, tag)
  class(AssemblyOperator), asynchronous, intent(inout) :: this
  type(MeshPartition),           intent(in)    :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:), intent(in)    :: v     !< mesh variable
  integer,                       intent(in)    :: tag   !< message tag

  call this % ToGhost_Transfer(mesh, v, tag)

end subroutine StartAssembly_S

!-------------------------------------------------------------------------------
!>

subroutine FinishAssembly_S(this, mesh, v)
  class(AssemblyOperator),       intent(inout) :: this
  type(MeshPartition),           intent(in)    :: mesh  !< mesh partition
  real(RNP), dimension(:,:,:,:), intent(inout) :: v     !< mesh variable

  integer :: e, n

  n = size(v, 1)

  call this % ToGhost_Merge(v, alpha=ZERO, beta=ONE)

  !$omp do
  do e = 1, mesh%ne

    ! faces ....................................................................

    associate(face => mesh % element(e)%face)

      if (face(1) % neighbor > 0) then
        v(1, :, :, e) = v(1, :, :, e) + v(n, :, :, face(1) % neighbor)
      end if

      if (face(2) % neighbor > 0) then
        v(n, :, :, e) = v(n, :, :, e) + v(1, :, :, face(2) % neighbor)
      end if

      if (face(3) % neighbor > 0) then
        v(:, 1, :, e) = v(:, 1, :, e) + v(:, n, :, face(3) % neighbor)
      end if

      if (face(4) % neighbor > 0) then
        v(:, n, :, e) = v(:, n, :, e) + v(:, 1, :, face(4) % neighbor)
      end if

      if (face(5) % neighbor > 0) then
        v(:, :, 1, e) = v(:, :, 1, e) + v(:, :, n, face(5) % neighbor)
      end if

      if (face(6) % neighbor > 0) then
        v(:, :, n, e) = v(:, :, n, e) + v(:, :, 1, face(6) % neighbor)
      end if

    end associate

    ! edges ....................................................................

    associate(edge => mesh % element(e)%edge)

      if (edge( 1) % neighbor > 0) then
        v(:, 1, 1, e) = v(:, 1, 1, e) + v(:, n, n, edge( 1) % neighbor)
      end if

      if (edge( 2) % neighbor > 0) then
        v(:, n, 1, e) = v(:, n, 1, e) + v(:, 1, n, edge( 2) % neighbor)
      end if

      if (edge( 3) % neighbor > 0) then
        v(:, 1, n, e) = v(:, 1, n, e) + v(:, n, 1, edge( 3) % neighbor)
      end if

      if (edge( 4) % neighbor > 0) then
        v(:, n, n, e) = v(:, n, n, e) + v(:, 1, 1, edge( 4) % neighbor)
      end if

      if (edge( 5) % neighbor > 0) then
        v(1, :, 1, e) = v(1, :, 1, e) + v(n, :, n, edge( 5) % neighbor)
      end if

      if (edge( 6) % neighbor > 0) then
        v(n, :, 1, e) = v(n, :, 1, e) + v(1, :, n, edge( 6) % neighbor)
      end if

      if (edge( 7) % neighbor > 0) then
        v(1, :, n, e) = v(1, :, n, e) + v(n, :, 1, edge( 7) % neighbor)
      end if

      if (edge( 8) % neighbor > 0) then
        v(n, :, n, e) = v(n, :, n, e) + v(1, :, 1, edge( 8) % neighbor)
      end if

      if (edge( 9) % neighbor > 0) then
        v(1, 1, :, e) = v(1, 1, :, e) + v(n, n, :, edge( 9) % neighbor)
      end if

      if (edge(10) % neighbor > 0) then
        v(n, 1, :, e) = v(n, 1, :, e) + v(1, n, :, edge(10) % neighbor)
      end if

      if (edge(11) % neighbor > 0) then
        v(1, n, :, e) = v(1, n, :, e) + v(n, 1, :, edge(11) % neighbor)
      end if

      if (edge(12) % neighbor > 0) then
        v(n, n, :, e) = v(n, n, :, e) + v(1, 1, :, edge(12) % neighbor)
      end if

    end associate

    ! vertices .................................................................

    associate(vertex => mesh % element(e)%vertex)

      if (vertex(1) % neighbor > 0) then
        v(1, 1, 1, e) = v(1, 1, 1, e) + v(n, n, n, vertex(1) % neighbor)
      end if

      if (vertex(2) % neighbor > 0) then
        v(n, 1, 1, e) = v(n, 1, 1, e) + v(1, n, n, vertex(2) % neighbor)
      end if

      if (vertex(3) % neighbor > 0) then
        v(1, n, 1, e) = v(1, n, 1, e) + v(n, 1, n, vertex(3) % neighbor)
      end if

      if (vertex(4) % neighbor > 0) then
        v(n, n, 1, e) = v(n, n, 1, e) + v(1, 1, n, vertex(4) % neighbor)
      end if

      if (vertex(5) % neighbor > 0) then
        v(1, 1, n, e) = v(1, 1, n, e) + v(n, n, 1, vertex(5) % neighbor)
      end if

      if (vertex(6) % neighbor > 0) then
        v(n, 1, n, e) = v(n, 1, n, e) + v(1, n, 1, vertex(6) % neighbor)
      end if

      if (vertex(7) % neighbor > 0) then
        v(1, n, n, e) = v(1, n, n, e) + v(n, 1, 1, vertex(7) % neighbor)
      end if

      if (vertex(8) % neighbor > 0) then
        v(n, n, n, e) = v(n, n, n, e) + v(1, 1, 1, vertex(8) % neighbor)
      end if

    end associate

  end do

end subroutine FinishAssembly_S

!===============================================================================

end module CART__Assembly_Operator
