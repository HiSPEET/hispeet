!> summary:  3D multilevel boundary variable
!> author:   Joerg Stiller
!> date:     2024/08/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Boundary_Variable__3D
  use Boundary_Variable__3D
  use ML__Mesh_Operators__3D
  implicit none
  private

  public :: ML_BoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> Array comprising variables on all boundaries on one mesh level

  type BoundaryVariableArray_3D
    type(BoundaryVariable_3D), allocatable :: var(:) !< variables per boundary
  end type BoundaryVariableArray_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel boundary variable

  type ML_BoundaryVariable_3D
    type(BoundaryVariableArray_3D), allocatable :: level(:)
      !< boundary variables per level
  contains
    procedure :: Init_ML_BoundaryVariable_3D
  end type ML_BoundaryVariable_3D

  ! constructor
  interface ML_BoundaryVariable_3D
    procedure New_ML_BoundaryVariable_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New multilevel boundary variable with nc components

  function New_ML_BoundaryVariable_3D(ml_op, nc) result(this)
    class(ML_MeshOperators_3D), intent(in) :: ml_op
    integer,                    intent(in) :: nc
    type(ML_BoundaryVariable_3D) :: this

    call Init_ML_BoundaryVariable_3D(this, ml_op, nc)

  end function New_ML_BoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> Initialize multilevel boundary variable with nc components

  subroutine Init_ML_BoundaryVariable_3D(this, ml_op, nc)
    class(ML_BoundaryVariable_3D), intent(inout) :: this
    class(ML_MeshOperators_3D),    intent(in)    :: ml_op
    integer,                       intent(in)    :: nc

    integer :: b, l

    allocate(this % level( size(ml_op % sem) ))

    do l = 1, size(this%level)
      allocate(this % level(l) % var( ml_op % sem(l) % mesh % n_bound ))
      do b = 1, ml_op % sem(l) % mesh % n_bound
        this % level(l) % var(b) = BoundaryVariable_3D                       &
                                       ( ml_op % sem(l) % mesh % boundary(b) &
                                       , ml_op % sem(l) % std_op % po        &
                                       , nc                                  )
      end do
    end do

  end subroutine Init_ML_BoundaryVariable_3D

  !=============================================================================

end module ML__Boundary_Variable__3D
