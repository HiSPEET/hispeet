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
    procedure :: Init => Init_ML_BoundaryVariable_3D
  end type ML_BoundaryVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Initialize multilevel boundary variable with nc components
  !>
  !> The optional argument `l_top` can be passed to restrict the top level to a
  !> value lower than `size(ml_op%sem)`

  subroutine Init_ML_BoundaryVariable_3D(this, ml_op, nc, l_top)
    class(ML_BoundaryVariable_3D), intent(inout) :: this
    class(ML_MeshOperators_3D),    intent(in)    :: ml_op
    integer,                       intent(in)    :: nc
    integer,             optional, intent(in)    :: l_top

    integer :: b, l, l_top_

    if (allocated(this%level)) deallocate(this%level)

    if (present(l_top)) then
      l_top_ = min(size(ml_op%sem), l_top)
    else
      l_top_ = size(ml_op%sem)
    end if

    allocate(this%level( l_top_ ))

    do l = 1, size(this%level)
      allocate(this % level(l) % var( ml_op % sem(l) % mesh % n_bound ))
      do b = 1, ml_op % sem(l) % mesh % n_bound
        call this % level(l) % var(b) %                    &
                 Init( ml_op % sem(l) % mesh % boundary(b) &
                     , ml_op % sem(l) % std_op % po        &
                     , nc                                  )
      end do
    end do

  end subroutine Init_ML_BoundaryVariable_3D

  !=============================================================================

end module ML__Boundary_Variable__3D
