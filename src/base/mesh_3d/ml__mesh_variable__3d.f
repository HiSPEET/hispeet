!> summary:  3D multilevel mesh variable
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Variable__3D
  use Kind_Parameters
  use ML__Mesh_Operators__3D
  implicit none
  private

  public :: ML_MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D mesh variable

  type MeshVariable_3D
    real(RNP), allocatable :: var(:,:,:,:,:) ! (0:po,0:po,0:po,1:ne,1:nc)
  end type MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh variable

  type ML_MeshVariable_3D
    type(MeshVariable_3D), allocatable :: level(:)
  contains
    procedure :: Init_ML_MeshVariable_3D
  end type ML_MeshVariable_3D

  ! constructor
  interface ML_MeshVariable_3D
    procedure New_ML_MeshVariable_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New multilevel variable with nc components from multilevel SE mesh

  function New_ML_MeshVariable_3D(ml_op, nc) result(this)
    class(ML_MeshOperators_3D), intent(in) :: ml_op
    integer,                    intent(in) :: nc
    type(ML_MeshVariable_3D) :: this

    call Init_ML_MeshVariable_3D(this, ml_op, nc)

  end function New_ML_MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> Initialize multilevel variable with nc components from multilevel SE mesh

  subroutine Init_ML_MeshVariable_3D(this, ml_op, nc)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_MeshOperators_3D), intent(in)    :: ml_op
    integer,                    intent(in)    :: nc

    integer :: l, ne, po

    allocate(this % level( size(ml_op%sem) ))

    do l = 1, size(this%level)
      po = ml_op % sem(l) % std_op % po
      ne = ml_op % sem(l) % mesh % n_elem
      allocate(this%level(l)%var(0:po,0:po,0:po,ne,nc))
    end do

  end subroutine Init_ML_MeshVariable_3D

  !=============================================================================

end module ML__Mesh_Variable__3D
