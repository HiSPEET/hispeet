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
    real(RNP), contiguous, pointer :: val(:,:,:,:,:) => null()
    logical, private :: is_original = .false.
  contains
    final :: Delete_MeshVariable_3D
  end type MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh variable

  type ML_MeshVariable_3D
    type(MeshVariable_3D), allocatable :: level(:) !< variable per level
    character(len=:),      allocatable :: name(:)  !< component names
  contains
    procedure :: Init_ML_MeshVariable_3D
    procedure :: ExportVTK
  end type ML_MeshVariable_3D

  ! constructor
  interface ML_MeshVariable_3D
    procedure New_ML_MeshVariable_3D
  end interface

  !=============================================================================
  ! module procedures

  interface

    !---------------------------------------------------------------------------
    !> Export to VTK

    module subroutine ExportVTK(this, ml_op, file, mode)
      class(ML_MeshVariable_3D),  intent(in) :: this  !< multilevel variable
      class(ML_MeshOperators_3D), intent(in) :: ml_op !< spectral element OPs
      character(len=*),           intent(in) :: file  !< export file base name
      integer,                    intent(in) :: mode  !< export mode {1,2,3}
    end subroutine ExportVTK

  end interface

contains

  !-----------------------------------------------------------------------------
  !> New multilevel variable with nc components from multilevel SE mesh

  function New_ML_MeshVariable_3D(ml_op, nc, name) result(this)
    class(ML_MeshOperators_3D), intent(in) :: ml_op
    integer,                    intent(in) :: nc
    character(len=*), optional, intent(in) :: name(nc)
    type(ML_MeshVariable_3D) :: this

    call Init_ML_MeshVariable_3D(this, ml_op, nc, name)

  end function New_ML_MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> Initialize multilevel variable with nc components from multilevel SE mesh

  subroutine Init_ML_MeshVariable_3D(this, ml_op, nc, name)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_MeshOperators_3D), intent(in)    :: ml_op
    integer,                    intent(in)    :: nc
    character(len=*), optional, intent(in)    :: name(nc)

    integer :: i, l, ne, po

    allocate(this % level( size(ml_op%sem) ))

    do l = 1, size(this%level)
      po = ml_op % sem(l) % std_op % po
      ne = ml_op % sem(l) % mesh % n_elem
      allocate(this%level(l)%val(0:po,0:po,0:po,ne,nc))
    end do
    this % level % is_original = .true.

    if (present(name)) then
      this % name = name
    else
      l = 2 + int(log10(dble(nc)))
      allocate(character(len=l) :: this % name(nc))
      do i = 1, nc
        write(this%name(i), '(A,I0)') 'v', i
      end do
    end if

  end subroutine Init_ML_MeshVariable_3D

  !=============================================================================
  ! Finalization

  !-----------------------------------------------------------------------------
  !>  Finalization of MeshVariable_3D

  subroutine Delete_MeshVariable_3D(this)
    type(MeshVariable_3D), intent(inout) :: this

    if (this%is_original .and. associated(this%val)) then
      deallocate(this%val)
    end if

    this % val => null()
    this % is_original = .false.

  end subroutine Delete_MeshVariable_3D

  !=============================================================================

end module ML__Mesh_Variable__3D
