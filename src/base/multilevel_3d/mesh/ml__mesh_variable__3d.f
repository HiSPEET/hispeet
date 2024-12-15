!> summary:  3D multilevel mesh variable
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Variable__3D
  use Kind_Parameters
  use Mesh_Variable__3D
  use ML__Mesh_Operators__3D
  implicit none
  private

  public :: ML_MeshVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D multilevel mesh variable

  type ML_MeshVariable_3D
    type(MeshVariable_3D), allocatable :: level(:) !< variable per level
    character(len=:),      allocatable :: name(:)  !< component names
  contains
    procedure :: Init => Init_ML_MeshVariable_3D
    procedure :: GetSlice
    procedure :: ExportVTK
  end type ML_MeshVariable_3D

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
  !> Initialize multilevel variable with nc components from multilevel SE mesh

  subroutine Init_ML_MeshVariable_3D(this, ml_op, nc, name)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_MeshOperators_3D), intent(in)    :: ml_op
    integer,                    intent(in)    :: nc
    character(len=*), optional, intent(in)    :: name(nc)

    integer :: i, l

    allocate(this % level( size(ml_op%sem) ))

    do l = 1, size(this%level)
      call this % level(l) % Init( mesh = ml_op % sem(l) % mesh        &
                                 , po   = ml_op % sem(l) % std_op % po &
                                 , nc   = nc                           )
    end do

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

  !-----------------------------------------------------------------------------
  !> Create a new multilevel mesh variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % val` if `copy` is false
  !> or absent. Otherwise they are stored in fresh memory, i.e. `slice % val`.
  !> In an OpenMP parallel section the routine is executed only by the master
  !> thread.

  subroutine GetSlice(this, slice, first, last, copy)
    class(ML_MeshVariable_3D), intent(in)    :: this
    class(ML_MeshVariable_3D), intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    integer :: l

    !$omp master
    allocate(slice % level(size(this % level)))
    slice % name = this % name(first:last)

    do l = 1, size(this%level)
      call this % level(l) % GetSlice(slice % level(l), first, last, copy)
    end do
    !$omp end master

  end subroutine GetSlice

  !=============================================================================

end module ML__Mesh_Variable__3D
