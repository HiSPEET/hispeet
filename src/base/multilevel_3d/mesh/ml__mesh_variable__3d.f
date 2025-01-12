!> summary:  3D multilevel mesh variable
!> author:   Joerg Stiller
!> date:     2024/06/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Mesh_Variable__3D
  use Kind_Parameters
  use Data_Exchange__3D
  use Mesh_Variable__3D
  use ML__Mesh__3D
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
    generic   :: Init => Init_ML_MeshVariable_3D__M, Init_ML_MeshVariable_3D__O
    procedure :: Init_ML_MeshVariable_3D__M, Init_ML_MeshVariable_3D__O
    procedure :: GetSlice
    procedure :: ExportVTK
    procedure :: FitAdapt
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

    !---------------------------------------------------------------------------
    !> Fit to adapted mesh
    module subroutine FitAdapt(this, ml_op, x_plan)
      class(ML_MeshVariable_3D),  intent(inout) :: this
      class(ML_MeshOperators_3D), intent(in)    :: ml_op     !< adapted operators
      class(DataExchangePlan_3D), intent(in)    :: x_plan(:) !< reassignment plan
    end subroutine FitAdapt

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Initialize variable from multilevel mesh
  !>
  !> Provides a variable with fixed polynomial order `po` and `nc` components
  !> fitting the given multilevel mesh

  subroutine Init_ML_MeshVariable_3D__M(this, ml_mesh, po, nc, name)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_Mesh_3D),          intent(in)    :: ml_mesh
    integer,                    intent(in)    :: po  !< polynomial order
    integer,                    intent(in)    :: nc  !< number of components
    character(len=*), optional, intent(in)    :: name(nc)

    integer :: i, l

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%name))  deallocate(this%name)

    allocate(this % level( size(ml_mesh%mesh) ))

    do l = 1, size(this%level)
      call this % level(l) % Init(ml_mesh % mesh(l), po, nc)
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

  end subroutine Init_ML_MeshVariable_3D__M

  !-----------------------------------------------------------------------------
  !> Initialize variable from multilevel mesh operators
  !>
  !> Provides a variable with `nc` components fitting the mesh and polynomial
  !> order of the given multilevel operators

  subroutine Init_ML_MeshVariable_3D__O(this, ml_op, nc, name)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_MeshOperators_3D), intent(in)    :: ml_op
    integer,                    intent(in)    :: nc
    character(len=*), optional, intent(in)    :: name(nc)

    integer :: i, l

    if (allocated(this%level)) deallocate(this%level)
    if (allocated(this%name )) deallocate(this%name)

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

  end subroutine Init_ML_MeshVariable_3D__O

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

    if (allocated(slice%level)) then
      deallocate(slice%level)
    end if

    allocate(slice % level(size(this % level)))
    slice % name = this % name(first:last)

    do l = 1, size(this%level)
      call this % level(l) % GetSlice(slice % level(l), first, last, copy)
    end do

  end subroutine GetSlice

  !=============================================================================

end module ML__Mesh_Variable__3D
