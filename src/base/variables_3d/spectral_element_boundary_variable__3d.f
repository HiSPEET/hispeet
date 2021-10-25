!> summary:  3D spectral element boundary variable
!> author:   Joerg Stiller
!> date:     2021/10/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Boundary_Variable__3D
  use Kind_Parameters, only: RNP
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: SpectralElementBoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> Type for keeping the spectral element data associated with one boundary

  type SpectralElementBoundaryData_3D
    real(RNP), contiguous,  pointer :: val(:,:,:,:) !< accessable values
    real(RNP), allocatable, private :: mem(:,:,:,:) !< memory allocated to val
  end type SpectralElementBoundaryData_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable
  !>
  !> Base type for keeping a set of variables on the boundaries of a spectral
  !> element mesh. The mesh, standard operators and metric coefficients are
  !> provided by the pointer `sem`, which can be shared with other entities.
  !> The variables are accessed through
  !>
  !>     bnd(1:nb) % val(0:po,0:po,1:nf,1:nc),
  !>
  !> where
  !>
  !>   - `nb` is the number of boundaries,
  !>   - `po` the polynomial order,
  !>   - `nf` the number of faces, and
  !>   - `nc` the number of components.
  !>
  !> The first three of these dimensions correspond to the spectral element
  !> mesh `sem`, which implies that `nb` is generally different for each
  !> boundary, whereas `po` is constant.
  !> In the current implementation, `nc` is identical for all boundaries.
  !>
  !> Depending on the creation of the boundary variable, the values `val` in
  !> `bnd` can be stored in its own `mem` component or refer to the `mem`
  !> component of another instance.

  type SpectralElementBoundaryVariable_3D
    class(SpectralElementMesh_3D), pointer :: sem
    class(SpectralElementBoundaryData_3D), allocatable :: bnd(:)
  contains
    procedure :: Init_SpectralElementBoundaryVariable_3D => Init_SEBV_nc0
    procedure :: GetSlice
  end type SpectralElementBoundaryVariable_3D

  ! constructor
  interface SpectralElementBoundaryVariable_3D
    module procedure New_SEBV_nc0
  end interface

contains

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable constructor

  function New_SEBV_nc0(sem, nc) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    integer, intent(in) :: nc  !< number of components
    type(SpectralElementBoundaryVariable_3D) :: this

    call Init_SEBV_nc0(this, sem, nc)

  end function New_SEBV_nc0

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable initialization

  subroutine Init_SEBV_nc0(this, sem, nc)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    integer, intent(in) :: nc  !< number of components

    integer :: b

    associate(mesh => sem % mesh, po => sem % std_op % po)

      if (allocated(this % bnd)) then
        deallocate(this % bnd)
      end if

      this % sem => sem

      allocate(SpectralElementBoundaryData_3D :: this % bnd(mesh % n_bound))

      do b = 1, mesh % n_bound
        allocate(this % bnd(b) % mem(0:po, 0:po, mesh%boundary(b)%n_face, nc))
        this % bnd(b) % val(0:,0:,1:,1:) => this % bnd(b) % mem
      end do

    end associate

  end subroutine Init_SEBV_nc0

  !-----------------------------------------------------------------------------
  !> Create a new spectral element boundary variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this%bnd%val` if `copy` is false
  !> or absent. Otherwise they are stored in fresh memory, i.e. `slice%bnd%mem`.
  !> In an OpenMP parallel section the routine is executed only by the master
  !> thread.

  subroutine GetSlice(this, slice, first, last, copy)
    class(SpectralElementBoundaryVariable_3D), target, intent(in)    :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    integer :: b, c1, c2, nb
    logical :: copy_

    !$omp master

    ! preliminaries ...........................................................

    nb = this % sem % mesh % n_bound

    c1 = first
    c2 = last
    do b = 1, nb
      c1 = max(min(c1, size(this%bnd(b)%val, 4)), 1)
      c2 = max(min(c2, size(this%bnd(b)%val, 4)), 0)
    end do

    if (c1 /= first .or. c2 /= last) then
      call Error( 'GetSlice', 'section exceeds bounds',    &
                  'Spectral_Element_Boundary_Variable__3D' )
    end if

    if (allocated(slice % bnd)) then
      if (size(slice % bnd) /= nb) deallocate(slice % bnd)
    end if

    if (.not. allocated(slice % bnd)) then
      allocate(slice % bnd(nb))
    end if

    if (present(copy)) then
      copy_ = copy
    else
      copy_ = .false.
    end if

    ! slice ....................................................................

    do b = 1, nb

      slice % bnd(b) % val => null()

      if (copy_) then
        slice % bnd(b) % mem = this % bnd(b) % val(:,:,:,c1:c2)
        slice % bnd(b) % val(0:,0:,1:,1:) => slice % bnd(b) % mem
      end if

      if (.not. associated(this % bnd(b) % val)) then
        slice % bnd(b) % val(0:,0:,1:,1:) => this % bnd(b) % val(:,:,:,c1:c2)
        if (allocated(slice % bnd(b) % mem)) then
          deallocate(slice % bnd(b) % mem)
        end if
      end if

    end do

    !$omp end master

  end subroutine GetSlice

  !=============================================================================

end module Spectral_Element_Boundary_Variable__3D
