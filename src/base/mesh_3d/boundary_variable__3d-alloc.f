!> summary:  3D boundary variable
!> author:   Joerg Stiller
!> date:     2022/10/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   -  works
!>   -  based on the original variant with unsafe constructors removed
!===============================================================================

module Boundary_Variable__3D
  use Kind_Parameters,   only: RNP
  use Execution_Control, only: Error
  use Mesh_Boundary__3D
  use Spectral_Element_Mesh__3D
  implicit none
  private

  public :: BoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D boundary variable
  !>
  !> Base type for keeping a set of variables on the boundaries of a spectral
  !> element mesh. The boundary values are accessed via
  !>
  !>       bv(b) % val(0:po,0:po,1:nf,1:nc),
  !>
  !> where
  !>
  !>   - `b ` is the boundary ID,
  !>   - `po` the polynomial order,
  !>   - `nf` the number of faces, and
  !>   - `nc` the number of components.
  !>
  !> Note that `po` is identical for all boundaries, while `nf` is not.
  !> The number of components `nc` can be equal or different.
  !>
  !> If the variable is initialized with the constructor, new memory is
  !> allocated in `mem` for storing the values.
  !> Variables generated via `GetSlice`are linke to the source by default
  !> and may get lost if the latter is deleted. Passing `copy = T` causes
  !> the values to copied into fresh memory which removes this risk.

  type BoundaryVariable_3D

    integer :: po = 0 !< polynomial order
    integer :: nc = 0 !< number of components

    class(MeshBoundary_3D), pointer :: boundary => null()
    real(RNP), contiguous,  pointer :: val(:,:,:,:) => null() !< value access
    real(RNP), allocatable, private :: mem(:,:,:,:) !< memory allocated to val

  contains

    procedure :: Init => Init_BoundaryVariable_3D
    procedure :: GetSlice

    generic :: Extract => Extract_A, Extract_S
    procedure, private :: Extract_A, Extract_S

    procedure :: ExtractNormalComponent

    generic :: MergeWithTraceVariable => MergeWithTraceVar_S, MergeWithTraceVar_A
    procedure, private :: MergeWithTraceVar_S, MergeWithTraceVar_A

  end type BoundaryVariable_3D

contains

  !=============================================================================
  ! Initializers

  !-----------------------------------------------------------------------------
  !> 3D boundary variable initialization

  subroutine Init_BoundaryVariable_3D(this, boundary, po, nc)
    class(BoundaryVariable_3D), target, intent(inout) :: this
    class(MeshBoundary_3D), target, intent(in) :: boundary
    integer, intent(in) :: po
    integer, intent(in) :: nc

    if (allocated(this % mem)) deallocate(this % mem)
    allocate(this % mem(0:po, 0:po, boundary%n_face, nc))

    this % po = po
    this % nc = nc

    this % boundary         => boundary
    this % val(0:,0:,1:,1:) => this % mem

  end subroutine Init_BoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> Create a new boundary variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % val` if `copy` is false
  !> or absent. Otherwise they are stored in fresh memory, i.e. `slice % mem`.
  !> In an OpenMP parallel section the routine is executed only by the master
  !> thread.

  subroutine GetSlice(this, slice, first, last, copy)
    class(BoundaryVariable_3D), target, intent(in) :: this
    class(BoundaryVariable_3D), target, intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    logical :: copy_

    !$omp master

    slice % po = this % po
    slice % nc = 1 + last - first
    slice % boundary => this % boundary

    if (present(copy)) then
      copy_ = copy
    else
      copy_ = .false.
    end if

    if (copy_) then
      slice % mem = this % val(:,:,:,first:last)
      slice % val(0:,0:,1:,1:) => slice % mem
    else
      slice % val(0:,0:,1:,1:) => this % val(:,:,:,first:last)
      if (allocated(slice % mem)) deallocate(slice % mem)
    end if

    !$omp end master

  end subroutine GetSlice

  !=============================================================================
  ! Extraction

  !-----------------------------------------------------------------------------
  !> Extract boundary variable from array-valued mesh variable
  !>
  !> `this` must be properly initialized on input!

  subroutine Extract_A(this, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), contiguous, target, intent(in) :: v(0:,0:,0:,:,:)

    integer :: c, e, f, i, j, k, m

    !$omp master
    if (ubound(v,1) /= this%po .or. size(v,5) /= this%nc) then
      call Error('Extract_A', 'arguments do not match', 'Boundary_Variable__3D')
    end if
    !$omp end master

    associate(vb => this % val, po => this%po, nc => this%nc)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        select case(m)

        case(1,2)
          i = (m - 1) * po
          do c = 1, nc
            do k = 0, po
            do j = 0, po
              vb(j,k,f,c) = v(i,j,k,e,c)
            end do
            end do
          end do

        case(3,4)
          j = (m - 3) * po
          do c = 1, nc
            do k = 0, po
            do i = 0, po
              vb(i,k,f,c) = v(i,j,k,e,c)
            end do
            end do
          end do

        case(5,6)
          k = (m - 5) * po
          do c = 1, nc
            do j = 0, po
            do i = 0, po
              vb(i,j,f,c) = v(i,j,k,e,c)
            end do
            end do
          end do

        end select

      end do
    end associate

  end subroutine Extract_A

  !-----------------------------------------------------------------------------
  !> Extract boundary variable from scalar mesh variable
  !>
  !> `this` must be properly initialized on input!

  subroutine Extract_S(this, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), contiguous, target, intent(in) :: v(0:,0:,0:,:)

    integer :: e, f, i, j, k, m

    !$omp master
    if (ubound(v,1) /= this%po) then
      call Error('Extract_S', 'arguments do not match', 'Boundary_Variable__3D')
    end if
    !$omp end master

    associate(vb => this % val, po => this%po)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        select case(m)

        case(1,2)
          i = (m - 1) * po
          do k = 0, po
          do j = 0, po
            vb(j,k,f,1) = v(i,j,k,e)
          end do
          end do

        case(3,4)
          j = (m - 3) * po
          do k = 0, po
          do i = 0, po
            vb(i,k,f,1) = v(i,j,k,e)
          end do
          end do

        case(5,6)
          k = (m - 5) * po
          do j = 0, po
          do i = 0, po
            vb(i,j,f,1) = v(i,j,k,e)
          end do
          end do

        end select

      end do
    end associate

  end subroutine Extract_S

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a spectral-element vector
  !>
  !> `this` must be properly initialized on input!

  subroutine ExtractNormalComponent(this, sem, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), contiguous, target, intent(in) :: v(0:,0:,0:,:,:)

    !$omp master
    if (this%po  /= ubound(v,1) .or. this%nc /= 1 .or. size(v,5) /= 3) then
      call Error( 'ExtractNormalComponent' &
                , 'arguments do not match' &
                , 'Boundary_Variable__3D'  )
    end if
    !$omp end master

    if (sem % mesh % regular) then
      call ExtractNormalComponent_R(this, v)
    else
      call ExtractNormalComponent_D(this, sem, v)
    end if

  end subroutine ExtractNormalComponent

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a vector: regular mesh

  subroutine ExtractNormalComponent_R(this, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), contiguous, target, intent(in) :: v(0:,0:,0:,:,:)

    integer :: e, f, i, j, k, m, n

    associate(vb_n => this % val, po => this%po)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        select case(m)

        case(1,2)
          i = (m - 1) * po
          n = (m - 1) * 2 - 1
          do k = 0, po
          do j = 0, po
            vb_n(j,k,f,1) = n * v(i,j,k,e,1)
          end do
          end do

        case(3,4)
          j = (m - 3) * po
          n = (m - 3) * 2 - 1
          do k = 0, po
          do i = 0, po
            vb_n(i,k,f,1) = n * v(i,j,k,e,2)
          end do
          end do

        case(5,6)
          k = (m - 5) * po
          n = (m - 5) * 2 - 1
          do j = 0, po
          do i = 0, po
            vb_n(i,j,f,1) = n * v(i,j,k,e,3)
          end do
          end do

        end select

      end do
    end associate

  end subroutine ExtractNormalComponent_R

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a vector: deformed mesh

  subroutine ExtractNormalComponent_D(this, sem, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), contiguous, target, intent(in) :: v(0:,0:,0:,:,:)

    integer :: e, f, i, j, k, m

    associate(vb_n => this%val, n => sem%metrics % n, po => this%po)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        select case(m)

        case(1,2)
          i = (m - 1) * po
          do k = 0, po
          do j = 0, po
            vb_n(j,k,f,1) = n(j,k,m,e,1) * v(i,j,k,e,1) &
                          + n(j,k,m,e,2) * v(i,j,k,e,2) &
                          + n(j,k,m,e,3) * v(i,j,k,e,3)
          end do
          end do

        case(3,4)
          j = (m - 3) * po
          do k = 0, po
          do i = 0, po
            vb_n(i,k,f,1) = n(i,k,m,e,1) * v(i,j,k,e,1) &
                          + n(i,k,m,e,2) * v(i,j,k,e,2) &
                          + n(i,k,m,e,3) * v(i,j,k,e,3)
          end do
          end do

        case(5,6)
          k = (m - 5) * po
          do j = 0, po
          do i = 0, po
            vb_n(i,j,f,1) = n(i,j,m,e,1) * v(i,j,k,e,1) &
                          + n(i,j,m,e,2) * v(i,j,k,e,2) &
                          + n(i,j,m,e,3) * v(i,j,k,e,3)
          end do
          end do

        end select

      end do
    end associate

  end subroutine ExtractNormalComponent_D

  !=============================================================================
  ! Merge boundary values with element-face variable

  !-----------------------------------------------------------------------------
  !> Merge first component of boundary variable with scalar element-face variable

  subroutine MergeWithTraceVar_S(this, cb, ct, vt)
    class(BoundaryVariable_3D), intent(in) :: this
    real(RNP), intent(in)    :: cb !< coefficient of boundary variable
    real(RNP), intent(in)    :: ct !< coefficient of trace variable
    real(RNP), intent(inout) :: vt(:,:,:,:)

    integer :: f, e, m

    !$omp master
    if (this%po + 1 /= size(vt,1)) then
      call Error( 'MergeWithTraceVar_S'    &
                , 'arguments do not match' &
                , 'Boundary_Variable__3D'  )
    end if
    !$omp end master

    !$omp do
    do f = 1, this % boundary % n_face

      e = this % boundary % face(f) % element_id
      m = this % boundary % face(f) % element_face

      vt(:,:,m,e) = ct * vt(:,:,m,e) + cb * this % val(:,:,f,1)

    end do

  end subroutine MergeWithTraceVar_S

  !-----------------------------------------------------------------------------
  !> Merge boundary variable with matching array-valued element-face variable

  subroutine MergeWithTraceVar_A(this, cb, ct, vt)
    class(BoundaryVariable_3D), intent(in) :: this
    real(RNP), intent(in)    :: cb !< coefficient of boundary variable
    real(RNP), intent(in)    :: ct !< coefficient of trace variable
    real(RNP), intent(inout) :: vt(:,:,:,:,:)

    integer :: c, e, f, m

    !$omp master
    if (this%po + 1 /= size(vt,1) .or. this%nc /= size(vt,5)) then
      call Error( 'MergeWithTraceVar_A'    &
                , 'arguments do not match' &
                , 'Boundary_Variable__3D'  )
    end if
    !$omp end master

    !$omp do
    do f = 1, this % boundary % n_face

      e = this % boundary % face(f) % element_id
      m = this % boundary % face(f) % element_face

      do c = 1, this%nc
        vt(:,:,m,e,c) = ct * vt(:,:,m,e,c) + cb * this % val(:,:,f,c)
      end do

    end do

  end subroutine MergeWithTraceVar_A

  !=============================================================================

end module Boundary_Variable__3D
