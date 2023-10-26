!> summary:  3D boundary variable
!> author:   Joerg Stiller
!> date:     2022/10/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Boundary_Variable__3D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
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

    procedure :: Create => Create_BoundaryVariable_3D
    procedure :: GetSlice

    generic :: Extract => Extract_A, Extract_S
    procedure, private :: Extract_A, Extract_S

    generic :: Merge   => Merge_A, Merge_S
    procedure, private :: Merge_A, Merge_S

    procedure :: ExtractNormalComponent
    procedure :: MergeNormalComponent

    generic :: MergeTrace => MergeTrace_S, MergeTrace_A
    procedure, private ::    MergeTrace_S, MergeTrace_A

    procedure :: ExtractNormalTrace
    procedure :: MergeNormalTrace

  end type BoundaryVariable_3D

contains

  !=============================================================================
  ! Initializers

  !-----------------------------------------------------------------------------
  !> 3D boundary variable initialization

  subroutine Create_BoundaryVariable_3D(this, boundary, po, nc)
    class(BoundaryVariable_3D), target, intent(inout) :: this
    class(MeshBoundary_3D), target, intent(in) :: boundary
    integer, intent(in) :: po
    integer, intent(in) :: nc

    if (allocated(this % mem)) deallocate(this % mem)
    allocate(this % mem(0:po, 0:po, boundary%n_face, nc), source = ZERO)

    this % po = po
    this % nc = nc

    this % boundary         => boundary
    this % val(0:,0:,1:,1:) => this % mem

  end subroutine Create_BoundaryVariable_3D

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
  ! Extraction from and merging with mesh variables

  !-----------------------------------------------------------------------------
  !> Extract boundary variable from array-valued mesh variable
  !>
  !> `this` must be properly initialized on input!

  subroutine Extract_A(this, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)

    call Merge_A(this, ZERO, ONE, v)

  end subroutine Extract_A

  !-----------------------------------------------------------------------------
  !> Extract boundary variable from scalar mesh variable
  !>
  !> `this` must be properly initialized on input!

  subroutine Extract_S(this, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:)

    call Merge_S(this, ZERO, ONE, v)

  end subroutine Extract_S

  !-----------------------------------------------------------------------------
  !> Merge boundary variable with array-valued mesh variable

  subroutine Merge_A(this, cb, cv, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: cv !< coefficient of mesh variable
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)

    integer :: e, f, i, j, k, l, m

    !$omp master
    if (ubound(v,1) /= this%po .or. size(v,5) /= this%nc) then
      call Error('Merge_A', 'argument mismatch', 'Boundary_Variable__3D')
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
          do l = 1, nc
            do k = 0, po
            do j = 0, po
              vb(j,k,f,l) = cb * vb(j,k,f,l) + cv * v(i,j,k,e,l)
            end do
            end do
          end do

        case(3,4)
          j = (m - 3) * po
          do l = 1, nc
            do k = 0, po
            do i = 0, po
              vb(i,k,f,l) = cb * vb(i,k,f,l) + cv * v(i,j,k,e,l)
            end do
            end do
          end do

        case(5,6)
          k = (m - 5) * po
          do l = 1, nc
            do j = 0, po
            do i = 0, po
              vb(i,j,f,l) = cb * vb(i,j,f,l) + cv * v(i,j,k,e,l)
            end do
            end do
          end do

        end select

      end do
    end associate

  end subroutine Merge_A

  !-----------------------------------------------------------------------------
  !> Merge boundary variable with scalar mesh variable

  subroutine Merge_S(this, cb, cv, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: cv !< coefficient of mesh variable
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:)

    integer :: e, f, i, j, k, m

    !$omp master
    if (ubound(v,1) /= this%po) then
      call Error('Merge_S', 'argument mismatch', 'Boundary_Variable__3D')
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
            vb(j,k,f,1) = cb * vb(j,k,f,1) + cv * v(i,j,k,e)
          end do
          end do

        case(3,4)
          j = (m - 3) * po
          do k = 0, po
          do i = 0, po
            vb(i,k,f,1) = cb * vb(i,k,f,1) + cv * v(i,j,k,e)
          end do
          end do

        case(5,6)
          k = (m - 5) * po
          do j = 0, po
          do i = 0, po
            vb(i,j,f,1) = cb * vb(i,j,f,1) + cv * v(i,j,k,e)
          end do
          end do

        end select

      end do
    end associate

  end subroutine Merge_S

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a spectral-element vector
  !>
  !> `this` must be properly initialized on input!

  subroutine ExtractNormalComponent(this, sem, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)

    call MergeNormalComponent(this, sem, ZERO, ONE, v)

  end subroutine ExtractNormalComponent

  !-----------------------------------------------------------------------------
  !> Merge with normal component of a vector-valued mesh variable

  subroutine MergeNormalComponent(this, sem, cb, cv, v)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: cv !< coefficient of mesh variable
    real(RNP), contiguous, intent(in) :: v(0:,0:,0:,:,:)

    integer :: e, f, i, j, k, m

    !$omp master
    if (this%po /= ubound(v,1) .or. this%nc /= 1 .or. size(v,5) < 3) then
      call Error( 'MergeNormalComponent'  &
                , 'argument mismatch'     &
                , 'Boundary_Variable__3D' )
    end if
    !$omp end master

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
            vb_n(j,k,f,1) = cb * vb_n(j,k,f,1)                 &
                          + cv * ( n(j,k,m,e,1) * v(i,j,k,e,1) &
                                 + n(j,k,m,e,2) * v(i,j,k,e,2) &
                                 + n(j,k,m,e,3) * v(i,j,k,e,3) )
          end do
          end do

        case(3,4)
          j = (m - 3) * po
          do k = 0, po
          do i = 0, po
            vb_n(i,k,f,1) = cb * vb_n(i,k,f,1)                 &
                          + cv * ( n(i,k,m,e,1) * v(i,j,k,e,1) &
                                 + n(i,k,m,e,2) * v(i,j,k,e,2) &
                                 + n(i,k,m,e,3) * v(i,j,k,e,3) )
          end do
          end do

        case(5,6)
          k = (m - 5) * po
          do j = 0, po
          do i = 0, po
            vb_n(i,j,f,1) = cb * vb_n(i,j,f,1)                 &
                          + cv * ( n(i,j,m,e,1) * v(i,j,k,e,1) &
                                 + n(i,j,m,e,2) * v(i,j,k,e,2) &
                                 + n(i,j,m,e,3) * v(i,j,k,e,3) )
          end do
          end do

        end select

      end do
    end associate

  end subroutine MergeNormalComponent

  !=============================================================================
  ! Merge boundary values with trace variable

  !-----------------------------------------------------------------------------
  !> Merge scalar trace variable with first component of boundary variable

  subroutine MergeTrace_S(this, cb, ct, vt)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: ct !< coefficient of trace variable
    real(RNP), intent(in) :: vt(0:,0:,:,:) !< scalar trace variable

    integer :: f, e, m

    !$omp master
    if (this%po  /= ubound(vt,1)) then
      call Error('MergeTrace_S', 'arguments mismatch', 'Boundary_Variable__3D')
    end if
    !$omp end master

    associate(vb => this % val)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        vb(:,:,f,1) = ct * vt(:,:,m,e) + cb * vb(:,:,f,1)

      end do

    end associate

  end subroutine MergeTrace_S

  !-----------------------------------------------------------------------------
  !> Merge array-valued trace variable with boundary variable

  subroutine MergeTrace_A(this, cb, ct, vt)
    class(BoundaryVariable_3D), intent(inout) :: this
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: ct !< coefficient of trace variable
    real(RNP), intent(in) :: vt(0:,0:,:,:,:) !< array-valued trace variable

    integer :: c, e, f, m

    !$omp master
    if (this%po  /= ubound(vt,1) .or. this%nc /= size(vt,5)) then
      call Error('MergeTrace_A', 'arguments mismatch', 'Boundary_Variable__3D')
    end if
    !$omp end master

    associate(vb => this % val)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        do c = 1, this%nc
          vb(:,:,f,c) = ct * vt(:,:,m,e,c) + cb * vb(:,:,f,c)
        end do

      end do

    end associate

  end subroutine MergeTrace_A

  !-----------------------------------------------------------------------------
  !> Extract normal component of vector trace variable

  subroutine ExtractNormalTrace(this, sem, vt)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), intent(in) :: vt(0:,0:,:,:,:) !< vector trace variable

    call MergeNormalTrace(this, sem, ZERO, ONE, vt)

  end subroutine ExtractNormalTrace

  !-----------------------------------------------------------------------------
  !> Merge normal component of vector trace variable

  subroutine MergeNormalTrace(this, sem, cb, ct, vt)
    class(BoundaryVariable_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), intent(in) :: sem
    real(RNP), intent(in) :: cb !< coefficient of boundary variable
    real(RNP), intent(in) :: ct !< coefficient of trace variable
    real(RNP), intent(in) :: vt(0:,0:,:,:,:) !< vector trace variable

    integer :: e, f, m

    !$omp master
    if (this%po  /= ubound(vt,1) .or. this%nc /= 1 .or. size(vt,5) < 3) then
      call Error( 'MergeNormalTrace'      &
                , 'argument mismatch'     &
                , 'Boundary_Variable__3D' )
    end if
    !$omp end master

    associate(vb => this%val, n => sem%metrics%n)

      !$omp do
      do f = 1, this % boundary % n_face

        e = this % boundary % face(f) % element_id
        m = this % boundary % face(f) % element_face

        vb(:,:,f,1) = cb * vb(:,:,f,1)                    &
                    + ct * ( n(:,:,m,e,1) * vt(:,:,m,e,1) &
                           + n(:,:,m,e,2) * vt(:,:,m,e,2) &
                           + n(:,:,m,e,3) * vt(:,:,m,e,3) )
      end do

    end associate

  end subroutine MergeNormalTrace

  !=============================================================================

end module Boundary_Variable__3D
