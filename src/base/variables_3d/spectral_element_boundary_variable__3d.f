!> summary:  3D spectral element boundary variable
!> author:   Joerg Stiller
!> date:     2021/10/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Element_Boundary_Variable__3D
  use Kind_Parameters, only: RNP
  use Mesh_Boundary__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Variable__3D
  use Spectral_Element_Vector__3D
  implicit none
  private

  public :: SpectralElementBoundaryVariable_3D

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable
  !>
  !> Base type for keeping a set of variables on the boundaries of a spectral
  !> element mesh. The mesh, standard operators and metric coefficients are
  !> provided by the pointer `sem`, which can be shared with other entities.
  !> The constructor and the `GetSlice` procedure are elemental and can be used
  !> to create scalar variable on one boundary or an array boundary variables,
  !> e.g.
  !>
  !>       class(SpectralElementMesh_3d) :: sem
  !>       type(SpectralElementBoundaryVariable_3D), allocatable :: sebv(:)
  !>       type(SpectralElementBoundaryVariable_3D) :: sebs
  !>
  !>       associate(boundary => sem % mesh boundary)
  !>
  !>         ! create an array with 4 components for all boundaries
  !>         sebv = SpectralElementBoundaryVariable_3D(sem, boundary, nc = 4)
  !>
  !>         ! create a slice containing components 1:2 for boundary b
  !>         call sebv(3) % GetSlice(sebs, first=1, last=2)
  !>
  !>       end associate
  !>
  !> The values of an array defined for all boundaries are accessed via
  !>
  !>       sebv(b) % val(0:po,0:po,1:nf,1:nc),
  !>
  !> where
  !>
  !>   - `b ` is the boundary ID,
  !>   - `po` the polynomial order,
  !>   - `nf` the number of faces, and
  !>   - `nc` the number of components.
  !>
  !> Note that `po` is identical for all boundaries, while `nf` is not.
  !> The number of components `nc` can be equal or different, depending
  !> on the creation of the variable.
  !>
  !> If the variable is initialized with the constructor, new memory is
  !> allocated in `mem` for storing the values.
  !> Variables generated via `GetSlice`are linke to the source by default
  !> and may get lost if the latter is deleted. Passing `copy = T` causes
  !> the values to copied into fresh memory which removes this risk.

  type SpectralElementBoundaryVariable_3D
    class(SpectralElementMesh_3D), pointer :: sem !< spectral element mesh
    real(RNP), contiguous,  pointer :: val(:,:,:,:) => null() !< value access
    real(RNP), allocatable, private :: mem(:,:,:,:) !< memory allocated to val
    integer :: bid = -1 !< associated mesh boundary identifier
  contains
    procedure :: Init_SpectralElementBoundaryVariable_3D => Init_SEBV
    procedure :: GetSlice
    procedure :: Extract
    procedure :: ExtractNormalComponent
    generic   :: CopyToElementFaceVariable => CopyToEFV_S, CopyToEFV_A
    procedure, private :: CopyToEFV_S, CopyToEFV_A
  end type SpectralElementBoundaryVariable_3D

  ! constructor
  interface SpectralElementBoundaryVariable_3D
    module procedure New_Scratch
    module procedure New_Slice
  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> New 3D spectral element boundary variable generated from scratch

  impure elemental function New_Scratch(sem, boundary, nc) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(MeshBoundary_3D), intent(in) :: boundary
    integer, intent(in) :: nc
    type(SpectralElementBoundaryVariable_3D) :: this

    call Init_SEBV(this, sem, boundary, nc)

  end function New_Scratch

  !-----------------------------------------------------------------------------
  !> New 3D spectral element boundary variable generated from slice

  impure elemental function New_Slice(sebv, first, last, copy) result(this)
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: sebv
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]
    type(SpectralElementBoundaryVariable_3D) :: this

    call sebv % GetSlice(this, first, last, copy)

  end function New_Slice

  !=============================================================================
  ! Initializers

  !-----------------------------------------------------------------------------
  !> 3D spectral element boundary variable initialization

  impure elemental subroutine Init_SEBV(this, sem, boundary, nc)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(MeshBoundary_3D), intent(in) :: boundary
    integer, intent(in) :: nc

    integer :: po, nf

    po = sem % std_op % po
    nf = boundary % n_face

    if (allocated(this % mem)) deallocate(this % mem)
    allocate(this % mem(0:po, 0:po, nf, nc))

    this % bid = boundary % id
    this % val(0:,0:,1:,1:) => this % mem
    this % sem              => sem

  end subroutine Init_SEBV

  !-----------------------------------------------------------------------------
  !> Create a new spectral element boundary variable as a slice of the given one
  !>
  !> The values of the new variable refer to `this % val` if `copy` is false
  !> or absent. Otherwise they are stored in fresh memory, i.e. `slice % mem`.
  !> In an OpenMP parallel section the routine is executed only by the master
  !> thread.

  impure elemental subroutine GetSlice(this, slice, first, last, copy)
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: slice
    integer,           intent(in) :: first !< first component of slice
    integer,           intent(in) :: last  !< last component of slice
    logical, optional, intent(in) :: copy  !< copy into fresh memory [F]

    logical :: copy_

    !$omp master

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

    slice % bid = this % bid

    !$omp end master

  end subroutine GetSlice

  !=============================================================================
  ! Extraction

  !-----------------------------------------------------------------------------
  !> Extract boundary variable from spectral-element variable
  !>
  !> `this` must be properly initialized on input!

  impure elemental subroutine Extract(this, sev, boundary)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementVariable_3D), intent(in) :: sev
    class(MeshBoundary_3D), intent(in) :: boundary

    integer :: c, e, f, i, j, k, m, nc, po

    associate(v => sev % val, vb => this % val)

      po = ubound(v,1)
      nc = ubound(v,5)

      !$omp do
      do f = 1, boundary % n_face

        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

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

  end subroutine Extract

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a spectral-element vector
  !>
  !> `this` must be properly initialized on input!

  impure elemental subroutine ExtractNormalComponent(this, sev, boundary)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementVector_3D), intent(in) :: sev
    class(MeshBoundary_3D), intent(in) :: boundary

    if (this % sem % mesh % regular) then
      call ExtractNormalComponent_R(this, sev, boundary)
    else
      call ExtractNormalComponent_D(this, sev, boundary)
    end if

  end subroutine ExtractNormalComponent

  !-----------------------------------------------------------------------------
  !> Extract the normal component of a spectral-element vector: regular mesh

  impure elemental subroutine ExtractNormalComponent_R(this, sev, boundary)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementVector_3D), intent(in) :: sev
    class(MeshBoundary_3D), intent(in) :: boundary

    integer :: e, f, i, j, k, m, n, po

    associate(v => sev % val, vb_n => this % val)

      po = ubound(v,1)

      !$omp do
      do f = 1, boundary % n_face

        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

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
  !> Extract the normal component of a spectral-element vector: deformed mesh

  impure elemental subroutine ExtractNormalComponent_D(this, sev, boundary)
    class(SpectralElementBoundaryVariable_3D), target, intent(inout) :: this
    class(SpectralElementVector_3D), intent(in) :: sev
    class(MeshBoundary_3D), intent(in) :: boundary

    integer :: e, f, i, j, k, m, po

    associate( n    => sev  % sem % metrics % n &
             , v    => sev  % val               &
             , vb_n => this % val               )

      po = ubound(v,1)

      !$omp do
      do f = 1, boundary % n_face

        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

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
  ! Copy boundary values to element-face variable

   !-----------------------------------------------------------------------------
   !> Copy first component of boundary variable to scalar element-face variable

  subroutine CopyToEFV_S(this, v)
    class(SpectralElementBoundaryVariable_3D), intent(in) :: this
    real(RNP), intent(inout) :: v(:,:,:,:)

    integer :: f, e, m

    associate(boundary => this % sem % mesh % boundary(this % bid))

      !$omp do
      do f = 1, boundary % n_face

        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

        v(:,:,m,e) = this % val(:,:,f,1)

      end do

    end associate

  end subroutine CopyToEFV_S

  !-----------------------------------------------------------------------------
  !> Copy boundary variable to matching array-valued element-face variable

  subroutine CopyToEFV_A(this, v)
    class(SpectralElementBoundaryVariable_3D), intent(in) :: this
    real(RNP), intent(inout) :: v(:,:,:,:,:)

    integer :: c, e, f, m, nc

    nc = size(this % val, 4)

    associate(boundary => this % sem % mesh % boundary(this % bid))

      !$omp do
      do f = 1, boundary % n_face

        e = boundary % face(f) % element_id
        m = boundary % face(f) % element_face

        do c = 1, nc
          v(:,:,m,e,c) = this % val(:,:,f,c)
        end do

      end do

    end associate

  end subroutine CopyToEFV_A

  !=============================================================================

end module Spectral_Element_Boundary_Variable__3D
