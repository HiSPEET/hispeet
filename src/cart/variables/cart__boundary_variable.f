!> summary:  Boundary variables on Equidistant Cartesian meshes
!> author:   Joerg Stiller
!> date:     2017/03/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Boundary variables on Equidistant Cartesian meshes
!===============================================================================

module CART__Boundary_Variable

  use Kind_Parameters, only: RNP
  use Execution_Control, only: Error
  use CART__Mesh_Partition

  implicit none
  private

  public :: BoundaryVariable

  !-----------------------------------------------------------------------------
  !> Type for handling boundary variables
  !>
  !> This type provides a structure for handling variables defined on a mesh
  !> boundary and related conditions. The variable values are stored in
  !> `val(0:po,0:po,1:nf,1:nc)`, where `po` defines the polynomial order,
  !> `nf` the number of faces in the corresponding mesh boundary, and
  !> `nc` the number of components. `bc` provides the boundary conditions
  !> for each component, or empty characters if no condition applies.
  !>
  !> Typically one generates an array of boundary variables, one for each
  !> boundary. This is achieved using the inbuilt `New` or `Extract` procedures.
  !> For example, given the local `mesh` partiton and the corresponding variable
  !> `u(0:po,0:po,0:po,1:ne,1:nc)`
  !>
  !>     type(BoundaryVariable) :: bv_u(mesh%nb)   ! array of boundary variables
  !>     character :: bc(mesh%nb) = 'D'            ! all Dirichlet conditions
  !>     do i = 1, mesh%nb
  !>       call bv_u(i) % Extract(mesh, u, i, bc(i))
  !>     end do
  !>
  !> The variable values `val` are private, but can be accessed via type-bound
  !> procedure as follows (using `bv_u` as defined in the above example):
  !>
  !>     real(RNP), contiguous, pointer :: s(:,:,:)   ! scalar variable
  !>     real(RNP), contiguous, pointer :: v(:,:,:,:) ! multi-component variable
  !>
  !>     ! for boundary i, assign component c to s(0:po,0:po,1:nf)
  !>     s => bv_u(i) % Component(c)
  !>
  !>     ! assign components c1:c2 to v(0:po,0:po,1:nf,1:1+c2-c1)
  !>     v => bv_u(i) % Components(c1, c2)
  !>
  !>     ! assign all components to v(0:po,0:po,1:nf,1:nc)
  !>     v => bv_u(i) % Components()
  !>
  !> Using `GetHandle` it is possible to define a new boundary variable as a
  !> subset of an existing one, e.g.
  !>
  !>     type(BoundaryVariable) :: bv_u(nb)
  !>     type(BoundaryVariable) :: bv_v(nb)
  !>
  !>     ! generate bv_v from components c1:c2
  !>     call bv_u % GetHandle(c1, c2, bv_v)
  !>
  !> Note that `GetHandle`is elemental: In the above example, `bv_v(i)` is
  !> generated from `bv_u` for `i = 1:nb`. The arguments `c1` and `c2` defining
  !> the range of components can be scalar or arrays of dimension `nb`.

  type BoundaryVariable
    private

    ! public entities ..........................................................

    integer, public :: po = -1  !< polynomial order
    integer, public :: nf = -1  !< number of faces
    integer, public :: nc = -1  !< number of components

    character, pointer, public :: bc(:) => null() !< boundary conditions

    ! private components .......................................................

    real(RNP), pointer, contiguous :: val(:,:,:,:) => null() !< values

    logical :: is_original = .false. !< indicates original boundary variable

  contains

    generic,   public  :: New        =>  New_S, New_A, New_O
    procedure, private :: New_S
    procedure, private :: New_A
    procedure, private :: New_O

    generic,   public  :: Extract    =>  Extract_S, Extract_A
    procedure, private :: Extract_S
    procedure, private :: Extract_A

    procedure, public  :: ExtractNormalComponent

    procedure, public  :: BoundaryCondition
    procedure, public  :: Component
    procedure, public  :: Components

    procedure, public  :: GetHandle

    final :: Delete_BoundaryVariable

  end type BoundaryVariable

contains

!===============================================================================
! New

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- scalar with BC

subroutine New_S(this, mesh, po, b, bc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: b     !< boundary ID
  character,               intent(in)    :: bc    !< boundary conditions

  call New_X(this, mesh, po, 1, b, bc)

end subroutine New_S

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- array with BC

subroutine New_A(this, mesh, po, b, bc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: b     !< boundary ID
  character,               intent(in)    :: bc(:) !< boundary conditions

  call New_X(this, mesh, po, size(bc), b, bc)

end subroutine New_A

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- no BC

subroutine New_O(this, mesh, po, b, nc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: nc    !< number of components
  integer,                 intent(in)    :: b     !< boundary ID

  call New_X(this, mesh, po, nc, b)

end subroutine New_O

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- eXplicit

subroutine New_X(this, mesh, po, nc, b, bc)
  class(BoundaryVariable), intent(inout) :: this   !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh   !< mesh partition
  integer,                 intent(in)    :: po     !< polynomial order
  integer,                 intent(in)    :: nc     !< number of components
  integer,                 intent(in)    :: b      !< boundary ID
  character,     optional, intent(in)    :: bc(nc) !< boundary conditions

  integer :: nf

  nf = mesh % boundary(b) % nf

  ! (re)allocate storage, if required
  if (this%po /= po .or. this%nf /= nf .or. this%nc /= nc) then
    if (this%is_original) then
      deallocate(this%bc, this%val)
    end if
    allocate(this%bc(nc), source = ' ')
    allocate(this%val(0:po, 0:po, nf, nc))
  end if

  if (present(bc)) this%bc = bc

  this % po = ubound(this % val, 1)
  this % nf = ubound(this % val, 3)
  this % nc = ubound(this % val, 4)

  this % is_original = .true.

end subroutine New_X

!===============================================================================
! Extract

!-------------------------------------------------------------------------------
!> Extract boundary variable from mesh variable -- scalar

subroutine Extract_S(this, mesh, v, b, bc)
  class(BoundaryVariable), intent(inout) :: this          !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh          !< mesh partition
  real(RNP),               intent(in)    :: v(0:,0:,0:,:) !< mesh variable
  integer,                 intent(in)    :: b             !< boundary ID
  character,     optional, intent(in)    :: bc            !< boundary condition

  call Extract_X(this, mesh, ubound(v,1), size(v,4), 1, v, b, bc)

end subroutine Extract_S

!-------------------------------------------------------------------------------
!> Extract boundary variable from mesh variable -- array

subroutine Extract_A(this, mesh, v, b, bc)
  class(BoundaryVariable), intent(inout) :: this            !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh            !< mesh partition
  real(RNP),               intent(in)    :: v(0:,0:,0:,:,:) !< mesh variable
  integer,                 intent(in)    :: b               !< boundary ID
  character,     optional, intent(in)    :: bc(:)           !< boundary cond.

  call Extract_X(this, mesh, ubound(v,1), size(v,4), size(v,5), v, b, bc)

end subroutine Extract_A

!-------------------------------------------------------------------------------
!> Extract boundary variable from mesh variable -- eXplicit shape

subroutine Extract_X(this, mesh, po, ne, nc, v, b, bc)
  class(BoundaryVariable), intent(inout) :: this   !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh   !< mesh partition
  integer,   intent(in) :: po                      !< polynomial order
  integer,   intent(in) :: ne                      !< number of elements
  integer,   intent(in) :: nc                      !< number of components
  real(RNP), intent(in) :: v(0:po,0:po,0:po,ne,nc) !< mesh variable
  integer,   intent(in) :: b                       !< boundary ID
  character, optional, intent(in) :: bc(nc)        !< boundary conditions

  integer :: c, e, f, i, j, k

  if (present(bc)) then
    call this % New(mesh, po, b, bc)
  else
    call this % New(mesh, po, b, nc)
  end if

  associate(boundary => mesh%boundary(b), vb => this%val)

    do f = 1, boundary%nf

      e = boundary % face(f) % mesh_element % id

      select case(boundary % face(f) % mesh_element % face)
      case(1)
        do c = 1, nc
        do k = 0, po
        do j = 0, po
          vb(j,k,f,c) = v(0,j,k,e,c)
        end do
        end do
        end do
      case(2)
        do c = 1, nc
        do k = 0, po
        do j = 0, po
          vb(j,k,f,c) = v(po,j,k,e,c)
        end do
        end do
        end do
      case(3)
        do c = 1, nc
        do i = 0, po
        do k = 0, po
          vb(i,k,f,c) = v(i,0,k,e,c)
        end do
        end do
        end do
      case(4)
        do c = 1, nc
        do i = 0, po
        do k = 0, po
          vb(i,k,f,c) = v(i,po,k,e,c)
        end do
        end do
        end do
      case(5)
        do c = 1, nc
        do i = 0, po
        do j = 0, po
          vb(i,j,f,c) = v(i,j,0,e,c)
        end do
        end do
        end do
      case(6)
        do c = 1, nc
        do i = 0, po
        do j = 0, po
          vb(i,j,f,c) = v(i,j,po,e,c)
        end do
        end do
        end do
      end select

    end do

  end associate

end subroutine Extract_X

!===============================================================================
! ExtractNormalComponent

!-------------------------------------------------------------------------------
!> Extract normal component of a vector field to boundary variable

subroutine ExtractNormalComponent(this, mesh, v, b, bc)
  class(BoundaryVariable), intent(inout) :: this            !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh            !< mesh partition
  real(RNP),               intent(in)    :: v(0:,0:,0:,:,:) !< mesh variable
  integer,                 intent(in)    :: b               !< boundary ID
  character,     optional, intent(in)    :: bc              !< boundary condition

  real(RNP), pointer :: vb_n(:,:,:)
  integer :: e, f, i, j, k, po

  po = ubound(v,1)

  if (size(v,5) /= 3) then
    call Error( 'ExtractNormalComponent', &
                'v must be a 3-vector',   &
                'CART__Boundary_Variable' )
  end if

  if (present(bc)) then
    call this % New(mesh, po, b, bc)
  else
    call this % New(mesh, po, b, bc = ' ')
  end if

  vb_n => this % Component(1)

  associate(boundary => mesh%boundary(b))

    do f = 1, boundary%nf

      e = boundary % face(f) % mesh_element % id

      select case(boundary % face(f) % mesh_element % face)
      case(1)
        do k = 0, po
        do j = 0, po
          vb_n(j,k,f) = -v(0,j,k,e,1)
        end do
        end do
      case(2)
        do k = 0, po
        do j = 0, po
          vb_n(j,k,f) =  v(po,j,k,e,1)
        end do
        end do
      case(3)
        do i = 0, po
        do k = 0, po
          vb_n(i,k,f) = -v(i,0,k,e,2)
        end do
        end do
      case(4)
        do i = 0, po
        do k = 0, po
          vb_n(i,k,f) =  v(i,po,k,e,2)
        end do
        end do
      case(5)
        do i = 0, po
        do j = 0, po
          vb_n(i,j,f) = -v(i,j,0,e,3)
        end do
        end do
      case(6)
        do i = 0, po
        do j = 0, po
          vb_n(i,j,f) =  v(i,j,po,e,3)
        end do
        end do
      end select

    end do

  end associate

end subroutine ExtractNormalComponent

!===============================================================================
! BoundaryCondition

!-------------------------------------------------------------------------------
!> Returns the boundary condition for a selected component

character function BoundaryCondition(this, c) result(bc)
  class(BoundaryVariable), intent(in) :: this
  integer,       optional, intent(in) :: c    !< component [1]

  if (present(c)) then
    bc = this % bc(c)
  else
    bc = this % bc(1)
  end if

end function BoundaryCondition

!===============================================================================
! Component

!-------------------------------------------------------------------------------
!> Provides a pointer to the values of a single component

function Component(this, c) result(vb)
  class(BoundaryVariable), intent(in) :: this
  integer,  optional, intent(in) :: c         !< component [1]
  real(RNP), contiguous, pointer :: vb(:,:,:) !< boundary-face variable

  if (present(c)) then
    vb(0:,0:,1:) => this%val(:,:,:,c)
  else
    vb(0:,0:,1:) => this%val(:,:,:,1)
  end if

end function Component

!-------------------------------------------------------------------------------
!> Provides a pointer to the values of a multiple components
!>
!> The optional argumenta `c1` and `c2` define the range of target components.
!> If omitted, the corresponding upper or lower bound is assumed.

function Components(this, c1, c2) result(vb)
  class(BoundaryVariable), intent(in)  :: this
  integer,  optional, intent(in) :: c1          !< first component [1]
  integer,  optional, intent(in) :: c2          !< last component [this%nc]
  real(RNP), contiguous, pointer :: vb(:,:,:,:) !< boundary-face variable

  integer :: k1, k2

  if (present(c1)) then
    k1 = max(c1, 1)
  else
    k1 = 1
  end if

  if (present(c2)) then
    k2 = min(c2, this%nc)
  else
    k2 = this%nc
  end if

  vb(0:,0:,1:,1:) => this%val(:,:,:,k1:k2)

end function Components

!===============================================================================
! GetHandle

impure elemental subroutine GetHandle(this, c1, c2, handle)
  class(BoundaryVariable), intent(in)  :: this
  integer,                 intent(in)  :: c1     !< first component
  integer,                 intent(in)  :: c2     !< last component
  type (BoundaryVariable), intent(out) :: handle !< subset of this

  integer :: k1, k2

  k1 = max(c1, 1)
  k2 = min(c2, this%nc)

  handle % po = this % po
  handle % nf = this % nf
  handle % nc = max(k2 - k1 + 1, 0)

  handle % bc => this % bc(k1:k2)

  handle % val(0:,0:,1:,1:) => this % val(:,:,:,k1:k2)

end subroutine GetHandle

!===============================================================================
! Finalization

!-------------------------------------------------------------------------------
!> Finalize BoundaryVariable

subroutine Delete_BoundaryVariable(this)
  type(BoundaryVariable), intent(inout) :: this

  if (this % is_original) then
    deallocate(this%bc, this%val)
  end if

  this % po = -1
  this % nf = -1
  this % nc = -1

  this % bc  => null()
  this % val => null()

  this % is_original = .false.

end subroutine Delete_BoundaryVariable

!===============================================================================

end module CART__Boundary_Variable
