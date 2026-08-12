!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Boundary variables on Equidistant Cartesian meshes
!> author:   Joerg Stiller
!> date:     2017/03/13
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
  !> boundary. This is achieved using the overloaded constructors or
  !> `Extract` procedures. For example, given the local `mesh` partiton
  !> and the corresponding variable `u(0:po,0:po,0:po,1:ne,1:nc)`
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

  type BoundaryVariable

    integer, public :: po = -1  !< polynomial order
    integer, public :: nf = -1  !< number of faces
    integer, public :: nc = -1  !< number of components

    character, allocatable, public :: bc(:)        !< boundary conditions
    real(RNP), allocatable, public :: val(:,:,:,:) !< values

  contains

    generic,   public  :: Init_BoundaryVariable => Init_S, Init_A, Init_O
    procedure, private :: Init_S
    procedure, private :: Init_A
    procedure, private :: Init_O

    generic,   public  :: Extract => Extract_S, Extract_A
    procedure, private :: Extract_S
    procedure, private :: Extract_A

    procedure, public  :: ExtractNormalComponent

    generic,   public  :: BoundaryCondition => BoundaryCondition_1, &
                                               BoundaryCondition_S
    procedure, private :: BoundaryCondition_1
    procedure, private :: BoundaryCondition_S

    procedure, public  :: Component
    procedure, public  :: Components

  end type BoundaryVariable

  ! constructor interface
  interface BoundaryVariable
    module procedure New_S
    module procedure New_A
    module procedure New_O
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New boundary variable -- scalar with BC

type(BoundaryVariable) function New_S(mesh, po, b, bc) result(this)
  class(MeshPartition), intent(in)    :: mesh  !< mesh partition
  integer,              intent(in)    :: po    !< polynomial order
  integer,              intent(in)    :: b     !< boundary ID
  character,            intent(in)    :: bc    !< boundary conditions

  call Init_S(this, mesh, po, b, bc)

end function New_S

!-------------------------------------------------------------------------------
!> New boundary variable -- array with BC

type(BoundaryVariable) function New_A(mesh, po, b, bc) result(this)
  class(MeshPartition), intent(in)    :: mesh  !< mesh partition
  integer,              intent(in)    :: po    !< polynomial order
  integer,              intent(in)    :: b     !< boundary ID
  character,            intent(in)    :: bc(:) !< boundary conditions

  call Init_A(this, mesh, po, b, bc)

end function New_A

!-------------------------------------------------------------------------------
!> New boundary variable -- no BC

type(BoundaryVariable) function New_O(mesh, po, b, nc) result(this)
  class(MeshPartition), intent(in)    :: mesh  !< mesh partition
  integer,              intent(in)    :: po    !< polynomial order
  integer,              intent(in)    :: nc    !< number of components
  integer,              intent(in)    :: b     !< boundary ID

  call Init_O(this, mesh, po, b, nc)

end function New_O

!===============================================================================
! Initialization

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- scalar with BC

subroutine Init_S(this, mesh, po, b, bc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: b     !< boundary ID
  character,               intent(in)    :: bc    !< boundary conditions

  character :: bc_(1)

  bc_ = bc
  call Init_X(this, mesh, po, 1, b, bc_)

end subroutine Init_S

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- array with BC

subroutine Init_A(this, mesh, po, b, bc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: b     !< boundary ID
  character,               intent(in)    :: bc(:) !< boundary conditions

  call Init_X(this, mesh, po, size(bc), b, bc)

end subroutine Init_A

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- no BC

subroutine Init_O(this, mesh, po, b, nc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: nc    !< number of components
  integer,                 intent(in)    :: b     !< boundary ID

  call Init_X(this, mesh, po, nc, b)

end subroutine Init_O

!-------------------------------------------------------------------------------
!> Generates a new boundary variable -- eXplicit
!>
!> `bc` is deliberately defined with assumed shape to allow for non-contiguous
!> arguments constructed from structure components

subroutine Init_X(this, mesh, po, nc, b, bc)
  class(BoundaryVariable), intent(inout) :: this  !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                 intent(in)    :: po    !< polynomial order
  integer,                 intent(in)    :: nc    !< number of components
  integer,                 intent(in)    :: b     !< boundary ID
  character,     optional, intent(in)    :: bc(:) !< boundary conditions (1:nc)

  integer :: nf

  nf = mesh % boundary(b) % nf

  ! (re)allocate storage, if required
  if (this%po /= po .or. this%nf /= nf .or. this%nc /= nc) then
    if (allocated(this % bc )) deallocate(this%bc)
    if (allocated(this % val)) deallocate(this%val)
    allocate(this % bc(nc), source = ' ')
    allocate(this % val(0:po, 0:po, nf, nc))
  end if

  if (present(bc)) this%bc = bc

  this % po = ubound(this % val, 1)
  this % nf = ubound(this % val, 3)
  this % nc = ubound(this % val, 4)

end subroutine Init_X

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

  character :: bc_(1)

  if (present(bc)) then
    bc_ = bc
  else
    bc_ = ' '
  end if

  call Extract_X(this, mesh, ubound(v,1), size(v,4), 1, v, b, bc_)

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
!>
!> `bc` is deliberately defined with assumed shape to allow for non-contiguous
!> arguments constructed from structure components

subroutine Extract_X(this, mesh, po, ne, nc, v, b, bc)
  class(BoundaryVariable), intent(inout) :: this   !< boundary variable
  class(MeshPartition),    intent(in)    :: mesh   !< mesh partition
  integer,   intent(in) :: po                      !< polynomial order
  integer,   intent(in) :: ne                      !< number of elements
  integer,   intent(in) :: nc                      !< number of components
  real(RNP), intent(in) :: v(0:po,0:po,0:po,ne,nc) !< mesh variable
  integer,   intent(in) :: b                       !< boundary ID
  character, optional, intent(in) :: bc(:)         !< boundary conditions

  integer :: c, e, f, i, j, k

  if (present(bc)) then
    call this % Init_BoundaryVariable(mesh, po, b, bc)
  else
    call this % Init_BoundaryVariable(mesh, po, b, nc)
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
    call this % Init_BoundaryVariable(mesh, po, b, bc)
  else
    call this % Init_BoundaryVariable(mesh, po, b, bc = ' ')
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
!> Returns the boundary condition of first component

impure elemental function BoundaryCondition_1(this) result(bc)
  class(BoundaryVariable), intent(in) :: this
  character :: bc   !< BC type of component 1

  if (allocated(this % bc)) then
    bc = this % bc(1)
  else
    bc = ''
  end if

end function BoundaryCondition_1

!-------------------------------------------------------------------------------
!> Returns the boundary condition for a selected component

impure elemental function BoundaryCondition_S(this, c) result(bc)
  class(BoundaryVariable), intent(in) :: this
  integer,  intent(in) :: c    !< component
  character            :: bc   !< BC type of component c

  if (allocated(this % bc)) then
    bc = this % bc(c)
  else
    bc = ''
  end if

end function BoundaryCondition_S

!===============================================================================
! Component

!-------------------------------------------------------------------------------
!> Provides a pointer to the values of a single component

function Component(this, c) result(vb)
  class(BoundaryVariable), target, intent(in) :: this
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
!> The optional arguments `c1` and `c2` define the range of target components.
!> If omitted, the corresponding lower or upper bound is assumed.

function Components(this, c1, c2) result(vb)
  class(BoundaryVariable), target, intent(in) :: this
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

end module CART__Boundary_Variable
