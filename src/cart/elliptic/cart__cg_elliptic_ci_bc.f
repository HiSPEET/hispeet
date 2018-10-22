!> summary:  Boundary treatment for CG elliptic operator: constant isotropic
!> author:   Joerg Stiller
!> date:     2018/10/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Boundary treatment for CG elliptic operator
!===============================================================================

module CART__CG_Elliptic_CI_BC

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error

  use CART__CG_Element_Operators
  use CART__Mesh_Partition
  use CART__Boundary_Variable

  implicit none
  private

  public :: ApplyBoundaryConditions

contains

!-------------------------------------------------------------------------------
!> Injects Dirichlet conditions to solution and Neumann contributions to RHS

subroutine ApplyBoundaryConditions(mesh, eop, nu, bv, u, f)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh  !< mesh partition
  class(CG_ElementOperators3D), intent(in)    :: eop   !< CG element operators
  real(RNP),                    intent(in)    :: nu    !< diffusivity
  type(BoundaryVariable),       intent(in)    :: bv(:) !< boundary conditions
  real(RNP),                    intent(inout) :: u     !< solution
  real(RNP),                    intent(inout) :: f     !< RHS

  dimension :: u(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne)
  dimension :: f(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne)

  ! local variables ............................................................

  real(RNP), pointer, contiguous :: ub(:,:,:), dn_u(:,:,:)
  real(RNP), allocatable :: nu_Mf(:,:)
  real(RNP) :: cf(3)
  integer   :: b, e, i, j, k, l, s

  associate(po => eop%po, Ms => eop%w, dx => eop%dx)

    ! initialization ...........................................................

    ! metric factors
    cf(1) = dx(2)*dx(3) / 4
    cf(2) = dx(3)*dx(1) / 4
    cf(3) = dx(1)*dx(2) / 4

    ! standard face mass matrix scaled with diffusivity
    allocate(nu_Mf(0:po,0:po))
    do j = 0, po
    do i = 0, po
      nu_Mf(i,j) = nu * Ms(i) * Ms(j)
    end do
    end do

    Boundaries: do b = 1, size(bv)

      Boundary_Faces: associate(face => mesh%boundary(b)%face)

        select case(bv(b) % BoundaryCondition())

        case('D')

          ! Dirichlet BC .......................................................

          ub => bv(b) % Component()

          do l = 1, size(face)
            e = face(l) % mesh_element % id
            s = face(l) % mesh_element % face

            select case(s)

            case(1) ! direction 1: west
              do k = 0, po
              do j = 0, po
                u(0,j,k,e) = ub(j,k,l)
              end do
              end do

            case(2) ! direction 1: east
              do k = 0, po
              do j = 0, po
                u(po,j,k,e) = ub(j,k,l)
              end do
              end do

            case(3) ! direction 2: south
              do k = 0, po
              do i = 0, po
                u(i,0,k,e) = ub(bi,k,l)
              end do
              end do

            case(4) ! direction 2: north
              do k = 0, po
              do i = 0, po
                u(i,po,k,e) = ub(i,k,l)
              end do
              end do

            case(5) ! direction 3: bottom
              do j = 0, po
              do i = 0, po
                u(i,j,0,e) = ub(i,j,l)
              end do
              end do

            case(6) ! direction 3: top
              do j = 0, po
              do i = 0, po
                u(i,j,po,e) = ub(i,j,l)
              end do
              end do

            end select
          end do

        case('N')

          ! Neumann BC .........................................................

          dn_u => bv(b) % Component()

          do l = 1, size(face)
            e = face(l) % mesh_element % id
            s = face(l) % mesh_element % face

            select case(s)

            case(1) ! direction 1: west
              do k = 0, po
              do j = 0, po
                f(0,j,k,e) = f(0,j,k,e) + cf(1) * nu_Mf(j,k) * dn_u(j,k,l)
              end do
              end do

            case(2) ! direction 1: east
              do k = 0, po
              do j = 0, po
                f(po,j,k,e) = f(po,j,k,e) + cf(1) * nu_Mf(j,k) * dn_u(j,k,l)
              end do
              end do

            case(3) ! direction 2: south
              do k = 0, po
              do i = 0, po
                f(i,0,k,e) = f(i,0,k,e) + cf(2) * nu_Mf(i,k) * dn_u(i,k,l)
              end do
              end do

            case(4) ! direction 2: north
              do k = 0, po
              do i = 0, po
                f(i,po,k,e) = f(i,po,k,e) + cf(2) * nu_Mf(i,k) * dn_u(i,k,l)
              end do
              end do

            case(5) ! direction 3: bottom
              do j = 0, po
              do i = 0, po
                f(i,j,0,e) = f(i,j,0,e) + cf(3) * nu_Mf(i,j) * dn_u(i,j,l)
              end do
              end do

            case(6) ! direction 3: top
              do j = 0, po
              do i = 0, po
                f(i,j,po,e) = f(i,j,po,e) + cf(3) * nu_Mf(i,j) * dn_u(i,j,l)
              end do
              end do

            end select
          end do

        end select

      end associate Boundary_Faces
    end do Boundaries
  end associate

end subroutine ApplyBoundaryConditions

!===============================================================================

end module CART__CG_Elliptic_CI_BC
