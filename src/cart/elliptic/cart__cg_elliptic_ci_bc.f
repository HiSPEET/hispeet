!> summary:  Boundary treatment for CG elliptic operator: constant isotropic
!> author:   Joerg Stiller
!> date:     2018/10/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   * not threadsafe, i.e. not parallelizable with OpenMP, since more than
!>     one element face may belong to the same boundary
!> @endnote
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

  public :: InjectDirichletConditions
  public :: ZeroDirichletEntries
  public :: ApplyBoundaryConditions

contains

!-------------------------------------------------------------------------------
!> Injects Dirichlet conditions into given mesh variable

subroutine InjectDirichletConditions(mesh, bv, u)
  class(MeshPartition),   intent(in)    :: mesh          !< mesh partition
  type(BoundaryVariable), intent(in)    :: bv(:)         !< boundary conditions
  real(RNP),              intent(inout) :: u(0:,0:,0:,:) !< mesh variable

  real(RNP), pointer, contiguous :: ub(:,:,:)
  integer :: b, e, i, j, k, l, s, po

  !$omp single

  po = ubound(u,1)

  Boundaries: do b = 1, size(bv)

    Boundary_Faces: associate(face => mesh%boundary(b)%face)

      select case(bv(b) % BoundaryCondition())

      case('D')

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
              u(i,0,k,e) = ub(i,k,l)
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

      end select

    end associate Boundary_Faces
  end do Boundaries

  !$omp end single

end subroutine InjectDirichletConditions

!-------------------------------------------------------------------------------
!> Sets the given mesh variable to zero in all Dirichlet points

subroutine ZeroDirichletEntries(mesh, bc, f)
  class(MeshPartition), intent(in)    :: mesh          !< mesh partition
  character,            intent(in)    :: bc(:)         !< BC types
  real(RNP),            intent(inout) :: f(0:,0:,0:,:) !< mesh variable

  integer :: b, e, i, j, k, l, s, po

  !$omp single

  po = ubound(f,1)

  Boundaries: do b = 1, size(bc)

    Boundary_Faces: associate(face => mesh%boundary(b)%face)

      if (bc(b) == 'D') then

        do l = 1, size(face)
          e = face(l) % mesh_element % id
          s = face(l) % mesh_element % face

          select case(s)

          case(1) ! direction 1: west
            do k = 0, po
            do j = 0, po
              f(0,j,k,e) = 0
            end do
            end do

          case(2) ! direction 1: east
            do k = 0, po
            do j = 0, po
              f(po,j,k,e) = 0
            end do
            end do

          case(3) ! direction 2: south
            do k = 0, po
            do i = 0, po
              f(i,0,k,e) = 0
            end do
            end do

          case(4) ! direction 2: north
            do k = 0, po
            do i = 0, po
              f(i,po,k,e) = 0
            end do
            end do

          case(5) ! direction 3: bottom
            do j = 0, po
            do i = 0, po
              f(i,j,0,e) = 0
            end do
            end do

          case(6) ! direction 3: top
            do j = 0, po
            do i = 0, po
              f(i,j,po,e) = 0
            end do
            end do

          end select
        end do

      end if

    end associate Boundary_Faces
  end do Boundaries

  !$omp end single

end subroutine ZeroDirichletEntries

!-------------------------------------------------------------------------------
!> Applies boundary conditions to RHS

subroutine ApplyBoundaryConditions(mesh, eop, nu, bv, f)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh  !< mesh partition
  class(CG_ElementOperators3D), intent(in)    :: eop   !< CG element operators
  real(RNP),                    intent(in)    :: nu    !< diffusivity
  type(BoundaryVariable),       intent(in)    :: bv(:) !< boundary conditions
  real(RNP),                    intent(inout) :: f     !< RHS

  dimension :: f(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne)

  ! local variables ............................................................

  real(RNP), pointer, contiguous :: dn_u(:,:,:)
  real(RNP), allocatable :: Mf(:,:)
  real(RNP) :: cf(3)
  integer   :: b, e, i, j, k, l, s

  !$omp single

  associate(po => eop%po, Ms => eop%w, dx => eop%dx)

    ! initialization ...........................................................

    ! metric factors
    cf(1) = dx(2)*dx(3) / 4
    cf(2) = dx(3)*dx(1) / 4
    cf(3) = dx(1)*dx(2) / 4

    ! standard face mass matrix
    allocate(Mf(0:po,0:po))
    do j = 0, po
    do i = 0, po
      Mf(i,j) = Ms(i) * Ms(j)
    end do
    end do

    Boundaries: do b = 1, size(bv)

      Boundary_Faces: associate(face => mesh%boundary(b)%face)

        select case(bv(b) % BoundaryCondition())

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
                f(0,j,k,e) = f(0,j,k,e) + cf(1) * nu * Mf(j,k) * dn_u(j,k,l)
              end do
              end do

            case(2) ! direction 1: east
              do k = 0, po
              do j = 0, po
                f(po,j,k,e) = f(po,j,k,e) + cf(1) * nu * Mf(j,k) * dn_u(j,k,l)
              end do
              end do

            case(3) ! direction 2: south
              do k = 0, po
              do i = 0, po
                f(i,0,k,e) = f(i,0,k,e) + cf(2) * nu * Mf(i,k) * dn_u(i,k,l)
              end do
              end do

            case(4) ! direction 2: north
              do k = 0, po
              do i = 0, po
                f(i,po,k,e) = f(i,po,k,e) + cf(2) * nu * Mf(i,k) * dn_u(i,k,l)
              end do
              end do

            case(5) ! direction 3: bottom
              do j = 0, po
              do i = 0, po
                f(i,j,0,e) = f(i,j,0,e) + cf(3) * nu * Mf(i,j) * dn_u(i,j,l)
              end do
              end do

            case(6) ! direction 3: top
              do j = 0, po
              do i = 0, po
                f(i,j,po,e) = f(i,j,po,e) + cf(3) * nu * Mf(i,j) * dn_u(i,j,l)
              end do
              end do

            end select
          end do

        end select

      end associate Boundary_Faces
    end do Boundaries
  end associate

  ! f = 0 in all Dirichlet points ..............................................

  call ZeroDirichletEntries(mesh, bv%BoundaryCondition(), f)

  !$omp end single

end subroutine ApplyBoundaryConditions

!===============================================================================

end module CART__CG_Elliptic_CI_BC
