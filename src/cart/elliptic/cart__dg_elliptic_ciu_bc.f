!> summary:  Boundary treatment for DG elliptic operator,
!>           constant isotropic + unstructured
!> author:   Joerg Stiller
!> date:     2018/09/16
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   * not threadsafe, i.e. not parallelizable with OpenMP, since more than
!>     one element face may belong to the same boundary
!> @endnote
!===============================================================================

module CART__DG_Elliptic_CIU_BC

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error

  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__Boundary_Variable

  implicit none
  private

  public :: ApplyBoundaryConditions

contains

!-------------------------------------------------------------------------------
!> Adds boundary contributions of the right hand side

subroutine ApplyBoundaryConditions(mesh, eop, nu, bv, f)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh  !< mesh partition
  class(DG_ElementOperators3D), intent(in)    :: eop   !< DG element operators
  real(RNP),                    intent(in)    :: nu    !< diffusivity
  type(BoundaryVariable),       intent(in)    :: bv(:) !< boundary conditions
  real(RNP),                    intent(inout) :: f     !< RHS

  dimension :: f(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne)

  ! local variables ............................................................

  real(RNP), pointer, contiguous :: u(:,:,:), dn_u(:,:,:)
  real(RNP), allocatable :: delta_0(:), delta_P(:), cd(:,:), nu_Mf(:,:)
  real(RNP) :: cx(3), cf(3)
  integer   :: b, e, i, j, k, l, s

  !$omp single

  associate(po => eop%po, Ms => eop%w, Ds => eop%D, mu => eop%mu)

    ! initialization ...........................................................

    ! metric factors
    cx(:) = eop%dx / 2
    cf(1) = cx(2)*cx(3)
    cf(2) = cx(3)*cx(1)
    cf(3) = cx(1)*cx(2)

    ! delta_i0
    allocate(delta_0(0:po), source = ZERO)
    delta_0(0) = ONE

    ! delta_iP
    allocate(delta_P(0:po), source = ZERO)
    delta_P(po) = ONE

    ! standard face mass matrix scaled with diffusivity
    allocate(nu_Mf(0:po,0:po))
    do j = 0, po
    do i = 0, po
      nu_Mf(i,j) = nu * Ms(i) * Ms(j)
    end do
    end do

    ! Dirichlet coefficients
    allocate(cd(0:po,6))
    cd(:,1) = cf(1) * ( Ds( 0,:)/cx(1) + 2*mu(1) * delta_0(:))
    cd(:,2) = cf(1) * (-Ds(po,:)/cx(1) + 2*mu(1) * delta_P(:))
    cd(:,3) = cf(2) * ( Ds( 0,:)/cx(2) + 2*mu(2) * delta_0(:))
    cd(:,4) = cf(2) * (-Ds(po,:)/cx(2) + 2*mu(2) * delta_P(:))
    cd(:,5) = cf(3) * ( Ds( 0,:)/cx(3) + 2*mu(3) * delta_0(:))
    cd(:,6) = cf(3) * (-Ds(po,:)/cx(3) + 2*mu(3) * delta_P(:))

    Boundaries: do b = 1, size(bv)

      Boundary_Faces: associate(face => mesh%boundary(b)%face)

        select case(bv(b) % BoundaryCondition())

        case('D')

          ! Dirichlet BC .......................................................

          u => bv(b) % Component()

          do l = 1, size(face)
            e = face(l) % mesh_element % id
            s = face(l) % mesh_element % face

            select case(s)

            case(1,2) ! direction 1
              do k = 0, po
              do j = 0, po
              do i = 0, po
                f(i,j,k,e) = f(i,j,k,e) + cd(i,s) * nu_Mf(j,k) * u(j,k,l)
              end do
              end do
              end do

            case(3,4) ! direction 2
              do k = 0, po
              do j = 0, po
              do i = 0, po
                f(i,j,k,e) = f(i,j,k,e) + cd(j,s) * nu_Mf(i,k) * u(i,k,l)
              end do
              end do
              end do

            case(5,6) ! direction 3
              do k = 0, po
              do j = 0, po
              do i = 0, po
                f(i,j,k,e) = f(i,j,k,e) + cd(k,s) * nu_Mf(i,j) * u(i,j,l)
              end do
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

  !$omp end single

end subroutine ApplyBoundaryConditions

!===============================================================================

end module CART__DG_Elliptic_CIU_BC
