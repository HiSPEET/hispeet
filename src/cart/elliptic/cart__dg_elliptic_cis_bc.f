!> summary:  Boundary treatment for DG elliptic operator,
!>           constant isotropic + structured
!> author:   Joerg Stiller
!> date:     2016/12/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Boundary treatment for DG elliptic operator
!===============================================================================

module CART__DG_Elliptic_CIS_BC

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error

  use CART__DG_Element_Operators
  use CART__Mesh_Partition
  use CART__Boundary_Variable

  implicit none
  private

  public :: Apply_BC_to_RHS

contains

!-------------------------------------------------------------------------------
!> Adds boundary contributions of the right hand side

subroutine Apply_BC_to_RHS(mesh, eop, nu, bv, f)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)    :: mesh  !< mesh partition
  class(DG_ElementOperators3D), intent(in)    :: eop   !< DG element operators
  real(RNP),                    intent(in)    :: nu    !< diffusivity
  type(BoundaryVariable),       intent(in)    :: bv(6) !< boundary conditions
  real(RNP),                    intent(inout) :: f     !< RHS

  dimension :: f(0:eop%po, 0:eop%po, 0:eop%po, mesh%ne1, mesh%ne2, mesh%ne3)

  target :: bv

  ! local variables ............................................................

  real(RNP) :: cx(3), cd, cf, cp
  real(RNP), allocatable :: delta_0(:), delta_P(:)
  real(RNP), pointer, contiguous :: u(:,:,:), dn_u(:,:,:)

  integer :: po, ne1, ne2, ne3
  integer :: i, j, k, l, m, n, q

  ! initialization .............................................................

  ! dimensions
  po  = eop  % po
  ne1 = mesh % ne1
  ne2 = mesh % ne2
  ne3 = mesh % ne3

  ! metric coefficients
  cx = eop%dx / 2

  ! delta_i0
  allocate(delta_0(0:po), source = ZERO)
  delta_0(0) = ONE

  ! delta_iP
  allocate(delta_P(0:po), source = ZERO)
  delta_P(po) = ONE

  associate(Ms => eop%w, Ds => eop%D, mu => eop%mu)

    ! metric coefficients for west and east boundaries .........................

    cf = cx(2) * cx(3)
    cd = cf / cx(1)
    cp = cf * 2 * mu(1)

    ! west boundary ............................................................

    WEST: if (mesh%boundary(1)%nf > 0) then

      select case( bv(1) % BoundaryCondition() )

      case('D') ! Dirichlet

        u => bv(1) % Component()

        q = 0
        do n = 1, ne3
        do m = 1, ne2

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,1,m,n) = f(i,j,k,1,m,n)                      &
                           + Ms(j) * Ms(k) * nu                  &
                             * (cd * Ds(0,i) + cp * delta_0(i))  &
                             * u(j,k,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(1) % Component()

        q = 0
        do n = 1, ne3
        do m = 1, ne2

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
            f(0,j,k,1,m,n) = f(0,j,k,1,m,n) &
                           + cf * Ms(j) * Ms(k) * nu * dn_u(j,k,q)
          end do
          end do
        end do
        end do

      end select

    end if WEST

    ! east boundary ............................................................

    EAST: if (mesh%boundary(2)%nf > 0) then

      select case( bv(2) % BoundaryCondition() )

      case('D') ! Dirichlet

        u => bv(2) % Component()

        q = 0
        do n = 1, ne3
        do m = 1, ne2

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,ne1,m,n) = f(i,j,k,ne1,m,n)                      &
                             + Ms(j) * Ms(k) * nu                    &
                               * (-cd * Ds(po,i) + cp * delta_P(i))  &
                               * u(j,k,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(2) % Component()

        q = 0
        do n = 1, ne3
        do m = 1, ne2

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
            f(po,j,k,ne1,m,n) = f(po,j,k,ne1,m,n) &
                              + cf * Ms(j) * Ms(k) * nu * dn_u(j,k,q)
          end do
          end do
        end do
        end do

      end select

    end if EAST

    ! metric coefficients for south and north boundaries .......................

    cf = cx(3) * cx(1)
    cd = cf / cx(2)
    cp = cf * 2 * mu(2)

    ! south boundary ...........................................................

    SOUTH: if (mesh%boundary(3)%nf > 0) then

      select case( bv(3) % BoundaryCondition() )

      case('D') ! Dirichlet

        u => bv(3) % Component()

        q = 0
        do n = 1, ne3
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,l,1,n) = f(i,j,k,l,1,n)                      &
                           + Ms(i) * Ms(k) * nu                  &
                             * (cd * Ds(0,j) + cp * delta_0(j))  &
                             * u(i,k,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(3) % Component()

        q = 0
        do n = 1, ne3
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do i = 0, po
            f(i,0,k,l,1,n) = f(i,0,k,l,1,n) &
                           + cf * Ms(i) * Ms(k) * nu * dn_u(i,k,q)
          end do
          end do
        end do
        end do

      end select

    end if SOUTH

    ! north boundary ...........................................................

    NORTH: if (mesh%boundary(4)%nf > 0) then

      select case( bv(4) % BoundaryCondition() )

      case('D') ! Dirichlet

        u => bv(4) % Component()

        q = 0
        do n = 1, ne3
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,l,ne2,n) = f(i,j,k,l,ne2,n)                      &
                             + Ms(i) * Ms(k) * nu                    &
                               * (-cd * Ds(po,j) + cp * delta_P(j))  &
                               * u(i,k,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(4) % Component()

        q = 0
        do n = 1, ne3
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do i = 0, po
            f(i,po,k,l,ne2,n) = f(i,po,k,l,ne2,n) &
                              + cf * Ms(i) * Ms(k) * nu * dn_u(i,k,q)
          end do
          end do
        end do
        end do

      end select

    end if NORTH

    ! metric coefficients for bottom and top boundaries ........................

    cf = cx(1) * cx(2)
    cd = cf / cx(3)
    cp = cf * 2 * mu(3)

    ! bottom boundary ..........................................................

    BOTTOM: if (mesh%boundary(5)%nf > 0) then

      select case( bv(5) % BoundaryCondition() )

      case('D') ! Dirichlet

        u => bv(5) % Component()

        q = 0
        do m = 1, ne2
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,l,m,1) = f(i,j,k,l,m,1)                      &
                           + Ms(i) * Ms(j) * nu                  &
                             * (cd * Ds(0,k) + cp * delta_0(k))  &
                             * u(i,j,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(5) % Component()

        q = 0
        do m = 1, ne2
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do j = 0, po
          do i = 0, po
            f(i,j,0,l,m,1) = f(i,j,0,l,m,1) &
                           + cf * Ms(i) * Ms(j) * nu * dn_u(i,j,q)
          end do
          end do
        end do
        end do

      end select

    end if BOTTOM

    ! top boundary .............................................................

    TOP: if (mesh%boundary(6)%nf > 0) then

      select case( bv(6) % BoundaryCondition() )

      ! Dirichlet
      case('D')

        u => bv(6) % Component()

        q = 0
        do m = 1, ne2
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do k = 0, po
          do j = 0, po
          do i = 0, po
            f(i,j,k,l,m,ne3) = f(i,j,k,l,m,ne3)                      &
                             + Ms(i) * Ms(j) * nu                    &
                               * (-cd * Ds(po,k) + cp * delta_P(k))  &
                               * u(i,j,q)
          end do
          end do
          end do
        end do
        end do

      case('N') ! Neumann

        dn_u => bv(6) % Component()

        q = 0
        do m = 1, ne2
        do l = 1, ne1

          ! increase face counter
          q = q + 1

          do j = 0, po
          do i = 0, po
            f(i,j,po,l,m,ne3) = f(i,j,po,l,m,ne3) &
                              + cf * Ms(i) * Ms(j) * nu * dn_u(i,j,q)
          end do
          end do
        end do
        end do

      end select

    end if TOP

  end associate

end subroutine Apply_BC_to_RHS

!===============================================================================

end module CART__DG_Elliptic_CIS_BC
