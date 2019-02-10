!> summary:  IP/DG boundary contribution to RHS
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IP/DG boundary contribution to RHS
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_BCtoRHS
  implicit none

contains

!-------------------------------------------------------------------------------
!> Adds the boundary contributions of the right hand side

module subroutine BcToRHS(this, bv, f)
  class(EllipticOperator3D_IP), intent(in)    :: this
  type(BoundaryVariable),       intent(in)    :: bv(:)      !< BC
  real(RNP),                    intent(inout) :: f(:,:,:,:) !< RHS

  if (allocated(this % nu_ci)) then
    call BcToRHS_CI(this, bv, f)
  else if (allocated(this % nu_vi)) then
    call BcToRHS_VI(this, bv, f)
  end if

end subroutine BcToRHS

!-------------------------------------------------------------------------------
!> Adds the boundary contributions of the right hand side: constant isotropic

subroutine BcToRHS_CI(this, bv, f)
  class(EllipticOperator3D_IP), intent(in)    :: this
  type(BoundaryVariable),       intent(in)    :: bv(:)         !< BC
  real(RNP),                    intent(inout) :: f(0:,0:,0:,:) !< RHS

  ! local variables ............................................................

  real(RNP), pointer, contiguous :: u(:,:,:), dn_u(:,:,:)
  real(RNP), allocatable :: delta_0(:), delta_P(:), cd(:,:), Mf(:,:)
  real(RNP) :: cx(3), cf(3), mu(3)
  integer   :: b, e, i, j, k, l, s

  select type(eop => this % eop)
  class is (IP_ElementOperators1D)

    associate( P  => eop % po , nu => this % nu_ci,     &
               Ms => eop % w  , dx => this % mesh % dx, &
               Ds => eop % D  , bc => this % bc         )

      ! initialization .........................................................

      ! metric factors
      cx(:) = dx / 2
      cf(1) = cx(2)*cx(3)
      cf(2) = cx(3)*cx(1)
      cf(3) = cx(1)*cx(2)

      ! delta function
      allocate(delta_0(0:P), source = [ ONE, (ZERO, i=1,P) ])
      allocate(delta_P(0:P), source = [ (ZERO, i=1,P), ONE ])

      ! penalties
      mu(1) = eop % PenaltyFactor(dx(1))
      mu(2) = eop % PenaltyFactor(dx(2))
      mu(3) = eop % PenaltyFactor(dx(3))

      ! standard face mass matrix scaled with diffusivity
      allocate(Mf(0:P,0:P))
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      ! Dirichlet coefficients
      allocate(cd(0:P,6))
      cd(:,1) = cf(1) * ( Ds(0,:)/cx(1) + 2*mu(1) * delta_0(:))
      cd(:,2) = cf(1) * (-Ds(P,:)/cx(1) + 2*mu(1) * delta_P(:))
      cd(:,3) = cf(2) * ( Ds(0,:)/cx(2) + 2*mu(2) * delta_0(:))
      cd(:,4) = cf(2) * (-Ds(P,:)/cx(2) + 2*mu(2) * delta_P(:))
      cd(:,5) = cf(3) * ( Ds(0,:)/cx(3) + 2*mu(3) * delta_0(:))
      cd(:,6) = cf(3) * (-Ds(P,:)/cx(3) + 2*mu(3) * delta_P(:))

      Boundaries: do b = 1, size(bv)

        Boundary_Faces: associate(face => this%mesh%boundary(b)%face)

          select case(bc(b))

          case('D')

            ! Dirichlet BC .....................................................

            u => bv(b) % Component()

            do l = 1, size(face)
              e = face(l) % mesh_element % id
              s = face(l) % mesh_element % face

              select case(s)

              case(1,2) ! direction 1
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(i,s) * Mf(j,k) * nu * u(j,k,l)
                end do
                end do
                end do

              case(3,4) ! direction 2
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(j,s) * Mf(i,k) * nu * u(i,k,l)
                end do
                end do
                end do

              case(5,6) ! direction 3
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) + cd(k,s) * Mf(i,j) * nu * u(i,j,l)
                end do
                end do
                end do

              end select
            end do

          case('N')

            ! Neumann BC .......................................................

            dn_u => bv(b) % Component()

            do l = 1, size(face)
              e = face(l) % mesh_element % id
              s = face(l) % mesh_element % face

              select case(s)

              case(1) ! direction 1: west
                do k = 0, P
                do j = 0, P
                  f(0,j,k,e) = f(0,j,k,e) + cf(1) * Mf(j,k) * nu * dn_u(j,k,l)
                end do
                end do

              case(2) ! direction 1: east
                do k = 0, P
                do j = 0, P
                  f(P,j,k,e) = f(P,j,k,e) + cf(1) * Mf(j,k) * nu * dn_u(j,k,l)
                end do
                end do

              case(3) ! direction 2: south
                do k = 0, P
                do i = 0, P
                  f(i,0,k,e) = f(i,0,k,e) + cf(2) * Mf(i,k) * nu * dn_u(i,k,l)
                end do
                end do

              case(4) ! direction 2: north
                do k = 0, P
                do i = 0, P
                  f(i,P,k,e) = f(i,P,k,e) + cf(2) * Mf(i,k) * nu * dn_u(i,k,l)
                end do
                end do

              case(5) ! direction 3: bottom
                do j = 0, P
                do i = 0, P
                  f(i,j,0,e) = f(i,j,0,e) + cf(3) * Mf(i,j) * nu * dn_u(i,j,l)
                end do
                end do

              case(6) ! direction 3: top
                do j = 0, P
                do i = 0, P
                  f(i,j,P,e) = f(i,j,P,e) + cf(3) * Mf(i,j) * nu * dn_u(i,j,l)
                end do
                end do

              end select
            end do

          end select

        end associate Boundary_Faces
      end do Boundaries

    end associate
  end select

end subroutine BCtoRHS_CI

!-------------------------------------------------------------------------------
!> Adds the boundary contributions of the right hand side: variable isotropic

subroutine BcToRHS_VI(this, bv, f)
  class(EllipticOperator3D_IP), intent(in)    :: this
  type(BoundaryVariable),       intent(in)    :: bv(:)         !< BC
  real(RNP),                    intent(inout) :: f(0:,0:,0:,:) !< RHS

  ! local variables ............................................................

  real(RNP), pointer, contiguous :: u(:,:,:), dn_u(:,:,:)
  real(RNP), allocatable :: delta_0(:), delta_P(:), cd(:,:), Mf(:,:)
  real(RNP) :: cx(3), cf(3), mu(3)
  integer   :: b, e, i, j, k, l, s

  select type(eop => this % eop)
  class is (IP_ElementOperators1D)

    associate( P  => this % eop % po , nu => this % nu_vi,     &
               Ms => this % eop % w  , dx => this % mesh % dx, &
               Ds => this % eop % D  , bc => this % bc         )

      ! initialization .........................................................

      ! metric factors
      cx(:) = dx / 2
      cf(1) = cx(2)*cx(3)
      cf(2) = cx(3)*cx(1)
      cf(3) = cx(1)*cx(2)

      ! delta function
      allocate(delta_0(0:P), source = [ ONE, (ZERO, i=1,P) ])
      allocate(delta_P(0:P), source = [ (ZERO, i=1,P), ONE ])

      ! penalties
      mu(1) = eop % PenaltyFactor(dx(1))
      mu(2) = eop % PenaltyFactor(dx(2))
      mu(3) = eop % PenaltyFactor(dx(3))

      ! standard face mass matrix scaled with diffusivity
      allocate(Mf(0:P,0:P))
      do j = 0, P
      do i = 0, P
        Mf(i,j) = Ms(i) * Ms(j)
      end do
      end do

      ! Dirichlet coefficients
      allocate(cd(0:P,6))
      cd(:,1) = cf(1) * ( Ds(0,:)/cx(1) + 2*mu(1) * delta_0(:))
      cd(:,2) = cf(1) * (-Ds(P,:)/cx(1) + 2*mu(1) * delta_P(:))
      cd(:,3) = cf(2) * ( Ds(0,:)/cx(2) + 2*mu(2) * delta_0(:))
      cd(:,4) = cf(2) * (-Ds(P,:)/cx(2) + 2*mu(2) * delta_P(:))
      cd(:,5) = cf(3) * ( Ds(0,:)/cx(3) + 2*mu(3) * delta_0(:))
      cd(:,6) = cf(3) * (-Ds(P,:)/cx(3) + 2*mu(3) * delta_P(:))

      Boundaries: do b = 1, size(bv)

        Boundary_Faces: associate(face => this%mesh%boundary(b)%face)

          select case(bc(b))

          case('D')

            ! Dirichlet BC .....................................................

            u => bv(b) % Component()

            do l = 1, size(face)
              e = face(l) % mesh_element % id
              s = face(l) % mesh_element % face

              select case(s)

              case(1) ! direction 1: west
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(i,s) * Mf(j,k) * nu(0,j,k,e) * u(j,k,l)
                end do
                end do
                end do

              case(2) ! direction 1: east
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(i,s) * Mf(j,k) * nu(P,j,k,e) * u(j,k,l)
                end do
                end do
                end do

              case(3) ! direction 2: south
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(j,s) * Mf(i,k) * nu(i,0,k,e) * u(i,k,l)
                end do
                end do
                end do

              case(4) ! direction 2: north
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(j,s) * Mf(i,k) * nu(i,P,k,e) * u(i,k,l)
                end do
                end do
                end do

              case(5) ! direction 3: bottom
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(k,s) * Mf(i,j) * nu(i,j,0,e) * u(i,j,l)
                end do
                end do
                end do

              case(6) ! direction 3: top
                do k = 0, P
                do j = 0, P
                do i = 0, P
                  f(i,j,k,e) = f(i,j,k,e) &
                             + cd(k,s) * Mf(i,j) * nu(i,j,P,e) * u(i,j,l)
                end do
                end do
                end do

              end select
            end do

          case('N')

            ! Neumann BC .......................................................

            dn_u => bv(b) % Component()

            do l = 1, size(face)
              e = face(l) % mesh_element % id
              s = face(l) % mesh_element % face

              select case(s)

              case(1) ! direction 1: west
                do k = 0, P
                do j = 0, P
                  f(0,j,k,e) = f(0,j,k,e) &
                             + cf(1) * Mf(j,k) * nu(0,j,k,e) * dn_u(j,k,l)
                end do
                end do

              case(2) ! direction 1: east
                do k = 0, P
                do j = 0, P
                  f(P,j,k,e) = f(P,j,k,e) &
                             + cf(1) * Mf(j,k) * nu(P,j,k,e) * dn_u(j,k,l)
                end do
                end do

              case(3) ! direction 2: south
                do k = 0, P
                do i = 0, P
                  f(i,0,k,e) = f(i,0,k,e) &
                             + cf(2) * Mf(i,k) * nu(i,0,k,e) * dn_u(i,k,l)
                end do
                end do

              case(4) ! direction 2: north
                do k = 0, P
                do i = 0, P
                  f(i,P,k,e) = f(i,P,k,e) &
                             + cf(2) * Mf(i,k) * nu(i,P,k,e) * dn_u(i,k,l)
                end do
                end do

              case(5) ! direction 3: bottom
                do j = 0, P
                do i = 0, P
                  f(i,j,0,e) = f(i,j,0,e) &
                             + cf(3) * Mf(i,j) * nu(i,j,0,e) * dn_u(i,j,l)
                end do
                end do

              case(6) ! direction 3: top
                do j = 0, P
                do i = 0, P
                  f(i,j,P,e) = f(i,j,P,e) &
                             + cf(3) * Mf(i,j) * nu(i,j,P,e) * dn_u(i,j,l)
                end do
                end do

              end select
            end do

          end select

        end associate Boundary_Faces
      end do Boundaries

    end associate
  end select

end subroutine BCtoRHS_VI

!===============================================================================

end submodule MP_BCtoRHS
