!> summary:  Base class for compressible Navier-Stokes problems
!> author:   Joerg Stiller, Benedikt Wex
!> date:     2023/11/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - investigate, select and implement simple + robust entropy fix
!>   - develop IP-DG formulation of viscous terms
!>   - develop IP-DG formulation of streamline diffusion (SD) terms
!>   - implement viscous term
!>   - implement SD terms
!>   - upcoming:
!>       * boundary conditions
!>       * residual
!>       * implicit solver
!>       * example problems
!>           + acoustic wave
!>           + shock tube
!>           + Shu-Osher test case
!===============================================================================

module CL__Problem__CNS__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO, HALF, THIRD, TWO
  use Execution_Control
  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D compressible Navier-Stokes problems
  !>
  !> TBD

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_CNS_1D

    real(RNP) :: r_gas   !< specific gas constant
    real(RNP) :: gamma   !< ratio of specific heats
    real(RNP) :: c_p     !< specific heat for const pressure
    real(RNP) :: c_v     !< specific heat for const volume
    real(RNP) :: eta     !< dynamic viscosity
    real(RNP) :: prandtl !< Prandtl
    real(RNP) :: lambda  !< thermal conductivity

  contains

    procedure :: GetConvectionTerm
   !#! Not jet Implemented
   ! procedure :: GetMaxVelocity
   ! procedure :: GetMaxDiffusivity

    ! CNS specific procedures ..................................................

    ! transformations
    procedure :: ConservativeToPrimitive
    procedure :: PrimitiveToConservative
    procedure :: ConservativeToCharacteristic
    procedure :: CharacteristicToConservative

    ! convection
    procedure :: ConvectiveFlux
    procedure :: ConvectiveJacobian
    procedure :: ConvectiveEigensystem
    procedure :: NumericalConvectiveFlux

    ! diffusion
    !#! the following routines will be provided in a separate submodule
    ! procedure :: PhysicalDiffusivity
    ! procedure :: StreamlineDiffusivity
    ! procedure :: GetDiffusionTerm
    ! procedure :: GetSDTerm
    ! procedure :: DiffusionSolver

  end type CL_Problem_CNS_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

  subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator
    real(RNP),                intent(in)  :: bv (:,:)
    real(RNP), contiguous,    intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,    intent(out) :: r_c(0:,:,:)

    real(RNP), allocatable :: MD_t(:,:), f_q(:,:), u_q(:,:), h_c(:,:)
    real(RNP), allocatable :: u_l(:,:), u_r(:,:)
    logical :: use_interpolation
    integer :: e, i, j, j0

    associate( eop      => cl_operator % eop      &
             , qop      => cl_operator % qop      &
             , iop      => cl_operator % iop_uq   &
             , po       => cl_operator % eop % po &
             , qo       => cl_operator % qop % po &
             , ne       => cl_operator % ne       &
             , activity => cl_operator % activity )

      ! initialization .........................................................

      !$omp master
      allocate(MD_t(0:po, 0:qo), f_q(0:qo,3), u_q(0:qo,3), h_c(0:ne,3))

      use_interpolation = qo /= po

      if (use_interpolation) then
        j0 = lbound(iop%A,1)
        do i = 0, po
        do j = 0, qo
          MD_t(i,j) = qop%w(j) * dot_product(iop%A(j0+j,:), eop%D(:,i))
        end do
        end do
      else
        do i = 0, po
        do j = 0, po
          MD_t(i,j) = eop%w(j) * eop%D(j,i)
        end do
        end do
      end if

      ! element integrals ......................................................

      !$omp do
      do e = 1, ne
        if (activity(e) > 0) then

          ! compute fluxes in quadrature points
          if (use_interpolation) then
            u_q(:,1) = matmul(iop%A, u(:,e,1))
            u_q(:,2) = matmul(iop%A, u(:,e,2))
            u_q(:,3) = matmul(iop%A, u(:,e,3))
          else
            u_q(:,1) = u(:,e,1)
            u_q(:,2) = u(:,e,2)
            u_q(:,3) = u(:,e,3)
          end if

          ! compute flux
          do i = 0, po
            f_q(i,:) = ConvectiveFlux(this, u_q(i,:))
          end do

          ! apply mass-weighted transposed diff-matrix
          r_c(:,e,1) = matmul(MD_t, f_q(:,1))
          r_c(:,e,2) = matmul(MD_t, f_q(:,2))
          r_c(:,e,3) = matmul(MD_t, f_q(:,3))

        else
          r_c(:,e,:) = 0
        end if
      end do

      ! element boundary fluxes ................................................

      allocate(u_l(0:ne,3), u_r(0:ne,3))

      u_l(1:ne  ,:) = u(po,1:ne,:)
      u_r(0:ne-1,:) = u( 0,1:ne,:)

      ! left boundary
      select case(this % bc(1))
        case('P')
          u_l(0,:) = u(po,ne,:)
        case default
          u_l(0,:) = u(0,1,:)
      end select

      ! right boundary
      select case(this % bc(2))
        case('P')
          u_r(ne,:) = u(0,1,:)
        case default
          u_r(ne,:) = u(po,ne,:)
      end select

      do e = 0, ne
        h_c(e,:) = NumericalConvectiveFlux(this, u_l(e,:), u_r(e,:))
      end do

      do e = 1, ne
        if (activity(e) > 0) then
          r_c( 0,e,:) = r_c( 0,e,:) + h_c(e-1,:)
          r_c(po,e,:) = r_c(po,e,:) - h_c(e,:)
        end if
      end do

      ! finalization ...........................................................

      deallocate(MD_t, f_q, u_q, h_c)
      !$omp end master

    end associate

  end subroutine GetConvectionTerm

  !=============================================================================
  ! CNS specific routines

  !-----------------------------------------------------------------------------
  !> Convert conservative variables into primitive ones

  pure subroutine ConservativeToPrimitive(this, u, v, T, p, s)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: u(:) !< conservative variables u(3)
    real(RNP),  optional,     intent(out) :: v    !< velocity
    real(RNP),  optional,     intent(out) :: T    !< temperature
    real(RNP),  optional,     intent(out) :: p    !< pressure
    real(RNP),  optional,     intent(out) :: s    !< specific entropy

    real(RNP) ::  T_, p_

    p_ = (this%gamma - 1) * (u(3) - HALF*u(2)*(u(2)/u(1)))
    T_ = p_ / (u(1)*this%r_gas)

    if (present(v))    v   = u(2)/u(1)
    if (present(p))    p   = p_
    if (present(T))    T   = T_
    if (present(s))    s   = this%c_p * log(T_) - this%r_gas * log(p_)

  end subroutine ConservativeToPrimitive

  !-----------------------------------------------------------------------------
  !> Converts primitive variables into conservative ones

  pure subroutine PrimitiveToConservative(this, v, T, p, u)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: v
    real(RNP),                intent(in)  :: T
    real(RNP),                intent(in)  :: p
    real(RNP),                intent(out) :: u(:)

    u = (p/(this%r_gas*T)) * [ONE, v, this%c_v * T + HALF*v*v]

  end subroutine PrimitiveToConservative

  !-----------------------------------------------------------------------------
  !> Converts conservative variables into characteristic ones

  pure subroutine ConservativeToCharacteristic(this, u, z)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: u(:)
    real(RNP),                intent(out) :: z(:)

    real(RNP) :: v, p, T, a

    v    = u(2)/u(1)
    p    = (this%gamma - 1) * (u(3) - HALF*u(2)*v)
    T    = p/(u(1)*this%r_gas)
    a    = sqrt(this%gamma * this%r_gas * T)

    z(1) = (TWO*a)/(this%gamma - 1) - v
    z(2) = this%c_p * log(T) - this%r_gas * log(p)
    z(3) = (TWO*a)/(this%gamma - 1) + v

  end subroutine ConservativeToCharacteristic

  !-----------------------------------------------------------------------------
  !> Converts conservative variable into conservative ones

  pure subroutine CharacteristicToConservative(this, z, u)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: z(:)
    real(RNP),                intent(out) :: u(:)

    real(RNP) :: v, T, p

    v = (z(3) - z(1)) * HALF
    T = (this%gamma - 1)**2 / (16 * this%r_gas * this%gamma) * (z(1) + z(3))**2
    p = exp(-z(2)/this%r_gas) * T**(this%gamma / (this%gamma - 1))

    u = p / (this%r_gas * T) * [ONE, v, this%c_v*T + HALF*v*v]

  end subroutine CharacteristicToConservative

  !-----------------------------------------------------------------------------
  !> Calculates the convective flux with given conservative variables

  pure function ConvectiveFlux(this, u) result(f_c)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:)
    real(RNP) :: f_c(3)

    real(RNP) ::  v, p

    v = u(2) / u(1)
    p = (this%gamma - 1) * (u(3) - u(2)*u(2) / (2*u(1)))
    f_c = v * u + [ ZERO, p, v*p ]

  end function ConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Calculates the Jacobian ot the Convective Term

  pure function ConvectiveJacobian(this, u) result(A)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:)
    real(RNP) :: A(3,3)

    real(RNP) :: v, e_k, h_t

    associate(gamma => this % gamma)

      v   = u(2)/u(1)
      e_k = HALF*v*v
      h_t = gamma*u(3)/u(1) - (gamma-1)*e_k

      A(1,:) = [ZERO                  	 ,  ONE                   ,  ZERO     ]
      A(2,:) = [(gamma-3)*e_k            , -(gamma-3)*v           ,  gamma-1  ]
      A(3,:) = [-v*(h_t - (gamma-1)*e_k) ,  h_t - 2*(gamma-1)*e_k ,  gamma*v  ]

    end associate

  end function ConvectiveJacobian

  !-----------------------------------------------------------------------------
  !> Calculates the Eigendecomposition of the Convective Jacobian

  pure subroutine ConvectiveEigensystem(this, u, Lambda, R, L)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: u(:)
    real(RNP), optional,      intent(out) :: Lambda(:)
    real(RNP), optional,      intent(out) :: R(:,:)
    real(RNP), optional,      intent(out) :: L(:,:)

    real(RNP) :: v, rho, p, aa, a, T, h, c1

    associate(gamma => this%gamma)

      v = u(2)/u(1)
      rho = u(1)
      p = (gamma - 1) * (u(3) - HALF*u(2)*v)
      aa = gamma * (p/rho)
      a = sqrt(aa)
      T = p/(u(1)*this%r_gas)
      h = this%c_p * T

      if(present(lambda)) then
        lambda = [v-a , v , v+a ]
      end if

      if(present(R)) then
        R(1,1:3) = [ONE     , ONE      , ONE     ]
        R(2,1:3) = [v-a     , v        , v+a     ]
        R(3,1:3) = [h-(v*a) , HALF*v*v , h+(v*a) ]
      end if

      if(present(L)) then
        c1 = a/(gamma-1)
        L(1,1:3) = [((v*v/2)+(v*c1)) , (-v-c1) ,  ONE ]
        L(2,1:3) = [(2*h-v*v)        , 2*v     , -TWO ]
        L(3,1:3) = [((v*v/2)-(v*c1)) , (-v+c1) ,  ONE ]
        L = HALF/h * L
      end if
    end associate

  end subroutine ConvectiveEigensystem

  !-----------------------------------------------------------------------------
  !> Calculates the Numerical Roe-Flux without entropy fix

  pure function NumericalConvectiveFlux(this, u_l, u_r) result(h_c)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u_l(:)
    real(RNP), intent(in) :: u_r(:)
    real(RNP) :: h_c(3)

    real(RNP) :: u_roe(3), Lambda(3)
    real(RNP) :: abs_A(3,3), R(3,3), L(3,3)
    integer   :: i,j

    !JS! Entropie-Fix dann z.B. durch Anpassung der Eigenwerte

  	!Calculation of the Roe-Average Value
    u_roe = RoeAverage(this, u_l, u_r)

    !Calculation of the Absolute Jacobian
    call ConvectiveEigensystem(this, u_roe, Lambda, R, L)
    Lambda = abs(Lambda)
    do i = 1, 3
    do j = 1, 3
      abs_A(i,j) = sum(R(i,:) * Lambda(:) * L(:,j))
    end do
    end do

    ! Calculation of the Roe-Flux
    h_c = ( ConvectiveFlux(this, u_l) &
          + ConvectiveFlux(this, u_r) &
          + matmul(abs_A, u_l - u_r)  &
          ) * HALF

  end function NumericalConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Calculates the Roe average of the Conservative Variables at two Points

  pure function RoeAverage(this, u_l, u_r) result(u_m)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u_l(:)
    real(RNP), intent(in) :: u_r(:)
    real(RNP) :: u_m(3)

    real(RNP) ::  r_l, v_l, p_l, T_l, H_l
    real(RNP) ::  r_r, v_r, p_r, T_r, H_r
    real(RNP) ::  rho_m, v_m, p_m, T_m, H_m

    r_l = sqrt(u_l(1))
    r_r = sqrt(u_r(1))
    rho_m = r_l * r_r

    call ConservativeToPrimitive(this, u_l, v=v_l, T=T_l, p=p_l)
    call ConservativeToPrimitive(this, u_r, v=v_r, T=T_r, p=p_r)
    v_m = (r_l*v_l + r_r*v_r) / (r_l + r_r)

    H_l = this%c_p * T_l + HALF * v_l * v_l
    H_r = this%c_p * T_r + HALF * v_r * v_r
    H_m = (r_l*H_l + r_r*H_r) / (r_l + r_r)

    u_m(1) = rho_m
    u_m(2) = rho_m * v_m
    u_m(3) = rho_m * (H_m + (this%gamma -1) * HALF * v_m * v_m) / this%gamma

  end function RoeAverage
  !=============================================================================

end module CL__Problem__CNS__1D

