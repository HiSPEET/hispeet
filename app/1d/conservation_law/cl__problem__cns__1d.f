!> summary:  Base class for compressible Navier-Stokes problems
!> author:   Joerg Stiller, Benedikt Wex
!> date:     2023/11/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - investigate, select and implement simple + robust entropy fix
!>   - upcoming:
!>       * boundary conditions
!>           + in convective term
!>           + in diffusive/SD term
!>       * example problems
!>           + acoustic wave
!>           + shock tube
!>           + Shu-Osher test case
!===============================================================================

module CL__Problem__CNS__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO, HALF, THIRD, TWO
  use Execution_Control
  use Logging_Levels
  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_1D
  public :: CL_Problem_CNS_Options_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D compressible Navier-Stokes problems
  !>
  !> TBD

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_CNS_1D

    real(RNP) :: r_gas    !< specific gas constant
    real(RNP) :: gamma    !< ratio of specific heats
    real(RNP) :: c_v      !< specific heat at constant volume
    real(RNP) :: c_p      !< specific heat at constant pressure
    real(RNP) :: eta      !< dynamic viscosity
    real(RNP) :: prandtl  !< Prandtl number
    real(RNP) :: lambda   !< thermal conductivity

    integer   :: n_krylov !< dimension of GMRES Krylov subspaces

  contains

    procedure :: HasDiffusion
    procedure :: GetConvectionTerm
    procedure :: GetDiffusionTerm
    procedure :: GetSDTerm
    procedure :: DiffusionSolver
    procedure :: GetMaxVelocity
    procedure :: GetMaxDiffusivity

    ! CNS specific procedures ..................................................

    ! initialization
    procedure :: Init_CL_Problem_CNS_1D

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
    procedure :: PhysicalDiffusivity
    procedure :: StreamlineDiffusivity
    procedure :: GetHybridDiffusionTerm

  end type CL_Problem_CNS_1D

  !-----------------------------------------------------------------------------
  !> Compressible Navier-Stokes options

  type CL_Problem_CNS_Options_1D

    ! fluid properties
    real(RNP) :: r_gas    = 287.280E+0_RNP !< specific gas constant
    real(RNP) :: gamma    =   1.400E+0_RNP !< ratio of specific heats
    real(RNP) :: eta      =   1.800E-5_RNP !< dynamic viscosity
    real(RNP) :: prandtl  =   0.691E+0_RNP !< Prandtl number

    ! solver options
    integer :: n_krylov = 5 !< dimension of GMRES Krylov subspaces

  end type CL_Problem_CNS_Options_1D

  !=============================================================================
  ! external module procedures

  interface

    !---------------------------------------------------------------------------
    !> Returns the physical diffusivity matrix

    pure module function PhysicalDiffusivity(this, u) result(A_pd)
      class(CL_Problem_CNS_1D), intent(in) :: this
      real(RNP), intent(in) :: u(:) !< conservative variables
      real(RNP) :: A_pd(3,3)
    end function PhysicalDiffusivity

    !---------------------------------------------------------------------------
    !> Returns the streamline diffusivity matrix

    pure module function StreamlineDiffusivity(this, theta, u) result(A_d)
      class(CL_Problem_CNS_1D), intent(in) :: this
      real(RNP), intent(in) :: u(:)  !< conservative variables
      real(RNP), intent(in) :: theta !< streamline-diffusion time scale
      real(RNP) :: A_d(3,3)
    end function StreamlineDiffusivity

    !---------------------------------------------------------------------------
    !> Unified CNS diffusion term combining physical and streamline contributions

    module subroutine GetHybridDiffusionTerm &
        (this, cl_operator, comp, theta, bv, u_0, u, r_d)

      class(CL_Problem_CNS_1D), intent(in)  :: this
      class(CL_Operator_1D),    intent(in)  :: cl_operator !< spatial operators
      character,                intent(in)  :: comp        !< composition flag
      real(RNP),                intent(in)  :: theta       !< SD time scale
      real(RNP), optional,      intent(in)  :: bv(:,:)     !< boundary values
      real(RNP), contiguous,    intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
      real(RNP), contiguous,    intent(in)  :: u  (0:,:,:) !< u(x,t)
      real(RNP), contiguous,    intent(out) :: r_d(0:,:,:) !< diffusion RHS

    end subroutine GetHybridDiffusionTerm

    !---------------------------------------------------------------------------
    !> Implicit CNS diffusion solver

    module subroutine DiffusionSolver &
        (this, cl_operator, dt, theta, bv, f, u_0, u, method, i_max, r_red, r_max)

      class(CL_Problem_CNS_1D), intent(in)    :: this
      class(CL_Operator_1D),    intent(in)    :: cl_operator
      real(RNP),                intent(in)    :: dt          !< ∆t = t - t₀
      real(RNP),                intent(in)    :: theta       !< SD time scale θ
      real(RNP),                intent(in)    :: bv (:,:)    !< boundary values
      real(RNP), contiguous,    intent(in)    :: f  (0:,:,:) !< sources
      real(RNP), contiguous,    intent(in)    :: u_0(0:,:,:) !< frozen solution
      real(RNP), contiguous,    intent(inout) :: u  (0:,:,:) !< approx solution
      integer,                  intent(in)    :: method      !< solution method
      integer,                  intent(in)    :: i_max       !< max num iterations
      real(RNP), optional,      intent(in)    :: r_red       !< residual reduction
      real(RNP), optional,      intent(in)    :: r_max       !< max residual

    end subroutine DiffusionSolver

  end interface

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the base type

  subroutine Init_CL_Problem_CNS_1D(this, opt)
    class(CL_Problem_CNS_1D),         intent(inout) :: this
    class(CL_Problem_CNS_Options_1D), intent(in)    :: opt

    this % nc = 3

    this % r_gas    = opt % r_gas
    this % gamma    = opt % gamma
    this % eta      = opt % eta
    this % prandtl  = opt % prandtl
    this % n_krylov = opt % n_krylov

    this % c_v      = this % r_gas / (this % gamma - 1)
    this % c_p      = this % c_v * this % gamma
    this % lambda   = this % c_p * this % eta / this % prandtl

  end subroutine Init_CL_Problem_CNS_1D

  !-----------------------------------------------------------------------------
  !> Query whether problem has non vanishing physical diffusion

  logical function HasDiffusion(this)
    class(CL_Problem_CNS_1D), intent(in) :: this

    HasDiffusion = this % eta > 0

  end function HasDiffusion

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
          do i = 0, qo
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
        case('S')
          u_l(0,:) = bv(:,1)
        case default
          u_l(0,:) = u(0,1,:)
      end select

      ! right boundary
      select case(this % bc(2))
        case('P')
          u_r(ne,:) = u(0,1,:)
        case('S')
          u_r(ne,:) = bv(:,2)
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

  !-----------------------------------------------------------------------------
  !> Diffusive contribution to RHS of DG-SEM formulation

  subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator
    real(RNP),                intent(in)  :: bv (:,:)    !< boundary values
    real(RNP), contiguous,    intent(in)  :: u  (0:,:,:) !< u(x,t)
    real(RNP), contiguous,    intent(out) :: r_d(0:,:,:) !< diffusion RHS

    character, parameter :: comp  = 'P'
    real(RNP), parameter :: theta =  0

    call GetHybridDiffusionTerm(this, cl_operator, comp, theta, bv, u, u, r_d)

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Streamline-diffusion contribution to RHS of DG-SEM formulation
  !>
  !> Evaluates the weak form of the streamline-diffusion operators for `u`
  !> using `u₀` for computing the streamline diffusivity.

  subroutine GetSDTerm(this, cl_operator, theta, bv, u_0, u, r_sd)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator
    real(RNP),                intent(in)  :: theta        !< SD time scale θ
    real(RNP),                intent(in)  :: bv  (:,:)    !< boundary values
    real(RNP), contiguous,    intent(in)  :: u_0 (0:,:,:) !< u₀(x,t)
    real(RNP), contiguous,    intent(in)  :: u   (0:,:,:) !< u(x,t)
    real(RNP), contiguous,    intent(out) :: r_sd(0:,:,:) !< SD-RHS

    character, parameter :: comp = 'S'

    call GetHybridDiffusionTerm(this, cl_operator, comp, theta, bv, u_0, u, r_sd)

  end subroutine GetSDTerm

  !-----------------------------------------------------------------------------
  !> Provides the maximum velocity based on eigenvalues of advective Jacobian

  subroutine GetMaxVelocity(this, cl_operator, u, v_max)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator
    real(RNP), contiguous,    intent(in)  :: u(0:,:,:) !< solution variable
    real(RNP),                intent(out) :: v_max     !< maximum velocity

    real(RNP) :: a, v, T
    integer   :: e, i

    v_max = 0

    do e = 1, cl_operator % ne
      if (cl_operator % activity(e) > 0) then
        do i = 0, cl_operator % eop % po
          call ConservativeToPrimitive(this, u(i,e,:), v, T)
          a = sqrt(this%gamma * this%r_gas * T)
          v_max = max(v_max, abs(v) + a)
        end do
      end if
    end do

  end subroutine GetMaxVelocity

  !-----------------------------------------------------------------------------
  !> Provides the maximum diffusivity

  subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator
    real(RNP), contiguous,    intent(in)  :: u(0:,:,:)  !< solution variable
    real(RNP),                intent(out) :: nu_max     !< maximum velocity

    real(RNP) :: c_eta
    integer   :: e

    nu_max = 0

    if (this%eta <= 0) return

    c_eta = this%eta * max(4*THIRD, 1/this%prandtl)

    do e = 1, cl_operator % ne
      if (cl_operator % activity(e) > 0) then
        nu_max = max(nu_max, c_eta / minval(u(:,e,1)))
      end if
    end do

  end subroutine GetMaxDiffusivity

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

    p_ = (this%gamma - 1) * (u(3) - u(2)*u(2) / (2*u(1)))
    T_ = p_ / (u(1) * this%r_gas)

    if (present(v)) v = u(2) / u(1)
    if (present(p)) p = p_
    if (present(T)) T = T_
    if (present(s)) s = this%c_p * log(T_) - this%r_gas * log(p_)

  end subroutine ConservativeToPrimitive

  !-----------------------------------------------------------------------------
  !> Converts primitive variables into conservative ones

  pure subroutine PrimitiveToConservative(this, v, T, p, u)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: v
    real(RNP),                intent(in)  :: T
    real(RNP),                intent(in)  :: p
    real(RNP),                intent(out) :: u(:)

    u = p / (this%r_gas*T) * [ONE, v, this%c_v * T + HALF*v*v]

  end subroutine PrimitiveToConservative

  !-----------------------------------------------------------------------------
  !> Converts conservative variables into characteristic ones

  pure subroutine ConservativeToCharacteristic(this, u, z)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: u(:)
    real(RNP),                intent(out) :: z(:)

    real(RNP) :: v, p, T, c

    v    = u(2) / u(1)
    p    = (this%gamma - 1) * (u(3) - HALF*v*u(2))
    T    = p / (this%r_gas * u(1))
    c    = 2 / (this%gamma - 1) * sqrt(this%gamma * this%r_gas * T)

    z(1) = c - v
    z(2) = this%c_p * log(T) - this%r_gas * log(p)
    z(3) = c + v

  end subroutine ConservativeToCharacteristic

  !-----------------------------------------------------------------------------
  !> Converts characteristic variable into conservative ones

  pure subroutine CharacteristicToConservative(this, z, u)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    real(RNP),                intent(in)  :: z(:)
    real(RNP),                intent(out) :: u(:)

    real(RNP) :: v, T, p

    associate(r_gas => this%r_gas, gamma => this%gamma)
      v = (z(3) - z(1)) * HALF
      T = (gamma - 1)**2 / (16 * gamma * r_gas) * (z(1) + z(3))**2
      p = T**(gamma /(gamma - 1)) * exp(-z(2)/r_gas)
      u = p / (r_gas * T) * [ONE, v, this%c_v*T + HALF*v*v]
    end associate

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

      A(1,:) = [ZERO                     ,  ONE                   ,  ZERO     ]
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

    real(RNP) :: a, aa, ek, h, p, rho, T, v, z2

    associate(gamma => this%gamma)

      v   = u(2)/u(1)
      ek  = HALF * v**2
      rho = u(1)
      p   = (gamma - 1) * (u(3) - HALF*u(2)*v)
      aa  = gamma * (p/rho)
      a   = sqrt(aa)
      T   = p/(u(1)*this%r_gas)
      h   = this%c_p * T

      if(present(lambda)) then
        lambda = [v-a , v , v+a ]
      end if

      if(present(R)) then
        R(1,1:3) = [ ONE     , ONE , ONE     ]
        R(2,1:3) = [ v-a     , v   , v+a     ]
        R(3,1:3) = [ h-(v*a) , ek  , h+(v*a) ]
      end if

      if(present(L)) then
        z2 = a/(gamma-1)
        L(1,1:3) = [ ek  + v*z2 , -v - z2 ,  ONE ]
        L(2,1:3) = [ 2*h - v*v  , 2*v     , -TWO ]
        L(3,1:3) = [ ek  - v*z2 , -v + z2 ,  ONE ]
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

    ! calculation of the Roe-average
    u_roe = RoeAverage(this, u_l, u_r)

    ! calculation of the absolute Jacobian
    call ConvectiveEigensystem(this, u_roe, Lambda, R, L)
    Lambda = abs(Lambda)
    do i = 1, 3
    do j = 1, 3
      abs_A(i,j) = sum(R(i,:) * Lambda(:) * L(:,j))
    end do
    end do

    ! calculation of the Roe-Flux
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

    real(RNP) ::  r_l, v_l, T_l, H_l
    real(RNP) ::  r_r, v_r, T_r, H_r
    real(RNP) ::  rho_m, v_m, H_m

    r_l = sqrt(u_l(1))
    r_r = sqrt(u_r(1))
    rho_m = r_l * r_r

    call ConservativeToPrimitive(this, u_l, v=v_l, T=T_l)
    call ConservativeToPrimitive(this, u_r, v=v_r, T=T_r)
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
