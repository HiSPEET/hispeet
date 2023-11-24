!> summary:  Base class for compressible Navier-Stokes problems
!> author:   Joerg Stiller, Benedikt Wex
!> date:     2023/11/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - add routines for transformation of variables
!>   - add routines for computing the advection term
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

! comments

!#!  subject of ongoing work
!?!  subject of upcoming work
!x!  obsolete or being replaced

module CL__Problem__CNS__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF
  use Execution_Control
  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D compressible Navier-Stokes problems
  !>
!#!!>       ∂u/∂t + ∂(u²/2)/∂x = ν∂²u/∂x² + f_s(x,t)
!#!!>
!#!!> with constant diffusivity ν. Available boundary conditions are
!#!!>
!#!!>   - Dirichlet with `bc = 'D'` and `bv = u`
!#!!>   - Neumann   with `bc = 'N'` and `bv = ∂u/∂x`
!#!!>   - Periodic  with `bc = 'P'`

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_CNS_1D

    real(RNP) :: r_gas   !< specific gas constant
    real(RNP) :: gamma   !< ratio of specific heats
    real(RNP) :: c_p     !< specific heat for const pressure
    real(RNP) :: c_v     !< specific heat for const volume
    real(RNP) :: eta     !< dynamic viscosity
    real(RNP) :: prandtl !< Prandtl

  contains

    procedure :: GetConvectionTerm
    procedure :: GetDiffusionTerm
    procedure :: GetSDTerm
    procedure :: DiffusionSolver
    procedure :: GetMaxVelocity
    procedure :: GetMaxDiffusivity

    ! CNS specific procedures ..................................................

!#! - small routines: declare pure and keep in this module to allow for inlining
!#! - adapt routines from `flow_1d_equations.f`
!#!     * convert to HiSPEET style and naming
!#!     * use HiSPEET Constants, e.d. I → ONE, O → ONE etc.
!#!     * do not change procedure type without reason (function or subroutine)
!#! - need to rethink the order of characteristic variables
!#!     * option 1: keep as is, order differs from theory
!#!     * option 2: reorder to comply with theory, must be done everywhere,
!#!       e.g., also in eigensystem

    ! transformations
    procedure :: ConservativeToPrimitive      !#! adapt old routines with
    procedure :: PrimitiveToConservative      !#! same name
    procedure :: ConservativeToCharacteristic !#!
    procedure :: CharacteristicToConservative !#!

    ! convection
    procedure :: ConvectiveFlux          !#! ← AdvectiveFluxVector
    procedure :: ConvectiveJacobian      !#! ← AdvectiveFluxJacobian
    procedure :: ConvectiveEigensystem   !#! ← Eigensystem
    procedure :: NumericalConvectiveFlux !#! ← NumericalAdvectiveFlux
                                         !#!   supplement with entropy fix

    !#! - equivalents of `RoeAverage`, `AbsJacobiMatrix` may be also required

    ! ...

  end type CL_Problem_CNS_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

!#! - adapt using CNS specific routines
!#! - keep in mind that in CNS `u` has 3 components instead of one

  subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_c(0:,:,:)

    real(RNP), allocatable :: MD_t(:,:), f_q(:), u_q(:), h_c(:)
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
      allocate(h_c(0:ne))

      allocate(MD_t(0:po, 0:qo), f_q(0:qo), u_q(0:qo))

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
            u_q = matmul(iop%A, u(:,e,1))
          else
            u_q = u(:,e,1)
          end if
          f_q = ConvectiveFlux(u_q)

          ! apply mass-weighted transposed diff-matrix
          r_c(:,e,1) = matmul(MD_t, f_q)

        else
          r_c(:,e,1) = 0
        end if
      end do

      ! element boundary fluxes ................................................

!#! the following call must be put into a loop since "NumericalConvectiveFlux"
!#! cannot be made elemental

      call GetNumericalConvectiveFlux(this, bv, u(:,:,1), h_c)

      do e = 1, ne
        if (activity(e) > 0) then
          r_c( 0,e,1) = r_c( 0,e,1) + h_c(e-1)
          r_c(po,e,1) = r_c(po,e,1) - h_c(e)
        end if
      end do

      ! finalization ...........................................................

      deallocate(h_c)
      !$omp end master

    end associate

  end subroutine GetConvectionTerm

!x! the following routines will be replaced
!x!
!x!  !-----------------------------------------------------------------------------
!x!  !> Convective flux fc(u)
!x!
!x!  elemental function ConvectiveFlux(u) result(f_c)
!x!    real(RNP), intent(in) :: u
!x!    real(RNP) :: f_c
!x!
!x!    f_c = HALF * u ** 2
!x!
!x!  end function ConvectiveFlux
!x!
!x!  !-----------------------------------------------------------------------------
!x!  !> Get numerical convective flux
!x!
!x!  subroutine GetNumericalConvectiveFlux(this, bv, u, h_c)
!x!    class(CL_Problem_CNS_1D), intent(in)  :: this
!x!    real(RNP),                    intent(in)  :: bv(:,:)
!x!    real(RNP), contiguous,        intent(in)  :: u(0:,:)
!x!    real(RNP), contiguous,        intent(out) :: h_c(0:)
!x!
!x!    real(RNP), allocatable :: ul(:), ur(:)
!x!
!x!    integer :: po, ne
!x!
!x!    po = ubound(u,1)
!x!    ne = ubound(u,2)
!x!
!x!    allocate(ul(0:ne), ur(0:ne))
!x!
!x!    ! inner traces
!x!    ul(1:ne)   = u(po, 1:ne)
!x!    ur(0:ne-1) = u( 0, 1:ne)
!x!
!x!    ! left boundary
!x!    select case(this % bc(1))
!x!    case('D')
!x!      ul(0) = 2 * bv(1,1) - u(0,1)
!x!    case('P')
!x!      ul(0) = u(po,ne)
!x!    case default
!x!      ul(0) = u(0,1)
!x!    end select
!x!
!x!    ! right boundary
!x!    select case(this % bc(2))
!x!    case('D')
!x!      ur(ne) = 2 * bv(1,2) - u(po,ne)
!x!    case('P')
!x!      ur(ne) = u(0,1)
!x!    case default
!x!      ur(ne) = u(po,ne)
!x!    end select
!x!
!x!    h_c = RiemannFlux(ul, ur)
!x!
!x!  end subroutine GetNumericalConvectiveFlux
!x!
!x!  !-----------------------------------------------------------------------------
!x!  !> Numerical convective flux hc(ul,ur) based on Riemann solver
!x!
!x!  elemental function RiemannFlux(ul, ur) result(h_c)
!x!    real(RNP), intent(in) :: ul
!x!    real(RNP), intent(in) :: ur
!x!    real(RNP) :: h_c
!x!
!x!    if (ul < ZERO .and. ur > ZERO) then
!x!      h_c = ZERO
!x!    else if (ul + ur > ZERO) then
!x!      h_c = HALF * ul ** 2
!x!    else
!x!      h_c = HALF * ur ** 2
!x!    end if
!x!
!x!  end function RiemannFlux

!?!  !-----------------------------------------------------------------------------
!?!  !> Diffusive contribution to RHS of DG-SEM formulation
!?!
!?!
!?!  subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
!?!    class(CL_Problem_CNS_1D), intent(in)  :: this
!?!    class(CL_Operator_1D),        intent(in)  :: cl_operator
!?!    real(RNP),                    intent(in)  :: bv (:,:)
!?!    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
!?!    real(RNP), contiguous,        intent(out) :: r_d(0:,:,:)
!?!
!?!    real(RNP), allocatable, save :: f(:,:)
!?!    logical,   allocatable, save :: mask(:)
!?!    character :: elliptic_bc(2)
!?!    real(RNP) :: elliptic_bv(2)
!?!
!?!    associate( po => cl_operator % eop % po &
!?!             , ne => cl_operator % ne       )
!?!
!?!      !$omp master
!?!      allocate(mask(ne), source = cl_operator % activity > 0)
!?!      allocate(f(0:po,ne), source = ZERO)
!?!
!?!      elliptic_bc = this % bc(:)(1:1)
!?!      elliptic_bv = bv(1,:)
!?!
!?!      call cl_operator % elliptic_op % Residual( elliptic_bc               &
!?!                                               , elliptic_bv               &
!?!                                               , dx     = cl_operator % dx &
!?!                                               , lambda = ZERO             &
!?!                                               , nu     = this % nu        &
!?!                                               , f      = f                &
!?!                                               , u      = u  (:,:,1)       &
!?!                                               , r      = r_d(:,:,1)       &
!?!                                               , mask   = mask             )
!?!
!?!      deallocate(f, mask)
!?!      !$omp end master
!?!
!?!    end associate
!?!
!?!  end subroutine GetDiffusionTerm
!?!
!?!  !-----------------------------------------------------------------------------
!?!  !> Streamline-diffusion contribution to RHS of DG-SEM formulation
!?!
!?!  subroutine GetSDTerm(this, cl_operator, tau, bv, u_0, u, r_sd)
!?!    class(CL_Problem_CNS_1D), intent(in)  :: this
!?!    class(CL_Operator_1D),        intent(in)  :: cl_operator
!?!    real(RNP),                    intent(in)  :: tau
!?!    real(RNP),                    intent(in)  :: bv (:,:)
!?!    real(RNP), contiguous,        intent(in)  :: u_0 (0:,:,:)
!?!    real(RNP), contiguous,        intent(in)  :: u   (0:,:,:)
!?!    real(RNP), contiguous,        intent(out) :: r_sd(0:,:,:)
!?!
!?!    real(RNP), allocatable, save :: f(:,:), nu_sd(:,:)
!?!    logical,   allocatable, save :: mask(:)
!?!    character :: elliptic_bc(2)
!?!    real(RNP) :: elliptic_bv(2)
!?!
!?!    associate( po => cl_operator % eop % po &
!?!             , ne => cl_operator % ne       )
!?!
!?!      !$omp master
!?!      allocate(mask(ne), source = cl_operator % activity > 0)
!?!      allocate(f(0:po,ne), source = ZERO)
!?!      allocate(nu_sd(0:po,ne))
!?!      call GetStreamlineDiffusivity(this, tau, u_0(:,:,1), nu_sd)
!?!
!?!      elliptic_bc = this % bc(:)(1:1)
!?!      elliptic_bv = bv(1,:)
!?!
!?!      call cl_operator % elliptic_op % Residual( elliptic_bc               &
!?!                                               , elliptic_bv               &
!?!                                               , dx     = cl_operator % dx &
!?!                                               , lambda = ZERO             &
!?!                                               , nu     = nu_sd            &
!?!                                               , f      = f                &
!?!                                               , u      = u   (:,:,1)      &
!?!                                               , r      = r_sd(:,:,1)      &
!?!                                               , mask   = mask             )
!?!
!?!      deallocate(f, mask, nu_sd)
!?!      !$omp end master
!?!
!?!    end associate
!?!
!?!  end subroutine GetSDTerm
!?!
!?!
!?!  !---------------------------------------------------------------------------
!?!  !> Implicit diffusion solver for Burgers problems
!?!  !>
!?!  !> Implicit method for solving or relaxing the diffusion subproblem
!?!  !>
!?!  !>       u = f + ∆t [r_d(bv,u) + r_ds(bv,τ,u₀,u)]
!?!  !>
!?!  !> The streamline-diffusion term `r_ds` is evaluated with `u₀` and included
!?!  !> only if τ > 0.
!?!  !> At present, the following solution methods are available:
!?!  !>
!?!  !> 1. Direct hybrid solver
!?!  !>      - requires `cl_operator%eop%hybrid = .true.`
!?!  !>      - not suited for cases with frozen elements
!?!  !>      - not suited for streamline-diffusion unless f_c'(u) is constant
!?!  !>
!?!  !> 2. Conjugate gradient method
!?!  !>      - basic solver
!?!  !>      - bad smoother
!?!  !>
!?!  !> 3. Schwarz method
!?!  !>      - requires proper initialization of `cl_operator%eop%schwarz`
!?!  !>      - bad solver
!?!  !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`
!?!  !>
!?!  !> 4. Inexact preconditioned conjugate gradient method (IPCG)
!?!  !>      - requires proper initialization of `cl_operator%eop%schwarz`
!?!  !>      - best iterative solver when used with no overlap, i.e.,
!?!  !>        `schwarz % delta = -1`, `schwarz % no_min = -1`
!?!  !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`
!?!
!?!  subroutine DiffusionSolver( this, cl_operator, dt, tau, bv, f, u_0, u &
!?!                            , method, i_max, r_red, r_max               )
!?!    class(CL_Problem_CNS_1D), intent(in) :: this
!?!    class(CL_Operator_1D), intent(in)    :: cl_operator
!?!    real(RNP),             intent(in)    :: dt          !< ∆t = t - t₀
!?!    real(RNP),             intent(in)    :: tau         !< SD time scale τ
!?!    real(RNP),             intent(in)    :: bv (:,:)    !< boundary values
!?!    real(RNP), contiguous, intent(in)    :: f  (0:,:,:) !< sources
!?!    real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< frozen solution
!?!    real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< approx solution
!?!    integer,               intent(in)    :: method      !< solution method
!?!    integer,               intent(in)    :: i_max       !< max num iterations
!?!    real(RNP), optional,   intent(in)    :: r_red       !< residual reduction
!?!    real(RNP), optional,   intent(in)    :: r_max       !< max residual
!?!
!?!    real(RNP), allocatable, save :: g(:,:), nu_tot(:,:)
!?!    logical,   allocatable, save :: mask(:)
!?!    character :: elliptic_bc(2)
!?!    real(RNP) :: elliptic_bv(2)
!?!    real(RNP) :: lambda
!?!    integer   :: e
!?!
!?!    associate( nu          => this % nu                 &
!?!             , po          => cl_operator % eop % po    &
!?!             , ne          => cl_operator % ne          &
!?!             , dx          => cl_operator % dx          &
!?!             , Me          => cl_operator % Me          &
!?!             , activity    => cl_operator % activity    &
!?!             , elliptic_op => cl_operator % elliptic_op )
!?!
!?!      !$omp master
!?!
!?!      ! elliptic boundary conditions and boundary values
!?!      elliptic_bc = this % bc(:)(1:1)
!?!      elliptic_bv = bv(1,:)
!?!
!?!      ! Helmholtz parameter
!?!      lambda = 1 / dt
!?!
!?!      ! activity mask
!?!      allocate(mask(ne), source = activity > 0)
!?!
!?!      ! sources
!?!      allocate(g(0:po,ne))
!?!      do e = 1, ne
!?!        if (mask(e)) then
!?!          g(:,e) = lambda * Me * f(:,e,1)
!?!        else
!?!          g(:,e) = 0
!?!        end if
!?!      end do
!?!
!?!      if (tau > ZERO) then
!?!
!?!        allocate(nu_tot(0:po,ne))
!?!        call GetStreamlineDiffusivity(this, tau, u_0(:,:,1), nu_tot)
!?!        nu_tot = nu_tot + this % nu
!?!
!?!        select case(method)
!?!        case(1)
!?!          call Error( 'DiffusionSolver'             &
!?!                    , 'method 1 not suited for ISD' &
!?!                    , 'CL__Problem__CNS__1D'    )
!?!        case(2)
!?!          call elliptic_op % CG_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        case(3)
!?!          call elliptic_op % Schwarz_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        case(4)
!?!          call elliptic_op % SchwarzPCG_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        end select
!?!
!?!      else
!?!
!?!        select case(method)
!?!        case(1)
!?!          call elliptic_op % HybridSolver &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) )
!?!        case(2)
!?!          call elliptic_op % CG_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        case(3)
!?!          call elliptic_op % Schwarz_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        case(4)
!?!          call elliptic_op % SchwarzPCG_Method &
!?!                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
!?!                  , i_max, r_red, r_max, mask )
!?!        end select
!?!
!?!      end if
!?!
!?!      deallocate(g, mask)
!?!      if (allocated(nu_tot)) deallocate(nu_tot)
!?!      !$omp end master
!?!
!?!    end associate
!?!
!?!  end subroutine DiffusionSolver
!?!
!?!  !---------------------------------------------------------------------------
!?!  !> Computation of streamline diffusivity
!?!
!?!  subroutine GetStreamlineDiffusivity(cl_problem, tau, u, nu_sd)
!?!    class(CL_Problem_CNS_1D), intent(in) :: cl_problem
!?!    real(RNP),             intent(in)  :: tau         !< SD time scale τ
!?!    real(RNP), contiguous, intent(in)  :: u    (0:,:) !< approx solution u
!?!    real(RNP), contiguous, intent(out) :: nu_sd(0:,:) !< approx solution u
!?!
!?!    integer :: e
!?!
!?!    nu_sd = tau/2 * u**2
!?!
!?!    select case(cl_problem % nu_sd_filter)
!?!    case(1)
!?!      ! element maximum
!?!      do e = 1, size(nu_sd,2)
!?!        nu_sd(:,e) = maxval(nu_sd(:,e))
!?!      end do
!?!    end select
!?!
!?!  end subroutine GetStreamlineDiffusivity
!?!
!?!  !-----------------------------------------------------------------------------
!?!  !> Provides the maximum velocity based on eigenvalues of advective Jacobian
!?!
!?!  subroutine GetMaxVelocity(this, cl_operator, u, v_max)
!?!    class(CL_Problem_CNS_1D), intent(in)  :: this
!?!    class(CL_Operator_1D),        intent(in)  :: cl_operator
!?!    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
!?!    real(RNP),                    intent(out) :: v_max     !< maximum velocity
!?!
!?!    integer :: e
!?!
!?!    v_max = 0
!?!
!?!    do e = 1, cl_operator % ne
!?!      if (cl_operator % activity(e) > 0) then
!?!        v_max = max(v_max, maxval(abs(u(:,e,1))))
!?!      end if
!?!    end do
!?!
!?!    ! avoid compiler warning
!?!    if (this % nc > 0) return
!?!
!?!  end subroutine GetMaxVelocity
!?!
!?!  !-----------------------------------------------------------------------------
!?!  !> Provides the maximum diffusivity
!?!
!?!  subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
!?!    class(CL_Problem_CNS_1D), intent(in)  :: this
!?!    class(CL_Operator_1D),        intent(in)  :: cl_operator
!?!    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
!?!    real(RNP),                    intent(out) :: nu_max     !< maximum velocity
!?!
!?!    nu_max = this % nu
!?!
!?!    ! avoid compiler warnings
!?!    if (this % nc > 0 .or. cl_operator % ne > 1 .or. size(u) > 0) return
!?!
!?!  end subroutine GetMaxDiffusivity

  !=============================================================================
  ! CNS specific routines

  !-----------------------------------------------------------------------------
  !>

  pure subroutine ConservativeToPrimitive(this, u, rho, v, T, p, s)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP),            intent(in)  :: u(3)
    real(RNP),  optional, intent(out) :: rho
    real(RNP),  optional, intent(out) :: v
    real(RNP),  optional, intent(out) :: T
    real(RNP),  optional, intent(out) :: p
    real(RNP),  optional, intent(out) :: s

  !#! ...
  !#! use `this % kappa` etc to access fluid properties

  end subroutine ConservativeToPrimitive

  !-----------------------------------------------------------------------------
  !>

  pure subroutine PrimitiveToConservative(this, v, T, p, u)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in)  :: v
    real(RNP), intent(in)  :: T
    real(RNP), intent(in)  :: p
    real(RNP), intent(out) :: u(3)

  !#! ...

 end subroutine PrimitiveToConservative

  !-----------------------------------------------------------------------------
  !>

  pure subroutine ConservativeToCharacteristic(this, u, z)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in)  :: u(3)
    real(RNP), intent(out) :: z(3)

    !#! ...

   end subroutine PrimitiveToConservative

  !-----------------------------------------------------------------------------
  !>

  pure subroutine CharacteristicToConservative(this, z, u)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in)  :: z(3)
    real(RNP), intent(out) :: u(3)

    !#! ...

  end subroutine CharacteristicToConservative

  !-----------------------------------------------------------------------------
  !>

  pure function ConvectiveFlux(this, u) result(f_c)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(3)
    real(RNP)             :: f_c(3)

    !#! ...

  end function AdvectiveFluxVector

  !-----------------------------------------------------------------------------
  !>

  pure function ConvectiveJacobian(this, u) result(A)
    real(RNP), intent(in) :: u(3)
    real(RNP)             :: A(3,3)

    !#! ...

  end function ConvectiveJacobian

  !-----------------------------------------------------------------------------
  !>

  pure subroutine ConvectiveEigensystem(this, u, Lambda, R, L)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in)            :: u(3)
    real(RNP), intent(out), optional :: Lambda(3)
    real(RNP), intent(out), optional :: R(3,3)
    real(RNP), intent(out), optional :: L(3,3)

    !#! ...

  end subroutine ConvectiveEigensystem

  !-----------------------------------------------------------------------------
  !>

  pure function NumericalConvectiveFlux(this, ul, ur) result(h_c)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: ul(3)
    real(RNP), intent(in) :: ur(3)
    real(RNP)             :: h_c(3)

    !#! ...

  end function NumericalConvectiveFlux

  !=============================================================================

end module CL__Problem__CNS__1D

