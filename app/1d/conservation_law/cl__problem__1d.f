!> summary:  Interface to 1D conservation problems discretized with DG-SEM
!> author:   Joerg Stiller
!> date:     2019/11/08, revised 2023/04/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__1D

  use Kind_Parameters, only: RNP, RDP
  use Constants,       only: ONE, HALF, THIRD
  use Eigenproblems

  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_1D
  public :: CL_Problem_Options_1D

  !-----------------------------------------------------------------------------
  !> Abstract type for defining and handling 1D conservation problems
  !>
  !> ## Boundary conditions
  !>
  !> The type of conditions at the left and right boundaries are encoded in
  !> component `bc`. The actual definition of boundary conditions is left to
  !> the derived problem class. Periodic 'P', static 'S' and do-nothing ' '
  !> conditions should be supported in all problem classes.
  !> Scalar conservation laws will typically include also Dirichlet 'D' and
  !> Neumann 'N' conditions.
  !> More complex BC can occur in systems, e.g. different types of inflow and
  !> outflow conditions. For example, a subsonic inlet with specified velocity,
  !> and temperature could be encoded as 'I:vT'.

  type, abstract :: CL_Problem_1D

    ! problem parameters .......................................................

    integer           :: nc    !< number of conservation variables
    real(RNP)         :: xb1   !< position of left boundary
    real(RNP)         :: xb2   !< position of right boundary
    character(len=80) :: bc(2) !< BC types at left and right boundaries

    ! control parameters .......................................................

    integer   :: ad_filter_method  !< artificial diffusivity filtering method:
                                   !! - `0` no filtering
                                   !! - `1` reduction to specified degree
                                   !! - `2` reduction to half degree
                                   !! - `3` diagonal maximum

    integer   :: ad_filter_basis   !< artificial diffusivity filtering basis:
                                   !! - `1` Legendre
                                   !! - `2` linear+bubble
                                   !! - `3` interpolation

    integer   :: ad_filter_degree  !< artificial diffusivity max degree:
                                   !! - with Legendre basis      ≥ 0
                                   !! - with linear+bubble basis ≥ 1

    integer   :: dc_diffusivity    !< discontinuity capturing (DC) diffusivity:
                                   !! - `0` none
                                   !! - `1` Persson & Preraire (2006)

    integer   :: dc_sensor_var     !< DC sensor variable
    real(RNP) :: dc_sensor_coeff   !< DC sensor coefficient
    real(RNP) :: dc_sensor_const   !< DC sensor constant
    real(RNP) :: dc_sensor_delta   !< DC sensor half width
    real(RNP) :: dc_scaling_coeff  !< DC max diffusivity scaling factor

    integer   :: limiting_method   !< limiting method:
                                   !! - `0` none
                                   !! - `1` momentum (Burbeau et al. 2001)

    integer   :: limiting_scope    !< limiting scope:
                                   !! - `1` step
                                   !! - `2` stages or substeps

    logical   :: limiting_initial  !< limit initial conditions

    integer   :: regularization    !< regularization method
                                   !!   - `0` none
                                   !!   - `1` strict
                                   !!   - `2` weak

    ! automatic parameters .....................................................

    logical :: has_exact_solution = .false.

  contains

    procedure :: Init_CL_Problem_1D

    procedure :: HasExactSolution
    procedure :: GetExactSolution
    procedure :: GetSources
    procedure :: GetTimeScales

    procedure :: MomentLimiter
    procedure :: RegularityFilter

    procedure(SetProblem            ), deferred :: SetProblem
    procedure(HasDiffusion          ), deferred :: HasDiffusion
    procedure(ConvectiveFlux        ), deferred :: ConvectiveFlux
    procedure(ConvectiveJacobian    ), deferred :: ConvectiveJacobian
    procedure(ConvectiveEigensystem ), deferred :: ConvectiveEigensystem
    procedure(GetInitialValues      ), deferred :: GetInitialValues
    procedure(GetBoundaryValues     ), deferred :: GetBoundaryValues
    procedure(GetConvectionTerm     ), deferred :: GetConvectionTerm
    procedure(GetHybridDiffusionTerm), deferred :: GetHybridDiffusionTerm
    procedure(GetDiffusionTerm      ), deferred :: GetDiffusionTerm
    procedure(GetSDTerm             ), deferred :: GetSDTerm
    procedure(DiffusionSolver       ), deferred :: DiffusionSolver
    procedure(GetMaxVelocity        ), deferred :: GetMaxVelocity
    procedure(GetMaxDiffusivity     ), deferred :: GetMaxDiffusivity

  end type CL_Problem_1D

  !-----------------------------------------------------------------------------
  !> 1D conservation problem options

  type CL_Problem_Options_1D

    real(RNP)         :: xb1   =  0    !< position of left boundary
    real(RNP)         :: xb2   =  1    !< position of right boundary
    character(len=80) :: bc(2) = 'P'   !< BC types at left and right boundaries

    integer   :: ad_filter_method =  0 !< SD diffusivity filtering method
    integer   :: ad_filter_basis  =  1 !< SD diffusivity filtering basis
    integer   :: ad_filter_degree =  0 !< SD diffusivity max degree

    integer   :: dc_diffusivity   =  0 !< DC diffusivity method
    integer   :: dc_sensor_var    =  1 !< DC sensor variable
    real(RNP) :: dc_sensor_coeff  = -4 !< DC sensor coefficient
    real(RNP) :: dc_sensor_const  = -4 !< DC sensor constant
    real(RNP) :: dc_sensor_delta  =  1 !< DC sensor half width
    real(RNP) :: dc_scaling_coeff =  1 !< DC max diffusivity scaling factor

    integer   :: limiting_method  =  0      !< limiting method
    integer   :: limiting_scope   =  0      !< limiting scope
    logical   :: limiting_initial = .false. !< limit initial conditions

    integer   :: regularization   = 0       !< regularization method

  end type CL_Problem_Options_1D

  !=============================================================================
  ! Procedures of type CL_Problem_1D defined in submodules

  interface

    !---------------------------------------------------------------------------
    !> Application of moment limiter

    module subroutine MomentLimiter(this, cl_operator, u)
      class(CL_Problem_1D),  intent(in)    :: this
      class(CL_Operator_1D), intent(in)    :: cl_operator
      real(RNP), contiguous, intent(inout) :: u(0:,:,:)
    end subroutine MomentLimiter

  end interface

  !=============================================================================
  ! Procedures that are defined in extensions of type CL_Problem_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Initialization of the flow problem

    subroutine SetProblem(this, file)
      import :: CL_Problem_1D
      class(CL_Problem_1D),       intent(inout) :: this
      character(len=*), optional, intent(in)    :: file !< (*.prm)
    end subroutine SetProblem

    !---------------------------------------------------------------------------
    !> Query whether problem has non vanishing physical diffusion

    logical function HasDiffusion(this)
      import :: CL_Problem_1D
      class(CL_Problem_1D), intent(in) :: this
    end function HasDiffusion

    !---------------------------------------------------------------------------
    !> Calculates the convective flux with given conservative variables

    pure function ConvectiveFlux(this, u) result(f_c)
      import :: CL_Problem_1D, RNP
      class(CL_Problem_1D), intent(in) :: this
      real(RNP), intent(in) :: u(:)
      real(RNP) :: f_c(this%nc)
    end function ConvectiveFlux

    !---------------------------------------------------------------------------
    !> Calculates the Jacobian ot the Convective Term

    pure function ConvectiveJacobian(this, u) result(A_c)
      import :: CL_Problem_1D, RNP
      class(CL_Problem_1D), intent(in) :: this
      real(RNP), intent(in) :: u(:)
      real(RNP) :: A_c(this%nc,this%nc)
    end function ConvectiveJacobian

    !---------------------------------------------------------------------------
    !> Calculates the eigendecomposition of the Convective Jacobian

    pure subroutine ConvectiveEigensystem(this, u, lambda, R, L)
      import :: CL_Problem_1D, RNP
      class(CL_Problem_1D), intent(in)  :: this
      real(RNP),                intent(in)  :: u(:)
      real(RNP), optional,      intent(out) :: lambda(:)
      real(RNP), optional,      intent(out) :: R(:,:)
      real(RNP), optional,      intent(out) :: L(:,:)
    end subroutine ConvectiveEigensystem

    !---------------------------------------------------------------------------
    !> Provides the initial values u(x,0)

    subroutine GetInitialValues(this, cl_operator, u)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:ne,1:nc)
    end subroutine GetInitialValues

    !---------------------------------------------------------------------------
    !> Provides the left and right boundary values for time t
    !>
    !> The boundary values are dimensioned as `bv(nc,2)`, where `nc = this%nc`
    !> is the number of solution components. The second index indicates the
    !> boundary, `1` for the left and `2`for the right.
    !> In general, the boundary values do not correspond to the conservation
    !> variables u. Their meaning depends on the actual problem and on the
    !> boundary type, which is given in `this%bc`.

    subroutine GetBoundaryValues(this, t, bv)
      import :: CL_Problem_1D, RNP
      class(CL_Problem_1D), intent(in)  :: this
      real(RNP),            intent(in)  :: t       !< time
      real(RNP),            intent(out) :: bv(:,:) !< boundary values
    end subroutine GetBoundaryValues

    !---------------------------------------------------------------------------
    !> Convective contribution to RHS of DG-SEM formulation

    subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP),             intent(in)  :: bv (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_c(0:,:,:) !< convective RHS
    end subroutine GetConvectionTerm

    !---------------------------------------------------------------------------
    !> Unified physical, artificial and streamline diffusion term

    subroutine GetHybridDiffusionTerm &
        (this, cl_operator, comp, theta, bv, u_0, u, r_d)

      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator !< spatial operators
      character(len=*),      intent(in)  :: comp        !< composition flag
      real(RNP),             intent(in)  :: theta       !< SD time scale
      real(RNP), optional,   intent(in)  :: bv(:,:)     !< boundary values
      real(RNP), contiguous, intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
      real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_d(0:,:,:) !< diffusion RHS

    end subroutine GetHybridDiffusionTerm

    !---------------------------------------------------------------------------
    !> Diffusive contribution to RHS of DG-SEM formulation

    subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP),             intent(in)  :: bv (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_d(0:,:,:) !< diffusion RHS
    end subroutine GetDiffusionTerm

    !---------------------------------------------------------------------------
    !> Streamline-diffusion contribution to RHS of DG-SEM formulation
    !>
    !> Evaluates the weak form of the streamline-diffusion operators for `u`
    !> using `u₀` for computing the streamline diffusivity.

    subroutine GetSDTerm(this, cl_operator, theta, bv, u_0, u, r_sd)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP),             intent(in)  :: theta        !< SD time scale θ
      real(RNP),             intent(in)  :: bv  (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)  :: u_0 (0:,:,:) !< u₀(x,t)
      real(RNP), contiguous, intent(in)  :: u   (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_sd(0:,:,:) !< SD-RHS
    end subroutine GetSDTerm

    !---------------------------------------------------------------------------
    !> Implicit diffusion solver
    !>
    !> Implicit method for solving or relaxing the diffusion subproblem
    !>
    !>       u = f + ∆t [r_d(bv,u) + r_ds(bv,θ,u₀,u)]
    !>
    !> The streamline-diffusion term `r_ds` is evaluated with `u₀` and included
    !> only if θ > 0.
    !> At present, the following solution methods are available:
    !>
    !> 1. Direct hybrid solver
    !>      - requires `cl_operator%eop%hybrid = .true.`
    !>      - not suited for cases with frozen elements
    !>      - not suited for streamline-diffusion unless f_c'(u) is constant
    !>
    !> 2. Conjugate gradient method
    !>      - basic solver
    !>      - bad smoother
    !>
    !> 3. Schwarz method
    !>      - requires proper initialization of `cl_operator%eop%schwarz`
    !>      - bad solver
    !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`
    !>
    !> 4. Inexact preconditioned conjugate gradient method (IPCG)
    !>      - requires proper initialization of `cl_operator%eop%schwarz`
    !>      - best iterative solver when used with no overlap, i.e.,
    !>        `schwarz % delta = -1`, `schwarz % no_min = -1`
    !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`

    subroutine DiffusionSolver( this, cl_operator, dt, theta, bv, f, u_0, u &
                              , method, i_max, r_red, r_max                 )

      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)    :: this
      class(CL_Operator_1D), intent(in)    :: cl_operator
      real(RNP),             intent(in)    :: dt          !< ∆t = t - t₀
      real(RNP),             intent(in)    :: theta       !< SD time scale θ
      real(RNP),             intent(in)    :: bv (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)    :: f  (0:,:,:) !< sources
      real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< frozen solution
      real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< approx solution
      integer,               intent(in)    :: method      !< solution method
      integer,               intent(in)    :: i_max       !< max num iterations
      real(RNP), optional,   intent(in)    :: r_red       !< residual reduction
      real(RNP), optional,   intent(in)    :: r_max       !< max residual

    end subroutine DiffusionSolver

    !---------------------------------------------------------------------------
    !> Provides the maximum velocity based on eigenvalues of advective Jacobian

    subroutine GetMaxVelocity(this, cl_operator, u, v_max)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(in)  :: u(0:,:,:) !< solution variable
      real(RNP),             intent(out) :: v_max     !< maximum velocity
    end subroutine GetMaxVelocity

    !---------------------------------------------------------------------------
    !> Provides the maximum diffusivity

    subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
      import :: CL_Problem_1D, CL_Operator_1D, RNP
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(in)  :: u(0:,:,:) !< solution variable
      real(RNP),             intent(out) :: nu_max    !< maximum diffusivity
    end subroutine GetMaxDiffusivity

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Initialization from options

  subroutine Init_CL_Problem_1D(this, opt)
    class(CL_Problem_1D),         intent(inout) :: this
    class(CL_Problem_Options_1D), intent(in)    :: opt

    this % xb1 = opt % xb1
    this % xb2 = opt % xb2
    this % bc  = opt % bc

    this % ad_filter_method = opt % ad_filter_method
    this % ad_filter_basis  = opt % ad_filter_basis
    this % ad_filter_degree = opt % ad_filter_degree

    this % dc_diffusivity   = opt % dc_diffusivity
    this % dc_sensor_var    = opt % dc_sensor_var
    this % dc_sensor_coeff  = opt % dc_sensor_coeff
    this % dc_sensor_const  = opt % dc_sensor_const
    this % dc_sensor_delta  = opt % dc_sensor_delta
    this % dc_scaling_coeff = opt % dc_scaling_coeff

    this % limiting_method  = opt % limiting_method
    this % limiting_scope   = opt % limiting_scope
    this % limiting_initial = opt % limiting_initial

    this % regularization   = opt % regularization

  end subroutine Init_CL_Problem_1D

  !-----------------------------------------------------------------------------
  !> Query whether an exact solution is available

  logical function HasExactSolution(this)
    class(CL_Problem_1D), intent(in) :: this

    HasExactSolution = this % has_exact_solution

  end function HasExactSolution

  !-----------------------------------------------------------------------------
  !> Dummy procedure for the exact solution u(x,t)

  subroutine GetExactSolution(this, cl_operator, t, u)
    class(CL_Problem_1D),  intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP),             intent(in)  :: t
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u(x,t) at mesh points

    u = 0

    ! just to avoid compiler warnings ;)
    if (this%nc * cl_operator%ne * t > 0) return

  end subroutine GetExactSolution

  !-----------------------------------------------------------------------------
  !> Dummy procedure for the sources f_s(x,t;u)
  !>
  !> The solution `u` is included to cover situations where independent source
  !> contributions are combined with solution components, as in the case of
  !> compressible flows.

  subroutine GetSources(this, cl_operator, t, u, f_s)
    class(CL_Problem_1D),  intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator !< provides x
    real(RNP),             intent(in)  :: t           !< time
    real(RNP), contiguous, intent(in)  :: u(0:,:,:)   !< u(x,t)
    real(RNP), contiguous, intent(out) :: f_s(0:,:,:) !< sources

    f_s = 0

    ! just to avoid compiler warnings ;)
    if (this%nc * cl_operator%ne * t * size(u) > 0) return

  end subroutine GetSources

  !-----------------------------------------------------------------------------
  !> Evaluation of convective and diffusive time scales
  !>
  !> The time scales are derived from the eigenvalues of the convection and
  !> diffusion operators for a single element with Dirichlet conditions on
  !> one or both sides, respectively. These eigenvalues are computed for the
  !> standard element with normalized coefficients and then scaled to the
  !> actual element size and maximum velocity or diffusivity. The convective
  !> and diffusive time scales are obtained from the minimum of the inverse
  !> magnitude of the respective eigenvalues.
  !>
  !> In case of no advection or diffusion, a negative value is returned for the
  !> corresponding time scale.

  subroutine GetTimeScales(this, cl_operator, u, tau_conv, tau_diff)
    class(CL_Problem_1D),  intent(in)  :: this        !< CL problem
    class(CL_Operator_1D), intent(in)  :: cl_operator !< CL operator
    real(RNP), contiguous, intent(in)  :: u(:,:,:)    !< conservative variables
    real(RNP),             intent(out) :: tau_conv    !< convective time scale
    real(RNP),             intent(out) :: tau_diff    !< diffusive  time scale

    real(RNP), parameter :: eps = epsilon(ONE)

    complex(RDP), allocatable :: lmb_conv(:)
    real(RDP),    allocatable :: lmb_diff(:), A_conv(:,:), A_diff(:,:)

    real(RNP) :: lmb_conv_max, lmb_diff_max, nu_max, v_max

    associate( eop => cl_operator % eop      &
             , po  => cl_operator % eop % po &
             , dx  => cl_operator % dx       )

      ! eigenvalues of model problems ..........................................

      select case(po)

      case(0:1)

        lmb_conv_max = HALF
        lmb_diff_max = ONE

      case default

        ! convection problem with unit velocity in standard element
        allocate(A_conv(po,po), lmb_conv(po))
        A_conv = eop % D(1:,1:)
        call SolveNonsymmetricEigenproblem(A_conv, lmb_conv)
        lmb_conv_max = real(maxval(abs(lmb_conv)), RNP)

        ! diffusion problem with Dirichlet conditions in standard element
        allocate(A_diff(po-1,po-1), lmb_diff(po-1))
        A_diff = eop % L(1:po-1,1:po-1)
        call SolveSymmetricEigenproblem(A_diff, lmb_diff)
        lmb_diff_max = real(maxval(abs(lmb_diff)), RNP)

      end select

      ! time scales ............................................................

      call this % GetMaxVelocity(cl_operator, u, v_max)
      if (v_max > 0) then
        tau_conv = dx / max(2 * lmb_conv_max * v_max, eps)
      else
        tau_conv = -1
      end if

      call this % GetMaxDiffusivity (cl_operator, u, nu_max)
      if (nu_max > 0) then
        tau_diff = dx**2 / max(4 * lmb_diff_max * nu_max, eps)
      else
        tau_diff = -1
      end if

    end associate

  end subroutine GetTimeScales

  !-----------------------------------------------------------------------------
  !> Regularity filter

  subroutine RegularityFilter(this, cl_operator, u)
    class(CL_Problem_1D),  intent(in)    :: this
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP), contiguous, intent(inout) :: u(0:,:,:)

    if (this%regularization == 0 .or. cl_operator%ne > 0 .or. size(u) > 0) then
      return
    end if

  end subroutine RegularityFilter

  !=============================================================================

end module CL__Problem__1D
