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

    ! solution parameters ......................................................

    integer :: nu_sd_filter = -1 !< streamline-diffusivity filtering mode

    ! control parameters .......................................................

    logical :: has_exact_solution = .false.

  contains

    procedure :: HasExactSolution
    procedure :: GetExactSolution
    procedure :: GetSources
    procedure :: GetTimeScales

    procedure(SetProblem       ), deferred :: SetProblem
    procedure(GetInitialValues ), deferred :: GetInitialValues
    procedure(GetBoundaryValues), deferred :: GetBoundaryValues
    procedure(GetConvectionTerm), deferred :: GetConvectionTerm
    procedure(GetDiffusionTerm ), deferred :: GetDiffusionTerm
    procedure(GetSDTerm        ), deferred :: GetSDTerm
    procedure(DiffusionSolver  ), deferred :: DiffusionSolver
    procedure(GetMaxVelocity   ), deferred :: GetMaxVelocity
    procedure(GetMaxDiffusivity), deferred :: GetMaxDiffusivity

  end type CL_Problem_1D

  !=============================================================================
  ! Procedures that are defined in extensions of type CL_Problem_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Initialization of the flow problem

    subroutine SetProblem(this, file)
      import
      class(CL_Problem_1D),       intent(inout) :: this
      character(len=*), optional, intent(in)    :: file !< (*.prm)
    end subroutine SetProblem

    !---------------------------------------------------------------------------
    !> Provides the initial values u(x,0)

    subroutine GetInitialValues(this, cl_operator, u)
      import
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:n1,1:nc)
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
      import
      class(CL_Problem_1D), intent(in)  :: this
      real(RNP),            intent(in)  :: t       !< time
      real(RNP),            intent(out) :: bv(:,:) !< boundary values
    end subroutine GetBoundaryValues

    !---------------------------------------------------------------------------
    !> Convective contribution to RHS of DG-SEM formulation

    subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
      import
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP),             intent(in)  :: bv (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_c(0:,:,:) !< convective RHS
    end subroutine GetConvectionTerm

    !---------------------------------------------------------------------------
    !> Diffusive contribution to RHS of DG-SEM formulation

    subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
      import
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
      import
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
      import
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
      import
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(in)  :: u(0:,:,:) !< solution variable
      real(RNP),             intent(out) :: v_max     !< maximum velocity
    end subroutine GetMaxVelocity

    !---------------------------------------------------------------------------
    !> Provides the maximum diffusivity

    subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
      import
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP), contiguous, intent(in)  :: u(0:,:,:) !< solution variable
      real(RNP),             intent(out) :: nu_max    !< maximum diffusivity
    end subroutine GetMaxDiffusivity

  end interface

contains

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

  !=============================================================================

end module CL__Problem__1D
