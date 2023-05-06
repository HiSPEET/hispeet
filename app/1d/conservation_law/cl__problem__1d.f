!> summary:  Interface to 1D conservation problems discretized with DG-SEM
!> author:   Joerg Stiller
!> date:     2019/11/08, revised 2023/04/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__1D

  use Kind_Parameters, only: RNP
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

    ! private control parameters ...............................................

    logical, private :: has_exact_solution = .false.
    logical, private :: has_sources        = .false.

  contains

    procedure :: GetExactSolution
    procedure :: GetSources

    procedure(SetProblem       ), deferred :: SetProblem
    procedure(GetInitialValues ), deferred :: GetInitialValues
    procedure(GetBoundaryValues), deferred :: GetBoundaryValues
    procedure(GetConvectionTerm), deferred :: GetConvectionTerm
    procedure(GetDiffusionTerm ), deferred :: GetDiffusionTerm
    procedure(GetSDTerm        ), deferred :: GetSDTerm
    procedure(DiffusionSolver  ), deferred :: DiffusionSolver

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

    subroutine GetSDTerm(this, cl_operator, tau, bv, u, r_sd)
      import
      class(CL_Problem_1D),  intent(in)  :: this
      class(CL_Operator_1D), intent(in)  :: cl_operator
      real(RNP),             intent(in)  :: tau          !< SD time scale τ
      real(RNP),             intent(in)  :: bv  (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)  :: u   (0:,:,:) !< u(x,t)
      real(RNP), contiguous, intent(out) :: r_sd(0:,:,:) !< SD-RHS
    end subroutine GetSDTerm

    !---------------------------------------------------------------------------
    !> Implicit diffusion solver
    !>
    !> Implicit method for solving or relaxing the diffusion subproblem
    !>
    !>       u = u₀ + ∆t [r_d(bv,u) + r_ds(bv,τ,u)]
    !>
    !> The streamline-diffusion term `r_ds` is included only if τ > 0.
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
    !>      - good smoother when used with overlap `schwarz%no ≈ po/4`
    !>
    !> 4. Inexact preconditioned conjugate gradient method (IPCG)
    !>      - requires proper initialization of `cl_operator%eop%schwarz`
    !>      - best iterative solver when used with no overlap, `schwarz%no = 0`
    !>      - good smoother when used with overlap `schwarz%no ≈ po/4`

    subroutine DiffusionSolver( this, cl_operator, dt, tau, bv, u_0, u &
                              , method, i_max, r_red, r_max            )
      import
      class(CL_Problem_1D),  intent(in)    :: this
      class(CL_Operator_1D), intent(in)    :: cl_operator
      real(RNP),             intent(in)    :: dt          !< ∆t = t - t₀
      real(RNP),             intent(in)    :: tau         !< SD time scale τ
      real(RNP),             intent(in)    :: bv (:,:)    !< boundary values
      real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< initial value u₀
      real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< approx solution u
      integer,               intent(in)    :: method      !< solution method
      integer,               intent(in)    :: i_max       !< max num iterations
      real(RNP), optional,   intent(in)    :: r_red       !< residual reduction
      real(RNP), optional,   intent(in)    :: r_max       !< max residual
    end subroutine DiffusionSolver

  end interface

contains

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

  !=============================================================================

end module CL__Problem__1D
