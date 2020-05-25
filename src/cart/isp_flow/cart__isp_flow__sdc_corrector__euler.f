!> summary:  SDC corrector based on IMEX Euler with dual splitting
!> author:   Joerg Stiller
!> date:     2020/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Corrector__Euler
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Array_Assignments

  use ISP_Flow_Problem

  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__Time_Derivative
  use CART__ISP_Flow__Time_Integrator

  use CART__ISP_Flow__SDC_Corrector

  implicit none
  private

  public :: SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Abstract type of a SDC corrector for incompressible flow

  type, extends(SDC_Corrector) :: SDC_Corrector_Euler
  contains
    procedure :: Init_SDC_Corrector_Euler
    procedure :: GetCorrectorRHS
    procedure :: CorrectionStep
  end type SDC_Corrector_Euler

  ! overloading the constructor
  interface SDC_Corrector_Euler
    module procedure New_SDC_Corrector_Euler
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Euler

  function New_SDC_Corrector_Euler(problem, flow_op) result(this)
    class(FlowProblem),   target,  intent(in)    :: problem !< flow problem
    class(FlowOperators), target,  intent(in)    :: flow_op !< flow operators
    type(SDC_Corrector_Euler) :: this

    call Init_SDC_Corrector_Euler(this, problem, flow_op)

  end function New_SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Euler object

  subroutine Init_SDC_Corrector_Euler(this, problem, flow_op)
    class(SDC_Corrector_Euler),   intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op !< flow operators

    call this % Init_SDC_Corrector(problem, flow_op)

  end subroutine Init_SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals
  !>
  !> For IMEX Euler with dual splitting:
  !>
  !>       F_ex = -∇⋅(v u) + ∇⋅[ν(∇v_ex)ᵀ]
  !>       F_im =  ∇⋅(ν ∇u)
  !>       F    =  F_ex + F_im - χ∇ν(∇·v) + f(x,t)
  !>
  !> Note
  !>  *  last terms in F_ex and F for velocity only

  subroutine GetCorrectorRHS(this, t, nu, u, F_ex, F_im, F)

    class(SDC_Corrector_Euler), intent(in) :: this
    real(RNP), intent(in)  :: t                !< time
    real(RNP), intent(in)  :: nu   (:,:,:,:,:) !< diffusivity
    real(RNP), intent(in)  :: u    (:,:,:,:,:) !< u
    real(RNP), intent(out) :: F_ex (:,:,:,:,:) !< explicit RHS for corrector
    real(RNP), intent(out) :: F_im (:,:,:,:,:) !< implicit RHS for corrector
    real(RNP), intent(out) :: F    (:,:,:,:,:) !< RHS for subintegrals

    optional :: nu

    real(RNP), allocatable, save :: F_d2(:,:,:,:,:)

    !$omp single
    allocate(F_d2, mold = u)
    !$omp end single

    call TimeDerivative( this % problem  &
                       , this % flow_op  &
                       , t               &
                       , u_c  = u        &
                       , u_d  = u        &
                       , nu   = nu       &
                       , F    = F        &
                       , F_c  = F_ex     & ! -∇⋅(v u)
                       , F_d1 = F_im     & !  ∇⋅(ν ∇u)
                       , F_d2 = F_d2     ) !  ∇⋅[ν(∇v_ex)ᵀ]

    call MergeArrays(ONE, F_ex, ONE, F_d2, multi=.true.)

    !$omp barrier
    !$omp master
    deallocate(F_d2)
    !$omp end master

  end subroutine GetCorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectionStep( this, t, dt               &
                           , F_ex_0_old, F_ex_0, F_ex  &
                           , F_im_0_old, F_im_0, F_im  &
                           , S, u_0, u, nu             )

    class(SDC_Corrector_Euler), intent(inout) :: this
    real(RNP), intent(inout) :: t                      !< time t₀ → t
    real(RNP), intent(in)    :: dt                     !< step size ∆t = t-t₀
    real(RNP), intent(in)    :: F_ex_0_old (:,:,:,:,:) !< F^ex (t₀)ᵏ⁻¹
    real(RNP), intent(in)    :: F_ex_0     (:,:,:,:,:) !< F^ex (t₀)ᵏ
    real(RNP), intent(inout) :: F_ex       (:,:,:,:,:) !< F^ex (t )ᵏ⁻¹
    real(RNP), intent(in)    :: F_im_0_old (:,:,:,:,:) !< F^im (t₀)ᵏ⁻¹
    real(RNP), intent(in)    :: F_im_0     (:,:,:,:,:) !< F^im (t₀)ᵏ
    real(RNP), intent(inout) :: F_im       (:,:,:,:,:) !< F^im (t )ᵏ⁻¹
    real(RNP), intent(in)    :: S          (:,:,:,:,:) !< S    (t₀)ᵏ⁻¹
    real(RNP), intent(in)    :: u_0        (:,:,:,:,:) !< u    (t₀)ᵏ
    real(RNP), intent(inout) :: u          (:,:,:,:,:) !< u    (t )ᵏ⁻¹ → (t)ᵏ
    real(RNP), intent(inout) :: nu         (:,:,:,:,:) !< v    (t )ᵏ⁻¹

    optional :: nu

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: w   (:,:,:,:,:) ! workspace for u
    real(RNP), allocatable, save :: dp  (:,:,:,:)   ! pressure correction

    real(RNP) :: t_0

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      t_0 = t
      t   = t + dt

      ! workspace
      !$omp single
      allocate(u_i, mold = u)
      allocate(w  , mold = u)
      allocate(dp , mold = p)
      !$omp end single

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      ! u' = u₀ + ∆t [F^ex(t₀)ᵏ - F^ex(t₀)ᵏ⁻¹] + Sᵏ⁻¹
      call SetArray(u_i, u_0, multi=.true.)
      call MergeArrays(ONE, u_i,  dt, F_ex_0    , multi=.true.)
      call MergeArrays(ONE, u_i, -dt, F_ex_0_old, multi=.true.)
      call MergeArrays(ONE, u_i, ONE, S         , multi=.true.)

      ! projection .............................................................

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! diffusion ...............................................................

      ! remove old diffusion termi i.e: u_i -= F_im ≡ ∇⋅[ν(t₀)∇u(t)]ᵏ⁻¹
      call MergeArrays(ONE, u_i, -dt, F_im, multi=.true.)

      ! solve u'''/∆t - ν(t₀)ᵏ u''' = u_i/∆t
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

      ! final projection .......................................................

      if (flow_op % control % div_final) then

        ! p = p" + ∆p  with ∆p such that  ∇²∆p = ∆t ∇⋅(v''' - v) = 0
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp master
      deallocate(u_i)
      deallocate(w  )
      deallocate(dp )
      !$omp end master

    end associate

    ! touch possibly unused arguments to suppress compiler warnings ;)
    if (size(F_im_0_old) > 0) return
    if (size(F_im_0    ) > 0) return
    if (size(F_ex_0_old) > 0) return
    if (size(F_ex_0    ) > 0) return
    if (size(F_ex      ) > 0) return

  end subroutine CorrectionStep

  !=============================================================================

end module CART__ISP_Flow__SDC_Corrector__Euler
