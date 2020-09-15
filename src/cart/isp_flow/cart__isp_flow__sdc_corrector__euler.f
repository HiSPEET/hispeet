!> summary:  SDC corrector based on IMEX Euler with dual splitting
!> author:   Joerg Stiller
!> date:     2020/05/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__SDC_Corrector__Euler
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Array_Assignments
  use XMPI

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
  public :: SDC_Corrector_Euler_Options

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC corrector for incompressible flow

  type, extends(SDC_Corrector) :: SDC_Corrector_Euler
  contains
    procedure :: Init_SDC_Corrector_Euler
    procedure :: Show => Show_SDC_Corrector_Euler
    procedure :: GetCorrectorRHS
    procedure :: CorrectionStep
  end type SDC_Corrector_Euler

  ! overloading the constructor
  interface SDC_Corrector_Euler
    module procedure New_SDC_Corrector_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler SDC-corrector options (none, so far)

  type, extends(SDC_Corrector_Options) :: SDC_Corrector_Euler_Options
  end type SDC_Corrector_Euler_Options

contains

  !=============================================================================
  ! SDC_Corrector_Euler: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Euler

  function New_SDC_Corrector_Euler(problem, flow_op, opt) result(this)
    class(FlowProblem),   target,       intent(in) :: problem !< flow problem
    class(FlowOperators), target,       intent(in) :: flow_op !< flow operators
    class(SDC_Corrector_Euler_Options), intent(in) :: opt     !< SDC options
    type(SDC_Corrector_Euler) :: this

    call Init_SDC_Corrector_Euler(this, problem, flow_op, opt)

  end function New_SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Euler object

  subroutine Init_SDC_Corrector_Euler(this, problem, flow_op, opt)
    class(SDC_Corrector_Euler),   intent(inout) :: this
    class(FlowProblem),   target, intent(in)    :: problem !< flow problem
    class(FlowOperators), target, intent(in)    :: flow_op !< flow operators
    class(SDC_Corrector_Euler_Options), intent(in) :: opt  !< SDC options

    ! intialize parent type
    call this % Init_SDC_Corrector(problem, flow_op, opt)
    this % name = 'IMEX Euler corrector'

  end subroutine Init_SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Euler settings

  subroutine Show_SDC_Corrector_Euler(this, unit)
    class(SDC_Corrector_Euler), intent(in) :: this
    integer,          optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_SDC_Corrector(unit)

    write(io,*)

  end subroutine Show_SDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals

  subroutine GetCorrectorRHS(this, t, nu_c, nu_s, u, F_ex, F_im, F)

    class(SDC_Corrector_Euler), intent(in) :: this
    real(RNP), intent(in)  :: t               !< time
    real(RNP), intent(in)  :: nu_c(:,:,:,:,:) !< variable ν for corrector RHS
    real(RNP), intent(in)  :: nu_s(:,:,:,:,:) !< variable ν for subintegral RHS
    real(RNP), intent(in)  :: u   (:,:,:,:,:) !< u
    real(RNP), intent(out) :: F_ex(:,:,:,:,:) !< explicit RHS for corrector
    real(RNP), intent(out) :: F_im(:,:,:,:,:) !< implicit RHS for corrector
    real(RNP), intent(out) :: F   (:,:,:,:,:) !< RHS for subintegrals

    optional :: nu_c, nu_s

    real(RNP), allocatable, save :: F_d2(:,:,:,:,:)
    real(RNP), allocatable, save :: F_d3(:,:,:,:,:)
    integer :: i

    !$omp single
    allocate(F_d2, mold = u)
    allocate(F_d3, mold = u)
    !$omp end single

    ! RHS for high-order subintegrals ..........................................

    call TimeDerivative( this % problem     &
                       , this % flow_op     &
                       , t                  &
                       , u    = u           &
                       , chi  = this % chi  &
                       , nu   = nu_s        &
                       , F    = F           &
                       , F_c  = F_ex        )

    ! RHS for low-order part of corrector ......................................

    call TimeDerivative( this % problem     &
                       , this % flow_op     &
                       , t                  &
                       , u    = u           &
                       , chi  = this % chi  &
                       , nu   = nu_c        &
                       , F_d1 = F_im        &
                       , F_d2 = F_d2        &
                       , F_d3 = F_d3        )

    do i = 1, 3
      call MergeArrays(ONE, F_ex(:,:,:,:,i), ONE, F_d2(:,:,:,:,i))
      if (this % chi /= 0 .and. this%cd3 /= ONE ) then
        call MergeArrays(ONE, F_ex(:,:,:,:,i), ONE-this%cd3, F_d3(:,:,:,:,i))
      end if
    end do

    !$omp single
    deallocate(F_d2, F_d3)
    !$omp end single

  end subroutine GetCorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectionStep( this, t, dt                   &
                           , F_ex_0_old, F_ex_0, F_ex_old  &
                           , F_im_0_old, F_im_0, F_im_old  &
                           , S, u_0, u, nu                 )

    class(SDC_Corrector_Euler), intent(inout) :: this
    real(RNP), intent(inout) :: t                      !< time t₀ → t
    real(RNP), intent(in)    :: dt                     !< step size ∆t = t-t₀
    real(RNP), intent(in)    :: F_ex_0_old (:,:,:,:,:) !< F^ex (t₀)ᵏ
    real(RNP), intent(in)    :: F_ex_0     (:,:,:,:,:) !< F^ex (t₀)ᵏ⁺¹
    real(RNP), intent(in)    :: F_ex_old   (:,:,:,:,:) !< F^ex (t )ᵏ    ! unused
    real(RNP), intent(in)    :: F_im_0_old (:,:,:,:,:) !< F^im (t₀)ᵏ    ! unused
    real(RNP), intent(in)    :: F_im_0     (:,:,:,:,:) !< F^im (t₀)ᵏ⁺¹  ! unused
    real(RNP), intent(in)    :: F_im_old   (:,:,:,:,:) !< F^im (t )ᵏ
    real(RNP), intent(in)    :: S          (:,:,:,:,:) !< S    (t₀)ᵏ
    real(RNP), intent(in)    :: u_0        (:,:,:,:,:) !< u    (t₀)ᵏ⁺¹
    real(RNP), intent(inout) :: u          (:,:,:,:,:) !< u    (t )ᵏ → (t)ᵏ⁺¹
    real(RNP), intent(inout) :: nu         (:,:,:,:,:) !< v    (t₀)ᵏ⁺¹  ! unused

    optional :: nu

    ! local variables  .........................................................

    real(RNP), allocatable, save :: u_i  (:,:,:,:,:) ! intermediate solution
    real(RNP), allocatable, save :: F_im (:,:,:,:,:) ! aproximate F_im
    real(RNP), allocatable, save :: F_d3 (:,:,:,:,:) ! aproximate F_d3
    real(RNP), allocatable, save :: w    (:,:,:,:,:) ! workspace for u
    real(RNP), allocatable, save :: dp   (:,:,:,:)   ! pressure correction

    real(RNP), save :: t_0
    integer :: i

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , cd3     => this % cd3              &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      !$omp single

      allocate(u_i , mold = u)
      allocate(F_im, mold = u)
      allocate(F_d3, mold = u)
      allocate(w   , mold = u)
      allocate(dp  , mold = p)

      t_0 = t
      t   = t + dt

      !$omp end single

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      ! compute F_im = F_d1(v₀ᵏ⁺¹,uᵏ) and F_d3(v₀ᵏ⁺¹,uᵏ)
      call TimeDerivative( this % problem     &
                         , this % flow_op     &
                         , t_0                &
                         , u    = u           &
                         , chi  = this % chi  &
                         , nu   = nu          &
                         , F_d1 = F_im        &
                         , F_d3 = F_d3        )

      ! u' = u₀ + Sᵏ + ∆t [ F^ex(t₀)ᵏ⁺¹    - F^ex(t₀)ᵏ
      !                   + F^im(v₀ᵏ⁺¹,uᵏ) - F^im(v₀ᵏ,uᵏ)
      !                   + F_d3(v₀ᵏ⁺¹,uᵏ) ]
      call SetArray(u_i, u_0, multi=.true.)
      do i = 1, problem % nc
        if (i == 4) cycle
        call MergeArrays(ONE, u_i(:,:,:,:,i), ONE, S         (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_ex_0    (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_ex_0_old(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i),  dt, F_im      (:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_im_old  (:,:,:,:,i))
        if (i < 4 .and. this % chi /= 0) then
          call MergeArrays(ONE, u_i(:,:,:,:,i), dt*this%cd3, F_d3(:,:,:,:,i))
        end if
      end do

      ! projection .............................................................

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! diffusion ...............................................................

      ! starting values and RHS for diffusion
      do i = 1, problem % nc
        if (i == 4) cycle
        call SetArray(u(:,:,:,:,i), u_i(:,:,:,:,i))
        call MergeArrays(ONE, u_i(:,:,:,:,i), -dt, F_im(:,:,:,:,i))
        if (i < 4 .and. this % cd3 /= 0) then
          call MergeArrays(ONE, u_i(:,:,:,:,i), -dt*cd3, F_d3(:,:,:,:,i))
        end if
      end do

      ! solve u'''/∆t - ν(t₀)ᵏ⁺¹ u''' = u_i/∆t
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

      ! final projection .......................................................

      if (this % project > 0) then

        ! p = p" + ∆p  with ∆p such that  ∇²∆p = ∆t ∇⋅(v''' - v) = 0
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp single
      deallocate(u_i )
      deallocate(F_im)
      deallocate(F_d3)
      deallocate(w   )
      deallocate(dp  )
      !$omp end single

    end associate

    ! touch possibly unused arguments to suppress compiler warnings ;)
    if (size(F_ex_0_old) > 0) return
    if (size(F_ex_0    ) > 0) return
    if (size(F_ex_old  ) > 0) return
    if (size(F_im_0_old) > 0) return
    if (size(F_im_0    ) > 0) return

  end subroutine CorrectionStep

  !=============================================================================

end module CART__ISP_Flow__SDC_Corrector__Euler
