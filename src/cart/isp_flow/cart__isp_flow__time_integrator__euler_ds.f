!> summary:  Euler method for incompressible flows with dual splitting
!> author:   Joerg Stiller
!> date:     2020/03/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Time_Integrator__Euler_DS
  use Kind_Parameters, only: RNP

  use CART__ISP_Flow__Time_Integrator
  use CART__ISP_Flow__Boundary_Values

  implicit none
  private

  public :: TimeIntegrator_EulerDS

  type, extends(TimeIntegrator) :: TimeIntegrator_EulerDS
  contains
    procedure :: TimeStep
  end type TimeIntegrator_EulerDS

contains

  !-----------------------------------------------------------------------------
  !>

  subroutine TimeStep(this, t, dt, u, nu)
    class(TimeIntegrator_EulerDS), intent(inout) :: this
    real(RNP),           intent(inout) :: t             !< time t₀ → t
    real(RNP),           intent(in)    :: dt            !< step size ∆t = t-t₀
    real(RNP), optional, intent(in)    :: nu(:,:,:,:,:) !< variable ν(x,t₀)
    real(RNP),           intent(inout) :: u (:,:,:,:,:) !< u(x,t₀) → u(x,t)

    ! local variables  .........................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: F_d1, F_d3
    real(RNP), dimension(:,:,:,:,:), allocatable, save :: u_i, w
    real(RNP), dimension(:,:,:,:),   allocatable, save :: dp

    associate( problem => this % problem          &
             , flow_op => this % flow_op          &
             , mesh    => this % flow_op % mesh   &
             , x       => this % flow_op % x      &
             , eop     => this % flow_op % eop_u  &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      t   = t + dt

      ! workspace
      !$omp single
      allocate(u_i , mold = u)
      allocate(F_d1, mold = u)
      allocate(F_d3, mold = u)
      allocate(w   , mold = u)
      allocate(dp  , mold = p)
      !$omp end single

      ! boundary conditions ....................................................

      call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

      ! explicit + extrapolated diffusive parts ................................

      ! u' = u₀ ≡ u(x,t₀)
      call SetArray(u_i, u, multi=.true.)

      ! w = -∇⋅v₀v₀ + ∇⋅ν₀(∇v₀)ᵀ - χ∇(ν₀∇⋅v₀)     for v
      ! w = -∇⋅v₀u₀                               for u \ (v,p)
       call TimeDerivative( problem, flow_op, t  &
                          , u_c  = u             &
                          , u_d  = u             &
                          , nu   = nu            &
                          , F    = w             &
                          , F_d1 = F_d1          &
                          , F_d3 = F_d3          &
                          )

      ! u' += ∆t w
      call MergeArrays(ONE, u_i, dt, w, multi=.true.)

      ! pressure, continuity and diffusion .....................................

      ! solve for p = p"
      call PressureSolver(problem, flow_op, dt, u_i, p, w)

      ! v" = v' - 1/∆t ∇p" - J(v")
      call ProjectionStep(problem, flow_op, dt, p, u_i, w)

      ! solve implicit diffusive part for u'''
      call MergeArrays(ONE, u_i, -dt, F_d1, multi=.true.)
      call MergeArrays(ONE, u_i, -dt, F_d3, multi=.true.)
      call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

      ! final projection .......................................................

      if (flow_op % control % div_final) then

        ! solve for p = p" + dp
        call SetArray(dp, ZERO)
        call PressureSolver(problem, flow_op, dt, u, dp, w)
        call MergeArrays(ONE, p, ONE, dp)

        ! v = v''' - 1/∆t ∇p - J(v)
        call ProjectionStep(problem, flow_op, dt, dp, u, w)

      end if

      ! clean-up ...............................................................

      !$omp barrier
      !$omp master
      if (allocated(F_d1)) deallocate(F_d1)
      if (allocated(F_d3)) deallocate(F_d3)
      if (allocated(u_i )) deallocate(u_i )
      if (allocated(w   )) deallocate(w   )
      if (allocated(dp  )) deallocate(dp  )
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CART__ISP_Flow__Time_Integrator__Euler_DS
