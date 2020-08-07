!> summary:  IMEX Euler with velocity or pressure correction scheme
!> author:   Joerg Stiller
!> date:     2018/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX Euler with velocity or pressure correction scheme
!>
!> @note
!>
!>
!>
!>
!> @endnote
!>
!> @todo
!>   *  clean-up
!>   *  eliminate repeated computation
!> @endtodo
!===============================================================================

module CART__ISP_Flow__Euler

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use TPO_sDDD
  use Array_Assignments
  use ISP_Flow_Problem
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__DG_Weak_Divergence
  use CART__DG_Weak_Gradient
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Boundary_Values
  use CART__ISP_Flow__Convection
  use CART__ISP_Flow__Diffusion
  use CART__ISP_Flow__Pressure
  use CART__ISP_Flow__Projection
  use CART__ISP_Flow__Time_Derivative
  implicit none
  private

  public :: EulerVC

contains

!-------------------------------------------------------------------------------
!> Velocity correction -- trying iteration
!>
!> At input, u is expected to provide at least a first-order approximation of
!> u(t+dt), e.g. u = u₀ = u(t).

subroutine EulerVC(problem, flow_op, t, dt, u_0, u, F, H, S)

  ! arguments ..................................................................

  class(FlowProblem),   intent(in)    :: problem !< flow problem
  class(FlowOperators), intent(inout) :: flow_op !< flow operators
  real(RNP),            intent(inout) :: t       !< time
  real(RNP),            intent(in)    :: dt      !< time step size
  real(RNP),            intent(in)    :: u_0     !< solution u(t)
  real(RNP),            intent(inout) :: u       !< solution u(t+dt)
  real(RNP),  optional, intent(out)   :: F       !< ∂u/∂t(t+dt)
  real(RNP),  optional, intent(inout) :: H       !< δu(t) → δu(t+dt)
  real(RNP),  optional, intent(in)    :: S       !< ∫∂u/dt over (t,t+dt)

  dimension :: u_0 (:,:,:,:,:)
  dimension :: u   (:,:,:,:,:)
  dimension :: F   (:,:,:,:,:)
  dimension :: H   (:,:,:,:,:)
  dimension :: S   (:,:,:,:,:)

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:,:,:), allocatable, save :: F_d1, F_d2, F_d3
  real(RNP), dimension(:,:,:,:,:), allocatable, save :: u_i, nu, w
  real(RNP), dimension(:,:,:,:),   allocatable, save :: dp

  real(RNP) :: t_0
  integer   :: np, ne, nc

  associate( mesh => flow_op % mesh  &
           , x    => flow_op % x     &
           , eop  => flow_op % eop_u &
           , p    => u(:,:,:,:,4)    )

    ! initialization ...........................................................

    t_0 = t
    t   = t + dt

    ! workspace
    !$omp single
    allocate(F_d1, mold = u)
    allocate(F_d2, mold = u)
    allocate(F_d3, mold = u)
    allocate(u_i , mold = u)
    allocate(w   , mold = u)
    if (problem % HasVariableProperties()) then
      allocate(nu, mold = u)
    end if
    allocate(dp, mold = p)
    !$omp end single

    np = size(u,1)
    ne = size(u,4)
    nc = size(u,5)

    ! boundary conditions ......................................................

    call GetBoundaryValues(problem, mesh, flow_op%bv_x, t, flow_op%bv_u)

    ! predictor ................................................................

    ! uᵢ' = uᵢ₋₁ ≡ u_0
    call SetArray(u_i, u_0, multi=.true.)

    ! SDC: uᵢ' += Sᵢᵏ⁻¹ - Hᵢᵏ⁻¹
    if (present(H) .and. present(S)) then
      call MergeArrays(ONE, u_i,  ONE, S, multi=.true.)
      call MergeArrays(ONE, u_i, -ONE, H, multi=.true.)
    end if

    ! nu = ν₀ = ν(tᵢ₋₁)
    if (problem % HasVariableProperties()) then
      call problem % GetDiffusivity(flow_op%x, t_0, u_0, nu)
    end if

    ! w = -∇⋅vᵢ₋₁vᵢ₋₁ + ∇⋅ν₀(∇vᵢ₋₁)ᵀ - χ∇(ν₀∇⋅vᵢ₋₁)     for v
    ! w = -∇⋅vᵢ₋₁uᵢ₋₁                                   for u \ (v,p)
     call TimeDerivative( problem, flow_op, t  &
                        , u_c  = u_0           &
                        , u_d  = u_0           &
                        , nu   = nu            &
                        , F    = w             &
                        , F_d1 = F_d1          &
                        , F_d2 = F_d2          &
                        , F_d3 = F_d3          &
                        )
     call MergeArrays(ONE, w, -ONE, F_d1, multi =.true.)

     ! w += ∇⋅ν₀∇uᵢᵏ⁻¹ ≡ F_d1, using current approximation of u
     call TimeDerivative( problem, flow_op, t  &
                        , u_d  = u             &
                        , nu   = nu            &
                        , F_d1 = F_d1          &
                        , F_d3 = F_d3          & !!! cut for published version
                        )
     call MergeArrays(ONE, w, ONE, F_d1, multi =.true.)


    ! uᵢ' += ∆t w
    call MergeArrays(ONE, u_i, dt, w, multi=.true.)

   ! pressure, continuity and diffusion .......................................

    ! solve for p = pᵢ"
    call PressureSolver(problem, flow_op, dt, u_i, p, w)

    ! vᵢ" = vᵢ' - 1/∆t ∇pᵢ" - J(vᵢ")
    call ProjectionStep(problem, flow_op, dt, p, u_i, w)

    ! solve for uᵢ'''
    call MergeArrays(ONE, u_i, -dt  , F_d1, multi=.true.)
    call MergeArrays(ONE, u_i, -dt/2, F_d3, multi=.true.) !
    call DiffusionStep(problem, flow_op, dt, f=u_i, u=u, w=w, nu=nu)

    !!! activate following command for published version
    !call TimeDerivative(problem, flow_op, t, u_d=u, nu=nu, F_d1 = F_d1)

    ! final projection .........................................................

    if (flow_op % control % div_final) then

      ! solve for pᵢᵏ
      call SetArray(dp, ZERO)
      call PressureSolver(problem, flow_op, dt, u, dp, w)
      call MergeArrays(ONE, p, ONE, dp)

      ! vᵢᵏ = vᵢ''' - 1/∆t ∇pᵢᵏ - J(vᵢᵏ)
      call ProjectionStep(problem, flow_op, dt, dp, u, w)

    end if

    ! SDC: low-order operator H ................................................

    if (present(H)) then

      ! convective part
      call TimeDerivative(problem, flow_op, t, u_c=u_0, F = H)

      ! diffusive part
      !!! deactivate following command for published version
      call TimeDerivative(problem, flow_op, t, u_d=u, nu=nu, F_d1=F_d1)
      call MergeArrays(ONE, H, ONE, F_d1, multi=.true.)
      call MergeArrays(ONE, H, ONE, F_d2, multi=.true.)
      call ScaleArray(H, dt, multi = .true.)

    end if

    ! SDC: time derivative .....................................................

    if (present(F)) then

      if (problem % HasVariableProperties()) then
        call problem % GetDiffusivity(flow_op%x, t, u, nu)
      end if

      call TimeDerivative(problem, flow_op, t, u, u, nu = nu, F = F)

    end if

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    if (allocated(F_d1)) deallocate(F_d1)
    if (allocated(F_d2)) deallocate(F_d2)
    if (allocated(F_d3)) deallocate(F_d3)
    if (allocated(u_i )) deallocate(u_i )
    if (allocated(nu  )) deallocate(nu  )
    if (allocated(w   )) deallocate(w   )
    if (allocated(dp  )) deallocate(dp  )
    !$omp end master

  end associate

end subroutine EulerVC

!===============================================================================

end module CART__ISP_Flow__Euler
