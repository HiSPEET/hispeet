!> summary:  Generic projection-diffusion step for incompressible flow
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - validate
!>   - revise interface
!>   - add standby mode
!>   - possibly find a better name
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_ProjectionStep
  use Array_Assignments
  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> extrapolation-projection-diffusion step for incompressible flow

  module subroutine ProjectionStep( this, tau, v_0, F_c, F_d, Q, bv_u, u &
                                  , i_max_p, i_max_v, r_red, r_max       )


    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    !< v₀, effective initial value of velocity
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    !< convective term at final time t
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    !< diffusion term at final time t
    real(RNP), contiguous, intent(in) :: Q(:,:,:,:,:)
    !< sources at time t, older contributions and possibly correction terms
    class(BoundaryVariable_3D), intent(inout) :: bv_u(:)
    !< boundary values at final time t
    !!   - for velocity, to be given in components 1-3 on input, and
    !!   - for pressure, returned in component 4 on output
    real(RNP), contiguous, intent(out) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure at final time u

    integer, intent(in) :: i_max_p !< max num iterations of pressure solver
    integer, intent(in) :: i_max_v !< max num iterations of diffusion solver
    real(RNP), optional, intent(in) :: r_red !< min L² residual reduction to reach
    real(RNP), optional, intent(in) :: r_max !< max L² residual allowed

    ! internal variables .......................................................

    real(RNP), allocatable, save :: pp (:,:,:,:)   ! outer pressure traces p⁺
    real(RNP), allocatable, save :: vm (:,:,:,:,:) ! inner velocity traces v⁻
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_vm

    type(BoundaryVariable_3D), allocatable, save :: bv_v(:)
    type(BoundaryVariable_3D), allocatable, save :: bv_p(:)

    integer :: np
    integer :: b, d, e

    associate( ins_op  => this % ins_op                   &
             , mesh    => this % ins_op % mesh            &
             , n_elem  => this % ins_op % mesh % n_elem   &
             , n_ghost => this % ins_op % mesh % n_ghost  &
             , v       => u(:,:,:,:,1:3)                  &
             , p       => u(:,:,:,:,4)                    )

      ! initialization .........................................................

      np = size(v,1)

      !$omp master

      allocate( pp (np, np,  6, n_elem             ), source = ZERO )
      allocate( vm (np, np,  6, n_elem + n_ghost, 3), source = ZERO )
      allocate( vp (np, np,  6, n_elem          , 3), source = ZERO )
      allocate( w  (np, np, np, n_elem          , 4), source = ZERO )

      buf_vm = ElementFaceTransferBuffer_3D(mesh, vm)

      ! handles for velocity and pressure boundary values, based on pointers
      allocate(bv_v( mesh%n_bound ))
      allocate(bv_p( mesh%n_bound ))
      do b = 1, mesh % n_bound
        bv_v(b) = BoundaryVariable_3D(bv_u(b), first=1, last=3)
        bv_p(b) = BoundaryVariable_3D(bv_u(b), first=4, last=4)
      end do
      !$omp end master
      !$omp barrier

      ! extrapolation step .....................................................

      ! Computes the intermediate velocity v' = v₀ + τ(Fc + Fd + Q) and extracts
      ! the inner traces v'⁻ needed in the next step.

      !$omp do collapse(2)
      do d = 1, 3
      do e = 1, mesh % n_elem

        ! intermediate velocity v = v'
        v(:,:,:,e,d) = v_0(:,:,:,e,d) + tau * ( F_c(:,:,:,e,d) &
                                              + F_d(:,:,:,e,d) &
                                              + Q  (:,:,:,e,d) )

        ! inner traces: vm = v'⁻
        vm(:,:,1,e,d) = v( 1,:,:,e,d)
        vm(:,:,2,e,d) = v(np,:,:,e,d)
        vm(:,:,3,e,d) = v(:, 1,:,e,d)
        vm(:,:,4,e,d) = v(:,np,:,e,d)
        vm(:,:,5,e,d) = v(:,:, 1,e,d)
        vm(:,:,6,e,d) = v(:,:,np,e,d)

      end do
      end do

      ! pressure computation ...................................................

      associate(div_v => w(:,:,:,:,4))

        ! generate outer traces of intermediate velocity
        call buf_vm % Transfer(mesh, vm, tag=100)       ! transfer v⁻ from masters
        call ins_op % SetVelocityBC(bv_v, vm)           ! set boundary values
        call buf_vm % Merge(vm)                         ! merge received traces
        call ConvertInnerToOuterTraces_3D(mesh, vm, vp) ! vm → vp = v⁺

        ! divergence of intermediate velocity
         call TPO_Div(ins_op % eop_v, ins_op % sem_v, v, vp, div_v)

        ! solve pressure equation
        call ins_op % PressureSolver( tau, bv_v, v, div_v, bv_p, p &
                                    , i_max_p, r_red, r_max        )

      end associate

      ! pressure correction ....................................................

      associate(grad_p => w(:,:,:,:,1:3))

        ! generate outer traces of pressure -- sufficient for Neumann BC
        call GetOuterTraces_3D(mesh, p, pp)

        ! pressure gradient
        call TPO_Grad(ins_op % eop_v, ins_op % sem_v, p, pp, grad_p)

        ! correction: v = v - τ∇p
        call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)

      end associate

      ! diffusive correction ...................................................

      associate(f => w(:,:,:,:,1:3))

        !$omp do collapse(2)
        do e = 1, mesh % n_elem
          do d = 1, 3
            f(:,:,:,e,d) = 1/tau * v(:,:,:,e,d) - F_d(:,:,:,e,d)
          end do
        end do

         call ins_op % DiffusionSolver(tau, f, bv_v, v, i_max_v, r_red, r_max)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w)
      deallocate(pp, vm, vp)
      deallocate(buf_vm)
      deallocate(bv_v, bv_p)
      !$omp end master

    end associate

  end subroutine ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
