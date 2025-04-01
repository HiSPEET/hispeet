!> summary:  Generic projection-diffusion step for incompressible flow
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - add standby mode
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_ProjectionStep
  use Array_Assignments
  use TPO__AAA__3D
  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> extrapolation-projection-diffusion step for incompressible flow

  module subroutine ProjectionStep( this, tau, t, v_0, F_c, F_d, Q &
                                  , bv_u, mu, nu, u, precon        )

    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), intent(in) :: t
    !< t, final time
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    !< v₀, effective initial value of velocity
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    !< convective term at final time t
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    !< diffusion term at final time t
    real(RNP), contiguous, intent(in) :: Q(:,:,:,:,:)
    !< sources at time t and further known terms
    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  [ v₁, v₂, v₃, - , -  ]
    !!   - Γᴼ :  [ - , - , - , p , ∆p ]
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure at final time u
    logical, optional, intent(in) :: precon
    !< if present, operate as preconditioner (regardless of the value)

    ! internal variables .......................................................

    real(RNP), allocatable, save :: pp (:,:,:,:)   ! outer pressure traces p⁺
    real(RNP), allocatable, save :: vm (:,:,:,:,:) ! inner velocity traces v⁻
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_vm

    type(BoundaryVariable_3D), allocatable, save :: bv_w(:), bv_p(:), bv_dp(:)

    real(RNP) :: r_red, r_max
    integer   :: i_max_p, i_max_v
    integer   :: b, d, e, np
    logical   :: update_viscosity

    associate( problem => this % problem                  &
             , ins_op  => this % ins_op                   &
             , sem_u   => this % ins_op % sem_u           &
             , mesh    => this % ins_op % mesh            &
             , n_elem  => this % ins_op % mesh % n_elem   &
             , n_ghost => this % ins_op % mesh % n_ghost  &
             , v       => u(:,:,:,:,1:3)                  &
             , p       => u(:,:,:,:,4)                    )

      ! initialization .........................................................

      if (present(precon)) then
        i_max_p = this % i_pre_p
        i_max_v = this % i_pre_v
        r_red   = this % r_pre_red
        update_viscosity = .false.
      else
        i_max_p = this % i_max_p
        i_max_v = this % i_max_v
        r_red   = this % r_red
        update_viscosity = present(nu) .and. problem % HasVariableProperties()
      end if
      r_max = this % r_max

      np = size(v,1)

      !$omp master

      allocate( pp (np, np,  6, n_elem             ), source = ZERO )
      allocate( vm (np, np,  6, n_elem + n_ghost, 3), source = ZERO )
      allocate( vp (np, np,  6, n_elem          , 3), source = ZERO )
      allocate( w  (np, np, np, n_elem          , 4), source = ZERO )

      buf_vm = ElementFaceTransferBuffer_3D(mesh, vm)

      ! provide handles for velocity and pressure boundary values
      allocate(bv_w ( mesh%n_bound ))
      allocate(bv_p ( mesh%n_bound ))
      allocate(bv_dp( mesh%n_bound ))
      do b = 1, mesh % n_bound
        ! copy bv_u to bv_w to keep the former unchanged
        call bv_u(b) % GetSlice(first=1, last=5, slice=bv_w(b), copy=.true.)
        ! generate pointer-based handles for pressure boundary values
        call bv_w(b) % GetSlice(first=4, last=4, slice=bv_p (b))
        call bv_w(b) % GetSlice(first=5, last=5, slice=bv_dp(b))
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

      ! optional filter
      if (allocated(this % A_fex)) then
        do d = 1, 3
          call SetArray(w(:,:,:,:,d), v(:,:,:,:,d))
          call TPO_AAA(this%A_fex, w(:,:,:,:,d), v(:,:,:,:,d))
        end do
      end if

      ! pressure computation ...................................................

      associate(div_v => w(:,:,:,:,4))

        ! outer traces of intermediate velocity using homogeneous Neumann BC
        call buf_vm % Transfer(mesh, vm, tag=100)
        call buf_vm % Merge(vm)
        call ConvertInnerToOuterTraces_3D(mesh, vm, vp)

        ! divergence of intermediate velocity
        call TPO_Div(ins_op % eop_u, ins_op % sem_u, v, vp, div_v)

        if (size(Q, 5) >= 4) then ! has additional RHS for mass conservation
          call MergeArrays(ONE, div_v, -ONE, Q(:,:,:,:,4))
        end if

        ! solve pressure equation
        call ins_op % PressureSolver( tau, bv_w, v, div_v, p &
                                    , i_max_p, r_red, r_max  &
                                    , precon = precon        )

      end associate

      ! pressure correction ....................................................

      associate(grad_p => w(:,:,:,:,1:3))

        ! generate outer traces of pressure -- sufficient for Neumann BC
        call GetOuterTraces_3D(mesh, p, pp)

        ! pressure gradient
        call TPO_Grad(ins_op % eop_u, ins_op % sem_u, p, pp, grad_p)

        ! correction: v = v - τ∇p
        call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)

      end associate

      ! diffusive correction ...................................................

      ! update boundary conditions
      do b = 1, mesh % n_bound
        select case(ins_op % bc_v(b))
        case('O')
          ! pᵇ = p - ∆pᵇ
          call bv_p(b) % Extract(p)
          call MergeArrays( ONE, bv_p (b) % val(:,:,:,1), &
                           -ONE, bv_dp(b) % val(:,:,:,1)  )
        end select
      end do

      associate(f => w(:,:,:,:,1:3))

        !$omp do collapse(2)
        do e = 1, mesh % n_elem
          do d = 1, 3
            f(:,:,:,e,d) = 1/tau * v(:,:,:,e,d) - F_d(:,:,:,e,d)
          end do
        end do

        if (update_viscosity) then
          call problem % GetViscosity(sem_u % metrics % x, t, u, nu)
        end if

        call ins_op % DiffusionSolver( tau, mu, nu, f, bv_w, v &
                                     , i_max_v, r_red, r_max   )

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w)
      deallocate(pp, vm, vp)
      deallocate(buf_vm)
      deallocate(bv_w, bv_p, bv_dp)
      !$omp end master

    end associate

  end subroutine ProjectionStep

  !=============================================================================

end submodule MP_ProjectionStep
