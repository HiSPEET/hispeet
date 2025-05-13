!> summary:  Generic projection-diffusion step for incompressible flow
!> author:   Joerg Stiller
!> date:     2022/09/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (INS__Operator__3D) MP_StokesProjection
  use Array_Assignments
  use TPO__AAA__3D
  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Projection-diffusion step for incompressible flow

  module subroutine StokesProjection(this, tau, f_d0, f, bv, mu, nu, u, precon)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: f_d0(:,:,:,:,:)
    !< approximate diffusion term
    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< RHS: f = v₀/τ + F_c + f_s + ...
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - , -  ]
    !!   - Γᴼ :  [ - , - , - , p , ∆p ]
    real(RNP), contiguous, optional, intent(in) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(in) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    logical, optional, intent(in) :: precon
    !< if present, operate as preconditioner (regardless of the value)

    ! internal variables .......................................................

    real(RNP), allocatable, save :: pp (:,:,:,:)   ! outer pressure traces p⁺
    real(RNP), allocatable, save :: vm (:,:,:,:,:) ! inner velocity traces v⁻
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_vm

    type(BoundaryVariable_3D), allocatable, save :: bv_w(:), bv_p(:), bv_dp(:)

    integer :: b, d, e, np

    associate( problem => this % problem          &
             , sem_u   => this %  sem_u           &
             , mesh    => this %  mesh            &
             , n_elem  => this %  mesh % n_elem   &
             , n_ghost => this %  mesh % n_ghost  &
             , v       => u(:,:,:,:,1:3)          &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

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
        ! copy bv to bv_w to keep the former unchanged
        call bv(b) % GetSlice(first=1, last=5, slice=bv_w(b), copy=.true.)
        ! generate pointer-based handles for pressure boundary values
        call bv_w(b) % GetSlice(first=4, last=4, slice=bv_p (b))
        call bv_w(b) % GetSlice(first=5, last=5, slice=bv_dp(b))
      end do
      !$omp end master
      !$omp barrier

      ! extrapolation step .....................................................

      ! Computes the intermediate velocity v' = τ(f + f_d0) and extracts
      ! the inner traces v'⁻ needed in the next step.

      !$omp do collapse(2)
      do d = 1, 3
      do e = 1, mesh % n_elem

        ! intermediate velocity v = v'
        v(:,:,:,e,d) =  tau * (f(:,:,:,e,d) + f_d0(:,:,:,e,d))

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

        ! outer traces of intermediate velocity using homogeneous Neumann BC
        call buf_vm % Transfer(mesh, vm, tag=100)
        call buf_vm % Merge(vm)
        call ConvertInnerToOuterTraces_3D(mesh, vm, vp)

        ! divergence of intermediate velocity
        call TPO_Div(this % eop_u, this % sem_u, v, vp, div_v)

        if (size(f, 5) >= 4) then ! has additional RHS for mass conservation
          call MergeArrays(ONE, div_v, -ONE, f(:,:,:,:,4))
        end if

        ! solve pressure equation
        call this % PressureSolver(tau, bv_w, v, div_v, p, precon)

      end associate

      ! pressure correction ....................................................

      associate(grad_p => w(:,:,:,:,1:3))

        ! generate outer traces of pressure -- sufficient for Neumann BC
        call GetOuterTraces_3D(mesh, p, pp)

        ! pressure gradient
        call TPO_Grad(this % eop_u, this % sem_u, p, pp, grad_p)

        ! correction: v = v - τ∇p
        call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)

      end associate

      ! diffusive correction ...................................................

      ! update boundary conditions
      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('O')
          ! pᵇ = p - ∆pᵇ
          call bv_p(b) % Extract(p)
          call MergeArrays( ONE, bv_p (b) % val(:,:,:,1), &
                           -ONE, bv_dp(b) % val(:,:,:,1)  )
        end select
      end do

      associate(q => w(:,:,:,:,1:3))

        !$omp do collapse(2)
        do e = 1, mesh % n_elem
          do d = 1, 3
            q(:,:,:,e,d) = 1/tau * v(:,:,:,e,d) - f_d0(:,:,:,e,d)
          end do
        end do

        call this % DiffusionSolver(tau, mu, nu, q, bv_w, v, precon)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w)
      deallocate(pp, vm, vp)
      deallocate(buf_vm)
      deallocate(bv_w, bv_p, bv_dp)
      !$omp end master

    end associate

  end subroutine StokesProjection

  !=============================================================================

end submodule MP_StokesProjection
