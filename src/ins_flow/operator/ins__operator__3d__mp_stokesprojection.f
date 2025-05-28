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
  !>
  !> On input `u` contains the approximate velocity `v` and pressure `p` in
  !> components 1:3 and 4, respectively.
  !> If `f_d0` is present, these values are overridden by an projection step
  !> based on the extrapolated velocity `v = τ(f + f_d0)`.
  !> If `f_d0` is absent, the projection is used to correct the given values.
  !> Finally, a diffusion step is performed to enforce the viscous terms and
  !> boundary conditions.
  !>
  !> The optional argument `precon` can be given to run the pressure and
  !> diffusion solvers in preconditioner mode.

  module subroutine StokesProjection(this, tau, f, bv, mu, nu, u, f_d0, precon)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< RHS: f = v₀/τ + F_c + f_s + ...
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]
    real(RNP), contiguous, optional, intent(in) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(in) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    real(RNP), contiguous, optional, intent(in) :: f_d0(:,:,:,:,:)
    !< approximate diffusion term
    logical, optional, intent(in) :: precon
    !< if present, operate as preconditioner (regardless of the value)

    ! internal variables .......................................................

    real(RNP), allocatable, save :: pp (:,:,:,:)   ! outer pressure traces p⁺
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work

    type(BoundaryVariable_3D), allocatable, save :: bv_w(:), bv_p(:), bv_dp(:)

    logical :: extrapolation
    integer :: b, d, e, np

    associate( problem => this % problem          &
             , sem_u   => this %  sem_u           &
             , mesh    => this %  mesh            &
             , n_elem  => this %  mesh % n_elem   &
             , n_ghost => this %  mesh % n_ghost  &
             , v       => u(:,:,:,:,1:3)          &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      extrapolation = present(f_d0)
      np = size(v,1)

      !$omp master

      allocate( pp (np, np,  6, n_elem   ), source = ZERO )
      allocate( vp (np, np,  6, n_elem, 3), source = ZERO )
      allocate( w  (np, np, np, n_elem, 4), source = ZERO )

      ! provide handles for velocity and pressure boundary values
      allocate(bv_w ( mesh%n_bound ))
      allocate(bv_p ( mesh%n_bound ))
      allocate(bv_dp( mesh%n_bound ))
      do b = 1, mesh % n_bound
        ! copy bv to bv_w to keep the former unchanged
        call bv(b) % GetSlice(first=1, last=4, slice=bv_w(b), copy=.true.)
        ! generate pointer-based handles for pressure boundary values
        call bv_w(b) % GetSlice(first=4, last=4, slice=bv_p (b))
        call bv_w(b) % GetSlice(first=3, last=3, slice=bv_dp(b))
      end do
      !$omp end master
      !$omp barrier

      ! extrapolation step .....................................................

      if (extrapolation) then
        !$omp do collapse(2)
        do d = 1, 3
        do e = 1, mesh % n_elem
          v(:,:,:,e,d) =  tau * (f(:,:,:,e,d) + f_d0(:,:,:,e,d))
        end do
        end do
      end if

      ! pressure correction.....................................................

      associate(div_v => w(:,:,:,:,1))

        ! divergence of approximate velocity
        call GetOuterTraces_3D(mesh, v, vp)
        call TPO_Div(this % eop_u, this % sem_u, v, vp, div_v)

        if (size(f, 5) >= 4) then ! has additional RHS for mass conservation
          call MergeArrays(ONE, div_v, -ONE, f(:,:,:,:,4))
        end if

        if (extrapolation) then
          associate(grad_p => w(:,:,:,:,1:3))
            ! compute pressure
            call this % PressureSolver(tau, bv_w, v, div_v, p, precon)
            ! compute pressure gradient
            call GetOuterTraces_3D(mesh, p, pp)
            call TPO_Grad(this % eop_u, this % sem_u, p, pp, grad_p)
            ! correct velocity: v = v - τ∇p
            call MergeArrays(ONE, v, -tau, grad_p, multi=.true.)
          end associate
        else
          associate(grad_dp => w(:,:,:,:,1:3), dp => w(:,:,:,:,4))
            ! compute pressure correction
            call SetArray(dp, ZERO)
            call this % PressureSolver(tau, bv_w, v, div_v, dp, precon)
            ! compute gradient of pressure correction δp
            call GetOuterTraces_3D(mesh, dp, pp)
            call TPO_Grad(this % eop_u, this % sem_u, dp, pp, grad_dp)
            ! correction: v = v - τ∇p, p = p + δp
            call MergeArrays(ONE, v, -tau, grad_dp, multi=.true.)
            call MergeArrays(ONE, p,  ONE, dp)
          end associate
        end if

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

        if (extrapolation) then
          !$omp do collapse(2)
          do d = 1, 3
          do e = 1, mesh % n_elem
            q(:,:,:,e,d) = 1/tau * v(:,:,:,e,d) - f_d0(:,:,:,e,d)
          end do
          end do
        else

          call GetOuterTraces_3D(mesh, p, pp)
          call TPO_Grad(this % eop_u, this % sem_u, p, pp, q)
          !$omp do collapse(2)
          do d = 1, 3
          do e = 1, mesh % n_elem
            q(:,:,:,e,d) = f(:,:,:,e,d) - q(:,:,:,e,d)
          end do
          end do
        end if

        call this % DiffusionSolver(tau, mu, nu, q, bv_w, v, precon)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w)
      deallocate(pp, vp)
      deallocate(bv_w, bv_p, bv_dp)
      !$omp end master

    end associate

  end subroutine StokesProjection

  !=============================================================================

end submodule MP_StokesProjection
