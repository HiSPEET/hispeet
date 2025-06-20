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

    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work
    real(RNP), allocatable, save :: wp (:,:,:,:,:) ! outer traces w⁺

    type(BoundaryVariable_3D), allocatable, save :: bv_w(:), bv_p(:), bv_dp(:)

    logical :: extrapolation, predictor
    integer :: b, c, e, na, ne, np

    associate( problem => this % problem          &
             , sem_u   => this %  sem_u           &
             , mesh    => this %  mesh            &
             , v       => u(:,:,:,:,1:3)          &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      predictor = present(f_d0)
      extrapolation = predictor .or. this%stokes_corrector(1:1) == 'X'

      np = size(v,1)
      na = mesh % n_elem_active
      ne = mesh % n_elem

      !$omp master

      allocate( wp (np, np,  6, ne, 6), source = ZERO )
      allocate( w  (np, np, np, ne, 5), source = ZERO )

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

      if (predictor) then

        !$omp do collapse(2)
        do c = 1, 3
        do e = 1, na
          v(:,:,:,e,c) =  tau * (f(:,:,:,e,c) + f_d0(:,:,:,e,c))
        end do
        end do

      else if (extrapolation) then
        ! corrector using extrapolation

        associate( f_d => w (:,:,:,:,1:3) &
                 , mm  => w (:,:,:,:, 4)  &
                 , vp  => wp(:,:,:,:,1:3) &
                 , sp  => wp(:,:,:,:,4:6) )

          if ((this % stokes_corrector(2:2) == 'R')) then
            ! using rotational form of the diffusion term
            call this % GetDiffusionTerm( mu, nu, v, vp, sp, f_d, bv_w &
                                        , xout = .true. , form = 2     )
          else
            ! using default form of the diffusion term
            call this % GetDiffusionTerm( mu, nu, v, vp, sp, f_d, bv_w &
                                        , xout = .true.                )
          end if

          call this % sem_u % Get_DG_DiagonalMassMatrix( mm )

          !$omp do
          do e = 1, na
            mm(:,:,:,e) = 1 / mm(:,:,:,e)
            do c = 1,3
              v(:,:,:,e,c) = tau * (f(:,:,:,e,c) + f_d(:,:,:,e,c) * mm(:,:,:,e))
            end do
          end do

        end associate
      end if

      ! projection step ........................................................

      ! sources
      associate( div_v => w (:,:,:,:, 4)  &
               , vp    => wp(:,:,:,:,1:3) )
        ! compute divergence of approximate velocity
        call GetOuterVectorTraces_3D(mesh, v, vp)
        call TPO_Div(this % eop_u, this % sem_u, v, vp, div_v)
        ! add additional sources
        if (size(f, 5) >= 4) then
          call MergeArrays(ONE, div_v, -ONE, f(:,:,:,:,4))
        end if
      end associate

      ! pressure and velocity correction
      if (extrapolation) then

        ! predictor or corrector with extrapolation
        associate( grad_p => w (:,:,:,:,1:3) &
                 , div_v  => w (:,:,:,:, 4)  &
                 , pp     => wp(:,:,:,:, 4 ) )
          ! compute pressure
          call this % PressureSolver(tau, bv_w, v, div_v, p, precon)
          ! compute pressure gradient
          call GetOuterTraces_3D(mesh, p, pp)
          call TPO_Grad(this % eop_u, this % sem_u, p, pp, grad_p)
          ! correct velocity: v = v - τ∇p
          do c = 1, 3
            call MergeArrays(ONE, v(:,:,:,:na,c), -tau, grad_p(:,:,:,:na,c))
          end do
        end associate

      else

        ! corrector without extrapolation
        associate( grad_q => w (:,:,:,:,1:3) &
                 , div_v  => w (:,:,:,:, 4)  &
                 , q      => w (:,:,:,:, 5 ) &
                 , qp     => wp(:,:,:,:, 4 ) )
          ! compute pressure correction δp = q
          call SetArray(q, ZERO)
          call this % PressureSolver(tau, bv_w, v, div_v, q, precon)
          ! compute gradient of pressure correction δp
          call GetOuterTraces_3D(mesh, q, qp)
          call TPO_Grad(this % eop_u, this % sem_u, q, qp, grad_q)
          ! correction: v = v - τ∇p, p = p + δp
          do c = 1, 3
            call MergeArrays(ONE, v(:,:,:,:na,c), -tau, grad_q(:,:,:,:na,c))
          end do
          call MergeArrays(ONE, p, ONE, q)
        end associate

      end if

      ! diffusion step .........................................................

      ! update outflow boundary conditions
      do b = 1, mesh % n_bound
        select case(problem % bc_v(b))
        case('O')
          ! pᵇ = p - ∆pᵇ
          call bv_p(b) % Extract(p)
          call MergeArrays( ONE, bv_p (b) % val(:,:,:,1), &
                           -ONE, bv_dp(b) % val(:,:,:,1)  )
        end select
      end do

      associate( q  => w (:,:,:,:,1:3) &
               , pp => wp(:,:,:,:, 4 ) )

        if (predictor) then

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, na
            q(:,:,:,e,c) = 1/tau * v(:,:,:,e,c) - f_d0(:,:,:,e,c)
          end do
          end do

        else

          ! q = ∇p
          call GetOuterTraces_3D(mesh, p, pp)
          call TPO_Grad(this % eop_u, this % sem_u, p, pp, q)

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, na
            ! q = f_m - ∇p
            q(:,:,:,e,c) = f(:,:,:,e,c) - q(:,:,:,e,c)
          end do
          end do

        end if

        call this % DiffusionSolver(tau, mu, nu, q, bv_w, v, precon)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w, wp)
      deallocate(bv_w, bv_p, bv_dp)
      !$omp end master

    end associate

  end subroutine StokesProjection

  !=============================================================================

end submodule MP_StokesProjection
