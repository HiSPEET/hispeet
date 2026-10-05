!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Generic projection-diffusion step for incompressible flow
!> author:   Joerg Stiller
!> date:     2022/09/24
!===============================================================================

submodule (INS__Operator__3D) MP_StokesProjection
  use Array_Assignments
  use TPO__AAA__3D
  use TPO__Div__3D
  use TPO__Grad__3D
  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

  ! enforce BC before computation of divergence (NOT RECOMMENDED)
  logical, parameter :: use_bc_for_div = .false.

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

  module subroutine StokesProjection(this, tau, mu, nu, bv, f_d0, f, u, precon)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, optional, intent(in) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(in) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]
    real(RNP), contiguous, optional, intent(in) :: f_d0(:,:,:,:,:)
    !< approximate diffusion term
    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< RHS: f = v₀/τ + F_c + f_s + ...
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    logical, optional, intent(in) :: precon
    !< if present, operate as preconditioner (regardless of the value)

    ! internal variables .......................................................

    real(RNP), allocatable, save :: w  (:,:,:,:,:) ! work
    real(RNP), allocatable, save :: tr (:,:,:,:,:) ! traces

    type(BoundaryVariable_3D), allocatable, save :: bv_w(:), bv_p(:), bv_dp(:)

    real(RNP), allocatable :: mm_inv(:,:,:)
    logical :: predictor
    integer :: b, c, e, na, ne, np

    associate( problem => this % problem          &
             , sem_u   => this %  sem_u           &
             , mesh    => this %  mesh            &
             , v       => u(:,:,:,:,1:3)          &
             , p       => u(:,:,:,:,4)            )

      ! initialization .........................................................

      predictor = present(f_d0)

      np = size(v,1)
      na = mesh % n_elem_active
      ne = mesh % n_elem

      !$omp master

      allocate( tr (np, np,  6, ne, 7), source = ZERO )
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

      associate( f_d    => w (:,:,:,:,1:3) &
               , q_d    => w (:,:,:,:,1:3) &
               , grad_p => w (:,:,:,:,1:3) &
               , div_v  => w (:,:,:,:, 4 ) &
               , mm     => w (:,:,:,:, 5 ) &
               , vp     => tr(:,:,:,:,1:3) &
               , sp     => tr(:,:,:,:,4:6) &
               , pp     => tr(:,:,:,:, 7 ) )

        ! extrapolation step ...................................................

        if (predictor) then

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, na
            v(:,:,:,e,c) =  tau * (f(:,:,:,e,c) + f_d0(:,:,:,e,c))
          end do
          end do

        else ! corrector

          call this % GetDiffusionTerm( mu, nu, bv_w, f_d, v, vp, sp &
                                      , xout = .true. , form = 2     )

          call this % sem_u % Get_DG_DiagonalMassMatrix( mm )

          !$omp do
          do e = 1, na
            mm_inv = 1 / mm(:,:,:,e)
            do c = 1, 3
              v(:,:,:,e,c) = tau * ( f(:,:,:,e,c) + mm_inv * f_d(:,:,:,e,c) )
            end do
          end do

        end if

        ! projection step ......................................................

        ! outer traces of extrapolated velocity
        call GetOuterVectorTraces_3D(mesh, v, vp)

        ! optionally enforce boundary conditions
        if (use_bc_for_div) then
          call this % ApplyEssentialBC(bv_w, vp, vp)
        end if

        ! divergence of approximate velocity
        call TPO_Div(this % eop_u, this % sem_u, v, vp, div_v)

        ! add additional sources
        if (size(f, 5) >= 4) then
          call MergeArrays(ONE, div_v, -ONE, f(:,:,:,:,4))
        end if

        ! compute pressure
        call this % PressureSolver(tau, bv_w, v, div_v, p, precon)

        ! compute pressure gradient
        call GetOuterTraces_3D(mesh, p, pp)
        call TPO_Grad(this % eop_u, this % sem_u, p, pp, grad_p)

        ! correct velocity: v = v - τ∇p
        do c = 1, 3
          call MergeArrays(ONE, v(:,:,:,:na,c), -tau, grad_p(:,:,:,:na,c))
        end do

        ! diffusion step .......................................................

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

        ! diffusion sources
        if (predictor) then

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, na
            q_d(:,:,:,e,c) = 1/tau * v(:,:,:,e,c) - f_d0(:,:,:,e,c)
          end do
          end do

        else

          !$omp do collapse(2)
          do c = 1, 3
          do e = 1, na
            ! q = f_m - ∇p
            q_d(:,:,:,e,c) = f(:,:,:,e,c) - grad_p(:,:,:,e,c)
          end do
          end do

        end if

        call this % DiffusionSolver(tau, mu, nu, bv_w, q_d, v, precon=precon)

      end associate

      ! cleanup ................................................................

      !$omp master
      deallocate(w, tr)
      deallocate(bv_w, bv_p, bv_dp)
      !$omp end master

    end associate

  end subroutine StokesProjection

  !=============================================================================

end submodule MP_StokesProjection
