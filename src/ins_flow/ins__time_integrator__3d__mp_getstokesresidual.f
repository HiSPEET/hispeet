!> summary:  Stokes residual for incompressible flow
!> author:   Joerg Stiller
!> date:     2024/11/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_GetStokesResidual
  use TPO__Div__3D,  only: TPO_Div
  use TPO__Grad__3D, only: TPO_Grad
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Stokes residual for incompressible flow

  module subroutine GetStokesResidual(this, tau, mm, f, bv_u, mu, nu, u, r)

    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this

    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: mm(:,:,:,:)
    !< M, diagonal mass matrix
    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< f, nodal RHS, including old values, sources, convection ...
    class(BoundaryVariable_3D), intent(in)  :: bv_u(:)
    !< boundary values
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(in) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure
    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< Stokes residual

    ! internal variables .......................................................

    real(RNP), allocatable, save :: up(:,:,:,:,:) ! outer traces
    real(RNP), allocatable, save :: w (:,:,:,:,:) ! work

    integer :: ne, np
    integer :: e

    ! initialization ...........................................................

    np = size(u,1)
    ne = size(u,4)

    !$omp master
    allocate(up (np, np,  6, ne, 4), source = ZERO)
    allocate(w  (np, np, np, ne, 4), source = ZERO)
    !$omp end master
    !$omp barrier

    !...........................................................................

    associate( ins_op => this % ins_op   &
             , f_m    => f (:,:,:,:,1:3) &
             , f_c    => f (:,:,:,:, 4 ) &
             , r_m    => r (:,:,:,:,1:3) &
             , r_c    => r (:,:,:,:, 4 ) &
             , v      => u (:,:,:,:,1:3) &
             , p      => u (:,:,:,:, 4 ) &
             , vp     => up(:,:,:,:,1:3) &
             , pp     => up(:,:,:,:, 4 ) &
             , grad_p => w (:,:,:,:,1:3) )

      ! pressure gradient and velocity divergence ..............................

      call GetOuterTraces_3D(ins_op % mesh, u(:,:,:,1:4), up)
      call TPO_Grad(ins_op % eop_v, ins_op % sem_v, p, pp, grad_p)
      call TPO_Div( ins_op % eop_v, ins_op % sem_v, v, vp, r_c)

      ! diffusion residual .....................................................

      if (present(mu) .and. present(nu)) then
        call ins_op % GetDiffusionResidual(tau, f, bv_u, v, r_m)
      else
        call ins_op % GetDiffusionResidual(tau, mu, nu, f, bv_u, v, r_m)
      end if

      ! complete residual ......................................................

      !$omp do
      do e = 1, ne
        r_m(:,:,:,e,1) = r_m(:,:,:,e,1) - mm * grad_p(:,:,:,e,1)
        r_m(:,:,:,e,2) = r_m(:,:,:,e,2) - mm * grad_p(:,:,:,e,2)
        r_m(:,:,:,e,3) = r_m(:,:,:,e,3) - mm * grad_p(:,:,:,e,3)
        r_c(:,:,:,e)   = mm * (f_c(:,:,:,e) - r_c(:,:,:,e))
      end do

    end associate

    !$omp master
    deallocate(up, w)
    !$omp end master

  end subroutine GetStokesResidual

  !=============================================================================

end submodule MP_GetStokesResidual

