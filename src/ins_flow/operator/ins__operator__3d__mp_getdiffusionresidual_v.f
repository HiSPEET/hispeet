!> summary:  Incompressible Navier-Stokes DG-SEM diffusion residual (DV)
!> author:   Joerg Stiller
!> date:     2024/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetDiffusionResidual_V
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion residual with constant viscosity
  !>
  !> Computes the DG-SEM residual of  of the viscous diffusion term including
  !> the implicit part of the discretized time derivative, i.e.,
  !>
  !>     r = F_d(v, vb, sb) - Mv/τ + Mf
  !>
  !> where `F_d` is the weak form of the diffusion term for the given velocity
  !> `v` and boundary values `vb`, `sb` obtained from the boundary variable
  !> `bv`, `M` is the diagonal mass matrix and `f` the nodal coefficients of
  !> the sources, which comprise the remaining coefficients of the momentum
  !> equation.

  module subroutine GetDiffusionResidual_V(this, tau, mu, nu, f, bv, v, r, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, intent(in) :: mu(:,:,:,:)
    !< kinematic bulk viscosity μ (np,np,np,ne)

    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    !< kinematic shear viscosity ν (np,np,np,ne)

    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< sources, f(np,np,np,ne,3)

    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values
    !!   - Γᴰ :  [ v₁, v₂, v₃, - ]
    !!   - Γᴼ :  [ - , - , ∆p, p ]

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< residual, r(np,np,np,ne,3)

    integer, optional, intent(in) :: form
    !< form of `∇⋅τ`: 0/1/2 ↔︎ default/diffusion/rotational [0]

    ! local variables ..........................................................

    real(RNP), allocatable, save :: mm(:,:,:,:)   ! diagonal mass matrix M
    real(RNP), allocatable, save :: vp(:,:,:,:,:) ! velocity traces v⁺
    real(RNP), allocatable, save :: sp(:,:,:,:,:) ! viscous flux traces s⁺

    real(RNP) :: lambda
    integer   :: np
    integer   :: d, e

    associate(mesh => this % sem_u % mesh)

      ! initialization .........................................................

      np = size(v,1)

      !$omp master
      allocate( mm (np, np, np, mesh%n_elem) )
      allocate( vp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      allocate( sp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      !$omp end master
      !$omp barrier

      call this % sem_u % Get_DG_DiagonalMassMatrix(mm)

      lambda = 1 / tau

      ! compute residual .......................................................

      call this % GetDiffusionTerm_V(mu, nu, v, vp, sp, r, bv, form=form)

      !$omp do collapse(2)
      do e = 1, mesh % n_elem
        do d = 1, 3
          r(:,:,:,e,d) = r(:,:,:,e,d) &
                       + mm(:,:,:,e) * (f(:,:,:,e,d) - lambda * v(:,:,:,e,d))
        end do
      end do

      ! cleanup ................................................................

      !$omp master
      deallocate(mm, vp, sp)
      !$omp end master

    end associate

  end subroutine GetDiffusionResidual_V

  !=============================================================================

end submodule MP_GetDiffusionResidual_V
