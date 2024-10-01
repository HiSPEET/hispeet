!> summary:  Homogeneous incompressible Navier-Stokes DG-SEM diffusion operator
!> author:   Joerg Stiller
!> date:     2024/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyDiffusionOperator_V
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Homogeneous diffusion operator with variable viscosity
  !>
  !> Computes the homogeneous DG-SEM viscous diffusion operator including the
  !> implicit part of the discretized time derivative, i.e.,
  !>
  !>     r = Mv/τ - Fd(v, vb=0, sb=0)
  !>
  !> where `Fd` is the weak form of the diffusion term for the given velocity
  !> `v` with zero boundary values `vb`, `sb` and `M` is diagonal mass matrix.

  module subroutine ApplyDiffusionOperator_V(this, tau, mu, nu, v, r, form)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, intent(in) :: mu(:,:,:,:)
    !< kinematic bulk viscosity μ (np,np,np,ne)

    real(RNP), contiguous, intent(in) :: nu(:,:,:,:)
    !< kinematic shear viscosity ν (np,np,np,ne)

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< result, r(np,np,np,ne,3)

    integer, optional, intent(in) :: form
    !< form of `∇⋅τ`: 0/1/2 ↔︎ default/diffusion/rotational [0]

    ! local variables ..........................................................

    real(RNP), allocatable, save :: mm(:,:,:,:)   ! diagonal mass matrix
    real(RNP), allocatable, save :: vp(:,:,:,:,:) ! velocity traces v⁺
    real(RNP), allocatable, save :: sp(:,:,:,:,:) ! viscous flux traces s⁺

    real(RNP) :: lambda
    integer   :: np
    integer   :: d, e

    associate(mesh => this % sem_v % mesh)

      ! initialization .........................................................

      np = size(v,1)

      !$omp master
      allocate( mm (np, np, np, mesh%n_elem) )
      allocate( vp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      allocate( sp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      !$omp end master
      !$omp barrier

      call this % sem_v % Get_DG_DiagonalMassMatrix(mm)

      lambda = 1 / tau

      ! computation ............................................................

      call this % GetDiffusionTerm_V(mu, nu, v, vp, sp, r, form=form)

      !$omp do collapse(2)
      do e = 1, mesh % n_elem
        do d = 1, 3
          r(:,:,:,e,d) = lambda * mm(:,:,:,e) * v(:,:,:,e,d) - r(:,:,:,e,d)
        end do
      end do

      ! cleanup ................................................................

      !$omp master
      deallocate(mm, vp, sp)
      !$omp end master

    end associate

  end subroutine ApplyDiffusionOperator_V

  !=============================================================================

end submodule MP_ApplyDiffusionOperator_V
