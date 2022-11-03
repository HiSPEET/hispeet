!> summary:  Homogeneous incompressible Navier-Stokes DG-SEM diffusion operator
!> author:   Joerg Stiller
!> date:     2022/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_ApplyDiffusionOperator_C
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Homogeneous diffusion operator with constant viscosity
  !>
  !> Computes the homogeneous DG-SEM viscous diffusion operator including the
  !>implicit part of the discretized time derivative, i.e.,
  !>
  !>     r = Mv/τ - Fd(v, vb=0, sb=0)
  !>
  !> where `Fd` is the weak form of the diffusion term for the given velocity
  !> `v` with zero boundary values `vb`, `sb` and `M` is diagonal mass matrix.

  module subroutine ApplyDiffusionOperator_C(this, tau, v, r)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< result, r(np,np,np,ne,3)

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

      ! compute residual .......................................................

      call this % GetDiffusionTerm(v, vp, sp, r)

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

  end subroutine ApplyDiffusionOperator_C

  !=============================================================================

end submodule MP_ApplyDiffusionOperator_C
