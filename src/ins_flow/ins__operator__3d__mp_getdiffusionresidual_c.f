!> summary:  Incompressible Navier-Stokes DG-SEM diffusion residual
!> author:   Joerg Stiller
!> date:     2022/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(INS__Operator__3D) MP_GetDiffusionResidual_C
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion residual with constant viscosity
  !>
  !> Computes the DG-SEM residual of  of the viscous diffusion term including
  !> the implicit part of the discretized time derivative, i.e.,
  !>
  !>     r = Fd(v, vb, sb) - λMv + Mf
  !>
  !> where `Fd` is the weak form of the diffusion term for the given velocity
  !> `v` and boundary values `vb`, `sb` extracted from the boundary variable
  !> `bv_v`, `M` is the diagonal mass matrix and `f` the nodal coefficients of
  !> the sources, which comprise the remaining coefficients of the momentum
  !> equation.
  !>
  !> Omitting `f` and `bv_v` results in the application of the homogeneous
  !> operator.

  module subroutine GetDiffusionResidual_C(this, lambda, nu, v, r, f, bv_v)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: lambda
    !< λ, coefficient of the linear term

    real(RNP), intent(in) :: nu
    !< ν, kinematic shear viscosity

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    real(RNP), contiguous, intent(out) :: r(:,:,:,:,:)
    !< residual, r(np,np,np,ne,3)

    real(RNP), contiguous, optional, intent(in) :: f(:,:,:,:,:)
    !< sources, f(np,np,np,ne,3)

    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv_v(:)
    !< velocity boundary values, bv_v(nb)

    ! local variables ..........................................................

    real(RNP), allocatable, save :: M(:,:,:,:)    ! diagonal mass matrix
    real(RNP), allocatable, save :: vp(:,:,:,:,:) ! velocity traces v⁺
    real(RNP), allocatable, save :: sp(:,:,:,:,:) ! viscous flux traces s⁺

    integer :: np
    integer :: d, e, i, j, k

    associate(mesh => this % sem_v % mesh)

      ! initialization .........................................................

      np = size(v,1)

      !$omp master
      allocate( M  (np, np, np, mesh%n_elem) )
      allocate( vp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      allocate( sp (np, np,  6, mesh%n_elem, 3), source = ZERO )
      !$omp end master
      !$omp barrier

      call this % sem_v % Get_DG_DiagonalMassMatrix(M)

      if (present(bv_v)) then
        call this % SetVelocityBC(bv_v, vp, sp)
      end if

      ! compute residual .......................................................

      call this % GetDiffusionTerm(nu, v, vp, sp, r)

      if (present(f)) then
        !$omp do
        do e = 1, mesh % n_elem
          do k = 1, np
          do j = 1, np
          do i = 1, np
            do d = 1, 3
              r(i,j,k,e,d) = r(i,j,k,e,d) &
                           + M(i,j,k,e) * (f(i,j,k,e,d) - lambda * v(i,j,k,e,d))
            end do
          end do
          end do
          end do
        end do
      else
        !$omp do
        do e = 1, mesh % n_elem
          do k = 1, np
          do j = 1, np
          do i = 1, np
            do d = 1, 3
              r(i,j,k,e,d) = r(i,j,k,e,d) - lambda * M(i,j,k,e) * v(i,j,k,e,d)
            end do
          end do
          end do
          end do
        end do
      end if

      ! cleanup ................................................................

      !$omp master
      deallocate(M, vp, sp)
      !$omp end master

    end associate

  end subroutine GetDiffusionResidual_C

  !=============================================================================

end submodule MP_GetDiffusionResidual_C
