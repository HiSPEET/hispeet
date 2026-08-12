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

!> summary:  Implementation of FGMRES for the incompressible Stokes problem
!> author:   Simon Ehrmanntraut, Joerg Stiller
!> date:     2024/11/12
!>
!> FGMRES: Van der Vorst, Fig. 6.4
!===============================================================================

submodule (INS__Operator__3D) MP_StokesFGMRES
  use Array_Assignments
  use Array_Reductions
  use Logging_Levels
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FGMRES for Stokes part using projection step as a preconditioner
  !>
  !> On input `u` contains the approximate velocity `v` and pressure `p` in
  !> components 1:3 and 4, respectively.
  !> If `f_d0` is present, these values are overridden in the initial projection
  !> step, based on the extrapolated velocity `v = τ(f + f_d0)`.
  !> If `f_d0` is absent, the projection is used to correct the given values.

  module subroutine StokesFGMRES(this, tau, f, bv, mu, nu, u, f_d0)

    ! arguments ................................................................

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes time integrator
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

    ! internal variables .......................................................

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3

    ! auxiliary variables for Krylov iteration
    real(RNP), allocatable, save :: h(:,:)         ! Hessenberg matrix
    real(RNP), allocatable, save :: v(:,:,:,:,:,:)
    real(RNP), allocatable, save :: y(:)
    real(RNP), allocatable, save :: z(:,:,:,:,:,:)
    real(RNP) :: beta

    ! auxiliary variables for solving the least-squares problem
    real(RNP), allocatable, save :: r(:,:) ! r       in VDV03
    real(RNP), allocatable, save :: b(:)   ! \hat{b} in VDV03
    real(RNP), allocatable, save :: c(:)   ! c       in VDV03
    real(RNP), allocatable, save :: s(:)   ! s       in VDV03
    real(RNP) :: delta, gamma, rho

    ! auxiliary variables for preconditioning
    type(BoundaryVariable_3D), allocatable, save :: bv_z(:)

    ! auxiliary variables for residual computation
    real(RNP), allocatable, save :: mm_inv(:,:,:,:)
    real(RNP), allocatable, save :: g(:,:,:,:,:)
    real(RNP), allocatable, save :: O(:,:,:,:,:)

    real(RNP), save :: r_term
    logical  , save :: converged

    logical :: check_convergence
    integer :: na, nb, ne, ni, np, po
    integer :: i, j, k

    associate(mesh => this % mesh)

      ! preliminaries ..........................................................

      ! number of Arnoldi iterations
      ni = this % k_max

      ! dimensions
      po = this % eop_u % po
      np = po + 1
      ne = mesh % n_elem
      na = mesh % n_elem_active
      nb = mesh % n_bound

      check_convergence = this % r_red > 0 .or. this % r_max > 0

      !$omp master

      allocate(h(ni+1,ni), r(ni,ni)         , source = ZERO)
      allocate(b(ni+1), c(ni), s(ni), y(ni) , source = ZERO)

      allocate(v(np,np,np,ne,4,ni+1)        , source = ZERO)
      allocate(z(np,np,np,ne,4,ni)          , source = ZERO)

      allocate(bv_z(nb))
      do i = 1, nb
        call bv(i) % GetClone(bv_z(i), copy = .true.)
      end do
      call bv_z % SetToZero()

      allocate(mm_inv(np,np,np,ne))
      allocate(g(np,np,np,ne,4), source=ZERO)
      allocate(O(np,np,np,ne,4), source=ZERO)

      !$omp end master
      !$omp barrier

      call this % sem_u % Get_DG_DiagonalMassMatrix(mm_inv)

      !$omp do
      do i = 1, ne
        mm_inv(:,:,:,i) = ONE / mm_inv(:,:,:,i)
      end do

      associate(v1 => v(:,:,:,:,:,1))

        ! initial approximation ................................................

        call StokesProjection(this, tau, f, bv, mu, nu, u, f_d0)

        ! initial residual, v₁ = f - Au ........................................

        call this % GetStokesResidual(tau, f, bv, mu, nu, u, v1)

        beta = sqrt(ScalarProduct(v1, v1, mesh%comm_parts))
        b(1) = beta

        ! convergence check ....................................................

        !$omp master
        if (check_convergence) then
          ! set terminal condition
          r_term = huge(r_term)
          if (this % r_max > 0) r_term = this % r_max
          if (this % r_red > 0) r_term = max(r_term, beta * this%r_red)
          converged = beta <= r_term
          call XMPI_Bcast(converged, 0, mesh%comm_parts)
        else
          converged = .false.
        end if

        if (mesh % part == 0 .and. log_level_outer_iteration > 0) then
          write(*,'(2X,A,I4,A,ES12.5)') &
              'INS FGMRES solver, iteration j =',0,': beta =', beta
        end if
        !$omp end master
        !$omp barrier

        if (converged) then
          ni = 0
        end if

        ! first Krylov vector ..................................................

        call ScaleArray(v1, 1/max(beta,eps), multi=.true.)

      end associate

      ARNOLDI_ITERATION: do j = 1, ni
        associate( vj => v(:,:,:,:,:,j)   &
                 , w  => v(:,:,:,:,:,j+1) &
                 , zj => z(:,:,:,:,:,j)   )

          ! preconditioning ....................................................

          !$omp do
          do i = 1, na
            g(:,:,:,i,1) = mm_inv(:,:,:,i) * vj(:,:,:,i,1)
            g(:,:,:,i,2) = mm_inv(:,:,:,i) * vj(:,:,:,i,2)
            g(:,:,:,i,3) = mm_inv(:,:,:,i) * vj(:,:,:,i,3)
            g(:,:,:,i,4) = mm_inv(:,:,:,i) * vj(:,:,:,i,4)
          end do

          ! projection with homogeneous BC and frozen viscosity: z(j) = K⁻¹v(j)
          call StokesProjection( this, tau, g, bv_z, mu, nu, zj &
                               , f_d0 = O, precon = .true.      )

          ! application of homogeneous operator: w = A v(j)
          call this % GetStokesResidual(tau, O, bv_z, mu, nu, zj, w)
          call ScaleArray(w, -ONE, multi=.true.)

          ! computation of new Krylov vector ...................................

          ! orthogonalization against old Krylov vectors
          do i = 1, j
            h(i,j) = ScalarProduct(v(:,:,:,:,:,i), w, mesh%comm_parts)
            call MergeArrays(ONE, w, -h(i,j), v(:,:,:,:,:,i), multi=.true.)
          end do

          ! normalization, v(j+1) = w / ‖w‖
          h(j+1,j) = sqrt(ScalarProduct(w, w, mesh%comm_parts))
          call ScaleArray(w, ONE/max(h(j+1,j),eps), multi=.true.)

          ! Givens rotation transforming h to upper triagonal matrix r .........

          r(1,j) = h(1,j)

          do i = 2, j
            gamma    =  c(i-1) * r(i-1,j) + s(i-1) * h(i,j)
            r(i,j)   = -s(i-1) * r(i-1,j) + c(i-1) * h(i,j)
            r(i-1,j) = gamma
          end do

          delta  =  max(sqrt(r(j,j)**2 + h(j+1,j)**2), eps)
          c(j)   =  r(j  ,j) / delta
          s(j)   =  h(j+1,j) / delta
          r(j,j) =  c(j) * r(j,j) + s(j) * h(j+1,j)
          b(j+1) = -s(j) * b(j)
          b(j)   =  c(j) * b(j)

          ! convergence test ...................................................

          ! residual norm if Krylov iterations were exited now
          rho = abs(b(j+1))

          if (mesh % part == 0 .and. log_level_outer_iteration > 0) then
            !$omp master
            write(*,'(2X,A,I4,A,ES12.5)') &
              'INS FGMRES solver, iteration j =',j,': rho  =', rho
            !$omp end master
          end if

          if (check_convergence) then
            !$omp master
            converged = beta <= r_term
            call XMPI_Bcast(converged, 0, mesh%comm_parts)
            !$omp end master
            !$omp barrier
          end if
          if (converged) exit

        end associate
      end do ARNOLDI_ITERATION

      ! improved solution ....................................................

      j = min(j, ni)

      if (j > 0) then

        ! solve least-squares problem for y using backward substitution
        y(j) = b(j) / r(j,j)
        do i = j-1, 1, -1
          y(i) = (b(i) - dot_product(r(i,i+1:j), y(i+1:j))) / r(i,i)
        end do

        ! improved approximate solution
        do i = 1, j
        do k = 1, 4
          call MergeArrays(ONE, u(:,:,:,1:na,k), y(i), z(:,:,:,1:na,k,i))
        end do
        end do

      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(b, c, h, r, s, y, v, z)
      deallocate(bv_z, mm_inv, g, O)
      !$omp end master

    end associate

  end subroutine StokesFGMRES

  !=============================================================================

end submodule MP_StokesFGMRES
