!> summary:  Implementation of FGMRES for the incompressible Stokes problem
!> author:   Simon Ehrmanntraut, Joerg Stiller
!> date:     2024/11/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> FGMRES: Van der Vorst, Fig. 6.4
!===============================================================================

submodule (INS__Time_Integrator__3D) MP_FGMRES_Step
  use Array_Assignments
  use Array_Reductions
  use Logging_Levels
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> FGMRES for Stokes part using projection step as a preconditioner

  module subroutine FGMRES_Step(this, tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u)

    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this
    !< incompressible Navier-Stokes time integrator
    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), intent(in) :: t
    !< t, final time
    real(RNP), contiguous, intent(in) :: v_0(:,:,:,:,:)
    !< v₀, effective initial value of velocity
    real(RNP), contiguous, intent(in) :: F_c(:,:,:,:,:)
    !< convective term at final time t
    real(RNP), contiguous, intent(in) :: F_d(:,:,:,:,:)
    !< diffusion term at final time t
    real(RNP), contiguous, intent(in) :: Q(:,:,:,:,:)
    !< sources at time t and further known terms
    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  [ v₁, v₂, v₃, - , -  ]
    !!   - Γᴼ :  [ - , - , - , p , ∆p ]
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure at final time u

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
    integer :: nb, ne, ni, np, po
    integer :: i, j

    associate( ins_op => this % ins_op                     &
             , mesh   => this % ins_op % mesh              &
             , comm   => this % ins_op % mesh % comm_parts )

      ! preliminaries ..........................................................

      ! number of Arnoldi iterations
      ni = this%i_krylov

      ! dimensions
      po = ins_op % eop_u % po
      np = po + 1
      nb = mesh % n_bound
      ne = mesh % n_elem

      check_convergence = this % r_red > 0 .or. this % r_max > 0

      !$omp master

      allocate(h(ni+1,ni), r(ni,ni)         , source = ZERO)
      allocate(b(ni+1), c(ni), s(ni), y(ni) , source = ZERO)

      allocate(v(np,np,np,ne,4,ni+1)        , source = ZERO)
      allocate(z(np,np,np,ne,4,ni)          , source = ZERO)

      allocate(bv_z(nb))
      do i = 1, nb
        call bv_u(i) % GetClone(bv_z(i), copy = .true.)
      end do
      call bv_z % SetToZero()

      allocate(mm_inv(np,np,np,ne))
      allocate(g(np,np,np,ne,4), source=ZERO)
      allocate(O(np,np,np,ne,4), source=ZERO)

      !$omp end master
      !$omp barrier

      call ins_op % sem_u % Get_DG_DiagonalMassMatrix(mm_inv)

      !$omp do
      do i = 1, ne
        mm_inv(:,:,:,i) = ONE / mm_inv(:,:,:,i)
        g(:,:,:,i,1) = Q(:,:,:,i,1) + 1/tau * v_0(:,:,:,i,1) + F_c(:,:,:,i,1)
        g(:,:,:,i,2) = Q(:,:,:,i,2) + 1/tau * v_0(:,:,:,i,2) + F_c(:,:,:,i,2)
        g(:,:,:,i,3) = Q(:,:,:,i,3) + 1/tau * v_0(:,:,:,i,3) + F_c(:,:,:,i,3)
        if (size(Q,5) >= 4) then
          g(:,:,:,i,4) = Q(:,:,:,i,4)
        end if
      end do

      !$omp master

      associate(v1 => v(:,:,:,:,:,1))

        ! initial approximation ................................................

        call this % ProjectionStep( tau, t, v_0, F_c, F_d, Q, bv_u, mu, nu, u &
                                  , this % i_max_p, this % i_max_v            )

        ! initial residual, v₁ = f - Au ........................................

        call ins_op % GetStokesResidual(tau, g, bv_u, mu, nu, u, v1)

        beta = sqrt(ScalarProduct(v1, v1, comm))
        b(1) = beta

        ! convergence check ....................................................

        !$omp master
        if (check_convergence) then
          ! set terminal condition
          r_term = huge(r_term)
          if (this % r_max > 0) r_term = this % r_max
          if (this % r_red > 0) r_term = max(r_term, beta * this%r_red)
          converged = beta <= r_term
          call XMPI_Bcast(converged, 0, comm)
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
          do i = 1, ne
            g(:,:,:,i,1) = mm_inv(:,:,:,i) * vj(:,:,:,i,1)
            g(:,:,:,i,2) = mm_inv(:,:,:,i) * vj(:,:,:,i,2)
            g(:,:,:,i,3) = mm_inv(:,:,:,i) * vj(:,:,:,i,3)
            g(:,:,:,i,4) = mm_inv(:,:,:,i) * vj(:,:,:,i,4)
          end do

          ! projection with homogeneous BC and frozen viscosity: z(j) = K⁻¹v(j)
          call this % ProjectionStep( tau, t, O, O, O, g, bv_z, mu, nu, zj &
                                    , this % i_pre_p, this % i_pre_v       &
                                    , freeze = .true.                      )

          ! application of homogeneous operator: w = A v(j)
          call ins_op % GetStokesResidual(tau, O, bv_z, mu, nu, zj, w)
          call ScaleArray(w, -ONE, multi=.true.)

          ! computation of new Krylov vector ...................................

          ! orthogonalization against old Krylov vectors
          do i = 1, j
            h(i,j) = ScalarProduct(v(:,:,:,:,:,i), w, comm)
            call MergeArrays(ONE, w, -h(i,j), v(:,:,:,:,:,i), multi=.true.)
          end do

          ! normalization, v(j+1) = w / ‖w‖
          h(j+1,j) = sqrt(ScalarProduct(w, w, comm))
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
            call XMPI_Bcast(converged, 0, comm)
            !$omp end master
            !$omp barrier
          end if

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
          call MergeArrays(ONE, u, y(i), z(:,:,:,:,:,i), multi=.true.)
        end do

      end if

      ! finalization ...........................................................

      !$omp master
      deallocate(b, c, h, r, s, y, v, z)
      deallocate(bv_z, mm_inv, g, O)
      !$omp end master

    end associate

  end subroutine FGMRES_Step

  !=============================================================================

end submodule MP_FGMRES_Step
