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
  implicit none

contains

  !-----------------------------------------------------------------------------
  !>

  module subroutine FGMRES_Step( this                                       &
                               , tau, t, v_0, F_c, F_d, Q , bv_u, mu, nu, u &
                               , n_krylov, i_pre_p, i_pre_v , r_red, r_max  )

    ! arguments ................................................................

    class(INS_TimeIntegrator_3D), intent(in) :: this

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
    class(BoundaryVariable_3D), intent(inout) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  [ v₁, v₂, v₃, - , -  ]
    !!   - Γᴼ :  [ - , - , - , p , ∆p ]
    real(RNP), contiguous, optional, intent(inout) :: mu(:,:,:,:)
    !< μ, kinematic bulk viscosity
    real(RNP), contiguous, optional, intent(inout) :: nu(:,:,:,:)
    !< ν, kinematic shear viscosity
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:)
    !< u = [v, p], velocity and pressure at final time u
    integer, intent(in) :: n_krylov
    !< max number of GMRES iterations (Krlyov subspace dimension)
    integer, intent(in) :: i_pre_p
    !< number of pressure iterations for preconditioning
    integer, intent(in) :: i_pre_v
    !< number of velocity iterations for preconditioning
    real(RNP), optional, intent(in) :: r_red
    !< min reduction of Euclidian residual norm
    real(RNP), optional, intent(in) :: r_max
    !< max Euclidian residual norm allowed

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
    real(RNP), allocatable, save :: mm(:,:,:,:)
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

      ! dimensions
      po = ins_op % eop_v % po
      np = po + 1
      nb = mesh % n_bound
      ne = mesh % n_elem
      ni = n_krylov

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

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

      allocate(mm(np,np,np,ne))
      allocate(g(np,np,np,ne,4))
      allocate(O(np,np,np,ne,4), source=ZERO)

      !$omp end master
      !$omp barrier

      call ins_op % sem_v % Get_DG_DiagonalMassMatrix(mm)

      !$omp do
      do i = 1, ne
        g(:,:,:,i,1:4) = Q  (:,:,:,i,1:4)              ! RHS
        g(:,:,:,i,1:3) = g  (:,:,:,i,1:3)           &  ! for
                       + v_0(:,:,:,i,1:3) * ONE/tau &  ! residual
                       + F_c(:,:,:,i,1:3)              ! computation
      end do

      !$omp master

      ! initial approximation ..................................................

      call this % ProjectionStep( tau, t, v_0, F_c, F_d, Q, bv_u &
                                , mu, nu, u, i_pre_p, i_pre_v    )

      ! initial residual, v₁ = f - Au ..........................................

      associate(v1 => v(:,:,:,:,:,1))

        call this % GetStokesResidual(tau, mm, g, bv_u, mu, nu, u, v1)
        beta = sqrt(ScalarProduct(v1, v1, comm))
        b(1) = beta

        ! convergence check ....................................................

        if (check_convergence) then

          !$omp master
          if (k == 1) then
            ! set terminal condition
            r_term = huge(r_term)
            if (present(r_max)) then
              if (r_max > 0) r_term = r_max
            end if
            if (present(r_red)) then
              if (r_red > 0) r_term = max(r_term, beta * r_red)
            end if
          end if
          converged = beta <= r_term
          call XMPI_Bcast(converged, 0, comm)
        else
          converged = .false.
        end if
        !$omp end master
        !$omp barrier

        if (converged) then
          ni = 0
        else
          call ScaleArray(v1, 1/max(beta,eps), multi=.true.)
        end if

      end associate

      KRYLOV_ITERATION: do j = 1, ni
        associate( vj => v(:,:,:,:,:,j)   &
                 , w  => v(:,:,:,:,:,j+1) &
                 , zj => z(:,:,:,:,:,j)   )

          ! preconditioning ....................................................

          ! projection with homogeneous BC and frozen viscosity: z(j) = K⁻¹v(j)
          call this % ProjectionStep( tau, t, O, O, O, vj, bv_z, mu, nu     &
                                    , zj, i_pre_p, i_pre_v, freeze = .true. )

          ! application of homogeneous operator: w = A v(j)
          call this % GetStokesResidual(tau, mm, O, bv_z, mu, nu, vj, w)
          call ScaleArray(w, -ONE, multi=.true.)

          ! computation of new Krylov vector ...................................

          ! orthogonalization against old Krylov vectors
          do i = 1, j
            h(i,j) = ScalarProduct(v(:,:,:,:,:,i), w, comm)
            call MergeArrays(ONE, w, -h(i,j), v(:,:,:,:,:,i), multi=.true.)
          end do

          ! normalization, v(j+1) = w / ‖w‖
          h(j+1,j) = sqrt(ScalarProduct(w, w, comm))
          call ScaleArray(w, 1//max(h(j+1,j),eps), multi=.true.)

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

          ! convergence test ....................................................

          ! residual norm if Krylov iterations were exited now
          rho = abs(b(j+1))

          if (mesh % part == 0) then
            !$omp master
            if (log_level_inner_iteration > 0) then
              write(*,'(2X,A,I4,A,ES12.5)') &
                'FGMRES iteration j =',j,': rho  =', rho
            end if
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
      end do KRYLOV_ITERATION

      ! improved solution ......................................................

      if (ni > 0) then

        j = min(j, ni)

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
      deallocate(bv_z, mm, g, O)
      !$omp end master

    end associate

  end subroutine FGMRES_Step

  !=============================================================================

end submodule MP_FGMRES_Step
