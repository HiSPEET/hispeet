!> summary:  Incompressible Navier-Stokes DG-SEM diffusion solver
!> author:   Joerg Stiller
!> date:     2022/09/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - variable viscosity
!>   - standby option?
!>   - consistent handling of empty partitions
!===============================================================================

submodule(INS__Operator__3D) MP_DiffusionSolver
! use Array_Assignments
! use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Diffusion solver --  so far with constant viscosity

  module subroutine DiffusionSolver(this, tau, f, bv_v, v, i_max, r_red, r_max, ni)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< sources, f(np,np,np,ne,3)

    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
    !< velocity boundary values, bv_v(nb)

    real(RNP), contiguous, intent(inout) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    integer,               intent(in)    :: i_max  !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max  !< max admissible residual
    integer,     optional, intent(out)   :: ni     !< executed num iterations

    ! internal variables .......................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: r, p, q, s, z
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, pq, rr
    logical   :: check_convergence
    integer   :: d, i, i_max_

    ! skip empty partition
    if (this % mesh % part < 0) return

    associate( mesh      => this % mesh      &
             , schwarz_v => this % schwarz_v )

      ! initialization .........................................................

      !$omp master
      allocate(r, mold = f)
      allocate(p, mold = f)
      allocate(q, mold = f)
      allocate(s, mold = f)
      allocate(z, mold = f)
      !$omp end master
      !$omp barrier

      check_convergence = present(r_red) .or. present(r_max)

      ! update Schwarz operators
      !!! TBD               !!!
      !!! lambda = 1 / tau  !!!

      ! initial residual
      call this % GetDiffusionResidual(tau, v, r, f, bv_v)

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r, r, mesh%comm_parts)
        !$omp master
        if (present(r_red)) then
          rr_term  = max(ZERO, sqrt(rr) * r_red)**2
        else
          rr_term = 0
        end if
        if (present(r_max)) then
          rr_term = max(rr_term, max(ZERO, r_max)**2)
        end if
        converged = rr <= rr_term
        call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
        !$omp end master
        !$omp barrier
      else
        !$omp master
        converged = .false.
        !$omp end master
        !$omp barrier
      end if

      if (converged) then
        i_max_ = 0
        i      = 0
      else
        i_max_ = i_max
      end if

      ! iteration ...............................................................

      do i = 1, i_max_

        ! Schwarz preconditioner, based on component-wise Helmholtz problems
        call SetArray(z, ZERO, multi = .true.)
        do d = 1, 3
          call schwarz_v(d) % Schwarz_Method( z(:,:,:,:,d)         &
                                            , r(:,:,:,:,d)         &
                                            , i_max   = 1          &
                                            , standby = i < i_max_ )
        end do

        ! set/update search vector
        if (i == 1) then
          call SetArray(p, z, multi = .true.)                 ! p = z
        else
          call SetArray(q, r, multi = .true.)                 ! q = r
          call MergeArrays(ONE, q, -ONE, s, multi = .true.)   ! q = r - s
          beta = ScalarProduct(q, z, mesh%comm_parts) / delta
          call MergeArrays(beta, p, ONE, z, multi = .true.)   ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r, multi = .true.)

        ! correction
        call this % GetDiffusionResidual(tau, p, q)  ! q = -Ap
        delta = ScalarProduct(r, z, mesh%comm_parts)
        pq    = ScalarProduct(p, q, mesh%comm_parts)
        alpha = delta / pq
        call MergeArrays(ONE, v, -alpha, p, multi = .true.)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          call this % GetDiffusionResidual(tau, v, r, f, bv_v)
        else
          call MergeArrays(ONE, r, -alpha, q, multi = .true.)
        end if

        if (check_convergence) then
          rr = ScalarProduct(r, r, mesh%comm_parts)
          !$omp master
          converged = rr <= rr_term
          call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
          !$omp end master
          !$omp barrier
        end if

        if (converged .or. i == i_max_) exit

      end do

      ! finalization ...........................................................

      if (present(ni)) ni = i

      !$omp master
      deallocate(p, q, r, s, z)
      !$omp end master

    end associate

  end subroutine DiffusionSolver

  !=============================================================================

end submodule MP_DiffusionSolver
