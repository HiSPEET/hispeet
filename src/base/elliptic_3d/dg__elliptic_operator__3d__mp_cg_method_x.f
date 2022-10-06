!> summary:  3D DG elliptic operator: Conjugate gradient method
!> author:   Joerg Stiller
!> date:     2022/10/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_CG_Method_X
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> CG method with constant or variable ν, modified for r = Au - f
  !>
  !> @remark
  !> One and only one of the parameters `nu_c` and `nu_v` is to be passed

  module subroutine CG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                               , i_max, r_red, r_max, ni            )

    class(DG_EllipticOperator_3D),             intent(in)    :: this
    real(RNP),                                 intent(in)    :: lambda
    real(RNP),                       optional, intent(in)    :: nu_c
    real(RNP), contiguous,           optional, intent(in)    :: nu_v(:,:,:,:)
    real(RNP), contiguous,                     intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv(:)
    integer,                                   intent(in)    :: i_max
    real(RNP),                       optional, intent(in)    :: r_red
    real(RNP),                       optional, intent(in)    :: r_max
    integer,                         optional, intent(out)   :: ni

    ! internal variables .......................................................

    real(RNP), dimension(:,:,:,:), allocatable, save :: r, p, q
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, pq, rr, rr_old
    logical   :: singular
    integer   :: i

    ! skip empty partition
    if (this % sem % mesh % part < 0) return

    ! initialization ...........................................................

    associate(mesh => this % sem % mesh)

      ! work space
      !$omp master
      allocate(r, mold = u)
      allocate(p, mold = u)
      allocate(q, mold = u)
      !$omp end master
      !$omp barrier

      singular = abs(lambda) < epsilon(ONE) .and. all(this%bc /= 'D')

      ! initial residual .......................................................

      if (present(nu_c)) then
        call this % Residual(lambda, nu_c, f, bv, u, r)
      else
        call this % Residual(lambda, nu_v, f, bv, u, r)
      end if

      if (singular) then
        call CalibrateArray(r, mesh%comm_parts)
      end if
      call SetArray(p, r)

      rr = ScalarProduct(r, r, mesh%comm_parts)

      !$omp single
      if (present(r_red)) then
        rr_term  = max(ZERO, sqrt(rr) * r_red)**2
        if (present(r_max)) then
          rr_term = max(rr_term, max(ZERO, r_max)**2)
        end if
      else
        rr_term = 0
      end if
      !$omp end single

      rr_old = 0

      ! iteration ..............................................................

      do i = 1, i_max

        ! termination check  . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! MPI master decides about termination
        !$omp master
        if (mesh%part == 0) then
          converged = rr <= rr_term
        end if
        call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
        !$omp end master
        !$omp barrier

        if (converged) exit

        rr_old = rr

        ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! operator application with no source and homogeneous BC
        if (present(nu_c)) then
          call this % Apply(lambda, nu_c, p, q)
        else
          call this % Apply(lambda, nu_v, p, q)
        end if

        ! correction
        pq = ScalarProduct(p, q, mesh%comm_parts)
        pq = sign(max(abs(pq),eps), pq)
        alpha = rr_old / pq
        call MergeArrays(ONE, u, alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          if (present(nu_c)) then
            call this % Residual(lambda, nu_c, f, bv, u, r)
          else
            call this % Residual(lambda, nu_v, f, bv, u, r)
          end if
          if (singular) then
            call CalibrateArray(r, mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
        end if

        rr = ScalarProduct(r, r, mesh%comm_parts)

        call  MergeArrays(rr/rr_old, p, ONE, r)

      end do

      if (present(ni)) ni = i - 1

      !$omp master
      deallocate(r, p, q)
      !$omp end master

    end associate

  end subroutine CG_Method_X

  !=============================================================================

end submodule MP_CG_Method_X
