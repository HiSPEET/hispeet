!> summary:  1D DG elliptic operator: Schwarz method
!> author:   Joerg Stiller
!> date:     2022/03/17
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__1D) MP_Schwarz_Method_RX
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Schwarz method with constant or variable ν
  !>
  !> @remark
  !> One and only one of the parameters `nu_c` and `nu_v` is to be passed

  module subroutine Schwarz_Method_RX( this, bc, bv, mask, dx    &
                                     , lambda, nu_c, nu_v , f, u &
                                     , i_max, r_red, r_max, ni   )

    class(DG_EllipticOperator_1D),   intent(in)    :: this
    character,                       intent(in)    :: bc(2)     !< BC {D,N,P}
    real(RNP),                       intent(in)    :: bv(2)     !< boundary vals
    logical,                         intent(in)    :: mask(:)   !< element mask
    real(RNP),                       intent(in)    :: dx
    real(RNP),                       intent(in)    :: lambda    !< λ
    real(RNP),             optional, intent(in)    :: nu_c      !< νᵖ+νˢ
    real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:) !< νᵖ
    real(RNP), contiguous,           intent(in)    :: f(:,:)
    real(RNP), contiguous,           intent(inout) :: u(:,:)
    integer,                         intent(in)    :: i_max
    real(RNP),             optional, intent(in)    :: r_red
    real(RNP),             optional, intent(in)    :: r_max
    integer,               optional, intent(out)   :: ni

    ! internal variables .......................................................

    integer,   allocatable, save :: cfg(:)    ! subdomain boundary configuration
    real(RNP), allocatable, save :: nu_avg(:) ! element averaged ν
    real(RNP), allocatable, save :: r (:,:)   ! residual
    real(RNP), allocatable, save :: rs(:,:)   ! subdomain residuals
    real(RNP), allocatable, save :: us(:,:)   ! subdomain corrections

    real(RNP), save :: rr_term, rr_red
    logical  , save :: converged

    real(RNP) :: rr
    integer   :: e, i, ne, np, ns
    logical   :: check_convergence, periodic

    associate(eop => this%eop, schwarz => this%schwarz)

      ! initialization .........................................................

      ne = size(mask)
      np = schwarz % po + 1
      ns = schwarz % no * 2 + np

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      periodic = all(bc == 'P')

      !$omp master

      ! termination condition
      if (check_convergence) then
        if (present(r_max)) then
          rr_term = r_max ** 2
        else
          rr_term = ZERO
        end if
        if (present(r_red)) then
          rr_red = r_red ** 2
        else
          rr_red = ZERO
        end if
      else
        converged = .false.
      end if

      allocate(cfg(ne), nu_avg(ne), r(np,ne), rs(ns,ne), us(ns,ne))

      !$omp end master
      !$omp barrier

      call schwarz % ConfigureSubdomains(bc, cfg, mask)

      ! element-averaged diffusivity
      if (present(nu_c)) then
        call SetArray(nu_avg, nu_c)
      else
        !$omp do
        do e = 1, ne
          if (mask(e)) then
            nu_avg(e) = HALF * dot_product(eop%w, nu_v(:,e))
          else
            nu_avg(e) = ONE
          end if
        end do
      end if

      ! Schwarz iterations .....................................................

      do i = 1, i_max

        ! residual
        if (present(nu_c)) then
          call this % Residual(bc, bv, dx, lambda, nu_c, f, u, r, mask)
        else
          call this % Residual(bc, bv, dx, lambda, nu_v, f, u, r, mask)
        end if

        ! termination check
        if (check_convergence) then
          rr = ScalarProduct(r, r)
          !$omp master
          if (i == 1) then
            rr_term  = max(rr_term, rr * rr_red)
          end if
          converged = rr <= rr_term
        end if
        !$omp end master
        !$omp barrier
        if (converged) exit

        ! Schwarz sweep
        call schwarz % RestrictResidual(periodic, r, rs)
        call schwarz % Apply(cfg, dx, lambda, nu_avg, rs, us)
        call schwarz % MergeCorrections(periodic, us, u)

      end do

      if (present(ni)) ni = min(i, i_max)

      ! cleanup ................................................................

      !$omp master
      deallocate(cfg, nu_avg, r, rs, us)
      !$omp end master

    end associate

  end subroutine Schwarz_Method_RX

  !=============================================================================

end submodule MP_Schwarz_Method_RX
