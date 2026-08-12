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

!> summary:  1D DG elliptic operator: IPCG with Schwarz preconditioner
!> author:   Joerg Stiller
!> date:     2022/03/17
!===============================================================================

submodule(DG__Elliptic_Operator__1D) MP_SchwarzPCG_Method_RX
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Schwarz-IPCG method with constant or variable ν
  !>
  !> @remark
  !> One and only one of the parameters `nu_c` and `nu_v` is to be passed

  module subroutine SchwarzPCG_Method_RX( this, bc, bv, mask, dx   &
                                        , lambda, nu_c, nu_v, f, u &
                                        , i_max, r_red, r_max, ni  )

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

    integer  , dimension(:)  , allocatable, save :: cfg
    real(RNP), dimension(:)  , allocatable, save :: nu_avg
    real(RNP), dimension(:,:), allocatable, save :: r, p, q, s, z
    real(RNP), dimension(:,:), allocatable, save :: rs, zs
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, rr
    logical   :: check_convergence, periodic, singular
    integer   :: ne, np, ns
    integer   :: e, i, i_max_

    associate(eop => this%eop, schwarz => this%schwarz)

      ! initialization .........................................................

      ne = size(mask)
      np = schwarz % po + 1
      ns = schwarz % no * 2 + np

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      periodic = all(bc == 'P')
      singular = abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')

      ! work space
      !$omp master
      allocate(cfg(ne))
      allocate(nu_avg(ne))
      allocate(r, mold = u)
      allocate(p, mold = u)
      allocate(q, mold = u)
      allocate(s, mold = u)
      allocate(z, mold = u)
      allocate(rs(ns,ne))
      allocate(zs(ns,ne))
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

      ! initial residual .......................................................

      ! r = f - Au
      if (present(nu_c)) then
        call this % Residual(bc, bv, dx, lambda, nu_c, f, u, r, mask)
      else
        call this % Residual(bc, bv, dx, lambda, nu_v, f, u, r, mask)
      end if
      if (singular) then
        call CalibrateArray(r)
      end if

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r, r)
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

        ! Schwarz preconditioner, z = (Aˢ)⁻¹ r
        call SetArray(z, ZERO)
        call schwarz % RestrictResidual(periodic, r, rs)
        call schwarz % Apply(cfg, dx, lambda, nu_avg, rs, zs)
        call schwarz % MergeCorrections(periodic, zs, z)

        ! set/update search vector
        if (i == 1) then
          if (singular) then
            call CalibrateArray(z)
          end if
          call SetArray(p, z)                   ! p = z
        else
          call SetArray(q, r)                   ! q = r
          call MergeArrays(ONE, q, -ONE, s)     ! q = r - s
          beta = ScalarProduct(q, z) / delta
          call MergeArrays(beta, p, ONE, z)     ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r)

        ! operator application with no source and homogeneous BC
        if (present(nu_c)) then
          call this % Apply(bc, dx, lambda, nu_c, p, q, mask)
        else
          call this % Apply(bc, dx, lambda, nu_v, p, q, mask)
        end if

        ! correction
        delta = ScalarProduct(r, z)
        alpha = delta / ScalarProduct(p, q)
        call MergeArrays(ONE, u, alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          if (present(nu_c)) then
            call this % Residual(bc, bv, dx, lambda, nu_c, f, u, r, mask)
          else
            call this % Residual(bc, bv, dx, lambda, nu_v, f, u, r, mask)
          end if
          if (singular) then
            call CalibrateArray(r)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
        end if

        if (check_convergence) then
          rr = ScalarProduct(r, r)
          !$omp master
          converged = rr <= rr_term
          !$omp end master
          !$omp barrier
        end if

        if (converged .or. i == i_max_) exit

      end do

      ! finalization ...........................................................

      if (present(ni)) ni = i

      !$omp master
      deallocate(cfg, nu_avg, r, p, q, s, z, rs, zs)
      !$omp end master

    end associate

  end subroutine SchwarzPCG_Method_RX

  !=============================================================================

end submodule MP_SchwarzPCG_Method_RX
