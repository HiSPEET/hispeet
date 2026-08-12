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

!> summary:  3D DG elliptic operator: Conjugate gradient method
!> author:   Joerg Stiller
!> date:     2022/10/03
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_CG_Method_X
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> CG method with constant or variable ν, modified for r = Au - f
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent
  !>
  !> @remark
  !> One and only one of the parameters `nu_c` and `nu_v` is to be passed

  module subroutine CG_Method_X( this, bc, lambda, nu_c, nu_v, u, f &
                               , bv, i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D),        intent(in)    :: this
    character,                            intent(in)    :: bc(:)
    real(RNP),                            intent(in)    :: lambda        !< λ
    real(RNP),                  optional, intent(in)    :: nu_c          !< νᵖ+νˢ
    real(RNP), contiguous,      optional, intent(in)    :: nu_v(:,:,:,:) !< νᵖ
    real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
    integer,                              intent(in)    :: i_max
    real(RNP),                  optional, intent(in)    :: r_red
    real(RNP),                  optional, intent(in)    :: r_max
    integer,                    optional, intent(out)   :: ni

    ! internal variables .......................................................

    real(RNP), dimension(:,:,:,:), allocatable, save :: r, p, q
    real(RNP), save :: rr_term
    logical  , save :: converged, singular, singular_loc

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, pq, rr, rr_old
    integer   :: i, na

    character(len=:), allocatable :: prefix
    logical :: logging

    ! skip empty partition
    if (this % sem % mesh % part < 0) then
      if (present(ni)) ni = -1
      return
    end if

    ! start logging ............................................................

    if (log_level > 0) then
      logging = this%sem%mesh%proc == 0 .or. log_level > 1
      prefix  = LoggingPrefix('CG_Method_X', this%sem%mesh%proc)
    else
      logging = .false.
    end if

    if (logging) then
      print '(2A)', prefix, 'start'
    end if

    ! ..........................................................................

    associate(mesh => this % sem % mesh)

      ! initialization .........................................................

      na = mesh % n_elem_active

      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !$omp master

      ! work space
      allocate(r, mold = u)
      allocate(p, mold = u)
      allocate(q, mold = u)

      singular = abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')
      if (singular .and. .not. mesh%is_root)  then
        singular_loc = mesh%n_elem_frozen == 0
        call XMPI_Allreduce(singular_loc, singular, MPI_LAND, mesh%comm_parts)
      end if

      !$omp end master
      !$omp barrier
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      ! initial residual .......................................................

      if (present(nu_c)) then
        call this % Residual(bc, lambda, nu_c, f, bv, u, r)
      else
        call this % Residual(bc, lambda, nu_v, f, bv, u, r)
      end if

      if (singular) then
        call CalibrateArray(r(:,:,:,:na), mesh%comm_parts)
      end if
      call SetArray(p, r)

      rr = ScalarProduct(r(:,:,:,:na), r(:,:,:,:na), mesh%comm_parts)

      !$omp single
      rr_term = 0
      if (present(r_red)) rr_term = rr * max(ZERO, r_red)**2
      if (present(r_max)) rr_term = max(rr_term, max(ZERO, r_max)**2)
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

        if (log_level_inner_iteration > 1 .and. mesh%part == 0) then
          print '(A,T25,A,I5,A,ES12.5)', &
                '#Elliptic:CG','>>>  i  =',i-1,',  |r| =', sqrt(rr)
        end if
        !$omp end master
        !$omp barrier

        if (converged) exit

        rr_old = rr

        ! next iteration . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! operator application with no source and homogeneous BC
        if (present(nu_c)) then
          call this % Apply(bc, lambda, nu_c, u=p, r=q)
        else
          call this % Apply(bc, lambda, nu_v, u=p, r=q)
        end if

        ! correction
        pq = ScalarProduct(p(:,:,:,:na), q(:,:,:,:na), mesh%comm_parts)
        pq = sign(max(abs(pq),eps), pq)
        alpha = rr_old / pq
        call MergeArrays(ONE, u(:,:,:,:na), alpha, p(:,:,:,:na))

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          if (present(nu_c)) then
            call this % Residual(bc, lambda, nu_c, f, bv, u, r)
          else
            call this % Residual(bc, lambda, nu_v, f, bv, u, r)
          end if
          if (singular) then
            call CalibrateArray(r(:,:,:,:na), mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r(:,:,:,:na), -alpha, q(:,:,:,:na))
        end if

        rr = ScalarProduct(r(:,:,:,:na), r(:,:,:,:na), mesh%comm_parts)

        call  MergeArrays(rr/rr_old, p(:,:,:,:na), ONE, r(:,:,:,:na))

      end do

      if (present(ni)) ni = min(i,i_max)

      !$omp master
      if (log_level_inner_iteration > 0 .and. mesh%part == 0) then
        print '(A,T25,A,I5,A,ES12.5)', &
              '#Elliptic:CG','>>>  ni =',min(i,i_max),',  |r| =', sqrt(rr)
      end if
      !$omp end master

      !$omp master
      deallocate(r, p, q)
      !$omp end master

    end associate

    ! exit logging .............................................................

    if (logging) then
      print '(2A)', prefix, 'exit'
    end if

  end subroutine CG_Method_X

  !=============================================================================

end submodule MP_CG_Method_X
