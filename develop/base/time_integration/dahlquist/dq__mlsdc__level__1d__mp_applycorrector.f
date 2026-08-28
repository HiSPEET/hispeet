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

!> summary:  Application of the SDC corrector
!> author:   Erik Pfister, Joerg Stiller
!> date:     2024/04/17
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_ApplyCorrector

  use Constants, only: ZERO

  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of the SDC corrector

  module subroutine ApplyCorrector(this, incremental, lambda, dt, G, u, n_sweep)
    class(DQ_MLSDC_Level_1D), intent(inout) :: this !inout
    logical,      intent(in)    :: incremental
    complex(RNP), intent(in)    :: lambda  !< lambda
    real(RNP),    intent(in)    :: dt      !< size of the time slice
    complex(RNP), intent(in)    :: G(0:,:) !< FAS RHS
    complex(RNP), intent(inout) :: u(0:,:) !< approximate solution
    integer,      intent(in)    :: n_sweep !< number of sweeps

    complex(RNP), allocatable, save :: F        (:) ! [Fᵢ]ᵏ
    complex(RNP), allocatable, save :: F_ex     (:) ! [F_exᵢ]ᵏ
    complex(RNP), allocatable, save :: F_im     (:) ! [F_imᵢ]ᵏ
    complex(RNP), allocatable, save :: F_ex_new (:) ! [F_exᵢ]ᵏ⁺¹
    complex(RNP), allocatable, save :: F_im_new (:) ! [F_imᵢ]ᵏ⁺¹

    real(RNP), allocatable :: t_sub(:), dt_sub(:)
    real(RNP) :: dt_step, t_0
    integer   :: k, m, n

    if (n_sweep < 1) return

    associate( n_sub    => this % sdc % n_sub &
             , n_time   => this % n_time      &
             , sdc      => this % sdc         )

      ! initialization .........................................................

      allocate(F        (0:n_sub))
      allocate(F_ex     (0:n_sub))
      allocate(F_im     (0:n_sub))
      allocate(F_ex_new (0:n_sub))
      allocate(F_im_new (0:n_sub))

      allocate(t_sub(0:n_sub), dt_sub(1:n_sub))

      t_0 = ZERO
      dt_step = dt / n_time

      ! time steps ............................................................

      Steps: do n = 1, n_time

        t_sub  = sdc % SubintervalPoints(t_0 + (n-1)*dt_step, dt_step)
        dt_sub = t_sub(1:n_sub) - t_sub(0:n_sub-1)

        ! prerequisites for SDC sweeps

        ! RHS for high-order quadrature
        F(0:n_sub) = lambda * u(0:n_sub,n)

        ! RHS for corrector
        do m = 0, n_sub
          call sdc % CorrectorRHS ( lambda    &
                                  , dt_sub(m) &
                                  , u(m,n)    &
                                  , F_ex(m)   &
                                  , F_im(m)   )

        end do

        F_ex_new(0) = F_ex(0)
        F_im_new(0) = F_im(0)

        ! SDC sweeps
        Sweeps: do k = 1, n_sweep

          do m = 1, n_sub
            select type (sdc)
            type is (DQ_MLSDC_Corrector_Euler)
              call sdc % CorrectorStepMLSDC( &
                  lambda, m, k, t_sub, u(:,n), F, F_ex, F_im &
                , F_ex_new, F_im_new, G(:,n), incremental    )
            type is (DQ_MLSDC_Corrector_ISD)
              call sdc % CorrectorStepMLSDC( &
                  lambda, m, k, t_sub, u(:,n), F, F_ex, F_im &
                , F_ex_new, F_im_new, G(:,n), incremental    )
            class default
              call Error('ApplyCorrector', 'unsupported MLSDC corrector')
            end select
          end do

          if (k == n_sweep) exit

          F(1:n_sub)    = lambda * u(1:n_sub,n)
          F_ex(1:n_sub) = F_ex_new(1:n_sub)
          F_im(1:n_sub) = F_im_new(1:n_sub)

        end do Sweeps

        ! set initial values for next step
        ! no assembly - doesn't work with RK
        if (n < n_time) then
          u(0,n+1) = u(n_sub,n)
        end if

      end do Steps

      ! clean-up ...............................................................

      deallocate(F, F_ex, F_im, F_ex_new, F_im_new)

    end associate

  end subroutine ApplyCorrector

 !=============================================================================

end submodule MP_ApplyCorrector
