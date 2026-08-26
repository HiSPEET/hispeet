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

!> summary:  Application of the SDC predictor
!> author:   Erik Pfister, Joerg Stiller
!> date:     2024/03/25
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_ApplyPredictor

  use Constants, only: ZERO
      
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of the SDC predictor

  module subroutine ApplyPredictor(this, lambda, dt, u)
    class(DQ_MLSDC_Level_1D), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda  !< lambda
    real(RNP)   , intent(in)    :: dt      !< size of the time slice
    complex(RNP), intent(inout) :: u(0:,:) !< approximate solution

    real(RNP), allocatable :: t_sub(:), dt_sub(:)
    real(RNP) :: t_0, dt_step
    integer   :: m, n

    associate( sdc    => this % sdc         &
             , n_sub  => this % sdc % n_sub &
             , n_time => this % n_time      )

      ! initialization .........................................................

      allocate(t_sub(0:n_sub), dt_sub(1:n_sub))

      t_0 = ZERO
      dt_step = dt / n_time

      ! time steps ............................................................

      do n = 1, n_time
        t_sub  = sdc % SubintervalPoints(t_0 + (n-1)*dt_step, dt_step)
        dt_sub = t_sub(1:n_sub) - t_sub(0:n_sub-1)

        ! sweep thru subintervals
        do m = 1, n_sub
          u(m,n) = u(m-1,n)
          call sdc % predictor % TimeStep( lambda = lambda    &
                                         , dt     = dt_sub(m) &
                                         , u      = u(m  ,n)  )
        end do

        ! set initial values for next step
        ! for now only one time step
        if (n < n_time) then
          u(0,n+1) = u(n_sub,n)
        end if

      end do

    end associate

  end subroutine ApplyPredictor

  !=============================================================================

end submodule MP_ApplyPredictor
