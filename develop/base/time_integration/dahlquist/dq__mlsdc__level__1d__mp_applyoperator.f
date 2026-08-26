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

!> summary:  Application of the operator for the time slice
!> author:   Erik Pfister, Jörg Stiller
!> date:     2024/04/18
!===============================================================================

submodule(DQ__MLSDC__Level__1D) MP_ApplyOperator
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of Operator v = L(u)

  module subroutine ApplyOperator(this, incremental, lambda, dt, u, r)
    class(DQ_MLSDC_Level_1D), intent(in) :: this
    logical     , intent(in)  :: incremental
    complex(RNP), intent(in)  :: lambda  !< lambda
    real(RNP)   , intent(in)  :: dt      !< size of the time slice
    complex(RNP), intent(in)  :: u(0:,:) !< approximate solution
    complex(RNP), intent(out) :: r(0:,:) !< u after operator applied
    
    real(RNP) :: dt_step
    integer   :: mt, nt, m, n, i

    mt = ubound(u, 1)
    nt = ubound(u, 2)

    dt_step = dt / nt

    r(0,:) = (0.0, 0.0)
    do n = 1, nt
    do m = 1, mt
      if (incremental) then
        r(m,n) = u(m,n) - u(m-1,n)
      else
        r(m,n) = u(m,n) - u(0,n)
      end if
    end do
    end do

    do n = 1, nt
    do i = 0, mt
      associate( w_nn => this % sdc % w_nn &
               , w_0n => this % sdc % w_0n )

        do m = 1, mt
          if (incremental) then
            r(m,n) = r(m,n) - dt_step * w_nn(m,i) * lambda * u(i,n)
          else
            r(m,n) = r(m,n) - dt_step * w_0n(m,i) * lambda * u(i,n)
          end if
        end do
      end associate
    end do
    end do

  end subroutine ApplyOperator

  !=============================================================================

end submodule MP_ApplyOperator
