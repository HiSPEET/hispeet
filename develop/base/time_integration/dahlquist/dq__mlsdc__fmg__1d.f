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

!> summary:  MLSDC FMG Start for Dahlquist problem
!> author:   Erik Pfister
!> date:     2024/04/28
!===============================================================================

module DQ__MLSDC__FMG__1D
  use Kind_Parameters
  use Constants
  use DQ__MLSDC__1D
  use DQ__MLSDC__Variable__1D
  use DQ__MLSDC__V_Cycle__1D
  implicit none
  private

  public :: DQ_MLSDC_FMG_1D

contains

  !-----------------------------------------------------------------------------
  !> Execution of a FMG start
  !> Variable names:
  !> I before a variable name indicates an interpolated variable
  !> Index f indicates the finer grid of the iteration
  !> Index c indicates the coarser grid of the iteration

  subroutine DQ_MLSDC_FMG_1D( &
      mlsdc, incremental, lambda, dt, n_s, n_cycle, n_coarse, u)
    type(DQ_MLSDC_1D), intent(inout) :: mlsdc
    logical,      intent(in)  :: incremental
    complex(RNP), intent(in)  :: lambda    !< lambda
    real(RNP),    intent(in)  :: dt        !< thickness of time slab
    integer,      intent(in)  :: n_s       !< number of smoothing sweeps
    integer,      intent(in)  :: n_cycle   !< cycles per added level
    integer,      intent(in)  :: n_coarse  !< number of coarse sweeps
    class(DQ_MLSDC_Variable_1D),   intent(inout) :: u !< approximate solution

    ! internal variables .......................................................
    integer :: l, l_top
    type(DQ_MLSDC_Variable_1D), allocatable, save :: G ! FAS RHS

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    G = DQ_MLSDC_Variable_1D(mlsdc)

    ! coarse to fine ...........................................................

    associate(Iu_c => u % level(2) % val &
             , u_c => u % level(1) % val &
             , G_c => G % level(1) % val )

      ! Maybe call cascade
      call mlsdc % level(1) % ApplyPredictor(lambda, dt, u_c)
      G_c = (0.0, 0.0)
      call mlsdc % level(1) % ApplyCorrector( &
          incremental, lambda, dt, G_c, u_c, n_s)
      call mlsdc % level(1) % Interpolate_CF(u_c, Iu_c, complete=.true.)

    end associate

    do l = 2, l_top-1
      associate( Iu_c => u % level(l+1) % val &
               , u_c  => u % level(l  ) % val &
               , G_c  => G % level(l  ) % val )

        call DQ_MLSDC_V_Cycle_1D( &
            mlsdc, incremental, lambda, dt, 1, 1, n_coarse, n_cycle, u, l)
        ! interpolate to finer grid
        call mlsdc % level(l) % Interpolate_CF(u_c, Iu_c, complete=.true.)

      end associate
    end do

  ! finalization .............................................................

  end subroutine DQ_MLSDC_FMG_1D

end module DQ__MLSDC__FMG__1D
