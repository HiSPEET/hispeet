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

!> summary:  MLSDC V-cycle for the Dahlquist problem
!> author:   Erik Pfister, Joerg Stiller
!> date:     2024/04/18
!===============================================================================

module DQ__MLSDC__V_Cycle__1D
  use Kind_Parameters
  use DQ__MLSDC__1D
  use DQ__MLSDC__Variable__1D
  implicit none
  private

  public :: DQ_MLSDC_V_Cycle_1D

contains

  !-----------------------------------------------------------------------------
  !> Execution of one MLSDC V-cycle
  !> Variable names:
  !> P before a variable name indicates a projected variable
  !> R before a variable name indicates a restricted variable
  !> I before a variable name indicates an interpolated variable
  !> Index f indicates the finer grid of the iteration
  !> Index c indicates the coarser grid of the iteration

  subroutine DQ_MLSDC_V_Cycle_1D( mlsdc, incremental, lambda, dt, n_s1, n_s2 &
                                , n_coarse, n_cycle, u, fmgl                )
    type(DQ_MLSDC_1D), intent(inout) :: mlsdc
    logical,      intent(in) :: incremental
    complex(RNP), intent(in) :: lambda   !< lambda
    real(RNP),    intent(in) :: dt       !< thickness of time slab
    integer,      intent(in) :: n_s1     !< number of pre-smoothing sweeps
    integer,      intent(in) :: n_s2     !< number of post-smoothing sweeps
    integer,      intent(in) :: n_coarse !< number of sweeps for coarse solution
    integer,      intent(in) :: n_cycle  !< number of cycles to perform
    class(DQ_MLSDC_Variable_1D), intent(inout) :: u    !< approximate solution
    integer, optional,           intent(in)    :: fmgl !< FMG top level

    ! internal variables .......................................................

    type(DQ_MLSDC_Variable_1D), allocatable, save :: G   ! FAS RHS
    type(DQ_MLSDC_Variable_1D), allocatable, save :: r   ! residual
    type(DQ_MLSDC_Variable_1D), allocatable, save :: v   ! auxiliary

    integer :: c, l, l_top

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    G   = DQ_MLSDC_Variable_1D(mlsdc)
    r   = DQ_MLSDC_Variable_1D(mlsdc)
    v   = DQ_MLSDC_Variable_1D(mlsdc)

    if (present(fmgl)) then
      l_top = fmgl
    end if

    do c = 1, n_cycle

      ! fine to coarse .........................................................

      ! set FAS RHS to zero for top level
      G % level(l_top) % val = (0.0, 0.0)

      do l = l_top, 2, -1
        associate( r_f  => r % level(l  ) % val &
                 , Rr_f => r % level(l-1) % val &
                 , G_f  => G % level(l  ) % val &
                 , G_c  => G % level(l-1) % val &
                 , Pu_f => v % level(l-1) % val &
                 , u_f  => u % level(l  ) % val )
    
          ! pre-smoothing
          call mlsdc % level(l) % ApplyCorrector( &
              incremental, lambda, dt, G_f, u_f, n_s1)

          ! get residual on fine grid
          ! r_f = g_f - L(u_f)
          call mlsdc % level(l) % GetResidual( &
              incremental, lambda, dt, G_f, u_f, r_f)
          ! restrict residual to coarse level
          call mlsdc % level(l-1) % Restrict_FC(r_f, Rr_f)
          
          ! restrict fine solution
          call mlsdc % level(l) % Project_FC(u_f, Pu_f)
          ! where regular refinement condition
          ! g_c = L(Pu_f)
          call mlsdc % level(l-1) % ApplyOperator( &
              incremental, lambda, dt, Pu_f, G_c)
          
          ! compute FAS RHS g (= f^ for Brandt)
          ! Gl. 8.5b in Multigrid techniques (1984, Achi Brandt)
          ! G_c = L(Pu_f)) + R(G_f - L(u_f))
          G_c = G_c + Rr_f

        end associate
      end do

      ! coarse solution .......................................................

      associate( G_c => G % level(1) % val &
               , u_c => u % level(1) % val )

        call mlsdc % level(1) % ApplyCorrector( &
            incremental, lambda, dt, G_c, u_c, n_coarse)

      end associate

      ! coarse to fine ........................................................

      do l = 2, l_top
        associate( G_f   => G % level(l  ) % val &
                 , v_cr  => r % level(l-1) % val &
                 , Iv_cr => r % level(l  ) % val &
                 , Pu_f  => v % level(l-1) % val &
                 , u_f   => u % level(l  ) % val &
                 , u_c   => u % level(l-1) % val )

          ! where regular refinement condition
          ! calculate correction v_cr
          v_cr = u_c - Pu_f

          ! interpolate to finer grid
          call mlsdc % level(l-1) % Interpolate_CF(v_cr, Iv_cr, complete=.true.)

          ! u_NEW
          u_f = u_f + Iv_cr

          if (l /= l_top .or. (l == l_top .and. c == n_cycle)) then
            call mlsdc % level(l) % ApplyCorrector( &
                incremental, lambda, dt, G_f, u_f, n_s2)
          end if

        end associate
      end do
      
    ! finalization .............................................................

    end do

    deallocate(G, r, v)

  end subroutine DQ_MLSDC_V_Cycle_1D

end module DQ__MLSDC__V_Cycle__1D
