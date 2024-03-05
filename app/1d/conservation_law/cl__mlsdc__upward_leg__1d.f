!> summary:  MLSDC Upward Leg for 1D conservation laws
!> author:   Erik Pfister
!> date:     2023/02/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__Upward_Leg__1D
  use Kind_Parameters
  use Constants
  use CL__MLSDC__1D
  use CL__MLSDC__Variable__1D
  implicit none
  private

  public :: CL_MLSDC_Upward_Leg_1D

contains

  !-----------------------------------------------------------------------------
  !> Execution of one MLSDC Upward Leg
  !> Variable names:
  !> I bevore a variable name indicates an interpolated variable
  !> Index f indicates the finer grid of the iteration
  !> Index c indicates the coarser grid of the iteration

  subroutine CL_MLSDC_Upward_Leg_1D(mlsdc, t, dt, n_s, u)
    class(CL_MLSDC_1D), intent(in) :: mlsdc
    real(RNP), intent(in) :: t   !< start time
    real(RNP), intent(in) :: dt  !< thickness of time slab
    integer,   intent(in) :: n_s !< number of smoothing sweeps 
    class(CL_MLSDC_Variable_1D), intent(inout) :: u   !< approximate solution

    ! internal variables .......................................................

    type(CL_MLSDC_Variable_1D), allocatable, save :: g ! FAS RHS
    integer :: l, l_top

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    g = CL_MLSDC_Variable_1D(mlsdc)

    ! coarse to fine ...........................................................

    do l = 2, l_top
      associate( Iu_c => u % level(l  ) % val &
               , u_c  => u % level(l-1) % val &
               , g_c  => g % level(l-1) % val )

        ! smoothing step
        g_c = ZERO
        call mlsdc % level(l-1) % ApplyPredictor(dt, t, u_c)
        !call mlsdc % level(l-1) % ApplyCorrector(dt, t, g_c, u_c, n_s) 

        ! interpolate to finer grid
        call mlsdc % level(l-1) % Interpolate_CF(u_c, Iu_c, complete=.true.)

      end associate
    end do
      
    ! finalization .............................................................

    ! deallocate(g, r, v)

  end subroutine CL_MLSDC_Upward_Leg_1D

end module CL__MLSDC__Upward_Leg__1D
