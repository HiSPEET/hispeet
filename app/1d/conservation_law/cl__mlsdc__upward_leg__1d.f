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

  subroutine CL_MLSDC_Upward_Leg_1D(mlsdc, t, dt, u)
    class(CL_MLSDC_1D), intent(in) :: mlsdc
    real(RNP), intent(in) :: t   !< start time
    real(RNP), intent(in) :: dt  !< thickness of time slab
    class(CL_MLSDC_Variable_1D), intent(inout) :: u   !< approximate solution

    ! internal variables .......................................................
    integer :: l, l_top

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    ! coarse to fine ...........................................................

    call mlsdc % level(1) % ApplyPredictor(dt, t, u % level(1) % val)

    do l = 1, l_top-1
      associate( Iu_c => u % level(l+1) % val &
               , u_c  => u % level(l  ) % val )

        ! predictor
        call mlsdc % level(l) % ApplyPredictor(dt, t, u_c)

        ! interpolate to finer grid
        call mlsdc % level(l) % Interpolate_CF(u_c, Iu_c, complete=.true.)

      end associate
    end do
      
    ! finalization .............................................................

    ! deallocate(g, r, v)

  end subroutine CL_MLSDC_Upward_Leg_1D

end module CL__MLSDC__Upward_Leg__1D
