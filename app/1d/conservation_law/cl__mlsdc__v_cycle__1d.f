!> summary:  MLSDC V-cycle for 1D conservation laws
!> author:   Erik Pfister, Joerg Stiller
!> date:     2023/11/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__V_Cycle__1D
  use Kind_Parameters
  use CL__MLSDC__1D
  use CL__MLSDC__Variable__1D
  implicit none
  private

  public :: CL_MLSDC_V_Cycle_1D

contains

  !-----------------------------------------------------------------------------
  !> Execution of one MLSDC V-cycle

  subroutine CL_MLSDC_V_Cycle_1D(mlsdc, t, dt, n_s1, n_s2, n_coarse, n_cycle, u)
    class(CL_MLSDC_1D), intent(in) :: mlsdc
    real(RNP), intent(in) :: t        !< start time
    real(RNP), intent(in) :: dt       !< thickness of time slab
    integer,   intent(in) :: n_s1     !< number of pre-smoothing sweeps
    integer,   intent(in) :: n_s2     !< number of post-smoothing sweeps
    integer,   intent(in) :: n_coarse !< number of sweeps for coarse solution
    integer,   intent(in) :: n_cycle  !< number of cycles tp perform
    class(CL_MLSDC_Variable_1D), intent(inout) :: u !< approximate solution

    ! internal variables .......................................................

    type(CL_MLSDC_Variable_1D), allocatable, save :: g ! FAS correction
    type(CL_MLSDC_Variable_1D), allocatable, save :: r ! residual
    type(CL_MLSDC_Variable_1D), allocatable, save :: v ! auxiliary

    integer :: l, l_top

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    ! allocate g, r, v

    ! fine to coarse ...........................................................

    do l = l_top, 2, -1
      !...
    end do

    ! coarse solution ..........................................................

    ! coarse to fine ...........................................................

    do l = 2, l_top
      !...
    end do

    ! finalization .............................................................

  end subroutine CL_MLSDC_V_Cycle_1D

  !=============================================================================

end module CL__MLSDC__V_Cycle__1D
