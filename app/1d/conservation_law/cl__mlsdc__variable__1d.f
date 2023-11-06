!> summary:  MultiLevel SDC variable for 1D conservation laws
!> author:   Joerg Stiller
!> date:     2023/10/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__Variable__1D
  use Kind_Parameters
  use Constants, only: ZERO
  use CL__MLSDC__1D

  implicit none
  private

  public :: CL_MLSDC_Variable_1D

  !-----------------------------------------------------------------------------
  !> 1D multilevel space-time variable

  type CL_MLSDC_Variable_1D
    type(LevelVariable), allocatable :: level(:)
  end type CL_MLSDC_Variable_1D

  ! constructor
  interface CL_MLSDC_Variable_1D
    procedure New_CL_MLSDC_Variable_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Structure representing one level of a multilevel space-time variable

  type LevelVariable
    real(RNP), allocatable :: val(:,:,:,:,:) !< variable values
  end type LevelVariable

contains

  !-----------------------------------------------------------------------------
  !> Construction of a new MLSDC variable
  !>
  !> The dimensions and, if not specified, the number of components are adopted
  !> from the given MLSDC data structure.

  function New_CL_MLSDC_Variable_1D(mlsdc, nc) result(var)
    class(CL_MLSDC_1D), intent(in) :: mlsdc !< MLSDC data structure
    integer,  optional, intent(in) :: nc    !< number of components
    type(CL_MLSDC_Variable_1D) :: var

    integer :: l, nc_, nl

    if (present(nc)) then
      nc_ = nc
    else
      nc_ = mlsdc % level(1) % cl_problem % nc
    end if

    nl = size(mlsdc % level)

    allocate(var % level(nl))

    do l = 1, nl
      associate( ps => mlsdc % level(l) % p_space &
               , ns => mlsdc % level(l) % n_space &
               , pt => mlsdc % level(l) % p_time  &
               , nt => mlsdc % level(l) % n_time  )
        allocate(var % level(l) % val(0:ps,1:ns,1:nc_,0:pt,1:nt), source = ZERO)
      end associate
    end do

  end function New_CL_MLSDC_Variable_1D

  !=============================================================================

end module CL__MLSDC__Variable__1D
