!> summary:  Legendre projection operator
!> author:   Joerg Stiller
!> date:     2019/02/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Legendre projection operator
!===============================================================================

module Legendre_Projection
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Gauss_Jacobi
  use Standard_Operators_1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Legendre projection operator

  type, public :: LegendreProjectionOperator

    integer :: no = -1 !< original DOF
    integer :: np = -1 !< projected DOF
    real(RNP), allocatable :: pop(:,:)  !< 1D projection operator

  contains

    procedure :: Init_LegendreProjectionOperator

  end type LegendreProjectionOperator

  ! constructor interface
  interface LegendreProjectionOperator
    module procedure New_LPO
  end interface

contains

!===============================================================================
! Constructors

type(LegendreProjectionOperator) New_LPO(eop, pp) result(this)
  class(StandardOperators1D), intent(in) :: eop  !< standard operators
  integer,                    intent(in) :: pp   !< order to project to

  call Init_LegendreProjectionOperator(this, eop, pp)

end function New_LPO

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization of projection operator

subroutine Init_LegendreProjectionOperator(this, eop, pp)
  class(LegendreProjectionOperator), intent(inout) :: this
  class(StandardOperators1D), intent(in) :: eop  !< standard operators
  integer,                    intent(in) :: pp   !< order to project to

  real(RNP), allocatable :: VL(:,:), VL_inv(:,:)
  integer :: i, j

  associate(po => eop%po)

    allocate(VL(0:po, 0:po))
    call eop % GetLegendreVDM(VL)

    allocate(VL_inv, mold=VL)
    call eop % GetInverseLegendreVDM(VL_inv)

    this % no = po + 1
    this % np = pp + 1

    allocate(this % pop(0:pp,0:po), source=ZERO)

    do i = 0, min(po,pp)
    do j = 0, po
      this % pop(i,j) = sum(VL(i,:) * VL_inv(:,j))
    end do
    end do

  end associate

end subroutine Init_LegendreProjectionOperator

!===============================================================================

end module Legendre_Projection