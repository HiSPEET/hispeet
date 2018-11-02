!> summary:  Embedded interpolation operator
!> author:   Joerg Stiller
!> date:     2018/03/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Embedded interpolation operator
!===============================================================================

module Embedded_Interpolation
  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation

  type, public :: InterpolationOperator

    integer :: no = -1 !< number of original points per direction
    integer :: ni = -1 !< number of interpolated points per direction
    real(RNP), allocatable :: iop(:,:)  !< 1D interpolation operator

  contains

    generic :: New  =>  New_IOP_P, New_IOP_X
    procedure, private :: New_IOP_P
    procedure, private :: New_IOP_X

    final :: Delete_InterpolationOperator

  end type InterpolationOperator

contains

!-------------------------------------------------------------------------------
!> Initialization using the polynomial order of the interpolant

subroutine New_IOP_P(this, po, pi)
  class(InterpolationOperator), intent(inout) :: this
  integer, intent(in) :: po  !< polynomial order of the original
  integer, intent(in) :: pi  !< polynomial order of the interpolant

  real(RNP) :: xo(0:po), xi(0:pi)
  integer   :: j, k

  ! clean-up ...................................................................

  call Delete_InterpolationOperator(this)

  ! prerequisites ..............................................................

  xo = GLL_Points(po)
  xi = GLL_Points(pi)

  ! set components .............................................................

  this % no = size(xo)
  this % ni = size(xi)

  allocate(this % iop(0:pi,0:po))

  do k = 0, po
  do j = 0, pi
    this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
  end do
  end do

end subroutine New_IOP_P

!-------------------------------------------------------------------------------
!> Initialization using the interpolation points

subroutine New_IOP_X(this, po, xi)
  class(InterpolationOperator), intent(inout) :: this
  integer,   intent(in) :: po     !< polynomial order of the original
  real(RNP), intent(in) :: xi(0:) !< interpolation points in [-1,1]

  real(RNP) :: xo(0:po)
  integer   :: j, k, pi

  ! clean-up ...................................................................

  call Delete_InterpolationOperator(this)

  ! prerequisites ..............................................................

  xo = GLL_Points(po)
  pi = ubound(xi,1)

  ! set components .............................................................

  this % no = size(xo)
  this % ni = size(xi)

  allocate(this % iop(0:pi,0:po))

  do k = 0, po
  do j = 0, pi
    this % iop(j,k) = GLL_Polynomial(k, xo, xi(j))
  end do
  end do

end subroutine New_IOP_X

!-------------------------------------------------------------------------------
!> Finalization of a InterpolationOperator object

subroutine Delete_InterpolationOperator(this)
  type(InterpolationOperator), intent(inout) :: this

  if (allocated(this % iop)) deallocate(this % iop)

end subroutine Delete_InterpolationOperator

!===============================================================================

end module Embedded_Interpolation
