!> summary:  Embedded interpolation operator
!> author:   Joerg Stiller
!> date:     2018/03/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Embedded_Interpolation__1D
  use Kind_Parameters, only: RNP
  use Gauss_Jacobi
  use Standard_Operators__1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation operator

  type, public :: EmbeddedInterpolation_1D

    integer :: no = -1 !< number of original points per direction
    integer :: ni = -1 !< number of interpolated points per direction
    real(RNP), allocatable :: A(:,:)  !< 1D interpolation operator

  contains

    generic :: Init_EmbeddedInterpolation_1D => Init_SX
    procedure, private :: Init_SX

  end type EmbeddedInterpolation_1D

  ! constructor interface
  interface EmbeddedInterpolation_1D
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New EmbeddedInterpolation_1D from 1D standard operators and interpolation points

type(EmbeddedInterpolation_1D) function New_SX(eop, xi) result(this)
  class(StandardOperators_1D), intent(in) :: eop    !< standard operators
  real(RNP),                  intent(in) :: xi(0:) !< points in [-1,1]

  call Init_SX(this, eop, xi)

end function New_SX

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization for given 1D standard operators and interpolation points

subroutine Init_SX(this, eop, xi)
  class(EmbeddedInterpolation_1D), intent(inout) :: this
  class(StandardOperators_1D),   intent(in)    :: eop    !< standard operators
  real(RNP),                    intent(in)    :: xi(0:) !< points in [-1,1]

  integer :: j, k, pi

  if (this % no > 0) call Delete_InterpolationOperator(this)

  associate(po => eop%po, xo => eop%x)

    pi = ubound(xi,1)

    this % no = size(xo)
    this % ni = size(xi)

    allocate(this % A(0:pi,0:po))

    select case(eop % basis)
    case('GL ') ! Gauss-Legendre
      do k = 0, po
      do j = 0, pi
        this % A(j,k) = GL_Polynomial(k, xo, xi(j))
      end do
      end do
    case('GRL') ! Gauss-Radau-Legendre
      do k = 0, po
      do j = 0, pi
        this % A(j,k) = GRL_Polynomial(k, xo, xi(j))
      end do
      end do
    case default
      do k = 0, po
      do j = 0, pi
        this % A(j,k) = GLL_Polynomial(k, xo, xi(j))
      end do
      end do
    end select

  end associate

end subroutine Init_SX

!-------------------------------------------------------------------------------
!> Delete EmbeddedInterpolation_1D object

subroutine Delete_InterpolationOperator(this)
  class(EmbeddedInterpolation_1D), intent(inout) :: this

  if (allocated(this % A)) deallocate(this % A)

end subroutine Delete_InterpolationOperator

!===============================================================================

end module Embedded_Interpolation__1D
