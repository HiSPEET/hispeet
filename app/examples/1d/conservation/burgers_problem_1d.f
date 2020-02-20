!> summary:  Interface to Burgers problem
!> author:   Joerg Stiller
!> date:     2019/11/23
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Burgers_Problem_1D

  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, HALF
  use Scalar_Problem_1D

  implicit none
  private

  public :: BurgersProblem1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D Burgers problems

  type, abstract, extends(ScalarProblem1D) :: BurgersProblem1D

  contains

    procedure :: ConvectiveFlux
    procedure :: NumericalConvectiveFlux

  end type BurgersProblem1D

contains

  !-----------------------------------------------------------------------------
  !> Convective flux of the Burgers problem

  function ConvectiveFlux(problem, u) result(fc)
    class(BurgersProblem1D), intent(in) :: problem
    real(RNP), contiguous,   intent(in) :: u(0:,:,:) !< u(x,t)
    real(RNP) :: fc(0:ubound(u,1), size(u,2), problem%nc)

    fc = HALF * u * u

  end function ConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Numerical flux of the Burgers problem (local Lax-Friedrichs)

  elemental function NumericalConvectiveFlux(problem, ul, ur) result(hc)
    class(BurgersProblem1D), intent(in) :: problem
    real(RNP),               intent(in) :: ul   !< left value
    real(RNP),               intent(in) :: ur   !< right value
    real(RNP)                           :: hc

    hc = ONE/4 * (ul**2 + ur**2) + HALF * max(abs(ul),abs(ur)) * (ul - ur)

    ! just to avoid compiler warnings ;)
    if (problem % ne > 0) return

  end function NumericalConvectiveFlux

  !=============================================================================

end module Burgers_Problem_1D
