!> summary:  Projection operator
!> author:   Joerg Stiller
!> date:     2019/02/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Projection operator
!===============================================================================

module Projection_Operator
  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Gauss_Jacobi
  use Standard_Operators_1D
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Embedded interpolation operator

  type, public :: ProjectionOperator

    integer :: nq = -1 !< number of quadrature points per direction
    integer :: np = -1 !< number of projection points per direction
    real(RNP), allocatable :: pop(:,:) !< 1D projection operator

  contains

    generic :: Init_ProjectionOperator => Init_SX
    procedure, private :: Init_SX

  end type ProjectionOperator

  ! constructor interface
  interface ProjectionOperator
    module procedure New_SX
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> New ProjectionOperator from 1D standard operators and interpolation points

type(ProjectionOperator) function New_SX(eop, xq, wq, dx) result(this)
  class(StandardOperators1D), intent(in) :: eop   !< standard operators
  real(RNP),                  intent(in) :: xq(:) !< quadrature points
  real(RNP),                  intent(in) :: wq(:) !< quadrature weights
  real(RNP),                  intent(in) :: dx    !< element length

  call Init_SX(this, eop, xq, wq, dx)

end function New_SX

!===============================================================================
! Type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine Init_SX(this, eop, xq, wq, dx)
  class(ProjectionOperator),  intent(inout) :: this
  class(StandardOperators1D), intent(in)    :: eop   !< standard operators
  real(RNP),                  intent(in)    :: xq(:) !< quadrature points
  real(RNP),                  intent(in)    :: wq(:) !< quadrature weights
  real(RNP),                  intent(in)    :: dx    !< element length

  integer :: j, k, po, nq

  if (this % nq > 0) call Delete_ProjectionOperator(this)

  associate(po => eop%po, xo => eop%x)

    nq = size(xq)

    this % nq = nq
    this % np = po + 1

    allocate(this % pop(0:po,1:nq))

    select case(eop % basis)
    case('GL ') ! Gauss-Legendre
      do k = 1, nq
      do j = 0, po
        this % pop(j,k) = GL_Polynomial(j, xo, xq(k)) * wq(k) * dx * HALF
      end do
      end do
    case('GRL') ! Gauss-Radau-Legendre
      do k = 1, nq
      do j = 0, po
        this % pop(j,k) = GRL_Polynomial(j, xo, xq(k)) * wq(k) * dx * HALF
      end do
      end do
    case default
      do k = 1, nq
      do j = 0, po
        this % pop(j,k) = GLL_Polynomial(j, xo, xq(k)) * wq(k) * dx * HALF
      end do
      end do
    end select

  end associate

end subroutine Init_SX

!-------------------------------------------------------------------------------
!> Delete ProjectionOperator object

subroutine Delete_ProjectionOperator(this)
  class(ProjectionOperator), intent(inout) :: this

  if (allocated(this % pop)) deallocate(this % pop)

end subroutine Delete_ProjectionOperator

!===============================================================================

end module Projection_Operator
