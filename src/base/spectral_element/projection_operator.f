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
  !> Element-based projection operators
  !>
  !> Provides the one-dimensional operators `A` and `MA` such that for a single
  !> element of length `dx`
  !>
  !>      (MA u)ᵢ = ∑ 𝓁ᵢ(ξ(:)) w(:) u(x(ξ(:))) ≈ ∫ 𝓁ᵢu dx
  !>
  !> is the mass-weighted projection of `u` using quadrature points `ξ(:)` and
  !> wights `w(:)` and
  !>
  !>      (A u)ᵢ = (M⁻¹ MA u)ᵢ ≈ ∫ 𝓁ᵢu dξ / ∫𝓁ᵢ𝓁ᵢ dξ
  !>
  !> is the L2 projection of `u`.

  type, public :: ProjectionOperator

    integer :: nq = -1 !< number of quadrature points per direction
    integer :: np = -1 !< number of projection points per direction
    real(RNP), allocatable :: A (:,:) !< 1D L2 projection operator
    real(RNP), allocatable :: MA(:,:) !< 1D mass-weighted projection operator

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

  real(RNP), allocatable :: VL(:,:)
  integer :: i, k, nq

  if (this % nq > 0) call Delete_ProjectionOperator(this)

  associate(po => eop%po, xo => eop%x)

    nq = size(xq)

    this % nq = nq
    this % np = po + 1

    allocate(this % A (0:po,1:nq))
    allocate(this % MA(0:po,1:nq))

    associate(A => this%A, MA => this%MA)

      ! mass-weighted L2 projection operator for standard element
      select case(eop % basis)
      case('GL ') ! Gauss-Legendre
        do k = 1, nq
        do i = 0, po
          MA(i,k) = GL_Polynomial(i, xo, xq(k)) * wq(k)
        end do
        end do
      case('GRL') ! Gauss-Radau-Legendre
        do k = 1, nq
        do i = 0, po
          MA(i,k) = GRL_Polynomial(i, xo, xq(k)) * wq(k)
        end do
        end do
      case default
        do k = 1, nq
        do i = 0, po
          MA(i,k) = GLL_Polynomial(i, xo, xq(k)) * wq(k)
        end do
        end do
      end select

      ! L2 projection operator
      allocate(VL(0:po,0:po))
      call eop % GetLegendreVDM(VL)
      A = matmul(transpose(VL), MA)
      do k = 1, nq
      do i = 0, po
        A(i,k) = (i + HALF) * A(i,k)
      end do
      end do
      A = matmul(VL, A)

      ! mass-weighted L2 projection operator for physical element
      MA = HALF * dx * MA

    end associate
  end associate

end subroutine Init_SX

!-------------------------------------------------------------------------------
!> Delete ProjectionOperator object

subroutine Delete_ProjectionOperator(this)
  class(ProjectionOperator), intent(inout) :: this

  if (allocated(this % MA)) deallocate(this % MA)

end subroutine Delete_ProjectionOperator

!===============================================================================

end module Projection_Operator
