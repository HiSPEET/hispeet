!> summary:  Projection operator
!> author:   Joerg Stiller
!> date:     2019/02/22, 2024/12/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Projection_Operator__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, TWO
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Standard_Element_Operators__1D
  implicit none
  private

  public :: ProjectionOperator_1D

  !-----------------------------------------------------------------------------
  !> Element-based projection operators
  !>
  !> Provides the one-dimensional operators `A` and `MA` such that
  !>
  !>      (MA f)ᵢ ≈ ∫ 𝓁ᵢ(ξ) f(ξ) dξ
  !>
  !> is the mass-weighted projection of the function f  in the standard element
  !> [-1,1] and
  !>
  !>      (A u)ᵢ = (M⁻¹ MA u)ᵢ ≈ ∫ 𝓁ᵢ(ξ) f(ξ) dξ / ∫𝓁ᵢ𝓁ᵢ dξ
  !>
  !> is the L2 projection of f.

  type ProjectionOperator_1D
    integer :: po_b = -1 !< polynomial order of the basis {𝓁ᵢ}
    integer :: po_f = -1 !< polynomial order of the projected function
    real(RNP), allocatable :: A (:,:) !< L2 projection operator (0:po_b,0:po_f)
    real(RNP), allocatable :: MA(:,:) !< mass-weighted operator (0:po_b,0:po_f)
  end type ProjectionOperator_1D

  ! constructor interface
  interface ProjectionOperator_1D
    module procedure New_ProjectionOperator_1D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> New ProjectionOperator_1D from 1D standard operators and interpolation

  function New_ProjectionOperator_1D(eop, x, nodes) result(this)
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP),              intent(in) :: x(0:) !< node coordinates
    character(*), optional, intent(in) :: nodes !< node type
    type(ProjectionOperator_1D) :: this

    call Init_ProjectionOperator_1D(this, eop, x, nodes)

  end function New_ProjectionOperator_1D

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine Init_ProjectionOperator_1D(this, eop, x, nodes)
    class(ProjectionOperator_1D), intent(inout) :: this
    class(StandardElementOperators_1D), intent(in) :: eop !< standard operators
    real(RNP),              intent(in) :: x(0:) !< node coordinates
    character(*), optional, intent(in) :: nodes !< node type

    real(RNP), allocatable :: x_g(:), w_g(:)
    real(RNP), allocatable :: x_q(:), w_q(:)
    real(RNP), allocatable :: C(:,:), VL_inv(:,:)
    real(RNP) :: z

    character(2) :: nodes_f
    integer :: po_b, po_f, po_q
    integer :: i, j, k, q

    ! initialization ...........................................................

    po_b = eop % po
    po_f = ubound(x, 1)
    po_q = max(po_b, po_f)

    if (present(nodes)) then
      nodes_f = nodes
    else
      nodes_f = ''
    end if

    this % po_b = po_b
    this % po_f = po_f

    allocate(this % A  (0:po_b, 0:po_f), source = ZERO)
    allocate(this % MA (0:po_b, 0:po_f), source = ZERO)

    ! Gauss points and weights to basis order
    allocate(x_g(0:po_b), source = GaussPoints(po_b))
    allocate(w_g(0:po_b), source = GaussWeights(x_g))

    ! Gauss points and weights for exact integration
    allocate(x_q(0:po_q), source = GaussPoints(po_q))
    allocate(w_q(0:po_q), source = GaussWeights(x_q))

    associate(A => this % A, MA => this % MA)

      ! intermediate projection to Gauss points ................................

      allocate(C(0:po_b,0:po_f), source = ZERO)

      do j = 0, po_f
      do q = 0, po_q
        select case(nodes_f)
        case('G')
          z = GaussPolynomial(j, x, x_q(q))
        case('L')
          z = LobattoPolynomial(j, x, x_q(q))
        case('RL','RR')
          z = RadauPolynomial(j, x, x_q(q))
        case default
          z = LagrangePolynomial(j, x, x_q(q))
        end select
        do k = 0, po_b
          C(k,j) = C(k,j) + w_q(q) / w_g(k) * z * GaussPolynomial(k, x_g, x_q(q))
        end do
      end do
      end do

      ! combine with interpolation to collocation points .......................

      do i = 0, po_b
      do k = 0, po_b
        z = GaussPolynomial(k, x_g, eop%x(i))
        do j = 0, po_f
          A(i,j) = A(i,j) + z * C(k,j)
        end do
      end do
      end do

      ! mass weighted projection operator ......................................

      ! inverse Vandermonde matrix V⁻¹
      allocate(VL_inv(0:po_b,0:po_b))
      call eop % Get_Inverse_Legendre_VDM(VL_inv)

      ! MA = V⁻ᵀ V⁻¹ A
      MA = matmul(VL_inv, A)
      do k = 0, po_f
      do i = 0, po_b
        MA(i,k) = TWO/(2*i + 1) * A(i,k)
      end do
      end do
      MA = matmul(transpose(VL_inv), MA)

    end associate

  end subroutine Init_ProjectionOperator_1D

  !=============================================================================

end module Projection_Operator__1D
