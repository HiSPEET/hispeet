!> summary:  Fine-to-coarse hp-projection operators
!> author:   Joerg Stiller
!> date:     2023/08/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module HP__Coarsening_Operator__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_Control

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> hp-coarsening operator
  !>
  !> Provides 1D operators to project data from a fine mesh to a coarser one.
  !> Pure p-coarsening as well as hp-coarsening are supported. In the latter
  !> case, a regular 1:2 refinement of the coarse element is assumed.
  !>
  !> For p-coarsening specify `mode = 1` and for hp-coarsening `mode = 2.
  !> To obtain the identity operator pass `mode = 0`.
  !>
  !> The fine and coarse variables are given as nodal values according to one
  !> of the following Lagrange bases:
  !>   - `E`  equidistant,
  !>   - `G`  Gauss,
  !>   - `L`  Lobatto,
  !>   - `RL` Radau left,
  !>   - `RR` Radau right.
  !>
  !> Two different projection methods are provided:
  !>   - `P`  L² projection,
  !>   - `I`  interpolation using the fine basis.
  !>
  !> In the case of hp-coarsening, the treatment of discontinuities between the
  !> two fine elements is controlled by the `smooth` parameter. The following
  !> choices are available:
  !>   - `0`  no discontinuity handling
  !>   - `1`  jump removal using antisymmetric linear correction
  !>   - `2`  jump removal by averaging interface coefficients
  !>
  !> Option `2` requires a boundary-interior decomposition and is therefore
  !> restricted to equidistant (`E`) or Lobatto (`L`) bases.

  type, public :: HP_CoarseningOperator_1D
    character(len=2) :: nodes          !< nodal basis type
    integer          :: po_f           !< polynomial order of fine mesh
    integer          :: po_c           !< polynomial order of coarse mesh
    integer          :: mode           !< coarsening mode
    character        :: method         !< projection method
    real(RNP), allocatable :: A(:,:,:) !< interpolation operator(s)
  end type HP_CoarseningOperator_1D

  ! constructor interface
  interface HP_CoarseningOperator_1D
    procedure New_HP_CoarseningOperator
  end interface

  !-----------------------------------------------------------------------------
  !> Options for initializing the hp-coarsening operator

  type, public :: HP_CoarseningOptions_1D
    character(len=2) :: nodes  = 'L' !< node set {'E','G','L','RL','RR'}
    integer          :: po_f   = -1  !< polynomial order of fine mesh
    integer          :: po_c   = -1  !< polynomial order of coarse mesh
    integer          :: mode   = -1  !< coarsening mode {0,1,2}
    character        :: method = 'I' !< L² projection 'P' or interpolation 'I'
    integer          :: smooth =  0  !< discontinuity handling {0,1,2}
  end type HP_CoarseningOptions_1D

contains

  !-----------------------------------------------------------------------------
  !> hp-coarsening operator constructor

  function New_HP_CoarseningOperator(opt) result(this)
    class(HP_CoarseningOptions_1D), intent(in) :: opt
    type(HP_CoarseningOperator_1D) :: this

    call Init_HP_CoarseningOperator(this, opt)

  end function New_HP_CoarseningOperator

  !-----------------------------------------------------------------------------
  !> Build hp-coarsening operator

  subroutine Init_HP_CoarseningOperator(this, opt)
    class(HP_CoarseningOperator_1D), intent(inout) :: this
    class(HP_CoarseningOptions_1D), intent(in) :: opt

    real(RNP), allocatable :: x_f(:), x_c(:), x_g(:), x_q(:), w_g(:), w_q(:)
    integer :: po_c, po_f, po_q, smooth

    ! initialization ...........................................................

    this % po_f   = opt % po_f
    this % po_c   = opt % po_c
    this % mode   = opt % mode
    this % nodes  = opt % nodes
    this % method = opt % method

    select case(this % nodes)
    case('E','L')
      smooth = max(min(opt % smooth, 2), 0)
    case default
      smooth = max(min(opt % smooth, 1), 0)
    end select

    if (allocated(this % A)) deallocate(this % A)

    ! shorthands
    po_f = this % po_f
    po_c = this % po_c
    po_q = max(po_c, po_f)

    select case(this % mode)
    case(0)
      return
    case(1)
      allocate(this % A(0:po_c,0:po_f,1), source = ZERO)
    case(2)
      allocate(this % A(0:po_c,0:po_f,2), source = ZERO)
    end select

    ! collocation points .......................................................

    select case(this % nodes)
    case('E')
      block
        integer :: i
        allocate(x_f(0:po_f), source = [(i * TWO/po_f - ONE, i = 0, po_f)])
        allocate(x_c(0:po_c), source = [(i * TWO/po_c - ONE, i = 0, po_c)])
      end block
    case('G')
      allocate(x_f(0:po_f), source = GaussPoints(po_f))
      allocate(x_c(0:po_c), source = GaussPoints(po_c))
    case('L')
      allocate(x_f(0:po_f), source = LobattoPoints(po_f))
      allocate(x_c(0:po_c), source = LobattoPoints(po_c))
    case('RL')
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .false.))
      allocate(x_c(0:po_c), source = RadauPoints(po_c, right = .false.))
    case('RR')
      allocate(x_f(0:po_f), source = RadauPoints(po_f, right = .true.))
      allocate(x_c(0:po_c), source = RadauPoints(po_c, right = .true.))
    end select

    ! quadrature points and weights ............................................

    if (this % mode > 0 .and. this % method == 'P') then

      ! Gauss points and weights to coarse element degree
      allocate(x_g(0:po_c), source = GaussPoints(po_c))
      allocate(w_g(0:po_c), source = GaussWeights(x_g))

      ! Gauss points and weights for exact integration
      allocate(x_q(0:po_q), source = GaussPoints(po_q))
      allocate(w_q(0:po_q), source = GaussWeights(x_q))

    end if

    ! projection operators .....................................................

    select case(this % mode)

    case(1)

      ! 1:1 projection
      select case(this % method)
      case('P')
        call Build_L2_Projection_Operator_1
      case('I')
        call Build_Interpolation_Operator_1
      end select

    case(2)

      ! 2:1 projection
      select case(this % method)
      case('P')
        call Build_L2_Projection_Operator_2
      case('I')
        call Build_Interpolation_Operator_2
      end select
      call Add_Smoothing

    end select

  contains

    !---------------------------------------------------------------------------
    !> Build 1:1 L² projection

    subroutine Build_L2_Projection_Operator_1

      real(RNP), allocatable :: C(:,:)
      real(RNP) :: z
      integer   :: i, j, k, q

      allocate(C(0:po_c,0:po_f), source = ZERO)

      ! intermediate projection to Gauss points ................................

      do j = 0, po_f
      do q = 0, po_q
        select case(this % nodes)
        case('E')
          z = LagrangePolynomial(j, x_f, x_q(q))
        case('G')
          z = GaussPolynomial(j, x_f, x_q(q))
        case('L')
          z = LobattoPolynomial(j, x_f, x_q(q))
        case('RL','RR')
          z = RadauPolynomial(j, x_f, x_q(q))
        end select
        do k = 0, po_c
          C(k,j) = C(k,j) + w_q(q) / w_g(k) * z * GaussPolynomial(k, x_g, x_q(q))
        end do
      end do
      end do

      ! combine with interpolation to collocation points .......................

      do i = 0, po_c
      do k = 0, po_c
        z = GaussPolynomial(k, x_g, x_c(i))
        do j = 0, po_f
          this % A(i,j,1) = this % A(i,j,1) + z * C(k,j)
        end do
      end do
      end do

    end subroutine Build_L2_Projection_Operator_1

    !---------------------------------------------------------------------------
    !> Build 1:1 projection based on embedded interpolation

    subroutine Build_Interpolation_Operator_1

      integer :: i, j

      select case(this % nodes)

      case('E') ! Nodal with equidistant spacing
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = LagrangePolynomial(j, x_f, x_c(i))
        end do
        end do

      case('G') ! Gauss
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = GaussPolynomial(j, x_f, x_c(i))
        end do
        end do

      case('L') ! Lobatto
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = LobattoPolynomial(j, x_f, x_c(i))
        end do
        end do

      case('RL','RR') ! Radau
        do i = 0, po_c
        do j = 0, po_f
          this % A(i,j,1) = RadauPolynomial(j, x_f, x_c(i))
        end do
        end do

      end select

    end subroutine Build_Interpolation_Operator_1

    !---------------------------------------------------------------------------
    !> Build 2:1 L² projection

    subroutine Build_L2_Projection_Operator_2

      real(RNP), allocatable :: C(:,:,:)
      real(RNP) :: w, z
      integer   :: i, j, k, q

      allocate(C(0:po_c,0:po_f,2), source = ZERO)

      ! intermediate projection to Gauss points ................................

      do j = 0, po_f
      do q = 0, po_q
        select case(this % nodes)
        case('E')
          z = LagrangePolynomial(j, x_f, x_q(q))
        case('G')
          z = GaussPolynomial(j, x_f, x_q(q))
        case('L')
          z = LobattoPolynomial(j, x_f, x_q(q))
        case('RL','RR') ! Radau
          z = RadauPolynomial(j, x_f, x_q(q))
        end select
        do k = 0, po_c
          w = w_q(q) / (2 * w_g(k)) * z
          C(k,j,1) = C(k,j,1) + w * GaussPolynomial(k, x_g, HALF*(x_q(q) - 1))
          C(k,j,2) = C(k,j,2) + w * GaussPolynomial(k, x_g, HALF*(x_q(q) + 1))
        end do
      end do
      end do

      ! combine with interpolation to collocation points .......................

      do i = 0, po_c
      do k = 0, po_c
        z = GaussPolynomial(k, x_g, x_c(i))
        do j = 0, po_f
          this % A(i,j,1) = this % A(i,j,1) + z * C(k,j,1)
          this % A(i,j,2) = this % A(i,j,2) + z * C(k,j,2)
        end do
      end do
      end do

    end subroutine Build_L2_Projection_Operator_2

    !---------------------------------------------------------------------------
    !> Build 2:1 projection based on embedded interpolation

    subroutine Build_Interpolation_Operator_2

      real(RNP) :: x_c2f(2)
      integer   :: o2_f, ph_c
      integer   :: i1, i2, j

      ph_c = po_c / 2
      o2_f = po_c - ph_c

      do i1 = 0, ph_c

        i2 = i1 + o2_f

        x_c2f(1) = 2*x_c(i1) + ONE ! coarse point mapped to left  element
        x_c2f(2) = 2*x_c(i2) - ONE ! coarse point mapped to right element

        select case(this % nodes)
        case('E') ! Nodal with equidistant spacing
          do j = 0, po_f
            this % A(i1,j,1) = LagrangePolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = LagrangePolynomial(j, x_f, x_c2f(2))
          end do
        case('G') ! Gauss
          do j = 0, po_f
            this % A(i1,j,1) = GaussPolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = GaussPolynomial(j, x_f, x_c2f(2))
          end do
        case('L') ! Lobatto
          do j = 0, po_f
            this % A(i1,j,1) = LobattoPolynomial(j, x_f, x_c2f(1))
            this % A(i2,j,2) = LobattoPolynomial(j, x_f, x_c2f(2))
          end do
        case('RL','RR') ! Radau: take care of asymmetry
          if (x_c2f(1) <= ONE) then
            do j = 0, po_f
              this % A(i1,j,1) = RadauPolynomial(j, x_f, x_c2f(1))
            end do
          end if
          if (x_c2f(2) >= -ONE) then
            do j = 0, po_f
              this % A(i2,j,2) = RadauPolynomial(j, x_f, x_c2f(2))
            end do
          end if
        end select

      end do

      select case(this % nodes)
      case('E','G','L')
        ! averaging at interface
        if (mod(po_c,2) == 0) then
          do j = 0, po_f
            this % A(ph_c,j,1) = HALF * this % A(ph_c,j,1)
            this % A(ph_c,j,2) = HALF * this % A(ph_c,j,2)
          end do
        end if
      end select

    end subroutine Build_Interpolation_Operator_2

    !---------------------------------------------------------------------------
    !> Build combined smoothing + projection operator

    subroutine Add_Smoothing

      real(RNP), allocatable :: B(:,:), L(:), R(:)
      real(RNP) :: AB
      integer :: i, j

      if (smooth == 0) return

      ! smothing operator ......................................................

      allocate(B(0:po_f,2), source = ZERO)

      select case(smooth)
      case(1)
        ! linear blending
        B(:,1) = -(ONE + x_f) / 4
        B(:,2) =  (ONE - x_f) / 4
      case(2)
        ! interface coefficient averaging
        B(po_f,1) = -HALF
        B( 0  ,2) =  HALF
      end select

      ! trace operators ........................................................

      allocate(L(0:po_f), R(0:po_f))

      do i = 0, po_f
        select case(this % nodes)
        case('E') ! Nodal with equidistant spacing
          L(i) = LagrangePolynomial(i, x_f, -ONE)
          R(i) = LagrangePolynomial(i, x_f,  ONE)
        case('G') ! Gauss
          L(i) = GaussPolynomial(i, x_f, -ONE)
          R(i) = GaussPolynomial(i, x_f,  ONE)
        case('L') ! Lobatto
          L(i) = LobattoPolynomial(i, x_f, -ONE)
          R(i) = LobattoPolynomial(i, x_f,  ONE)
        case('RL','RR') ! Radau: take care of asymmetry
          L(i) = RadauPolynomial(i, x_f, -ONE)
          R(i) = RadauPolynomial(i, x_f,  ONE)
        end select
      end do

      ! combined operator ......................................................

      associate(A => this % A)
        do i = 0, po_c
          AB =  sum(A(i,:,1) * B(:,1) + A(i,:,2) * B(:,2))
          do j = 0, po_f
            A(i,j,1) = A(i,j,1) + AB * R(j)
            A(i,j,2) = A(i,j,2) - AB * L(j)
          end do
        end do
      end associate

    end subroutine Add_Smoothing

    !---------------------------------------------------------------------------

  end subroutine Init_HP_CoarseningOperator

  !=============================================================================

end module HP__Coarsening_Operator__1D
