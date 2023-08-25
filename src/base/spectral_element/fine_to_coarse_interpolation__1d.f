!> summary:  Fine-to-coarse hp-interpolation operators
!> author:   Joerg Stiller
!> date:     2023/08/18
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Fine_To_Coarse_Interpolation__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_Control

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse hp-interpolation operator
  !>
  !> Provides 1D operators for interpolating data from a fine mesh to a coarser
  !> a coarser one. Pure p-coarsening as well as hp-coarsening are supported. In
  !> the latter case, a regular 1:2 refinement of the coarse element is assumed.
  !> When initializing set `mode` to 1 for p-coarsening and 2 for hp-coarsening.
  !> `mode = 0` indicates the identity, i.e., no interpolation is performed and,
  !> hence, no interpolation operator is provided.
  !> In the case of hp-coarsening, the treatment of discontinuities between the
  !> two fine elements is controlled by the `smoothing` parameter. The following
  !> choices are available:
  !>   - `0`  no discontinuity handling
  !>   - `1`  jump removal using antisymmetric linear correction
  !>   - `2`  jump removal by averaging interface coefficients
  !> Option `2` requires a boundary-interior decomposition and, hence, is
  !> available only with equidistant (`E`) or Lobatto (`L`) nodal bases.

  type, public :: FineToCoarseInterpolation_1D
    integer   :: po_f      = -1        !< polynomial order of fine mesh
    integer   :: po_c      = -1        !< polynomial order of coarse mesh
    integer   :: mode      = -1        !< coarsening mode {0,1,2}
    character :: basis     = 'L'       !< basis type {'E','G','L'}
    integer   :: smoothing =  0        !< discontinuity handling {0,1,2}
    real(RNP), allocatable :: A(:,:,:) !< interpolation operators
    real(RNP), allocatable :: B(:,:)   !< blending operators
  end type FineToCoarseInterpolation_1D

  ! constructor interface
  interface FineToCoarseInterpolation_1D
    procedure New_FineToCoarseInterpolation
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse interpolation constructor

  function New_FineToCoarseInterpolation &
               (po_f, po_c, mode, basis, smoothing) result(this)

    integer,             intent(in) :: po_f      !< polynomial order of child
    integer,             intent(in) :: po_c      !< polynomial order of parent
    integer,             intent(in) :: mode      !< coarsening mode
    character, optional, intent(in) :: basis     !< basis type  ['L']
    integer,   optional, intent(in) :: smoothing !< discontinuity handling [0]
    type(FineToCoarseInterpolation_1D) :: this

    call Init_FineToCoarseInterpolation(this, po_f, po_c, mode, basis, smoothing)

  end function New_FineToCoarseInterpolation

  !-----------------------------------------------------------------------------
  !> Build fine-to-coarse interpolation

  subroutine Init_FineToCoarseInterpolation &
                 (this, po_f, po_c, mode, basis, smoothing)

    class(FineToCoarseInterpolation_1D), intent(inout) :: this
    integer,             intent(in) :: po_f      !< polynomial order of child
    integer,             intent(in) :: po_c      !< polynomial order of parent
    integer,             intent(in) :: mode      !< coarsening mode
    character, optional, intent(in) :: basis     !< basis type  ['L']
    integer,   optional, intent(in) :: smoothing !< discontinuity handling [0]

    real(RNP), allocatable :: x_f(:), x_c(:)
    integer   :: i, i2, j, qc

    ! initialization ...........................................................

    this % po_f = po_f
    this % po_c = po_c
    this % mode = mode

    if (present(basis)) then
      if (scan(basis, 'EGL') > 0) then
        this % basis = basis
      else
        call Error( 'Init_FineToCoarseInterpolation'        &
                  , 'basis "' // basis // '" not supported' &
                  , 'Fine_To_Coarse_Interpolation__1D'      )
      end if
    end if

    if (present(smoothing)) then
      select case(this % basis)
      case('E','L')
        this % smoothing = max(min(smoothing, 2), 0)
      case default
        this % smoothing = max(min(smoothing, 1), 0)
      end select
    end if

    if (allocated(this % A)) deallocate(this % A)
    if (allocated(this % B)) deallocate(this % B)

    select case(this % mode)
    case(0)
      return
    case(1)
      allocate(this % A(0:po_c,0:po_f,1))
    case(2)
      allocate(this % A(0:po_c/2,0:po_f,2))
      if (this % smoothing == 1) then
        allocate(this % B(0:po_f,2))
      end if
    end select

    ! collocation points .......................................................

    select case(this % basis)
    case('E')
      allocate(x_f(0:po_f), source = [(i * TWO/po_f - ONE, i = 0, po_f)])
      allocate(x_c(0:po_c), source = [(i * TWO/po_c - ONE, i = 0, po_c)])
    case('G')
      allocate(x_f(0:po_f), source = GaussPoints(po_f))
      allocate(x_c(0:po_c), source = GaussPoints(po_c))
    case('L')
      allocate(x_f(0:po_f), source = LobattoPoints(po_f))
      allocate(x_c(0:po_c), source = LobattoPoints(po_c))
    end select

    ! interpolation matrices ...................................................

    select case(this % mode)

    case(1)
      ! 1:1 interpolation
      select case(this % basis)
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
      end select

    case(2)
      ! 2:1 interpolation
      qc = po_c/2
      i2 = qc + mod(po_c,2)
      select case(this % basis)
      case('E') ! Nodal with equidistant spacing
        do i = 0, qc
        do j = 0, po_f
          this % A(i,j,1) = LagrangePolynomial(j, x_f, 2*x_c(i   ) + ONE)
          this % A(i,j,2) = LagrangePolynomial(j, x_f, 2*x_c(i+i2) - ONE)
        end do
        end do
      case('G') ! Gauss
        do i = 0, qc
        do j = 0, po_f
          this % A(i,j,1) = GaussPolynomial(j, x_f, 2*x_c(i   ) + ONE)
          this % A(i,j,2) = GaussPolynomial(j, x_f, 2*x_c(i+i2) - ONE)
        end do
        end do
      case('L') ! Lobatto
        do i = 0, qc
        do j = 0, po_f
          this % A(i,j,1) = LobattoPolynomial(j, x_f, 2*x_c(i   ) + ONE)
          this % A(i,j,2) = LobattoPolynomial(j, x_f, 2*x_c(i+i2) - ONE)
        end do
        end do
      end select

      ! averaging at interface
      if (mod(po_c,2) == 0) then
        do j = 0, po_f
          this % A(qc,j,1) = HALF * this % A(qc,j,1)
          this % A( 0,j,2) = HALF * this % A( 0,j,2)
        end do
      end if

      ! smoothing
      if (this % smoothing == 1) then
        do i = 0, po_f
          this % B(i,1) = -(ONE + x_f(i)) / 4
          this % B(i,2) =  (ONE - x_f(i)) / 4
        end do
      end if

    end select

  end subroutine Init_FineToCoarseInterpolation

  !=============================================================================

end module Fine_To_Coarse_Interpolation__1D
