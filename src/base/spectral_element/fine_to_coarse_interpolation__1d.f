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
  !> In the case of hp-coarsening, the switch `smooth` can be used to activate
  !> or deactivate the removal of jumps between the 2 fine elements.

  type, public :: FineToCoarseInterpolation_1D
    integer   :: po_f   = -1           !< polynomial order of fine mesh
    integer   :: po_c   = -1           !< polynomial order of coarse mesh
    integer   :: mode   = -1           !< coarsening mode {0,1,2}
    character :: basis  = 'L'          !< basis type {'E','G','L'}
    logical   :: smooth = .true.       !< apply linear blending to remove jumps
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
               (po_f, po_c, mode, basis, smooth) result(this)

    integer,             intent(in) :: po_f   !< polynomial order of child
    integer,             intent(in) :: po_c   !< polynomial order of parent
    integer,             intent(in) :: mode   !< coarsening mode
    character, optional, intent(in) :: basis  !< basis type  ['L']
    logical,   optional, intent(in) :: smooth !< switch for blending  [T]
    type(FineToCoarseInterpolation_1D) :: this

    call Init_FineToCoarseInterpolation(this, po_f, po_c, mode, basis, smooth)

  end New_FineToCoarseInterpolation

  !-----------------------------------------------------------------------------
  !> Build fine-to-coarse interpolation

  subroutine Init_FineToCoarseInterpolation &
                 (this, po_f, po_c, mode, basis, smooth)

    class(FineToCoarseInterpolation_1D), intent(inout) :: this
    integer,             intent(in) :: po_f   !< polynomial order of child
    integer,             intent(in) :: po_c   !< polynomial order of parent
    integer,             intent(in) :: mode   !< coarsening mode
    character, optional, intent(in) :: basis  !< basis type  ['L']
    logical,   optional, intent(in) :: smooth !< switch for blending  [T]

    real(RNP), allocatable :: x_f(:), x_c(:)
    integer :: i, j, k, qc

    ! initialization ...........................................................

    this % po_f = po_f
    this % po_c = po_c
    this % mode = mode

    if (present(smooth)) then
      this % smooth = smooth
    end if

    if (present(basis)) then
      if (scan(basis, 'EGL') > 0) then
        this % basis = basis
      else
        call Error( 'Init_FineToCoarseInterpolation'        &
                  , 'basis "' // basis // '" not supported' &
                  , 'Fine_To_Coarse_Interpolation__1D'      )
      end if
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
      if (this % smooth) then
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
      do k = 1, 2
        x0 = 3 - 2 * k
        select case(this % basis)
        case('E') ! Nodal with equidistant spacing
          do i = 0, qc
          do j = 0, po_f
            this % A(i,j,k) = LagrangePolynomial(j, x_f, x0 + 2*x_c(i))
          end do
          end do
        case('G') ! Gauss
          do i = 0, qc
          do j = 0, po_f
            this % A(i,j,k) = GaussPolynomial(j, x_f, x0 + 2*x_c(i))
          end do
          end do
        case('L') ! Lobatto
          do i = 0, qc
          do j = 0, po_f
            this % A(i,j,k) = LobattoPolynomial(j, x_f, x0 + 2*x_c(i))
          end do
          end do
        end select
        if (this % smooth) then
          do i = 0, po_f
            this % B(i,1) = -(ONE + x_f(i)) / 4
            this % B(i,2) =  (ONE - x_f(i)) / 4
          end do
        end if
      end do

    end select

  end subroutine Init_FineToCoarseInterpolation

  !=============================================================================

end module Fine_To_Coarse_Interpolation__1D
