!> summary:  Coarse-to-fine hp-interpolation operators
!> author:   Joerg Stiller
!> date:     2023/08/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Coarse_To_Fine_Interpolation__1D
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use Execution_Control

  !-----------------------------------------------------------------------------
  !> Coarse-to-fine hp-interpolation operator
  !>
  !> Provides 1D operators for interpolating data from a fine mesh to a coarser
  !> a coarser one. Pure p-refinement as well as hp-refinement are supported. In
  !> the latter case, a regular 2:1 refinement of the coarse element is assumed.
  !> When initializing set `mode` to 1 for p-refinement and 2 for hp-refinement.
  !> `mode = 0` indicates the identity, i.e., no interpolation is performed and,
  !> hence, no interpolation operator is provided.

  type CoarseToFineInterpolation_1D
    integer   :: pc     = -1           !< polynomial order of coarse mesh
    integer   :: pf     = -1           !< polynomial order of fine mesh
    integer   :: mode   = -1           !< refinement mode {0,1,2}
    character :: basis  = 'L'          !< basis type {'E','G','L'}
    real(RNP), allocatable :: A(:,:,:) !< interpolation operators
  end type CoarseToFineInterpolation_1D

  ! constructor interface
  interface CoarseToFineInterpolation_1D
    procedure New_CoarseToFineInterpolation
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse interpolation constructor

  function New_CoarseToFineInterpolation(pc, pf, mode, basis) result(this)
    integer,             intent(in) :: pc    !< polynomial order of parent
    integer,             intent(in) :: pf    !< polynomial order of child
    integer,             intent(in) :: mode  !< coarsening mode
    character, optional, intent(in) :: basis !< basis type  ['L']
    type(CoarseToFineInterpolation_1D) :: this

    call Init_CoarseToFineInterpolation(this, pc, pf, mode, basis)

  end New_CoarseToFineInterpolation

  !-----------------------------------------------------------------------------
  !> Build fine-to-coarse interpolation

  subroutine Init_CoarseToFineInterpolation(this, pc, pf, mode, basis)
    class(CoarseToFineInterpolation_1D), intent(inout) :: this
    integer,             intent(in) :: pf    !< polynomial order of child
    integer,             intent(in) :: pc    !< polynomial order of parent
    integer,             intent(in) :: mode  !< coarsening mode
    character, optional, intent(in) :: basis !< basis type  ['L']

    real(RNP), allocatable :: xf(:), xc(:)
    integer :: i, j

    ! initialization ...........................................................

    this % pc = pc
    this % pf = pf
    this % mode = mode

    if (present(basis)) then
      if (scan(basis, 'EGL') > 0) then
        this % basis = basis
      else
        call Error( 'Init_CoarseToFineInterpolation'        &
                  , 'basis "' // basis // '" not supported' &
                  , 'Coarse_To_Fine_Interpolation__1D'      )
      end if
    end if

    if (allocated(this % A)) deallocate(this % A)

    select case(this % mode)
    case(0)
      return
    case(1:2)
      allocate(this % A(0:pf,0:pc,mode))
    end select

    ! collocation points .......................................................

    select case(this % basis)
    case('E')
      allocate(xc(0:pc), source = [(i * TWO/pc - ONE, i = 0, pc)])
      allocate(xf(0:pf), source = [(i * TWO/pf - ONE, i = 0, pf)])
    case('G')
      allocate(xc(0:pc), source = GaussPoints(pc))
      allocate(xf(0:pf), source = GaussPoints(pf))
    case('L')
      allocate(xc(0:pc), source = LobattoPoints(pc))
      allocate(xf(0:pf), source = LobattoPoints(pf))
    end select

    ! interpolation matrices ...................................................

    select case(this % mode)

    case(1)
      ! 1:1 interpolation
      select case(this % basis)
      case('E') ! Nodal with equidistant spacing
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = LagrangePolynomial(j, xc, xf(i))
        end do
        end do
      case('G') ! Gauss
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = GaussPolynomial(j, xc, xf(i))
        end do
        end do
      case('L') ! Lobatto
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = LobattoPolynomial(j, xc, xf(i))
        end do
        end do
      end select

    case(2)
      ! 1:2 interpolation
      select case(this % basis)
      case('E') ! Nodal with equidistant spacing
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = LagrangePolynomial(j, xc, HALF * (xf(i) - ONE))
          this % A(i,j,2) = LagrangePolynomial(j, xc, HALF * (xf(i) + ONE))
        end do
        end do
      case('G') ! Gauss
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = GaussPolynomial(j, xc, HALF * (xf(i) - ONE))
          this % A(i,j,2) = GaussPolynomial(j, xc, HALF * (xf(i) + ONE))
        end do
        end do
      case('L') ! Lobatto
        do i = 0, pf
        do j = 0, pc
          this % A(i,j,1) = LobattoPolynomial(j, xc, HALF * (xf(i) - ONE))
          this % A(i,j,2) = LobattoPolynomial(j, xc, HALF * (xf(i) + ONE))
        end do
        end do
      end select

    end select

  end subroutine Init_CoarseToFineInterpolation

  !=============================================================================

end module Coarse_To_Fine_Interpolation__1D
