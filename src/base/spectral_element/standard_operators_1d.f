!> summary:   Spectral element operators in the one-dimensional standard region
!> author:    Immo Huismann, Joerg Stiller
!> date:      2014/11/24
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Spectral element operators in the one-dimensional standard region
!===============================================================================

module Standard_Operators_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO
  use Execution_Control, only: Warning, Error
  use Gauss_Jacobi
  use Matrix_Operators,  only: Inverse
  implicit none
  private

  !-----------------------------------------------------------------------------
  !> Standard operators for one polynomial order.

  type, public :: StandardOperators1D
    private

    ! public components
    integer               , public :: po = -1  !< polynomial order
    real(RNP), allocatable, public :: x(:)     !< collocation points
    real(RNP), allocatable, public :: w(:)     !< quadrature weights
    real(RNP), allocatable, public :: D(:,:)   !< differention matrix
    real(RNP), allocatable, public :: L(:,:)   !< stiffness (Laplace) matrix

    ! private comoponents
    real(RNP), allocatable :: V(:,:)   !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: VI(:,:)  !< inverse Legendre-Vandermonde matrix

  contains

    generic,   public  :: New => New_StandardOperators1D
    procedure, private :: New_StandardOperators1D
    procedure, public  :: PolynomialOrder
    procedure, public  :: GetVandermondeMatrix
    procedure, public  :: GetInverseVandermondeMatrix

    final :: Delete_StandardOperators1D

  end type StandardOperators1D

contains

!-------------------------------------------------------------------------------
!> Returns the polymial order of operators, or -1 if none

pure integer function PolynomialOrder(this) result(po)
  class(StandardOperators1D), intent(in) :: this !< standard operators

  po = this % po

end function PolynomialOrder

!-------------------------------------------------------------------------------
!> Intializes the public components.
!>
!> The routine provides 1D standard operators for the following basis systems
!>
!>   *  Lagrange polynomials to Gauss-Legendre points (`basis = 'GL'`)
!>   *  Lagrange polynomials to Gauss-Lobatto-Legendre points (`basis = 'GLL'`)
!>   *  Lagrange polynomials to Gauss-Radau-Legendre points (`basis = 'GRL'`)
!>
!> GLL is the default, but po = 0 implies GL irrespective of the chosen basis.

subroutine New_StandardOperators1D(this, po, basis)
  class(StandardOperators1D), intent(inout) :: this  !< standard operators
  integer,                    intent(in)    :: po    !< polynomial order
  character(len=*), optional, intent(in)    :: basis !< point set {GL,GLL,GRL}

  integer :: i, j
  character(len=3) :: chosen_basis

  ! safeguard ..................................................................

  call Delete_StandardOperators1D(this)

  if (po < 0) then

    call Error('New_StandardOperators1D', 'Invalid polynomial order (po < 0)')

  else if (po == 0) then

    chosen_basis = 'GL'

  else

    if (present(basis)) then
      if (any(basis == [ 'GL ', 'GLL', 'GRL' ])) then
        chosen_basis = basis
      else
        call Error('New_StandardOperators1D', 'Invalid basis argument')
      end if
    else
      chosen_basis = 'GLL'
    end if

  end if

  ! operators available from Gauss-Jacobi module ...............................

  this%po = po

  if (po == 0) then

    ! 1-point Gauss
    allocate(this%x(0:0),       source = ZERO)
    allocate(this%w(0:0),       source = TWO)
    allocate(this%D(0:0,0:0),   source = ZERO)

  else if (chosen_basis == 'GL') then

    ! Gauss-Legendre
    allocate(this%x(0:po),      source = GL_Points(po))
    allocate(this%w(0:po),      source = GL_Weights(this%x))
    allocate(this%D(0:po,0:po), source = GL_DiffMatrix(this%x))

  else if (chosen_basis == 'GLL') then

    ! Gauss-Lobatto-Legendre points
    allocate(this%x(0:po),      source = GLL_Points(po))
    allocate(this%w(0:po),      source = GLL_Weights(this%x))
    allocate(this%D(0:po,0:po), source = GLL_DiffMatrix(this%x))

  else

    ! Gauss-Radau-Legendre
    allocate(this%x(0:po),      source = GRL_Points(po))
    allocate(this%w(0:po),      source = GRL_Weights(this%x))
    allocate(this%D(0:po,0:po), source = GRL_DiffMatrix(this%x))

  end if

  ! stiffness matrix ...........................................................

  allocate(this%L(0:po,0:po), source = ZERO)

  do i = 0, po
  do j = 0, po
    this%L(i,j) = this%L(i,j) + sum(this%w * this%D(:,i) * this%D(:,j))
  end do
  end do

end subroutine New_StandardOperators1D

!-------------------------------------------------------------------------------
!> Get the Legendre-Vandermonde matrix

subroutine GetVandermondeMatrix(sop, V)
  class(StandardOperators1D), intent(inout) :: sop !< standard operators
  real(RNP), intent(out) :: V(0:sop%po,0:sop%po)   !< Vandermonde matrix

  if (sop%po < 0) then
    call Error('GetVandermondeMatrix', 'sop must be initialized!')
  end if

  if (.not. allocated(sop%V)) then
    call InitVandermondeMatrix(sop)
  end if

  V = sop % V

end subroutine GetVandermondeMatrix

!-------------------------------------------------------------------------------
!> Get the inverse Legendre-Vandermonde matrix

subroutine GetInverseVandermondeMatrix(sop, VI)
  class(StandardOperators1D), intent(inout) :: sop !< standard operators
  real(RNP), intent(out) :: VI(0:sop%po,0:sop%po)  !< inverse Vandermonde matrix

  if (sop%po < 0) then
    call Error('GetVandermondeMatrix', 'sop must be initialized!')
  end if

  if (.not. allocated(sop%VI)) then
    call InitVandermondeMatrix(sop)
  end if

  VI = sop % VI

end subroutine GetInverseVandermondeMatrix

!-------------------------------------------------------------------------------
!> Intializes the Legendre-Vandermonde matrix and its inverse.

subroutine InitVandermondeMatrix(sop)
  class(StandardOperators1D), intent(inout) :: sop !< standard operators

  integer :: i, j, po

  ! safeguard ..................................................................

  if (allocated(sop%V )) deallocate(sop%V )
  if (allocated(sop%VI)) deallocate(sop%VI)

  ! prerequisites ..............................................................

  po = sop%po

  ! Vandermonde matrix .........................................................

  allocate(sop%V(0:po,0:po))
  do i = 0, po
  do j = 0, po
     sop%V(i,j) = JacobiPolynomial(a=ZERO, b=ZERO, n=j, x=sop%x(i))
  end do
  end do

  ! inverse Vandermonde matrix .................................................

  allocate(sop%VI(0:po,0:po), source=sop%V)
  sop%VI = Inverse(sop%V)

end subroutine InitVandermondeMatrix

!-------------------------------------------------------------------------------
!> Finalization

subroutine Delete_StandardOperators1D(this)
  type(StandardOperators1D), intent(inout) :: this  !< standard operators

  if(allocated(this%x))  deallocate(this%x)
  if(allocated(this%w )) deallocate(this%w )
  if(allocated(this%D )) deallocate(this%D )
  if(allocated(this%L )) deallocate(this%L )
  if(allocated(this%V )) deallocate(this%V )
  if(allocated(this%VI)) deallocate(this%VI)

end subroutine Delete_StandardOperators1D

!===============================================================================

end module Standard_Operators_1D
