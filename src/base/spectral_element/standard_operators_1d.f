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
    real(RNP), allocatable :: V(:,:)     !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: V_inv(:,:) !< inverse Legendre-Vandermonde matrix

  contains

    generic,   public  :: New => New_StandardOperators1D
    procedure, private :: New_StandardOperators1D
    procedure, public  :: PolynomialOrder
    procedure, public  :: InitVandermondeMatrix
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

subroutine New_StandardOperators1D(this, po, basis, no_vdm)

  !> standard operators that will be initialized
  class(StandardOperators1D), intent(inout) :: this

  !> polynomial order
  integer, intent(in) :: po

  !> point set for generating the Lagrange basis {GL,GLL,GRL} [GLL]
  character(len=*), optional, intent(in) :: basis

  !> switch to skip the generation of the Vandermonde matrix and its inverse
  logical, optional, intent(in) :: no_vdm

  character(len=3) :: chosen_basis
  logical :: build_vdm
  integer :: i, j

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

  ! Vandermonde matrix and its inverse .........................................

  if (present(no_vdm)) then
    build_vdm = .not. no_vdm
  else
    build_vdm = .true.
  end if

  if (build_vdm) then
    call InitVandermondeMatrix(this)
  end if

end subroutine New_StandardOperators1D

!-------------------------------------------------------------------------------
!> Intializes the Legendre-Vandermonde matrix and its inverse.

subroutine InitVandermondeMatrix(this)
  class(StandardOperators1D), intent(inout) :: this !< standard operators

  integer :: i, j, po

  ! safeguard ..................................................................

  if (allocated(this % V    )) deallocate(this % V    )
  if (allocated(this % V_inv)) deallocate(this % V_inv)

  ! prerequisites ..............................................................

  po = this % po

  ! Vandermonde matrix .........................................................

  allocate(this % V(0:po,0:po))
  do i = 0, po
  do j = 0, po
     this % V(i,j) = JacobiPolynomial(a=ZERO, b=ZERO, n=j, x=this%x(i))
  end do
  end do

  ! inverse Vandermonde matrix .................................................

  allocate(this % V_inv(0:po,0:po), source=this%V)
  this % V_inv = Inverse(this % V)

end subroutine InitVandermondeMatrix

!-------------------------------------------------------------------------------
!> Get the Legendre-Vandermonde matrix

subroutine GetVandermondeMatrix(this, V)
  class(StandardOperators1D), intent(in) :: this   !< standard operators
  real(RNP), intent(out) :: V(0:this%po,0:this%po) !< Vandermonde matrix

  if (.not. allocated(this % V)) then
    call Error( 'GetVandermondeMatrix'               &
              , 'Vandermonde matrix not initialized' &
              , 'Standard_Operators_1D'              )
  end if

  V = this % V

end subroutine GetVandermondeMatrix

!-------------------------------------------------------------------------------
!> Get the inverse Legendre-Vandermonde matrix

subroutine GetInverseVandermondeMatrix(this, V_inv)
  class(StandardOperators1D), intent(in) :: this       !< standard operators
  real(RNP), intent(out) :: V_inv(0:this%po,0:this%po) !< inverse VDM matrix

  if (.not. allocated(this%V_inv)) then
    call Error( 'GetInverseVandermondeMatrix'        &
              , 'Vandermonde matrix not initialized' &
              , 'Standard_Operators_1D'              )
  end if

  V_inv = this % V_inv

end subroutine GetInverseVandermondeMatrix

!-------------------------------------------------------------------------------
!> Finalization

subroutine Delete_StandardOperators1D(this)
  type(StandardOperators1D), intent(inout) :: this  !< standard operators

  if(allocated(this%x    )) deallocate(this%x    )
  if(allocated(this%w    )) deallocate(this%w    )
  if(allocated(this%D    )) deallocate(this%D    )
  if(allocated(this%L    )) deallocate(this%L    )
  if(allocated(this%V    )) deallocate(this%V    )
  if(allocated(this%V_inv)) deallocate(this%V_inv)

end subroutine Delete_StandardOperators1D

!===============================================================================

end module Standard_Operators_1D
