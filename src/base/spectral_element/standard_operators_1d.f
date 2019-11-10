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
  !> Standard operators for one polynomial order
  !>
  !> Supports the following types of nodal base functions
  !>
  !>   * Lagrange polynomials to Gauss-Legendre points:         `basis = 'GL'`
  !>   * Lagrange polynomials to Gauss-Lobatto-Legendre points: `basis = 'GLL'`
  !>   * Lagrange polynomials to Gauss-Radau-Legendre points:   `basis = 'GRL'`

  type, public :: StandardOperators1D
    private

    ! public components
    character(len=3)      , public :: basis    !< basis type
    integer               , public :: po = -1  !< polynomial order
    real(RNP), allocatable, public :: x(:)     !< collocation points
    real(RNP), allocatable, public :: w(:)     !< quadrature weights
    real(RNP), allocatable, public :: D(:,:)   !< differention matrix
    real(RNP), allocatable, public :: L(:,:)   !< stiffness (Laplace) matrix

    ! private components
    real(RNP), allocatable :: VL(:,:)     !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: VL_inv(:,:) !< inverse Legendre-Vandermonde matrix

    ! for now SVV parameters are not implemented as component of the standard
    ! operator, as this might not be necessary, when they are constant anyway
    ! real(RNP) :: amplitude_svv   !< epsilon, maximum amplitude of spectral
    !                              !< viscosity
    ! integer   :: cutoff_mode_svv !< M, SVV is applied for all modes larger than M

    real(RNP), allocatable :: D_svv(:,:)
    real(RNP), allocatable :: L_svv(:,:)

  contains

    procedure :: Init_StandardOperators1D
    procedure :: PolynomialOrder
    procedure :: InitLegendreVDM
    procedure :: HasLegendreVDM
    procedure :: GetLegendreVDM
    procedure :: GetInverseLegendreVDM

    procedure :: InitSVV
    procedure :: HasSVV
    procedure :: GetSVV_DiffMatrix
    procedure :: GetSVV_StiffnessMatrix

  end type StandardOperators1D

  ! Constructor interface
  interface StandardOperators1D
    module procedure New_StandardOperators1D
  end interface

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> Constructor for StandardOperators1D

  function New_StandardOperators1D(po, basis, no_vdm, svv) result(this)
    integer,                    intent(in) :: po     !< polynomial order
    character(len=*), optional, intent(in) :: basis  !< points {GL,GLL,GRL} [GLL]
    logical,          optional, intent(in) :: no_vdm !< skip Vandermonde matrix [F]
    logical,          optional, intent(in) :: svv    !< activate SVV model [F]

    type(StandardOperators1D) :: this

    call Init_StandardOperators1D(this, po, basis, no_vdm, svv)

  end function New_StandardOperators1D

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Returns the polymial order of operators, or -1 if none

  pure integer function PolynomialOrder(this) result(po)
    class(StandardOperators1D), intent(in) :: this !< standard operators

    po = this % po

  end function PolynomialOrder

  !-----------------------------------------------------------------------------
  !> Intializes the public components.
  !>
  !> The routine provides 1D standard operators for the chosen nodal basis.
  !> GLL is the default, except for `po = 0` which always implies GL.

  subroutine Init_StandardOperators1D(this, po, basis, no_vdm, svv)

    !> standard operators that will be initialized
    class(StandardOperators1D), intent(inout) :: this

    !> polynomial order
    integer, intent(in) :: po

    !> point set for generating the Lagrange basis {GL,GLL,GRL} [GLL]
    character(len=*), optional, intent(in) :: basis

    !> switch to skip the generation of the Vandermonde matrix and its inverse
    logical, optional, intent(in) :: no_vdm

    !> switch to activate the SVV model
    logical, optional, intent(in) :: svv

    logical :: build_vdm
    integer :: i, j

    ! safeguard ................................................................

    call Delete_StandardOperators1D(this)

    if (po < 0) then

      call Error('Init_StandardOperators1D', 'Invalid polynomial order (po < 0)')

    else if (po == 0) then

      this%basis = 'GL'

    else

      if (present(basis)) then
        if (any(basis == [ 'GL ', 'GLL', 'GRL' ])) then
          this%basis = basis
        else
          call Error('Init_StandardOperators1D', 'Invalid basis argument')
        end if
      else
        this%basis = 'GLL'
      end if

    end if

    ! operators available from Gauss-Jacobi module .............................

    this%po = po

    if (po == 0) then

      ! 1-point Gauss
      allocate(this%x(0:0),       source = ZERO)
      allocate(this%w(0:0),       source = TWO)
      allocate(this%D(0:0,0:0),   source = ZERO)

    else if (this%basis == 'GL') then

      ! Gauss-Legendre
      allocate(this%x(0:po),      source = GL_Points(po))
      allocate(this%w(0:po),      source = GL_Weights(this%x))
      allocate(this%D(0:po,0:po), source = GL_DiffMatrix(this%x))

    else if (this%basis == 'GLL') then

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

    ! stiffness matrix .........................................................

    allocate(this%L(0:po,0:po), source = ZERO)

    do i = 0, po
    do j = 0, po
      this%L(i,j) = this%L(i,j) + sum(this%w * this%D(:,i) * this%D(:,j))
    end do
    end do

    ! Vandermonde matrix and its inverse .......................................

    if (present(no_vdm)) then
      build_vdm = .not. no_vdm
    else
      build_vdm = .true.
    end if

    if (present(svv) .and. svv) then
      build_vdm = .true.
    end if

    if (build_vdm) then
      call InitLegendreVDM(this)
    end if

    ! SVV differentation and stiffness matrix ..................................

    if (present(svv) .and. svv) then
      call InitSVV(this)
    end if

  end subroutine Init_StandardOperators1D

  !-----------------------------------------------------------------------------
  !> Intializes the Legendre-Vandermonde matrix and its inverse.

  subroutine InitLegendreVDM(this)
    class(StandardOperators1D), intent(inout) :: this !< standard operators

    integer :: i, j, po

    ! safeguard ................................................................

    if (allocated(this % VL    )) deallocate(this % VL    )
    if (allocated(this % VL_inv)) deallocate(this % VL_inv)

    ! prerequisites ............................................................

    po = this % po

    ! Vandermonde matrix .......................................................

    allocate(this % VL(0:po,0:po))
    do i = 0, po
    do j = 0, po
      this % VL(i,j) = JacobiPolynomial(a=ZERO, b=ZERO, n=j, x=this%x(i))
    end do
    end do

    ! inverse Vandermonde matrix ...............................................

    allocate(this % VL_inv(0:po,0:po), source=this%VL)
    this % VL_inv = Inverse(this % VL)

  end subroutine InitLegendreVDM

  !-----------------------------------------------------------------------------
  !> Query if Legendre VDM is available

  logical function HasLegendreVDM(this) result(has)
    class(StandardOperators1D), intent(in) :: this    !< standard operators
    has = allocated(this % VL)
  end function HasLegendreVDM

  !-----------------------------------------------------------------------------
  !> Get the Legendre-Vandermonde matrix

  subroutine GetLegendreVDM(this, VL)
    class(StandardOperators1D), intent(in) :: this    !< standard operators
    real(RNP), intent(out) :: VL(0:this%po,0:this%po) !< Vandermonde matrix

    if (.not. allocated(this % VL)) then
      call Error( 'GetLegendreVDM'                     &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Operators_1D'              )
    end if

    VL = this % VL

  end subroutine GetLegendreVDM

  !-----------------------------------------------------------------------------
  !> Get the inverse Legendre-Vandermonde matrix

  subroutine GetInverseLegendreVDM(this, VL_inv)
    class(StandardOperators1D), intent(in) :: this        !< standard operators
    real(RNP), intent(out) :: VL_inv(0:this%po,0:this%po) !< inverse VDM matrix

    if (.not. allocated(this%VL_inv)) then
      call Error( 'GetInverseLegendreVDM'              &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Operators_1D'              )
    end if

    VL_inv = this % VL_inv

  end subroutine GetInverseLegendreVDM

  !-----------------------------------------------------------------------------
  !> Initializes the SVV differentitation and stiffness matrix D_SVV and L_SVV

  subroutine InitSVV(this)
    class(StandardOperators1D), intent(inout) :: this !< standard operators

    real(RNP), allocatable :: Q_modal(:,:) ! SVV operator in spectral space
    real(RNP), allocatable :: Q(:,:)       ! SVV operator in nodal space

    real(RNP) :: epsilon_svv               ! svv amplitude
    real(RNP) :: viscosity                 ! mocked viscosity, will be removed
                                           ! in future versions as it will be
                                           ! given as a parameter
    integer   :: cutoff_mode_svv           ! mode up until no viscosity is applied
    integer   :: i, j, k, po

    ! safeguard ................................................................

    if (allocated(this % D_SVV)) deallocate(this % D_SVV)
    if (allocated(this % L_SVV)) deallocate(this % L_SVV)

    ! prerequisites ............................................................

    po = this % po

    ! according to Xu04
    epsilon_svv     = ONE / po
    viscosity       = 1.0 ! to be removed
    cutoff_mode_svv = floor(po / TWO)

    ! SVV operator in modal space ..............................................

    allocate(Q_modal(0:po,0:po), source = ZERO)

    ! computes the square root of the altered viscous prefactor, only non-unity
    ! for modes larger than the cutoff mode number M
    do k = cutoff_mode_svv+1, po
      Q_modal(k,k) = sqrt(ONE + epsilon_svv / viscosity *                      &
                        exp(-(real(po-k,RNP)/real(cutoff_mode_svv-k,RNP))**2))
    end do

    ! SVV operator in nodal space ..............................................

    allocate(Q(0:po,0:po), source = ZERO)

    ! Q is transfered into physical space utilizing the Vandermonde matrix as
    ! passage matrix
    Q = matmul(this%VL_inv, matmul(Q_modal, this%VL))

    ! SVV differentiation matrix ...............................................

    allocate(this%D_SVV(0:po,0:po), source = ZERO)

    ! application of the SVV operator on the differentiation matrix
    this%D_SVV = matmul(Q, this%D)

    ! SVV stiffness matrix .....................................................

    allocate(this%L_SVV(0:po,0:po), source = ZERO)

    do i = 0, po
    do j = 0, po
      this%L_SVV(i,j) = this%L_SVV(i,j) +                                      &
                            sum(this%w * this%D_SVV(:,i) * this%D_SVV(:,j))
    end do
    end do

  end subroutine InitSVV

  !-----------------------------------------------------------------------------
  !> Query if SVV is used

  logical function HasSVV(this) result(has)
    class(StandardOperators1D), intent(in) :: this !< standard operators
    has = allocated(this % D_SVV)
  end function HasSVV

  !-----------------------------------------------------------------------------
  !> Get the SVV differentiation matrix D_SVV

  subroutine GetSVV_DiffMatrix(this, D_SVV)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: D_SVV(0:this%po,0:this%po) !< SVV diff matrix

    if (.not. allocated(this % D_SVV)) then
      call Error( 'GetSVV_DiffMatrix'                          &
                , 'SVV differentiation matrix not initialized' &
                , 'Standard_Operators_1D'                      )
    end if

    D_SVV = this % D_SVV

  end subroutine GetSVV_DiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV stiffness matrix L_SVV

  subroutine GetSVV_StiffnessMatrix(this, L_SVV)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: L_SVV(0:this%po,0:this%po) !< SVV stiffness matrix

    if (.not. allocated(this % L_SVV)) then
      call Error( 'GetSVV_StiffnessMatrix'                     &
                , 'SVV differentiation matrix not initialized' &
                , 'Standard_Operators_1D'                      )
    end if

    L_SVV = this % L_SVV

  end subroutine GetSVV_StiffnessMatrix

  !-----------------------------------------------------------------------------
  !> Finalization

  subroutine Delete_StandardOperators1D(this)
    type(StandardOperators1D), intent(inout) :: this  !< standard operators

    if(allocated(this%x     )) deallocate(this%x     )
    if(allocated(this%w     )) deallocate(this%w     )
    if(allocated(this%D     )) deallocate(this%D     )
    if(allocated(this%L     )) deallocate(this%L     )
    if(allocated(this%VL    )) deallocate(this%VL    )
    if(allocated(this%VL_inv)) deallocate(this%VL_inv)
    if(allocated(this%D_SVV )) deallocate(this%D_SVV )
    if(allocated(this%L_SVV )) deallocate(this%L_SVV )

  end subroutine Delete_StandardOperators1D

  !=============================================================================

end module Standard_Operators_1D
