!> summary:   Spectral element operators in the one-dimensional standard region
!> author:    Immo Huismann, Joerg Stiller, Gustav Tschirschnitz
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
    character(len=3)      , public :: basis   !< basis type
    integer               , public :: po = -1 !< polynomial order
    real(RNP), allocatable, public :: x(:)    !< collocation points
    real(RNP), allocatable, public :: w(:)    !< quadrature weights
    real(RNP), allocatable, public :: D(:,:)  !< differention matrix
    real(RNP), allocatable, public :: L(:,:)  !< stiffness (Laplace) matrix

    ! private components
    real(RNP), allocatable :: VL(:,:)         !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: VL_inv(:,:)     !< inverse Legendre-Vandermonde
                                              !! matrix

    real(RNP), allocatable :: D_root_svv(:,:) !< SVV root-based diff matrix: √Q D
    real(RNP), allocatable :: D_svv(:,:)      !< SVV differentiation matrix:  Q D
    real(RNP), allocatable :: L_svv(:,:)      !< SVV stiffness matrix:
                                              !! (√Q D)ᵀ M (√Q D)

  contains

    procedure :: Init_StandardOperators1D
    procedure :: PolynomialOrder

    procedure :: InitLegendreVDM
    procedure :: HasLegendreVDM
    procedure :: GetLegendreVDM
    procedure :: GetInverseLegendreVDM

    procedure :: InitSVV
    procedure :: HasSVV
    procedure :: GetSVV_RootDiffMatrix
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

  function New_StandardOperators1D(po, basis, no_vdm, svv, po_cut_svv)         &
     result(this)

    integer,                    intent(in) :: po    !< polynomial order
    character(len=*), optional, intent(in) :: basis !< points {GL,GLL,GRL} [GLL]

    logical, optional, intent(in) :: no_vdm     !< skip Vandermonde matrix [F]
    logical, optional, intent(in) :: svv        !< activate SVV model [F]
    integer, optional, intent(in) :: po_cut_svv !< cut-off PO for SVV [po/2]

    type(StandardOperators1D) :: this

    call Init_StandardOperators1D(this, po, basis, no_vdm, svv, po_cut_svv)

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

  subroutine Init_StandardOperators1D(this, po, basis, no_vdm, svv, po_cut_svv)

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

    !> cut-off polynomial degree for the SVV model
    integer, optional, intent(in) :: po_cut_svv

    logical :: build_vdm
    logical :: build_svv
    integer :: i, j, po_cut

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
      this%L(i,j) = sum(this%w * this%D(:,i) * this%D(:,j))
    end do
    end do

    ! check what else to build .................................................

    ! Vandermonde matrix
    if (present(no_vdm)) then
      build_vdm = .not. no_vdm
    else
      build_vdm = .true.
    end if

    ! SVV
    if (present(svv)) then
      build_svv = svv
    else
      build_svv = .false.
    end if

    ! require Vandermonde matrix when SVV is activated
    build_vdm = build_vdm .or. build_svv

    ! Vandermonde matrix and its inverse .......................................

    if (build_vdm) then
      call InitLegendreVDM(this)
    end if

    ! SVV differentation and stiffness matrix ..................................

    po_cut = floor(po / TWO) ! default value according to Xu04
    if (build_svv) then
      ! only if po_cut_svv is present and not the default value use it
      if (present(po_cut_svv)) then
        if (po_cut_svv /= -huge(1)) po_cut = po_cut_svv
      end if
      call InitSVV(this, po_cut)
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
  !> Initializes SVV operators

  subroutine InitSVV(this, po_cut)
    class(StandardOperators1D), intent(inout) :: this   !< standard operators
    integer,                    intent(in)    :: po_cut !< cut-off polynomial degree

    real(RNP), allocatable :: Q_hat(:,:) ! SVV filter coefficients

    integer :: i, j, k, po, po_init

    ! safeguard ................................................................

    if (allocated(this % D_root_svv)) deallocate(this % D_root_svv)
    if (allocated(this % D_svv     )) deallocate(this % D_svv     )
    if (allocated(this % L_svv     )) deallocate(this % L_svv     )

    ! prerequisites ............................................................

    po      = this % po
    po_init = abs(max(po_cut+1,0)) ! lowest order to which the filter is applied

    ! SVV filter coeffients ....................................................

    allocate(Q_hat(0:po,0:po), source = ZERO)

    ! computes the SVV filter coeffients based on the SVV kernel
    do k = po_init, po
      Q_hat(k,k) = exp(-(real(po-k,RNP) / real(po_cut-k,RNP))**2)
    end do

    ! SVV operators ............................................................

    allocate(this%D_root_svv (0:po,0:po), source = ZERO)
    allocate(this%D_svv      (0:po,0:po), source = ZERO)
    allocate(this%L_svv      (0:po,0:po), source = ZERO)

    associate( D          => this % D           &
             , VL         => this % VL          &
             , VL_inv     => this % VL_inv      &
             , D_root_svv => this % D_root_svv  &
             , D_svv      => this % D_svv       &
             , L_svv      => this % L_svv       )

      ! diff matrices
      D_root_svv = matmul( matmul( matmul( VL, sqrt(Q_hat) ), VL_inv ), D )
      D_svv      = matmul( matmul( matmul( VL,      Q_hat  ), VL_inv ), D )

      ! stiffness matrix
      do i = 0, po
      do j = 0, po
        L_svv(i,j) = sum( this%w * D_root_svv(:,i) * D_root_svv(:,j) )
      end do
      end do

    end associate

  end subroutine InitSVV

  !-----------------------------------------------------------------------------
  !> Query if SVV is used

  logical function HasSVV(this) result(has)
    class(StandardOperators1D), intent(in) :: this !< standard operators
    has = allocated(this % D_root_svv)
  end function HasSVV

  !-----------------------------------------------------------------------------
  !> Get the SVV root-based differentiation matrix `D_root_svv`

  subroutine GetSVV_RootDiffMatrix(this, D_root_svv)
    class(StandardOperators1D), intent(in) :: this            !< standard operators
    real(RNP), intent(out) :: D_root_svv(0:this%po,0:this%po) !< SVV diff matrix

    if (.not. allocated(this % D_root_svv)) then
      call Error( 'GetSVV_RootDiffMatrix'                                 &
                , 'SVV root-based differentiation matrix not initialized' &
                , 'Standard_Operators_1D'                                 )
    end if

    D_root_svv = this % D_root_svv

  end subroutine GetSVV_RootDiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV differentiation matrix `root-based`

  subroutine GetSVV_DiffMatrix(this, D_svv)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: D_svv(0:this%po,0:this%po) !< SVV diff matrix

    if (.not. allocated(this % D_svv)) then
      call Error( 'GetSVV_DiffMatrix'                          &
                , 'SVV differentiation matrix not initialized' &
                , 'Standard_Operators_1D'                      )
    end if

    D_svv = this % D_svv

  end subroutine GetSVV_DiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV stiffness matrix `L_svv`

  subroutine GetSVV_StiffnessMatrix(this, L_svv)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: L_svv(0:this%po,0:this%po) !< SVV stiffness matrix

    if (.not. allocated(this % L_svv)) then
      call Error( 'GetSVV_StiffnessMatrix'                     &
                , 'SVV differentiation matrix not initialized' &
                , 'Standard_Operators_1D'                      )
    end if

    L_svv = this % L_svv

  end subroutine GetSVV_StiffnessMatrix

  !-----------------------------------------------------------------------------
  !> Finalization

  subroutine Delete_StandardOperators1D(this)
    type(StandardOperators1D), intent(inout) :: this  !< standard operators

    if(allocated(this%x          )) deallocate(this%x          )
    if(allocated(this%w          )) deallocate(this%w          )
    if(allocated(this%D          )) deallocate(this%D          )
    if(allocated(this%L          )) deallocate(this%L          )
    if(allocated(this%VL         )) deallocate(this%VL         )
    if(allocated(this%VL_inv     )) deallocate(this%VL_inv     )
    if(allocated(this%D_root_svv )) deallocate(this%D_root_svv )
    if(allocated(this%D_svv      )) deallocate(this%D_svv      )
    if(allocated(this%L_svv      )) deallocate(this%L_svv      )

  end subroutine Delete_StandardOperators1D

  !=============================================================================

end module Standard_Operators_1D
