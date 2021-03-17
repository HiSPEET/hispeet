!> summary:   Spectral element operators in the one-dimensional standard region
!> author:    Immo Huismann, Joerg Stiller, Gustav Tschirschnitz
!> date:      2014/11/24
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Standard_Operators_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO
  use Execution_Control, only: Warning, Error
  use Gauss_Jacobi
  use Matrix_Operators,  only: Inverse
  use XMPI
  implicit none
  private

  public :: StandardOperators1D
  public :: StandardOperatorOptions1D

  !-----------------------------------------------------------------------------
  !> Standard operators for one polynomial order
  !>
  !> Supports the following types of nodal base functions
  !>
  !>   * Lagrange polynomials to Gauss-Legendre points:         `basis = 'GL'`
  !>   * Lagrange polynomials to Gauss-Lobatto-Legendre points: `basis = 'GLL'`
  !>   * Lagrange polynomials to Gauss-Radau-Legendre points:   `basis = 'GRL'`

  type StandardOperators1D
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

    procedure :: Init_Legendre_VDM
    procedure :: Has_Legendre_VDM
    procedure :: Get_Legendre_VDM
    procedure :: Get_Inverse_Legendre_VDM

    procedure :: Init_SVV
    procedure :: Has_SVV
    procedure :: Get_SVV_StandardRootDiffMatrix
    procedure :: Get_SVV_StandardDiffMatrix
    procedure :: Get_SVV_StandardStiffnessMatrix

  end type StandardOperators1D

  ! Constructor interface
  interface StandardOperators1D
    module procedure New_StandardOperators1D__f
    module procedure New_StandardOperators1D__b
  end interface

  !-----------------------------------------------------------------------------
  !> Options for StandardOperators1D

  type StandardOperatorOptions1D
    integer          :: po         = -1       !< polynomial order
    character(len=3) :: basis      = 'GLL'    !< basis type
    logical          :: no_vdm     = .false.  !< skip Vandermonde matrix
    logical          :: svv        = .false.  !< activate SVV model
    integer          :: po_cut_svv = -huge(1) !< cut-off PO for SVV
  contains
    procedure :: Bcast => Bcast_StandardOperatorOptions1D
  end type StandardOperatorOptions1D

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> Constructor for StandardOperators1D -- flat interface

  function New_StandardOperators1D__f(po, basis, no_vdm, svv, po_cut_svv)      &
     result(this)

    integer,                    intent(in) :: po    !< polynomial order
    character(len=*), optional, intent(in) :: basis !< points {GL,GLL,GRL} [GLL]

    logical, optional, intent(in) :: no_vdm     !< skip Vandermonde matrix [F]
    logical, optional, intent(in) :: svv        !< activate SVV model      [F]
    integer, optional, intent(in) :: po_cut_svv !< cut-off PO for SVV   [po/2]

    type(StandardOperators1D)       :: this
    type(StandardOperatorOptions1D) :: opt

    opt % po = po

    if (present(basis     )) opt % basis      = basis
    if (present(no_vdm    )) opt % no_vdm     = no_vdm
    if (present(svv       )) opt % svv        = svv
    if (present(po_cut_svv)) opt % po_cut_svv = po_cut_svv

    call Init_StandardOperators1D(this, opt)

  end function New_StandardOperators1D__f

  !-----------------------------------------------------------------------------
  !> Constructor for StandardOperators1D -- bundled arguments

  function New_StandardOperators1D__b(opt) result(this)
    class(StandardOperatorOptions1D), intent(in) :: opt
    type(StandardOperators1D) :: this

    call Init_StandardOperators1D(this, opt)

  end function New_StandardOperators1D__b

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

  subroutine Init_StandardOperators1D(this, opt)
    !> standard operators that will be initialized
    class(StandardOperators1D),       intent(inout) :: this
    !> options variable for the standard operator construction
    class(StandardOperatorOptions1D), intent(in)    :: opt

    integer :: i, j, po_cut

    associate( po          =>  opt % po          &
             , basis       =>  opt % basis       &
             , no_vdm      =>  opt % no_vdm      &
             , svv         =>  opt % svv         &
             , po_cut_svv  =>  opt % po_cut_svv  )

      ! safeguard ..............................................................

      call Delete_StandardOperators1D(this)

      if (po < 0) then
        call Error('Init_StandardOperators1D', 'Invalid order (po < 0)')
      else if (po == 0) then
        this%basis = 'GL'
      else if (any(basis == [ 'GL ', 'GLL', 'GRL' ])) then
        this%basis = basis
      else
        call Error('Init_StandardOperators1D', 'Invalid basis')
      end if

      ! operators available from Gauss-Jacobi module ...........................

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

      ! stiffness matrix .......................................................

      allocate(this%L(0:po,0:po), source = ZERO)

      do i = 0, po
      do j = 0, po
        this%L(i,j) = sum(this%w * this%D(:,i) * this%D(:,j))
      end do
      end do

      ! Vandermonde matrix and its inverse .....................................

      if (.not. no_vdm .or. svv) then
        call Init_Legendre_VDM(this)
      end if

      ! SVV differentiation and stiffness matrix ...............................

      if (svv) then
        if (po_cut_svv > -huge(1)) then
          po_cut = po_cut_svv
        else
          po_cut = floor(po / TWO) ! default value according to Xu04
        end if
        call Init_SVV(this, po_cut)
      end if

    end associate

  end subroutine Init_StandardOperators1D

  !-----------------------------------------------------------------------------
  !> Intializes the Legendre-Vandermonde matrix and its inverse.

  subroutine Init_Legendre_VDM(this)
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

  end subroutine Init_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Query if Legendre VDM is available

  logical function Has_Legendre_VDM(this) result(has)
    class(StandardOperators1D), intent(in) :: this    !< standard operators
    has = allocated(this % VL)
  end function Has_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Get the Legendre-Vandermonde matrix

  subroutine Get_Legendre_VDM(this, VL)
    class(StandardOperators1D), intent(in) :: this    !< standard operators
    real(RNP), intent(out) :: VL(0:this%po,0:this%po) !< Vandermonde matrix

    if (.not. allocated(this % VL)) then
      call Error( 'Get_Legendre_VDM'                   &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Operators_1D'              )
    end if

    VL = this % VL

  end subroutine Get_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Get the inverse Legendre-Vandermonde matrix

  subroutine Get_Inverse_Legendre_VDM(this, VL_inv)
    class(StandardOperators1D), intent(in) :: this        !< standard operators
    real(RNP), intent(out) :: VL_inv(0:this%po,0:this%po) !< inverse VDM matrix

    if (.not. allocated(this%VL_inv)) then
      call Error( 'Get_Inverse_Legendre_VDM'            &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Operators_1D'              )
    end if

    VL_inv = this % VL_inv

  end subroutine Get_Inverse_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Initializes SVV operators

  subroutine Init_SVV(this, po_cut)
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

  end subroutine Init_SVV

  !-----------------------------------------------------------------------------
  !> Query if SVV is used

  logical function Has_SVV(this) result(has)
    class(StandardOperators1D), intent(in) :: this !< standard operators
    has = allocated(this % D_root_svv)
  end function Has_SVV

  !-----------------------------------------------------------------------------
  !> Get the SVV standard differentiation matrix based on the root of SVV kernel

  subroutine Get_SVV_StandardRootDiffMatrix(this, D_root_svv)
    class(StandardOperators1D), intent(in) :: this            !< standard operators
    real(RNP), intent(out) :: D_root_svv(0:this%po,0:this%po) !< SVV diff matrix

    if (.not. allocated(this % D_root_svv)) then
      call Error( 'Get_SVV_StandardRootDiffMatrix' &
                , 'SVV not initialized'            &
                , 'Standard_Operators_1D'          )
    end if

    D_root_svv = this % D_root_svv

  end subroutine Get_SVV_StandardRootDiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV standard differentiation matrix

  subroutine Get_SVV_StandardDiffMatrix(this, D_svv)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: D_svv(0:this%po,0:this%po) !< SVV diff matrix

    if (.not. allocated(this % D_svv)) then
      call Error( 'Get_SVV_StandardDiffMatrix' &
                , 'SVV not initialized'        &
                , 'Standard_Operators_1D'      )
    end if

    D_svv = this % D_svv

  end subroutine Get_SVV_StandardDiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV standard stiffness matrix

  subroutine Get_SVV_StandardStiffnessMatrix(this, L_svv)
    class(StandardOperators1D), intent(in) :: this       !< standard operators
    real(RNP), intent(out) :: L_svv(0:this%po,0:this%po) !< SVV stiffness matrix

    if (.not. allocated(this % L_svv)) then
      call Error( 'Get_SVV_StandardStiffnessMatrix' &
                , 'SVV not initialized'             &
                , 'Standard_Operators_1D'           )
    end if

    L_svv = this % L_svv

  end subroutine Get_SVV_StandardStiffnessMatrix

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

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type StandardOperatorOptions1D

  subroutine Bcast_StandardOperatorOptions1D(this, root, comm)
    class(StandardOperatorOptions1D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast( this % po        , root, comm )
    call XMPI_Bcast( this % basis     , root, comm )
    call XMPI_Bcast( this % no_vdm    , root, comm )
    call XMPI_Bcast( this % svv       , root, comm )
    call XMPI_Bcast( this % po_cut_svv, root, comm )

  end subroutine Bcast_StandardOperatorOptions1D

  !=============================================================================

end module Standard_Operators_1D
