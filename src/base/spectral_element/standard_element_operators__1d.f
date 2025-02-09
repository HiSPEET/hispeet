!> summary:   Spectral element operators in the one-dimensional standard region
!> author:    Immo Huismann, Joerg Stiller, Gustav Tschirschnitz
!> date:      2014/11/24
!> license:   Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Standard_Element_Operators__1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, TWO
  use Execution_Control, only: Warning, Error
  use Gauss_Jacobi
  use Matrix_Operators,  only: Inverse
  use XMPI
  implicit none
  private

  public :: StandardElementOperators_1D
  public :: StandardElementOptions_1D

  !-----------------------------------------------------------------------------
  !> Standard operators for one polynomial order
  !>
  !> Supports the following types of nodal base functions
  !>
  !>   - Lagrange polynomials to Gauss-Legendre nodes:             nodes = 'G'
  !>   - Lagrange polynomials to Gauss-Lobatto-Legendre nodes:     nodes = 'L'
  !>   - Lagrange polynomials to left  Gauss-Radau-Legendre nodes: nodes = 'RL'
  !>   - Lagrange polynomials to right Gauss-Radau-Legendre nodes: nodes = 'RR'

  type StandardElementOperators_1D
    private

    ! public components
    character(len=2)      , public :: nodes   !< node set
    integer               , public :: po = -1 !< polynomial order
    real(RNP), allocatable, public :: x(:)    !< collocation points (nodes)
    real(RNP), allocatable, public :: w(:)    !< quadrature weights
    real(RNP), allocatable, public :: D(:,:)  !< differention matrix
    real(RNP), allocatable, public :: L(:,:)  !< stiffness (Laplace) matrix

    ! private components
    real(RNP), allocatable :: VL(:,:)     !< Legendre-Vandermonde matrix
    real(RNP), allocatable :: VL_inv(:,:) !< inverse Legendre-Vandermonde matrix
    real(RNP), allocatable :: VB(:,:)     !< Bubble-Vandermonde matrix
    real(RNP), allocatable :: VB_inv(:,:) !< inverse Bubble-Vandermonde matrix

    real(RNP), allocatable :: D_root_svv(:,:) !< SVV root-based diff matrix: √Q D
    real(RNP), allocatable :: D_svv(:,:)      !< SVV differentiation matrix:  Q D
    real(RNP), allocatable :: L_svv(:,:)      !< SVV stiffness matrix:
                                              !! (√Q D)ᵀ M (√Q D)

  contains

    procedure :: Init_StandardElementOperators_1D
    procedure :: PolynomialOrder

    procedure :: Init_Legendre_VDM
    procedure :: Has_Legendre_VDM
    procedure :: Get_Legendre_VDM
    procedure :: Get_Inverse_Legendre_VDM
    procedure :: Get_Legendre_CutoffFilter

    procedure :: Init_Bubble_VDM
    procedure :: Has_Bubble_VDM
    procedure :: Get_Bubble_VDM
    procedure :: Get_Inverse_Bubble_VDM
    procedure :: Get_Bubble_CutoffFilter

    procedure :: Init_SVV
    procedure :: Has_SVV
    procedure :: Get_SVV_StandardRootDiffMatrix
    procedure :: Get_SVV_StandardDiffMatrix
    procedure :: Get_SVV_StandardStiffnessMatrix

  end type StandardElementOperators_1D

  ! Constructor interface
  interface StandardElementOperators_1D
    module procedure New_StandardElementOperators_1D__f
    module procedure New_StandardElementOperators_1D__b
  end interface

  !-----------------------------------------------------------------------------
  !> Options for StandardElementOperators_1D

  type StandardElementOptions_1D
    integer      :: po         = -1       !< polynomial order
    character(2) :: nodes      = 'L'      !< node set
    logical      :: no_vdm     = .false.  !< skip Vandermonde matrix
    logical      :: svv        = .false.  !< activate SVV model
    integer      :: po_cut_svv = -huge(1) !< cut-off order for SVV
  contains
    procedure :: Bcast => Bcast_StandardOperatorOptions1D
  end type StandardElementOptions_1D

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> Constructor for StandardElementOperators_1D -- flat interface

  function New_StandardElementOperators_1D__f &
              (po, nodes, no_vdm, svv, po_cut_svv) result(this)

    integer,                intent(in) :: po     !< polynomial order
    character(*), optional, intent(in) :: nodes  !< nodes {G,L,RR,RL}        [L]
    logical     , optional, intent(in) :: no_vdm !< skip Vandermonde matrix  [F]
    logical     , optional, intent(in) :: svv    !< activate SVV model       [F]
    integer     , optional, intent(in) :: po_cut_svv !< SVV cut-off order [po/2]

    type(StandardElementOperators_1D) :: this
    type(StandardElementOptions_1D)   :: opt

    opt % po = po

    if (present(nodes     )) opt % nodes      = nodes
    if (present(no_vdm    )) opt % no_vdm     = no_vdm
    if (present(svv       )) opt % svv        = svv
    if (present(po_cut_svv)) opt % po_cut_svv = po_cut_svv

    call Init_StandardElementOperators_1D(this, opt)

  end function New_StandardElementOperators_1D__f

  !-----------------------------------------------------------------------------
  !> Constructor for StandardElementOperators_1D -- bundled arguments

  function New_StandardElementOperators_1D__b(opt) result(this)
    class(StandardElementOptions_1D), intent(in) :: opt
    type(StandardElementOperators_1D) :: this

    call Init_StandardElementOperators_1D(this, opt)

  end function New_StandardElementOperators_1D__b

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Returns the polymial order of operators, or -1 if none

  pure integer function PolynomialOrder(this) result(po)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators

    po = this % po

  end function PolynomialOrder

  !-----------------------------------------------------------------------------
  !> Intializes the public components.
  !>
  !> The routine provides 1D standard operators for the chosen nodal basis.
  !> Lobatto is the default, except for `po = 0` which always implies Gauss.

  subroutine Init_StandardElementOperators_1D(this, opt)
    class(StandardElementOperators_1D), intent(inout) :: this
      !< standard operators that will be initialized
    class(StandardElementOptions_1D), intent(in) :: opt
      !< options variable for the standard operator construction

    integer :: i, j, po_cut

    associate( po          =>  opt % po          &
             , nodes       =>  opt % nodes       &
             , no_vdm      =>  opt % no_vdm      &
             , svv         =>  opt % svv         &
             , po_cut_svv  =>  opt % po_cut_svv  )

      ! safeguard ..............................................................

      call Delete_StandardOperators1D(this)

      if (po < 0) then
        call Error('Init_StandardElementOperators_1D', 'Invalid order (po < 0)')
      else if (po == 0) then
        this%nodes = 'G'
      else
        select case(nodes)
        case('G','L','RL','RR')
          this%nodes = nodes
        case default
          call Error('Init_StandardElementOperators_1D', 'Invalid nodes')
        end select
      end if

      ! operators available from Gauss-Jacobi module ...........................

      this%po = po

      if (po == 0) then

        ! 1-point Gauss
        allocate(this%x(0:0),       source = ZERO)
        allocate(this%w(0:0),       source = TWO)
        allocate(this%D(0:0,0:0),   source = ZERO)

      else
        select case(this%nodes)

        case('G')
          ! Gauss-Legendre
          allocate(this%x(0:po),      source = GaussPoints(po))
          allocate(this%w(0:po),      source = GaussWeights(this%x))
          allocate(this%D(0:po,0:po), source = GaussDiffMatrix(this%x))

        case('L')
          ! Gauss-Lobatto-Legendre
          allocate(this%x(0:po),      source = LobattoPoints(po))
          allocate(this%w(0:po),      source = LobattoWeights(this%x))
          allocate(this%D(0:po,0:po), source = LobattoDiffMatrix(this%x))

        case('RL')
          ! left-sided Gauss-Radau-Legendre
          allocate(this%x(0:po),      source = RadauPoints(po, right = .false.))
          allocate(this%w(0:po),      source = RadauWeights(this%x))
          allocate(this%D(0:po,0:po), source = RadauDiffMatrix(this%x))

        case('RR')
          ! right-sided Gauss-Radau-Legendre
          allocate(this%x(0:po),      source = RadauPoints(po, right = .true.))
          allocate(this%w(0:po),      source = RadauWeights(this%x))
          allocate(this%D(0:po,0:po), source = RadauDiffMatrix(this%x))

        end select
      end if

      ! stiffness matrix .......................................................

      allocate(this%L(0:po,0:po), source = ZERO)

      do i = 0, po
      do j = 0, po
        this%L(i,j) = sum(this%w * this%D(:,i) * this%D(:,j))
      end do
      end do

      ! Vandermonde matrices and their inverse .................................

      if (.not. no_vdm .or. svv) then
        call Init_Legendre_VDM(this)
        call Init_Bubble_VDM(this)
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

  end subroutine Init_StandardElementOperators_1D

  !-----------------------------------------------------------------------------
  !> Intializes the Legendre-Vandermonde matrix and its inverse.

  subroutine Init_Legendre_VDM(this)
    class(StandardElementOperators_1D), intent(inout) :: this
      !< standard operators

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
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    has = allocated(this % VL)
  end function Has_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Get the Legendre-Vandermonde matrix

  subroutine Get_Legendre_VDM(this, VL)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: VL(0:this%po,0:this%po)      !< Vandermonde matrix

    if (.not. allocated(this % VL)) then
      call Error( 'Get_Legendre_VDM'                   &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Element_Operators__1D'     )
    end if

    VL = this % VL

  end subroutine Get_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Get the inverse Legendre-Vandermonde matrix

  subroutine Get_Inverse_Legendre_VDM(this, VL_inv)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: VL_inv(0:this%po,0:this%po)  !< inverse VDM matrix

    if (.not. allocated(this%VL_inv)) then
      call Error( 'Get_Inverse_Legendre_VDM'            &
                , 'Vandermonde matrix not initialized'  &
                , 'Standard_Element_Operators__1D'      )
    end if

    VL_inv = this % VL_inv

  end subroutine Get_Inverse_Legendre_VDM

  !-----------------------------------------------------------------------------
  !> Get a filter matrix based on spectral cutoff in Legendre space
  !>
  !> The returned filter matrix removes all Legendre modes with a degree greater
  !> than `po_cut`.

  subroutine Get_Legendre_CutoffFilter(this, po_cut, A)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    integer,   intent(in)  :: po_cut                 !< cutoff degree
    real(RNP), intent(out) :: A(0:this%po,0:this%po) !< filter matrix

    if (.not. this%Has_Legendre_VDM()) then
      call Error( 'Get_Legendre_CutoffFilter'           &
                , 'Vandermonde matrix not initialized'  &
                , 'Standard_Element_Operators__1D'              )
    end if

    if (po_cut >= 0 .and. po_cut <= this%po) then
      A = matmul(this % VL(:,0:po_cut), this % VL_inv(0:po_cut,:))
    else
      A = ZERO
    end if

  end subroutine Get_Legendre_CutoffFilter

  !-----------------------------------------------------------------------------
  !> Intializes the Vandermonde matrix for basis with linear and bubble modes.

  subroutine Init_Bubble_VDM(this)
    class(StandardElementOperators_1D), intent(inout) :: this
      !< standard operators

    real(RNP) :: x
    integer   :: i, j, po

    ! safeguard ................................................................

    if (allocated(this % VB    )) deallocate(this % VB    )
    if (allocated(this % VB_inv)) deallocate(this % VB_inv)

    ! prerequisites ............................................................

    po = this % po

    ! Vandermonde matrix .......................................................

    allocate(this % VB(0:po,0:po))
    associate(VB => this % VB)
      do i = 0, po
        x = this % x(i)
        VB(i,0) = (1 - x) / 2
        VB(i,1) = (1 + x) / 2
        do j = 2, po
          VB(i,j) = ((1-x)*(1+x)/4) * JacobiPolynomial(a=ONE, b=ONE, n=j-2, x=x)
        end do
      end do
    end associate

    ! inverse Vandermonde matrix ...............................................

    allocate(this % VB_inv(0:po,0:po), source=this%VB)
    this % VB_inv = Inverse(this % VB)

  end subroutine Init_Bubble_VDM

  !-----------------------------------------------------------------------------
  !> Query if Bubble VDM is available

  logical function Has_Bubble_VDM(this) result(has)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    has = allocated(this % VB)
  end function Has_Bubble_VDM

  !-----------------------------------------------------------------------------
  !> Get the Bubble-Vandermonde matrix

  subroutine Get_Bubble_VDM(this, VB)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: VB(0:this%po,0:this%po)      !< Vandermonde matrix

    if (.not. allocated(this % VB)) then
      call Error( 'Get_Bubble_VDM'                     &
                , 'Vandermonde matrix not initialized' &
                , 'Standard_Element_Operators__1D'             )
    end if

    VB = this % VB

  end subroutine Get_Bubble_VDM

  !-----------------------------------------------------------------------------
  !> Get the inverse Bubble-Vandermonde matrix

  subroutine Get_Inverse_Bubble_VDM(this, VB_inv)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: VB_inv(0:this%po,0:this%po ) !< inverse VDM matrix

    if (.not. allocated(this%VB_inv)) then
      call Error( 'Get_Inverse_Bubble_VDM'              &
                , 'Vandermonde matrix not initialized'  &
                , 'Standard_Element_Operators__1D'              )
    end if

    VB_inv = this % VB_inv

  end subroutine Get_Inverse_Bubble_VDM

  !-----------------------------------------------------------------------------
  !> Get a filter matrix based on spectral cutoff of bubble modes
  !>
  !> The returned filter matrix removes all modes with a degree greater than
  !> `po_cut`.

  subroutine Get_Bubble_CutoffFilter(this, po_cut, A)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    integer,   intent(in)  :: po_cut                 !< cutoff degree
    real(RNP), intent(out) :: A(0:this%po,0:this%po) !< filter matrix

    if (.not. this%Has_Bubble_VDM()) then
      call Error( 'Get_Bubble_CutoffFilter'             &
                , 'Vandermonde matrix not initialized'  &
                , 'Standard_Element_Operators__1D'              )
    end if

    if (po_cut >= 0 .and. po_cut <= this%po) then
      A = matmul(this % VB(:,0:po_cut), this % VB_inv(0:po_cut,:))
    else
      A = ZERO
    end if

  end subroutine Get_Bubble_CutoffFilter

  !-----------------------------------------------------------------------------
  !> Initializes SVV operators

  subroutine Init_SVV(this, po_cut)
    class(StandardElementOperators_1D), intent(inout) :: this
      !< standard operators
    integer, intent(in) :: po_cut
      !< cutoff polynomial degree

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
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    has = allocated(this % D_root_svv)
  end function Has_SVV

  !-----------------------------------------------------------------------------
  !> Get the SVV standard differentiation matrix based on the root of SVV kernel

  subroutine Get_SVV_StandardRootDiffMatrix(this, D_root_svv)
    class(StandardElementOperators_1D), intent(in) :: this
      !< standard operators
    real(RNP), intent(out) :: D_root_svv(0:this%po,0:this%po)
      !< SVV diff matrix

    if (.not. allocated(this % D_root_svv)) then
      call Error( 'Get_SVV_StandardRootDiffMatrix' &
                , 'SVV not initialized'            &
                , 'Standard_Element_Operators__1D'         )
    end if

    D_root_svv = this % D_root_svv

  end subroutine Get_SVV_StandardRootDiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV standard differentiation matrix

  subroutine Get_SVV_StandardDiffMatrix(this, D_svv)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: D_svv(0:this%po,0:this%po)   !< SVV diff matrix

    if (.not. allocated(this % D_svv)) then
      call Error( 'Get_SVV_StandardDiffMatrix' &
                , 'SVV not initialized'        &
                , 'Standard_Element_Operators__1D'     )
    end if

    D_svv = this % D_svv

  end subroutine Get_SVV_StandardDiffMatrix

  !-----------------------------------------------------------------------------
  !> Get the SVV standard stiffness matrix

  subroutine Get_SVV_StandardStiffnessMatrix(this, L_svv)
    class(StandardElementOperators_1D), intent(in) :: this !< standard operators
    real(RNP), intent(out) :: L_svv(0:this%po,0:this%po) !< SVV stiffness matrix

    if (.not. allocated(this % L_svv)) then
      call Error( 'Get_SVV_StandardStiffnessMatrix' &
                , 'SVV not initialized'             &
                , 'Standard_Element_Operators__1D'          )
    end if

    L_svv = this % L_svv

  end subroutine Get_SVV_StandardStiffnessMatrix

  !-----------------------------------------------------------------------------
  !> Finalization

  subroutine Delete_StandardOperators1D(this)
    type(StandardElementOperators_1D), intent(inout) :: this
      !< standard operators

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
  !> MPI_Bcast for objects of type StandardElementOptions_1D

  subroutine Bcast_StandardOperatorOptions1D(this, root, comm)
    class(StandardElementOptions_1D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast( this % po        , root, comm )
    call XMPI_Bcast( this % nodes     , root, comm )
    call XMPI_Bcast( this % no_vdm    , root, comm )
    call XMPI_Bcast( this % svv       , root, comm )
    call XMPI_Bcast( this % po_cut_svv, root, comm )

  end subroutine Bcast_StandardOperatorOptions1D

  !=============================================================================

end module Standard_Element_Operators__1D
