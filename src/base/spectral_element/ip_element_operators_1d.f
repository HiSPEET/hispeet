!> summary:  Element operators for IP/DG-SEM
!> author:   Joerg Stiller
!> date:     2016/03/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Element operators for IP/DG-SEM
!===============================================================================

module IP_Element_Operators_1D
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  use XMPI
  implicit none
  private

  public :: IP_ElementOperators1D
  public :: IP_ElementOptions1D

  !-----------------------------------------------------------------------------
  !> Element operators for the symmetric interior penalty IP/DG-SEM
  !>
  !> The operator accommodates classical IP as well as hybridizable IP-H. The
  !> particular method is selected during initialization and controlled by then
  !> component `hybrid`.
  !>
  !> In case of IP-H, the operator provides  the column matrix of generalized
  !> eigenvectors `S` and the diagonal matrix of eigenvalues `Λ = Lambda` such
  !> that
  !>
  !>     Sᵀ Lᵢᵢ S = Λ
  !>     Sᵀ Mᵢᵢ S = I
  !>
  !> where `Lᵢᵢ` and `Mᵢᵢ` the standard stiffness matrix and the standard
  !> diagonal mass matrix restricted to the interior points.

  type, extends(StandardOperators1D) :: IP_ElementOperators1D
    real(RNP) :: penalty = 2       !< penalty parameter > 1
    logical   :: hybrid  = .false. !< switch to hybridized method
  contains

    procedure :: Init_IP_ElementOperators1D

    generic :: PenaltyFactor => PenaltyFactor_NE, PenaltyFactor_EQ
    procedure, private :: PenaltyFactor_NE
    procedure, private :: PenaltyFactor_EQ

    procedure :: Get_DiffusionMatrix

    generic   :: Get_EllipticEigensystem => Get_EllipticEigensystem__w_svv, &
                                            Get_EllipticEigensystem__n_svv

    procedure, private :: Get_EllipticEigensystem__w_svv
    procedure, private :: Get_EllipticEigensystem__n_svv

    generic   :: Get_EllipticSuboperators => Get_EllipticSuboperators__w_svv, &
                                             Get_EllipticSuboperators__n_svv

    procedure, private :: Get_EllipticSuboperators__w_svv
    procedure, private :: Get_EllipticSuboperators__n_svv

  end type IP_ElementOperators1D

  ! constructor interface
  interface IP_ElementOperators1D
    module procedure New_IP_ElementOperators1D__f
    module procedure New_IP_ElementOperators1D__b
  end interface

  !-----------------------------------------------------------------------------
  !> Options for IP_ElementOperators1D

  type, extends(StandardOperatorOptions1D) :: IP_ElementOptions1D
    real(RNP) :: penalty    =  2       !< penalty parameter > 1
    logical   :: hybrid     = .false.  !< switch to hybridized method
  contains
    procedure :: Bcast => IP_ElementOptions1D_Bcast
  end type IP_ElementOptions1D

  ! constructor interface
  interface IP_ElementOptions1D
    module procedure New_IP_ElementOptions1D_o
  end interface

contains

!===============================================================================
! Constructors

!-------------------------------------------------------------------------------
!> Constructor for IP_ElementOperators1D -- flat interface

function New_IP_ElementOperators1D__f( po, penalty, hybrid     &
                                     , no_vdm, svv, po_cut_svv ) result(this)
  integer,             intent(in) :: po         !< polynomial order
  real(RNP), optional, intent(in) :: penalty    !< penalty parameter > 1     [2]
  logical,   optional, intent(in) :: hybrid     !< switch to hyb. method     [F]
  logical,   optional, intent(in) :: no_vdm     !< skip Vandermonde matrix   [F]
  logical,   optional, intent(in) :: svv        !< activate SVV model        [F]
  integer,   optional, intent(in) :: po_cut_svv !< cutoff PO for SVV model  [-∞]

  type(IP_ElementOperators1D) :: this
  type(IP_ElementOptions1D)   :: opt

  opt % po = po
  if (present(penalty   )) opt % penalty    = penalty
  if (present(hybrid    )) opt % hybrid     = hybrid
  if (present(no_vdm    )) opt % no_vdm     = no_vdm
  if (present(svv       )) opt % svv        = svv
  if (present(po_cut_svv)) opt % po_cut_svv = po_cut_svv

  call Init_IP_ElementOperators1D(this, opt)

end function New_IP_ElementOperators1D__f

!-------------------------------------------------------------------------------
!> Constructor for IP_ElementOperators1D -- bundled arguments

function New_IP_ElementOperators1D__b(opt) result(this)
  type(IP_ElementOptions1D), intent(in) :: opt
  type(IP_ElementOperators1D) :: this

  call Init_IP_ElementOperators1D(this, opt)

end function New_IP_ElementOperators1D__b

!===============================================================================
! IP_ElementOperators1D type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization

subroutine Init_IP_ElementOperators1D(this, opt)
  class(IP_ElementOperators1D), intent(inout) :: this
  class(IP_ElementOptions1D),   intent(in)    :: opt

  call this % Init_StandardOperators1D(opt)

  this % penalty = opt % penalty
  this % hybrid  = opt % hybrid

end subroutine Init_IP_ElementOperators1D

!-------------------------------------------------------------------------------
!> Penalty factor for non-equidistant spacing

real(RNP) function PenaltyFactor_NE(this, dx) result(mu)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in) :: dx(2)  !< element extensions

  mu = this%penalty/4 * this%po * (this%po + 1) * (1/dx(1) + 1/dx(2))

end function PenaltyFactor_NE

!-------------------------------------------------------------------------------
!> Penalty factor for equidistant spacing

real(RNP) function PenaltyFactor_EQ(this, dx) result(mu)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in) :: dx  !< element extension

  mu = this%penalty/4 * this%po * (this%po + 1) * 2/dx

end function PenaltyFactor_EQ

!-------------------------------------------------------------------------------
!> Returns the 1D element diffusion matrix for the interior penalty DGM, i.e.
!> the product of diffusivity and the element stiffness matrix.
!>
!>  C = ν L + νˢ Lˢ
!>
!> Here the element diffusion matrix `Ce` represents the nontrivial row entries
!> of the global diffusion matrix corresponding to the given element. It
!> must be dimensioned as `Ce(0:P,0:P,-1:1)`, where `P = this%po` is the
!> polynomial order. The third index refers to the
!>
!>   * preceding (-1),
!>   * current (0) and
!>   * succeeding (1) elements,
!>
!> with lengths `dx(-1:1)è respectively.
!>
!> The diffusion matrix is available in two forms
!>
!>   * `'primal'`: all numeric fluxes û are eliminated (default)
!>   * `'flux'`  : û is retained, all corresponding terms are removed from `Ce`
!>
!> Except for the single element case, i.e. `all(bc /= '')`, the flux form can
!> be activated by passing `form = 'flux'`.

subroutine Get_DiffusionMatrix(this, dx, bc, nu, nu_svv, Ce, form)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
  character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N','P'}
  real(RNP), intent(in)  :: nu            !< diffusivity
  real(RNP), intent(in)  :: nu_svv        !< SVV diffusivity [0]
  real(RNP), intent(out) :: Ce(0:,0:,-1:) !< regular stiffness matrix
  character(len=*), optional, intent(in) :: form !< operator form ['primal']

  logical   :: primal
  integer   :: P
  real(RNP) :: g(-1:1), mu_0, mu_P, c_0, c_P, h_0, h_P
  real(RNP) :: Cs(0:this%po, 0:this%po), Fs(0:this%po, 0:this%po) ! nu * D
  real(RNP), allocatable :: delta_0(:), delta_P(:)

  ! initialization .............................................................

  P = this % po
  g = ONE / dx

  mu_0 = this % PenaltyFactor(dx(-1:0))
  mu_P = this % PenaltyFactor(dx( 0:1))

  allocate(delta_0(0:P), source = ZERO)
  delta_0(0) = ONE

  allocate(delta_P(0:P), source = ZERO)
  delta_P(P) = ONE

  if (present(form)) then
    primal = form == 'primal' .or. all(bc /= ' ')
  else
    primal = .true.
  end if

  if (primal) then

    ! left boundary condition
    select case(bc(1))
    case('D')    ! Dirichlet
      c_0 = 2
    case('N')    ! Neumann
      c_0 = 0
    case default ! none
      c_0 = 1
    end select

    ! right boundary condition
    select case(bc(2))
    case('D')    ! Dirichlet
      c_P = 2
    case('N')    ! Neumann
      c_P = 0
    case default ! none
      c_P = 1
    end select

  else

    ! left boundary condition
    select case(bc(1))
    case('N')    ! Neumann
      c_0 = 0
    case default
      c_0 = 2
    end select

    ! right boundary condition
    select case(bc(2))
    case('N')    ! Neumann
      c_P = 0
    case default
      c_P = 2
    end select

  end if

  if (this % Has_SVV()) then
    call this % Get_SVV_StandardStiffnessMatrix(Cs)
    call this % Get_SVV_StandardDiffMatrix(Fs)
    Cs = nu_svv * Cs
    Fs = nu_svv * Fs
  else
    Cs = 0
    Fs = 0
  end if

  Cs = Cs + nu * this%L
  Fs = Fs + nu * this%D

  Ce = DiffusionMatrix(Cs, Fs)

contains

  function DiffusionMatrix(Cs, Fs) result (Ce)
    real(RNP), intent(in) :: Cs(0:P,0:P) !< 1D standard diffusion operator
    real(RNP), intent(in) :: Fs(0:P,0:P) !< 1D standard "flux" operator
    real(RNP) :: Ce(0:P,0:P,-1:1)
    integer   :: i, j

    ! contribution from preceding element (Le⁻) ................................

    if (scan(bc(1), 'DN') > 0) then
      Ce(:,:,-1) = 0
    else
      do j = 0, P
      do i = 0, P
        Ce(i,j,-1) = - g( 0) * Fs   (0,i) * delta_P(j)  &
                     + g(-1) * delta_0(i) * Fs   (P,j)  &
                     - mu_0  * delta_0(i) * delta_P(j)
      end do
      end do

      if (primal .and. this%hybrid) then
        h_0 = 1 / (dx(-1) * dx(0) * mu_0)
        do j = 0, P
        do i = 0, P
          Ce(i,j,-1) = Ce(i,j,-1) + h_0 * Fs(0,i) * Fs(P,j)
        end do
        end do
      end if
    end if

    ! own contribution (Le⁰) ...................................................

    do j = 0, P
    do i = 0, P
      Ce(i,j,0) = 2 * g(0) * Cs(i,j)                          &

                + c_0 * (   g(0) * Fs   (0,i) * delta_0(j)    &
                          + g(0) * delta_0(i) * Fs   (0,j)    &
                          + mu_0 * delta_0(i) * delta_0(j) )  &

                + c_P * ( - g(0) * Fs   (P,i) * delta_P(j)    &
                          - g(0) * delta_P(i) * Fs   (P,j)    &
                          + mu_P * delta_P(i) * delta_P(j) )
    end do
    end do

    if (primal .and. this%hybrid) then
      select case(bc(1))
      case(' ','P')
        h_0 = 1 / (dx(0) * dx(0) * mu_0)
        do j = 0, P
        do i = 0, P
          Ce(i,j,0) = Ce(i,j,0) - h_0 * Fs(0,i) * Fs(0,j)
        end do
        end do
      end select

      select case(bc(2))
      case(' ','P')
        h_P = 1 / (dx(0) * dx(0) * mu_P)
        do j = 0, P
        do i = 0, P
          Ce(i,j,0) = Ce(i,j,0) - h_P * Fs(P,i) * Fs(P,j)
        end do
        end do
      end select
    end if

    ! contribution from following element (Le⁺) ................................

    if (scan(bc(2), 'DN') > 0) then
      Ce(:,:, 1) = 0
    else
      do j = 0, P
      do i = 0, P
        Ce(i,j,1) =   g(0) * Fs   (P,i) * delta_0(j)  &
                    - g(1) * delta_P(i) * Fs   (0,j)  &
                    - mu_P * delta_P(i) * delta_0(j)
      end do
      end do

      if (primal .and. this%hybrid) then
        h_P = 1 / (dx(0) * dx(1) * mu_P)
        do j = 0, P
        do i = 0, P
          Ce(i,j,1) = Ce(i,j,1) + h_P * Fs(P,i) * Fs(0,j)
        end do
        end do
      end if
    end if

    ! special case: single periodic element ....................................

    if (all(bc == 'P')) then
      Ce(:,:, 0) = Ce(:,:,0) + Ce(:,:,-1) + Ce(:,:,1)
      Ce(:,:,-1) = 0
      Ce(:,:, 1) = 0
    end if

  end function DiffusionMatrix

end subroutine Get_DiffusionMatrix

!-------------------------------------------------------------------------------
!> Provides the generalized eigensystem for interior diffusion operator
!>
!> Returns the column matrix of generalized eigenvectors `S` and the diagonal
!> matrix of eigenvalues `Λ = Lambda` to the interior element stiffness matrix
!> `Cᵢᵢ` and diagonal mass matrix `Mᵢᵢ` of the hybridized element system such
!> that
!>
!>     Sᵀ Cᵢᵢ S = Λ
!>     Sᵀ Mᵢᵢ S = I
!>
!> The diffusion matrix comprises a regular part with diffusivity `ν` and
!> an SVV part with diffusivity `νˢ`
!>
!>     Cᵢᵢ = ν Lᵢᵢ + νˢ Lˢᵢᵢ
!>
!> `Lᵢᵢ` and `Lˢᵢᵢ` are the corresponding interior stiffness matrices

subroutine Get_EllipticEigensystem__w_svv(this, dx, bc, nu, nu_svv, S, Lambda)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)   !< element extensions
  character, intent(in)  :: bc(2)      !< boundary conditions {'','D','N','P'}
  real(RNP), intent(in)  :: nu         !< diffusivity
  real(RNP), intent(in)  :: nu_svv     !< SVV diffusivity [0]
  real(RNP), intent(out) :: S(0:,0:)   !< eigenvectors
  real(RNP), intent(out) :: Lambda(0:) !< eigenvalues

  real(RNP), allocatable :: Mii(:), Cii(:,:,:)

  if (.not. this%hybrid) then
    call Error( 'Get_EllipticEigensystem',             &
                'available only for hybridizable IP', &
                'IP_Element_Operators_1D'             )
  end if

  ! hybrid element operators
  allocate(Mii(0:this%po), Cii(0:this%po, 0:this%po, -1:1))
  Mii = dx(0)/2 * this % w
  call this % Get_DiffusionMatrix(dx, bc, nu, nu_svv, Cii, form='flux')

  ! solve eigenproblem
  call SolveGeneralizedEigenproblem(Cii(:,:,0), Mii, Lambda, S)

  ! singular case
  if (all(bc == 'N') .or. all(bc == 'P')) then
    Lambda(0) = 0
  end if

end subroutine Get_EllipticEigensystem__w_svv

!-------------------------------------------------------------------------------
!> Provides the generalized eigensystem for the interior diffusion operator
!> with no SVV

subroutine Get_EllipticEigensystem__n_svv(this, dx, bc, nu, S, Lambda)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)   !< element extensions
  character, intent(in)  :: bc(2)      !< boundary conditions {'','D','N','P'}
  real(RNP), intent(in)  :: nu         !< diffusivity
  real(RNP), intent(out) :: S(0:,0:)   !< eigenvectors
  real(RNP), intent(out) :: Lambda(0:) !< eigenvalues

  call Get_EllipticEigensystem__w_svv(this, dx, bc, nu, ZERO, S, Lambda)

end subroutine Get_EllipticEigensystem__n_svv

!-------------------------------------------------------------------------------
!> Computes operators for hybrid IP/DG-SEM diffusion problem including SVV
!>
!> To cope with the singular case (Neumann or periodic with c = 0),
!> we use the Moore-Penrose inverse, i.e. `Aii_inv = Ã⁺`

subroutine Get_EllipticSuboperators__w_svv( this, dx, bc, c, nu, nu_svv &
                                          , Aib, Aii_inv                )
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)       !< element extensions
  character, intent(in)  :: bc(2)          !< boundary conds {'','D','N','P'}
  real(RNP), intent(in)  :: c              !< coefficient of linear term
  real(RNP), intent(in)  :: nu             !< diffusivity
  real(RNP), intent(in)  :: nu_svv         !< SVV diffusivity [0]
  real(RNP), intent(out) :: Aib(0:,:)      !< interior-boundary part, Â(0:P,1:2)
  real(RNP), intent(out) :: Aii_inv(0:,0:) !< inv interior part, Ã⁺(0:P,0:P)

  real(RNP), allocatable :: S(:,:), Lambda(:), D_inv(:)
  real(RNP), allocatable :: delta_0(:), delta_P(:)
  real(RNP) :: Fs(0:this%po, 0:this%po), mu_0, mu_P, nu_total
  integer   :: P, i, j

  ! initialization .............................................................

  P = this % po

  mu_0 = this % PenaltyFactor(dx(-1:0))
  mu_P = this % PenaltyFactor(dx( 0:1))

  allocate(delta_0(0:P), source = ZERO)
  delta_0(0) = ONE

  allocate(delta_P(0:P), source = ZERO)
  delta_P(P) = ONE

  allocate(S(0:P,0:P), Lambda(0:P), D_inv(0:P))
  call this % Get_EllipticEigensystem(dx, bc, nu, nu_svv, S, Lambda)

  if (this % Has_SVV()) then
    call this % Get_SVV_StandardDiffMatrix(Fs)
    Fs = nu_svv * Fs
    nu_total = nu + nu_svv
  else
    Fs = 0
    nu_total = nu
  end if
  Fs = Fs + nu * this%D

  ! interior-boundary part .....................................................

  if (all(bc == 'P')) then
    Aib(:,1) =  0
    Aib(:,2) =  0
  else

  select case (bc(1))
  case('D','N')
    Aib(:,1) =  0
  case default ! interior or periodic
    Aib(:,1) = -2/dx(0) * Fs(0,:) - 2 * nu_total * mu_0 * delta_0
  end select

  select case (bc(2))
  case('D','N')
    Aib(:,2) =  0
  case default ! interior or periodic
    Aib(:,2) =  2/dx(0) * Fs(P,:) - 2 * nu_total * mu_P * delta_P
  end select

  end if

  ! inverse interior part ......................................................

  where(abs(c + Lambda) > epsilon(ONE))
    D_inv = ONE / (c + Lambda)
  elsewhere
    D_inv = ZERO
  end where

  do j = 0, P
  do i = 0, P
    Aii_inv(i,j) = sum(S(i,:) * D_inv * S(j,:))
  end do
  end do

end subroutine Get_EllipticSuboperators__w_svv

!-------------------------------------------------------------------------------
!> Computes operators for hybrid IP/DG-SEM diffusion problem without SVV

subroutine Get_EllipticSuboperators__n_svv(this, dx, bc, c, nu, Aib, Aii_inv)
  class(IP_ElementOperators1D), intent(in) :: this
  real(RNP), intent(in)  :: dx(-1:1)       !< element extensions
  character, intent(in)  :: bc(2)          !< boundary conds {'','D','N','P'}
  real(RNP), intent(in)  :: c              !< coefficient of linear term
  real(RNP), intent(in)  :: nu             !< diffusivity
  real(RNP), intent(out) :: Aib(0:,:)      !< interior-boundary part, Â(0:P,1:2)
  real(RNP), intent(out) :: Aii_inv(0:,0:) !< inv interior part, Ã⁺(0:P,0:P)

  call Get_EllipticSuboperators__w_svv(this, dx, bc, c, nu, ZERO, Aib, Aii_inv)

end subroutine Get_EllipticSuboperators__n_svv

!===============================================================================
! IP_ElementOptions1D constructors and type-bound procedures

!-------------------------------------------------------------------------------
!> IP_ElementOptions1D from given operators, optionally overriding the order

function New_IP_ElementOptions1D_o(eop, po) result(this)
  class(StandardOperators1D), intent(in) :: eop !< element operators
  integer,          optional, intent(in) :: po  !< polynomial order

  type(IP_ElementOptions1D) :: this

  if (present(po)) then
    this % po = po
  else
    this % po = eop % po
  end if

  select type(eop)
  class is(IP_ElementOperators1D)
    this % penalty = eop % penalty
    this % hybrid  = eop % hybrid
  end select

  this % no_vdm  = .not. ( eop % Has_Legendre_VDM() )

end function New_IP_ElementOptions1D_o

!-------------------------------------------------------------------------------
!> Extension of MPI_Bcast to objects of type SchwarzOptions3D

subroutine IP_ElementOptions1D_Bcast(this, root, comm)
  class(IP_ElementOptions1D), intent(inout) :: this
  integer,                    intent(in)    :: root !< rank of broadcast root
  type(MPI_Comm),             intent(in)    :: comm !< MPI communicator

  call XMPI_Bcast( this % po      , root, comm )
  call XMPI_Bcast( this % penalty , root, comm )
  call XMPI_Bcast( this % hybrid  , root, comm )
  call XMPI_Bcast( this % no_vdm  , root, comm )

end subroutine IP_ElementOptions1D_Bcast

!===============================================================================

end module IP_Element_Operators_1D
