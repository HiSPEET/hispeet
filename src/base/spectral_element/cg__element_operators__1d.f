!> summary:  Element operators for continuous Galerkin-SEM
!> author:   Joerg Stiller, Gustav Tschirschnitz
!> date:     2018/09/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CG__Element_Operators__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use Eigenproblems,   only: SolveGeneralizedEigenproblem
  use Standard_Operators__1D
  use XMPI
  implicit none
  private

  public :: CG_ElementOperators_1D
  public :: CG_ElementOptions_1D

  !-----------------------------------------------------------------------------
  !> Element operators for continuous Galerkin-SEM

  type, extends(StandardOperators_1D) :: CG_ElementOperators_1D
  contains

    procedure :: Init_CG_ElementOperators_1D
    procedure :: Get_StiffnessMatrix

    generic   :: Get_EllipticEigensystem => Get_EllipticEigensystem__w_svv, &
                                            Get_EllipticEigensystem__n_svv

    procedure, private :: Get_EllipticEigensystem__w_svv
    procedure, private :: Get_EllipticEigensystem__n_svv

    generic   :: Get_EllipticSuboperators => Get_EllipticSuboperators__w_svv, &
                                             Get_EllipticSuboperators__n_svv

    procedure, private :: Get_EllipticSuboperators__w_svv
    procedure, private :: Get_EllipticSuboperators__n_svv

  end type CG_ElementOperators_1D

  ! Constructor interface
  interface CG_ElementOperators_1D
    module procedure New_CG_ElementOperators_1D__f
    module procedure New_CG_ElementOperators_1D__b
  end interface

  !-----------------------------------------------------------------------------
  !> Options for CG_ElementOperators_1D
  !>
  !> @note: identical to StandardOperatorOptions_1D, so far

  type, extends(StandardOperatorOptions_1D) :: CG_ElementOptions_1D
  contains
    procedure :: Bcast => Bcast_CG_ElementOptions1D
  end type CG_ElementOptions_1D

contains

  !=============================================================================
  ! Constructor

  !-----------------------------------------------------------------------------
  !> Constructor for CG_ElementOperators_1D -- flat interface

  function New_CG_ElementOperators_1D__f(po, no_vdm, svv, po_cut_svv) &
      result(this)

    integer,             intent(in) :: po         !< polynomial order
    logical,   optional, intent(in) :: no_vdm     !< skip Vandermonde matrix [F]
    logical,   optional, intent(in) :: svv        !< activate SVV model      [F]
    integer,   optional, intent(in) :: po_cut_svv !< cut-off PO for SVV     [-∞]

    type(CG_ElementOperators_1D) :: this
    type(CG_ElementOptions_1D)   :: opt

    opt % po = po
    if (present(no_vdm    )) opt % no_vdm     = no_vdm
    if (present(svv       )) opt % svv        = svv
    if (present(po_cut_svv)) opt % po_cut_svv = po_cut_svv

    call Init_CG_ElementOperators_1D(this, opt)

  end function New_CG_ElementOperators_1D__f

  !-----------------------------------------------------------------------------
  !> Constructor for CG_ElementOperators_1D -- bundled arguments

  function New_CG_ElementOperators_1D__b(opt) result(this)
    type(CG_ElementOptions_1D), intent(in) :: opt
    type(CG_ElementOperators_1D) :: this

    call Init_CG_ElementOperators_1D(this, opt)

  end function New_CG_ElementOperators_1D__b

  !=============================================================================
  ! Type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization

  subroutine Init_CG_ElementOperators_1D(this, opt)
    class(CG_ElementOperators_1D), intent(inout) :: this
    class(CG_ElementOptions_1D),   intent(in)    :: opt

    call this % Init_StandardOperators_1D(opt)

  end subroutine Init_CG_ElementOperators_1D

  !-----------------------------------------------------------------------------
  !> Returns the regular and, optionally, SVV stiffness matrices for one element
  !>
  !> The element stiffness matrix `Le` represents the nontrivial row entries
  !> of the global stiffness matrix corresponding to a single element.
  !> It must be dimensioned as `Le(0:P,0:P,-1:1)`, where `P = this%po` is the
  !> polynomial order. The third index refers to the
  !>
  !>   * preceding (-1),
  !>   * current (0) and
  !>   * succeeding (1) elements,
  !>
  !> with lengths `dx(-1:1) respectively.
  !>
  !> The corresponding SVV stiffness matrix `Le_svv` is returned if requested
  !> If SVV is not initialized, `Le_svv` is set to zero.

  subroutine Get_StiffnessMatrix(this, dx, bc, Le, Le_svv)
    class(CG_ElementOperators_1D), intent(in) :: this
    real(RNP), intent(in)  :: dx(-1:1)      !< element extensions
    character, intent(in)  :: bc(2)         !< boundary conditions {'','D','N'}
    real(RNP), intent(out) :: Le(0:,0:,-1:) !< regular stiffness matrix
    real(RNP), optional, intent(out) :: Le_svv(0:,0:,-1:) !< SVV stiffness matrix

    integer   :: P
    real(RNP) :: g(-1:1), Ls(0:this%po, 0:this%po)
    real(RNP), allocatable :: delta_0(:), delta_P(:)

    P = this % po
    g = 2 / dx

    allocate(delta_0(0:P), source = ZERO)
    delta_0(0) = ONE

    allocate(delta_P(0:P), source = ZERO)
    delta_P(P) = ONE

    Le = StiffnessMatrix(this%L)

    if (present(Le_svv)) then
      if (this % Has_SVV()) then
        call this % Get_SVV_StandardStiffnessMatrix(Ls)
        Le_svv = StiffnessMatrix(Ls)
      else
        Le_svv = 0
      end if
    end if

  contains

    function StiffnessMatrix(Ls) result(Le)
      real(RNP), intent(in) :: Ls(0:P,0:P) !< 1D standard stiffness matrix
      real(RNP) :: Le(0:P,0:P,-1:1)
      integer   :: i, j

      ! contribution from preceding element (Le⁻) ..............................

      if (scan(bc(1), 'DN') > 0) then
        Le(:,:,-1) = 0
      else
        do j = 0, P
        do i = 0, P
          Le(i,j,-1) = g(-1) * delta_0(i) * Ls(P,j)
        end do
        end do
      end if

      ! own contribution (Le⁰) .................................................

      Le(:,:,0) = g(0) * Ls

      ! nullify Dirichlet entries
      if (bc(1) == 'D') then
        Le(0,:,0) = 0
        Le(:,0,0) = 0
      end if
      if (bc(2) == 'D') then
        Le(P,:,0) = 0
        Le(:,P,0) = 0
      end if

      ! contribution from following element (Le⁺) ..............................

      if (scan(bc(2), 'DN') > 0) then
        Le(:,:,1) = 0
      else
        do j = 0, P
        do i = 0, P
          Le(i,j,1) = g(1) * delta_P(i) * Ls(0,j)
        end do
        end do
      end if

    end function StiffnessMatrix

  end subroutine Get_StiffnessMatrix

  !-----------------------------------------------------------------------------
  !> Provides the generalized eigensystem for the interior diffusion operator
  !>
  !> Returns the column matrix of generalized eigenvectors `S` and the diagonal
  !> matrix of eigenvalues `Λ = Lambda` to the interior element diffusion matrix
  !> `Cᵢᵢ` and diagonal mass matrix `Mᵢᵢ` such that
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

  subroutine Get_EllipticEigensystem__w_svv(this, dx, nu, nu_svv, S, Lambda)
    class(CG_ElementOperators_1D), intent(in) :: this
    real(RNP), intent(in)  :: dx        !< element length
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(in)  :: nu_svv    !< SVV diffusivity [0]
    real(RNP), intent(out) :: S(:,:)    !< eigenvectors
    real(RNP), intent(out) :: Lambda(:) !< eigenvalues

    real(RNP), allocatable :: Mii(:), Cii(:,:), L_svv(:,:)
    integer :: np

    np = size(Lambda)
    if (np < 1) return

    allocate(Mii, source = dx/2 * this % w(1:np))
    allocate(Cii, source = nu * 2/dx * this % L(1:np,1:np))

    if (this % Has_SVV()) then
      allocate(L_svv(0:this%po, 0:this%po))
      call this % Get_SVV_StandardStiffnessMatrix(L_svv)
      Cii = Cii + nu_svv * 2/dx * L_svv(1:np,1:np)
    end if

    call SolveGeneralizedEigenproblem(Cii, Mii, Lambda, S)

  end subroutine Get_EllipticEigensystem__w_svv

  !-----------------------------------------------------------------------------
  !> Provides the generalized eigensystem for the interior diffusion operator
  !> with no SVV

  subroutine Get_EllipticEigensystem__n_svv(this, dx, nu, S, Lambda)
    class(CG_ElementOperators_1D), intent(in) :: this
    real(RNP), intent(in)  :: dx        !< element length
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(out) :: S(:,:)    !< eigenvectors
    real(RNP), intent(out) :: Lambda(:) !< eigenvalues

    call Get_EllipticEigensystem__w_svv(this, dx, nu, ZERO, S, Lambda)

  end subroutine Get_EllipticEigensystem__n_svv

  !-----------------------------------------------------------------------------
  !> Computes operators for condensed CG-SEM diffusion problem including SVV

  subroutine Get_EllipticSuboperators__w_svv( this, dx, c, nu, nu_svv  &
                                            , Aib, Abb, Aii_inv        )

    class(CG_ElementOperators_1D), intent(in) :: this
    real(RNP), intent(in)  :: dx           !< element length
    real(RNP), intent(in)  :: c            !< coefficient of linear term, c ≥ 0
    real(RNP), intent(in)  :: nu           !< regular diffusivity
    real(RNP), intent(in)  :: nu_svv       !< SVV diffusivity
    real(RNP), intent(out) :: Aib(:,:)     !< interior-boundary part, dim (po-1,2)
    real(RNP), intent(out) :: Abb(:,:)     !< boundary-boundary part, dim (2,2)
    real(RNP), intent(out) :: Aii_inv(:,:) !< Aᵢᵢ⁻¹, dimension (po-1,po-1)

    real(RNP), allocatable :: Cs(:,:), S(:,:), Lambda(:), D_inv(:)
    real(RNP) :: g0, g1
    integer   :: po, i, j, np

    po = this%po

    ! standard diffusion matrix comprising regular and SVV contributions
    allocate(Cs(0:po,0:po), source = ZERO)
    if (this % Has_SVV()) then
      call this % Get_SVV_StandardStiffnessMatrix(Cs)
    end if
    Cs = nu * this%L + nu_svv * Cs

    associate(Ms => this%w)

      np = po - 1

      allocate(S(np,np), Lambda(np), D_inv(np))
      call this % Get_EllipticEigensystem(dx, nu, nu_svv, S, Lambda)

      g0 = c * dx / 2
      g1 = 2 / dx

      do i = 1, np
        Aib(i,1)  =  g1 * Cs( 0,i)
        Aib(i,2)  =  g1 * Cs(po,i)
      end do

      Abb(1,1)  =  g0 * Ms( 0)  +  g1 * Cs( 0, 0)
      Abb(2,1)  =                  g1 * Cs(po, 0)
      Abb(1,2)  =                  g1 * Cs( 0,po)
      Abb(2,2)  =  g0 * Ms(po)  +  g1 * Cs(po,po)

      where (abs(c + Lambda) > epsilon(ONE))
        D_inv = ONE / (c + Lambda)
      elsewhere
        D_inv = ZERO
      end where

      do j = 1, np
      do i = 1, np
        Aii_inv(i,j) = sum(S(i,:) * D_inv * S(j,:))
      end do
      end do

    end associate

  end subroutine Get_EllipticSuboperators__w_svv

  !-----------------------------------------------------------------------------
  !> Computes operators for condensed CG-SEM diffusion problem without SVV

  subroutine Get_EllipticSuboperators__n_svv(this, dx, c, nu, Aib, Abb, Aii_inv)
    class(CG_ElementOperators_1D), intent(in) :: this
    real(RNP), intent(in)  :: dx           !< element length
    real(RNP), intent(in)  :: c            !< coefficient of linear term, c ≥ 0
    real(RNP), intent(in)  :: nu           !< regular diffusivity
    real(RNP), intent(out) :: Aib(:,:)     !< interior-boundary part, dim (po-1,2)
    real(RNP), intent(out) :: Abb(:,:)     !< boundary-boundary part, dim (2,2)
    real(RNP), intent(out) :: Aii_inv(:,:) !< Aᵢᵢ⁻¹, dimension (po-1,po-1)

    call Get_EllipticSuboperators__w_svv(this, dx, c, nu, ZERO, Aib, Abb, Aii_inv)

  end subroutine Get_EllipticSuboperators__n_svv

  !-----------------------------------------------------------------------------
  !> Extension of MPI_Bcast to objects of type CG_ElementOptions_1D

  subroutine Bcast_CG_ElementOptions1D(this, root, comm)
    class(CG_ElementOptions_1D), intent(inout) :: this
    integer,                     intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),              intent(in)    :: comm !< MPI communicator

    call this % StandardOperatorOptions_1D % Bcast(root, comm)

  end subroutine Bcast_CG_ElementOptions1D

  !=============================================================================

end module CG__Element_Operators__1D
