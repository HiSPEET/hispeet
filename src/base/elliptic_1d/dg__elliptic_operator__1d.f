!> summary:  1D elliptic operator
!> author:   Joerg Stiller
!> date:     2023/03/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Elliptic_Operator__1D
  use Kind_Parameters, only: RNP, RDP
  use Constants      , only: ZERO, ONE, HALF
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__1D

  implicit none
  private

  public :: DG_EllipticOperator_1D
  public :: DG_ElementOptions_1D
  public :: DG_SchwarzOptions_1D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 1D DG-SEM

  type DG_EllipticOperator_1D

    type(DG_ElementOperators_1D) :: eop
    type(DG_SchwarzOperator_1D)  :: schwarz
    real(RNP) :: r_nu_s !< ratio νˢ/(νᵖ+νˢ)

  contains

    procedure :: Init_DG_EllipticOperator_1D

    generic :: Apply              =>  Apply_RC, Apply_RV
    procedure, private ::             Apply_RC, Apply_RV

    generic :: Residual           =>  Residual_RC, Residual_RV
    procedure, private ::             Residual_RC, Residual_RV

    generic :: CG_Method          =>  CG_Method_RC, CG_Method_RV
    procedure, private ::             CG_Method_RC, CG_Method_RV

    generic :: Schwarz_Method     =>  Schwarz_Method_RC, Schwarz_Method_RV
    procedure, private ::             Schwarz_Method_RC, Schwarz_Method_RV

    generic :: SchwarzPCG_Method  =>  SchwarzPCG_Method_RC, SchwarzPCG_Method_RV
    procedure, private ::             SchwarzPCG_Method_RC, SchwarzPCG_Method_RV

    procedure :: HybridSolver     =>  HybridSolver_RC

    procedure :: PhysicalDiffusivity
    procedure :: SpectralDiffusivity

  end type DG_EllipticOperator_1D

  ! constructors
  interface DG_EllipticOperator_1D
    module procedure New_DG_EllipticOperator_1D
  end interface

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Evaluation with regular mesh and constant ν

    module subroutine Eval_RC(this, bc, mask, dx, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_1D),   intent(in)  :: this
      character,                       intent(in)  :: bc(2)   !< BC {D,N,P}
      logical,                         intent(in)  :: mask(:) !< element mask
      real(RNP),                       intent(in)  :: dx      !< ∆xᵉ
      real(RNP),                       intent(in)  :: lambda  !< λ
      real(RNP),                       intent(in)  :: nu      !< ν = νᵖ+νˢ
      real(RNP), contiguous,           intent(in)  :: u(0:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(0:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(0:,:) !< RHS
      real(RNP),             optional, intent(in)  :: bv(2)   !< boundary vals
    end subroutine Eval_RC

    !---------------------------------------------------------------------------
    !> Evaluation with regular mesh and variable ν

    module subroutine Eval_RV(this, bc, mask, dx, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_1D),   intent(in)  :: this
      character,                       intent(in)  :: bc(2)    !< BC {D,N,P}
      logical,                         intent(in)  :: mask(:)  !< element mask
      real(RNP),                       intent(in)  :: dx       !< ∆xᵉ
      real(RNP),                       intent(in)  :: lambda   !< λ
      real(RNP), contiguous,           intent(in)  :: nu(0:,:) !< ν = νᵖ
      real(RNP), contiguous,           intent(in)  :: u (0:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r (0:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f (0:,:) !< RHS
      real(RNP),             optional, intent(in)  :: bv(2)    !< boundary vals
    end subroutine Eval_RV

    !---------------------------------------------------------------------------
    !> Conjugate gradient method with either constant or variable ν

    module subroutine CG_Method_RX( this, bc, bv, mask, dx   &
                                  , lambda, nu_c, nu_v, f, u &
                                  , i_max, r_red, r_max, ni  )

      class(DG_EllipticOperator_1D),   intent(in)    :: this
      character,                       intent(in)    :: bc(2)     !< BC {D,N,P}
      real(RNP),                       intent(in)    :: bv(2)     !< BV, u or u'
      logical,                         intent(in)    :: mask(:)   !< elem mask
      real(RNP),                       intent(in)    :: dx        !< ∆xᵉ
      real(RNP),                       intent(in)    :: lambda    !< λ
      real(RNP),             optional, intent(in)    :: nu_c      !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:) !< νᵖ
      real(RNP), contiguous,           intent(in)    :: f(:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine CG_Method_RX

    !---------------------------------------------------------------------------
    !> Overlapping Schwarz method with either constant or variable ν

    module subroutine Schwarz_Method_RX( this, bc, bv, mask, dx   &
                                       , lambda, nu_c, nu_v, f, u &
                                       , i_max, r_red, r_max, ni  )

      class(DG_EllipticOperator_1D),   intent(in)    :: this
      character,                       intent(in)    :: bc(2)     !< BC {D,N,P}
      real(RNP),                       intent(in)    :: bv(2)     !< BV, u or u'
      logical,                         intent(in)    :: mask(:)   !< elem mask
      real(RNP),                       intent(in)    :: dx        !< ∆xᵉ
      real(RNP),                       intent(in)    :: lambda    !< λ
      real(RNP),             optional, intent(in)    :: nu_c      !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:) !< νᵖ
      real(RNP), contiguous,           intent(in)    :: f(:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine Schwarz_Method_RX

    !---------------------------------------------------------------------------
    !> Schwarz-preconditioned CG method with either constant or variable ν

    module subroutine SchwarzPCG_Method_RX( this, bc, bv, mask, dx   &
                                          , lambda, nu_c, nu_v, f, u &
                                          , i_max, r_red, r_max, ni  )

      class(DG_EllipticOperator_1D),   intent(in)    :: this
      character,                       intent(in)    :: bc(2)     !< BC {D,N,P}
      real(RNP),                       intent(in)    :: bv(2)     !< BV, u or u'
      logical,                         intent(in)    :: mask(:)   !< elem mask
      real(RNP),                       intent(in)    :: dx        !< ∆xᵉ
      real(RNP),                       intent(in)    :: lambda    !< λ
      real(RNP),             optional, intent(in)    :: nu_c      !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:) !< νᵖ
      real(RNP), contiguous,           intent(in)    :: f(:,:)
      real(RNP), contiguous,           intent(inout) :: u(:,:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine SchwarzPCG_Method_RX

    !---------------------------------------------------------------------------
    !> Direct elliptic solver based on hybridization including SVV

    module subroutine HybridSolver_RC( this, bc, bv, dx, lambda, nu, f, u &
                                     , standby                            )

      class(DG_EllipticOperator_1D), intent(in) :: this
      character,             intent(in)    :: bc(2)   !< BC {D,N,P}
      real(RNP),             intent(in)    :: bv(2)   !< boundary values, u or u'
      real(RNP),             intent(in)    :: dx      !< ∆xᵉ
      real(RNP),             intent(in)    :: lambda  !< λ
      real(RNP),             intent(in)    :: nu      !< νᵖ+νˢ
      real(RNP), contiguous, intent(inout) :: f(0:,:) !< source, being destroyed
      real(RNP), contiguous, intent(out)   :: u(0:,:) !< solution
      logical,     optional, intent(in)    :: standby !< keep suboperators [F]

    end subroutine HybridSolver_RC

  end interface

contains

  !=============================================================================
  ! Constructor and initialization

  !-----------------------------------------------------------------------------
  !> New diffusion operator

  function New_DG_EllipticOperator_1D(dg_opt, schwarz_opt, r_nu_s) result(this)
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    class(DG_SchwarzOptions_1D), intent(in) :: schwarz_opt
    real(RNP),         optional, intent(in) :: r_nu_s !< νˢ/(νᵖ+νˢ) [0]

    type(DG_EllipticOperator_1D) :: this

    call Init_DG_EllipticOperator_1D(this, dg_opt, schwarz_opt, r_nu_s)

  end function New_DG_EllipticOperator_1D

  !-----------------------------------------------------------------------------
  !> Initialization of the diffusion operator

  subroutine Init_DG_EllipticOperator_1D(this, dg_opt, schwarz_opt, r_nu_s)
    class(DG_EllipticOperator_1D), intent(inout) :: this
    class(DG_ElementOptions_1D),   intent(in)    :: dg_opt
    class(DG_SchwarzOptions_1D),   intent(in)    :: schwarz_opt
    real(RNP),           optional, intent(in)    :: r_nu_s !< νˢ/(νᵖ+νˢ) [0]

    this % eop = DG_ElementOperators_1D(dg_opt)

    if (present(r_nu_s)) then
      this % r_nu_s = r_nu_s
    else
      this % r_nu_s = 0
    end if

    this % schwarz = DG_SchwarzOperator_1D(schwarz_opt, this%eop, this%r_nu_s)

  end subroutine Init_DG_EllipticOperator_1D

  !=============================================================================
  ! Application of the homogeneous elliptic operator, r = Au with bv = 0

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator with constant diffusivity

  subroutine Apply_RC(this, bc, dx, lambda, nu, u, r, mask)
    class(DG_EllipticOperator_1D), intent(in)  :: this
    character,                     intent(in)  :: bc(2)   !< BC {'D','N','P'}
    real(RNP),                     intent(in)  :: dx      !< ∆xᵉ
    real(RNP),                     intent(in)  :: lambda  !< λ
    real(RNP),                     intent(in)  :: nu      !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: u(:,:)  !< operand
    real(RNP), contiguous,         intent(out) :: r(:,:)  !< result
    logical,     optional,         intent(in)  :: mask(:) !< element mask

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Eval_RC(this, bc, mask, dx, lambda, nu, u, r)
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call Eval_RC(this, bc, mask_, dx, lambda, nu, u, r)
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Apply_RC

  !-----------------------------------------------------------------------------
  !> Application of the elliptic operator with variable diffusivity

  subroutine Apply_RV(this, bc, dx, lambda, nu, u, r, mask)
    class(DG_EllipticOperator_1D), intent(in)  :: this
    character,                     intent(in)  :: bc(2)   !< BC {'D','N','P'}
    real(RNP),                     intent(in)  :: dx      !< ∆xᵉ
    real(RNP),                     intent(in)  :: lambda  !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: u (:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r (:,:) !< result
    logical,     optional,         intent(in)  :: mask(:) !< element mask

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Eval_RV(this, bc, mask, dx, lambda, nu, u, r)
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call Eval_RV(this, bc, mask_, dx, lambda, nu, u, r)
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Apply_RV

  !=============================================================================
  ! Residual of the elliptic equation, r = Au - f

  !-----------------------------------------------------------------------------
  !> Residual for constant diffusivity
  !>
  !> Boundary conditions and values
  !>   - Dirichlet:  `bc = 'D',  bv = u`
  !>   - Neumann:    `bc = 'N',  bv = q = ν∂u/∂x`
  !>
  !> Note that the flux `q` is aligned with the x-direction and not with the
  !> normal!

  subroutine Residual_RC(this, bc, bv, dx, lambda, nu, f, u, r, mask)
    class(DG_EllipticOperator_1D), intent(in)  :: this
    character,                     intent(in)  :: bc(2)  !< BC {'D','N','P'}
    real(RNP),                     intent(in)  :: bv(2)  !< boundary values
    real(RNP),                     intent(in)  :: dx     !< ∆xᵉ
    real(RNP),                     intent(in)  :: lambda !< λ
    real(RNP),                     intent(in)  :: nu     !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: f(:,:) !< RHS
    real(RNP), contiguous,         intent(in)  :: u(:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r(:,:) !< result
    logical,     optional,         intent(in)  :: mask(:) !< element mask

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Eval_RC(this, bc, mask, dx, lambda, nu, u, r, f, bv)
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call Eval_RC(this, bc, mask_, dx, lambda, nu, u, r, f, bv)
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Residual_RC

  !-----------------------------------------------------------------------------
  !> Residual for variable diffusivity
  !>
  !> Boundary conditions and values
  !>   - Dirichlet:  `bc = 'D',  bv = u`
  !>   - Neumann:    `bc = 'N',  bv = q = ν∂u/∂x`
  !>
  !> Note that the flux `q` is aligned with the x-direction and not with the
  !> normal!

  subroutine Residual_RV(this, bc, bv, dx, lambda, nu, f, u, r, mask)
    class(DG_EllipticOperator_1D), intent(in)  :: this
    character,                     intent(in)  :: bc(2)   !< BC {'D','N','P'}
    real(RNP),                     intent(in)  :: bv(2)   !< boundary values
    real(RNP),                     intent(in)  :: dx      !< ∆xᵉ
    real(RNP),                     intent(in)  :: lambda  !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: f (:,:) !< RHS
    real(RNP), contiguous,         intent(in)  :: u (:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r (:,:) !< result
    logical,     optional,         intent(in)  :: mask(:) !< element mask

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Eval_RV(this, bc, mask, dx, lambda, nu, u, r, f, bv)
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call Eval_RV(this, bc, mask_, dx, lambda, nu, u, r, f, bv)
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Residual_RV

  !=============================================================================
  ! Conjugate gradient method

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with constant ν

  subroutine CG_Method_RC( this, bc, bv, dx, lambda, nu, f, u &
                         , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP),             intent(in)    :: nu      !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(in)    :: f(:,:)  !< right hand side
    real(RNP), contiguous, intent(inout) :: u(:,:)  !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call CG_Method_RX( this, bc, bv, mask, dx, lambda, nu, null() &
                       , f, u, i_max, r_red, r_max, ni              )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      call CG_Method_RX( this, bc, bv, mask_, dx, lambda, nu, null() &
                       , f, u, i_max, r_red, r_max, ni               )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine CG_Method_RC

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with variable ν

  subroutine CG_Method_RV( this, bc, bv, dx, lambda, nu, f, u &
                         , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(in)    :: f (:,:) !< right hand side
    real(RNP), contiguous, intent(inout) :: u (:,:) !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call CG_Method_RX( this, bc, bv, mask, dx, lambda, null(), nu &
                       , f, u, i_max, r_red, r_max, ni              )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      call CG_Method_RX( this, bc, bv, mask_, dx, lambda, null(), nu &
                       , f, u, i_max, r_red, r_max, ni               )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine CG_Method_RV

  !=============================================================================
  ! Schwarz method

  !-----------------------------------------------------------------------------
  !> Overlapping Schwarz method with constant ν

  subroutine Schwarz_Method_RC( this, bc, bv, dx, lambda, nu, f, u &
                              , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP),             intent(in)    :: nu      !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(in)    :: f(:,:)  !< right hand side
    real(RNP), contiguous, intent(inout) :: u(:,:)  !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Schwarz_Method_RX( this, bc, bv, mask, dx, lambda, nu, null() &
                            , f, u, i_max, r_red, r_max, ni              )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      call Schwarz_Method_RX( this, bc, bv, mask_, dx, lambda, nu, null() &
                            , f, u, i_max, r_red, r_max, ni               )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Schwarz_Method_RC

  !-----------------------------------------------------------------------------
  !> Overlapping Schwarz method with variable ν

  subroutine Schwarz_Method_RV( this, bc, bv, dx, lambda, nu, f, u &
                              , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(in)    :: f (:,:) !< right hand side
    real(RNP), contiguous, intent(inout) :: u (:,:) !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call Schwarz_Method_RX( this, bc, bv, mask, dx, lambda, null(), nu &
                            , f, u, i_max, r_red, r_max, ni              )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      call Schwarz_Method_RX( this, bc, bv, mask_, dx, lambda, null(), nu &
                            , f, u, i_max, r_red, r_max, ni               )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine Schwarz_Method_RV

  !=============================================================================
  ! Schwarz-preconditioned conjugate gradient method

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with constant ν

  subroutine SchwarzPCG_Method_RC( this, bc, bv, dx, lambda, nu, f, u &
                                 , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP),             intent(in)    :: nu      !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(in)    :: f(:,:)  !< right hand side
    real(RNP), contiguous, intent(inout) :: u(:,:)  !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call SchwarzPCG_Method_RX( this, bc, bv, mask, dx   &
                               , lambda, nu, null(), f, u &
                               , i_max, r_red, r_max, ni  )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call SchwarzPCG_Method_RX( this, bc, bv, mask_, dx  &
                               , lambda, nu, null(), f, u &
                               , i_max, r_red, r_max, ni  )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine SchwarzPCG_Method_RC

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with variable ν

  subroutine SchwarzPCG_Method_RV( this, bc, bv, dx, lambda, nu, f, u &
                                 , i_max, r_red, r_max, mask, ni      )

    class(DG_EllipticOperator_1D), intent(in) :: this

    character,             intent(in)    :: bc(2)   !< BC {'D','N','P'}
    real(RNP),             intent(in)    :: bv(2)   !< boundary values
    real(RNP),             intent(in)    :: dx      !< ∆xᵉ
    real(RNP),             intent(in)    :: lambda  !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(in)    :: f (:,:) !< right hand side
    real(RNP), contiguous, intent(inout) :: u (:,:) !< approximate solution
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    logical,     optional, intent(in)    :: mask(:) !< element mask
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    logical, allocatable, save :: mask_(:)

    if (present(mask)) then
      call SchwarzPCG_Method_RX( this, bc, bv, mask, dx, lambda, null() &
                               , nu, f, u, i_max, r_red, r_max, ni      )
    else
      !$omp master
      allocate(mask_(size(u,2)), source = .true.)
      !$omp end master
      !$omp barrier
      call SchwarzPCG_Method_RX( this, bc, bv, mask_, dx, lambda, null() &
                               , nu, f, u, i_max, r_red, r_max, ni       )
      !$omp master
      deallocate(mask_)
      !$omp end master
    end if

  end subroutine SchwarzPCG_Method_RV

  !-----------------------------------------------------------------------------
  !> Physical diffusivity coefficient νᵖ

  pure real(RNP) function PhysicalDiffusivity(this, nu) result(nu_p)
    class(DG_EllipticOperator_1D), intent(in) :: this
    real(RNP), intent(in) :: nu !< ν = νᵖ + νˢ

    associate(r => this % r_nu_s)
      nu_p = (1 - r) * nu
    end associate

  end function PhysicalDiffusivity

  !-----------------------------------------------------------------------------
  !> Spectral diffusivity coefficient νˢ

  pure real(RNP) function SpectralDiffusivity(this, nu) result(nu_s)
    class(DG_EllipticOperator_1D), intent(in) :: this
    real(RNP), intent(in) :: nu !< ν = νᵖ + νˢ

    associate(r => this % r_nu_s)
      nu_s = r * nu
    end associate

  end function SpectralDiffusivity

  !=============================================================================

end module DG__Elliptic_Operator__1D
