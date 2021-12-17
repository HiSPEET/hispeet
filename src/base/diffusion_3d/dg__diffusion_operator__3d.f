!> summary:  Abstract 3D elliptic operator
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Diffusion_Operator__3D
  use Kind_Parameters, only: RNP, RDP, RSP
  use Constants      , only: ZERO, ONE, HALF
  use DG__Element_Operators__1D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM

  type DG_DiffusionOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    real(RNP) :: lambda = 0                  !< Helmholtz parameter
    real(RNP) :: nu_pc  = 0                  !< constant physical diffusivity
    real(RNP) :: nu_sc  = 0                  !< constant spectral diffusivity
    real(RNP), allocatable :: nu_pv(:,:,:,:) !< variable physical diffusivity
    real(RNP), allocatable :: nu_mf(:,:,:)   !< maximum  diffusivity on faces
    character, allocatable :: bc(:)          !< boundary conditions {P,D,N}
    type(DG_ElementOperators_1D) :: eop

contains

    generic   :: Init_DG_DiffusionOperator_3D  =>  Init_C0, Init_CC, Init_V
    generic   :: SetDiffusivity  =>  SetDiffusivity_C, SetDiffusivity_V
    procedure :: Apply
    procedure :: AddBC

    procedure, private :: Init_C0, Init_CC, Init_V
    procedure, private :: SetDiffusivity_C, SetDiffusivity_V

  end type DG_DiffusionOperator_3D

  ! constructors
  interface DG_DiffusionOperator_3D
    module procedure New_C0
    module procedure New_CC
    module procedure New_V
  end interface

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Set constant physical and spectral diffusivities

    module subroutine SetDiffusivity_C(this, nu_p, nu_s)
      class(DG_DiffusionOperator_3D), intent(inout) :: this
      real(RNP),           intent(in) :: nu_p !< physical diffusivity
      real(RNP), optional, intent(in) :: nu_s !< spectral diffusivity [0]
    end subroutine SetDiffusivity_C

    !---------------------------------------------------------------------------
    !> Set variable physical diffusivity

    module subroutine SetDiffusivity_V(this, nu_p)
      class(DG_DiffusionOperator_3D), intent(inout) :: this
      real(RNP), contiguous, intent(in) :: nu_p(:,:,:,:) !< physical diffusivity
    end subroutine SetDiffusivity_V

    !---------------------------------------------------------------------------
    !> Application of the diffusion operator

    module subroutine Apply(this, u, v, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      real(RNP), contiguous, intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous, intent(out) :: v(:,:,:,:) !< result
      real(RNP), contiguous, intent(in), optional :: f(:,:,:,:) !< RHS
    end subroutine Apply

    !---------------------------------------------------------------------------
    !> Addition of BC to the RHS for regular mesh and constant ν
    !>
    !> `bv` is a boundary variable which contains the Dirichlet or Neumann
    !> boundary values four each boundary. These values are applied to the
    !> right hand side `f` according boundary type specified in `this % bc`.

    module subroutine AddBC_RC(this, bv, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
    end subroutine AddBC_RC

    !-----------------------------------------------------------------------------
    !> Addition of BC to the RHS for regular mesh and variable ν

    module subroutine AddBC_RV(this, bv, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
    end subroutine AddBC_RV

    !---------------------------------------------------------------------------
    !> Addition of BC to the RHS for deformed mesh and constant ν

    module subroutine AddBC_DC(this, bv, f)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
      real(RNP), contiguous, intent(inout) :: f(:,:,:,:)
    end subroutine AddBC_DC

  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> New diffusion operator with constant physical diffusivity

  function New_C0(sem, dg_opt, lambda, nu_p, bc) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    type(DG_DiffusionOperator_3D) :: this

    call Init_CC(this, sem, dg_opt, lambda, nu_p, ZERO, bc)

  end function New_C0

  !-----------------------------------------------------------------------------
  !> New diffusion operator with constant physical and spectral diffusivities

  function New_CC(sem, dg_opt, lambda, nu_p, nu_s, bc) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    real(RNP), intent(in) :: nu_s   !< spectral diffusivity [0]
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}
    type(DG_DiffusionOperator_3D) :: this

    call Init_CC(this, sem, dg_opt, lambda, nu_p, nu_s, bc)

  end function New_CC

  !-----------------------------------------------------------------------------
  !> New diffusion operator with variable physical diffusivity

  function New_V(sem, dg_opt, lambda, nu_p, bc) result(this)
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP),             intent(in) :: lambda        !< Helmholtz parameter
    real(RNP), contiguous, intent(in) :: nu_p(:,:,:,:) !< variable physical ν
    character,             intent(in) :: bc(:)         !< BC {'D','N','P'}
    type(DG_DiffusionOperator_3D) :: this

    call Init_V(this, sem, dg_opt, lambda, nu_p, bc)

  end function New_V

  !=============================================================================
  ! Initialization procedures

  !-----------------------------------------------------------------------------
  !> Initialization with constant physical diffusivity

  subroutine Init_C0(this, sem, dg_opt, lambda, nu_p, bc)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}

    call Init_CC(this, sem, dg_opt, lambda, nu_p, ZERO, bc)

  end subroutine Init_C0

  !-----------------------------------------------------------------------------
  !> Initialization with constant physical and spectral diffusivities

  subroutine Init_CC(this, sem, dg_opt, lambda, nu_p, nu_s, bc)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p   !< physical diffusivity
    real(RNP), intent(in) :: nu_s   !< spectral diffusivity [0]
    character, intent(in) :: bc(:)  !< BC {'D','N','P'}

    this % sem     => sem
    this % eop     =  DG_ElementOperators_1D(dg_opt)
    this % lambda  =  lambda
    this % bc      =  bc

    call this % SetDiffusivity(nu_p, nu_s)

  end subroutine Init_CC

  !-----------------------------------------------------------------------------
  !> Initialization with variable physical diffusivity

  subroutine Init_V(this, sem, dg_opt, lambda, nu_p, bc)
    class(DG_DiffusionOperator_3D), intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D), intent(in) :: dg_opt
    real(RNP), intent(in) :: lambda        !< Helmholtz parameter
    real(RNP), intent(in) :: nu_p(:,:,:,:) !< variable physical diffusivity
    character, intent(in) :: bc(:)         !< BC {'D','N','P'}

    this % sem     => sem
    this % eop     =  DG_ElementOperators_1D(dg_opt)
    this % lambda  =  lambda
    this % bc      =  bc

    call this % SetDiffusivity(nu_p)

  end subroutine Init_V

  !-----------------------------------------------------------------------------
  !> Addition of boundary conditions to the right hande side
  !>
  !> `bv` is a boundary variable which contains the Dirichlet or Neumann
  !> boundary values four each boundary. These values are applied to the
  !> right hand side `f` according boundary type specified in `this % bc`.

  subroutine AddBC(this, bv, f)
    class(DG_DiffusionOperator_3D), intent(in) :: this
    class(SpectralElementBoundaryVariable_3D), target, intent(in) :: bv(:)
    real(RNP), contiguous, intent(inout) :: f(:,:,:,:)

    if (this % sem % mesh % regular) then
      if (allocated(this % nu_pv)) then
        ! regular variable
        call AddBC_RV(this, bv, f)
      else
        ! regular constant
        call AddBC_RC(this, bv, f)
      end if
    else
      ! deformed constant
      call AddBC_DC(this, bv, f)
    end if

  end subroutine AddBC

  !=============================================================================

end module DG__Diffusion_Operator__3D
