!> summary:  Abstract 3D elliptic operator
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DG__Diffusion_Operator__3D
  use Kind_Parameters, only: RNP
  use Constants      , only: ZERO, ONE, HALF
  use DG__Element_Operators__1D
  use Spectral_Element_Mesh__3D

  !-----------------------------------------------------------------------------
  !> Abstract type accommodating 3D elliptic operators

  type, abstract :: DG_DiffusionOperator_3D

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
!!  procedure :: BcToRHS
!!  procedure :: Residual

    procedure, private :: Init_C0, Init_CC, Init_V
    procedure, private :: SetDiffusivity_C, SetDiffusivity_V

  end type DG_DiffusionOperator_3D

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
      real(RNP), intent(in) :: nu_p(:,:,:,:) !< physical diffusivity
    end subroutine SetDiffusivity_V

    !---------------------------------------------------------------------------
    !> Application of the diffusion operator

    module subroutine Apply(this, u, v)
      class(DG_DiffusionOperator_3D), intent(in) :: this
      real(RNP), intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), intent(out) :: v(:,:,:,:) !< result
    end subroutine Apply

  end interface

contains

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

  !=============================================================================

end module DG__Diffusion_Operator__3D
