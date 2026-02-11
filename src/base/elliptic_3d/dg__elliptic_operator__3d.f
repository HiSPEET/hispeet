!> summary:  3D elliptic operator
!> author:   Joerg Stiller
!> date:     2021/08/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - The iterators (CG, Schwarz, SchwarzPCG) return no correct value for `ni`
!>     when called from empty partition
!===============================================================================

module DG__Elliptic_Operator__3D
  use Kind_Parameters, only: RNP, RDP, RSP
  use Constants      , only: ZERO, ONE, HALF
  use XMPI
  use Logging_Levels
  use Array_Assignments
  use Execution_Control
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D
  use Mesh_Element__3D
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D

  implicit none
  private

  public :: DG_EllipticOperator_3D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM
  !>
  !> When applied with local refinement, the component `interior_bc` defines the
  !> conditions at interior boundaries: a blank space yields a direct coupling
  !> to adjacent frozen elements, whereas `D` results in Dirichlet boundary
  !> conditions.

  type DG_EllipticOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    type(DG_ElementOperators_1D) :: eop
    type(DG_SchwarzOperator_3D)  :: schwarz
    character :: interior_bc = ' ' !< coupling with frozen elements {' ','D'}

  contains

    procedure :: Init_DG_EllipticOperator_3D

    generic :: Apply              =>  Apply_C, Apply_V
    procedure, private ::             Apply_C, Apply_V

    generic :: Residual           =>  Residual_C, Residual_V
    procedure, private ::             Residual_C, Residual_V

    generic :: CG_Method          =>  CG_Method_C, CG_Method_V
    procedure, private ::             CG_Method_C, CG_Method_V

    generic :: Schwarz_Method     =>  Schwarz_Method_C, Schwarz_Method_V
    procedure, private ::             Schwarz_Method_C, Schwarz_Method_V

    generic :: SchwarzPCG_Method  =>  SchwarzPCG_Method_C, SchwarzPCG_Method_V
    procedure, private ::             SchwarzPCG_Method_C, SchwarzPCG_Method_V

    ! procedures intended for internal use .....................................

    procedure :: EnforceBoundaryConditions

    generic   :: GetElementBoundaryFluxes => GetElementBoundaryFluxes_C, &
                                             GetElementBoundaryFluxes_V

    procedure, private :: GetElementBoundaryFluxes_C, &
                          GetElementBoundaryFluxes_V

    procedure :: CG_Method_X
    procedure :: Schwarz_Method_X
    procedure :: SchwarzPCG_Method_X

  end type DG_EllipticOperator_3D

  ! constructors
  interface DG_EllipticOperator_3D
    module procedure New_DG_EllipticOperator_3D
  end interface

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Evaluation with regular mesh and constant ν

    module subroutine Eval_RC(this, bc, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      character,                       intent(in)  :: bc(:)      !< {P,D,N}
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν = νᵖ + νˢ
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_RC

    !---------------------------------------------------------------------------
    !> Evaluation with regular mesh and variable ν

    module subroutine Eval_RV(this, bc, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      character,                       intent(in)  :: bc(:)       !< {P,D,N}
      real(RNP),                       intent(in)  :: lambda      !< λ
      real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
      real(RNP), contiguous,           intent(in)  :: u (:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r (:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f (:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_RV

    !---------------------------------------------------------------------------
    !> Evaluation with deformed mesh and constant ν

    module subroutine Eval_DC(this, bc, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      character,                       intent(in)  :: bc(:)      !< {P,D,N}
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν = νᵖ
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_DC

    !---------------------------------------------------------------------------
    !> Evaluation with deformed mesh and variable isotropic ν

    module subroutine Eval_DV(this, bc, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      character,                       intent(in)  :: bc(:)       !< {P,D,N}
      real(RNP),                       intent(in)  :: lambda      !< λ
      real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:)  !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:)  !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:)  !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_DV

    !---------------------------------------------------------------------------
    !> Conjugate gradient method with either constant or variable ν

    module subroutine CG_Method_X( this, bc, lambda, nu_c, nu_v, u, f &
                                 , bv, i_max, r_red, r_max, ni        )

      class(DG_EllipticOperator_3D),        intent(in)    :: this
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP),             optional,      intent(in)    :: nu_c
      real(RNP), contiguous, optional,      intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
      integer,                              intent(in)    :: i_max
      real(RNP),                  optional, intent(in)    :: r_red
      real(RNP),                  optional, intent(in)    :: r_max
      integer,                    optional, intent(out)   :: ni

    end subroutine CG_Method_X

    !---------------------------------------------------------------------------
    !> Overlapping Schwarz method with either constant or variable ν

    module subroutine Schwarz_Method_X( this, bc, lambda, nu_c, nu_v, u, f &
                                      , bv, i_max, r_red, r_max, ni        )

      class(DG_EllipticOperator_3D),        intent(in)    :: this
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP),                  optional, intent(in)    :: nu_c
      real(RNP), contiguous,      optional, intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
      integer,                              intent(in)    :: i_max
      real(RNP),                  optional, intent(in)    :: r_red
      real(RNP),                  optional, intent(in)    :: r_max
      integer,                    optional, intent(out)   :: ni

    end subroutine Schwarz_Method_X

    !---------------------------------------------------------------------------
    !> Schwarz-preconditioned CG method with either constant or variable ν

    module subroutine SchwarzPCG_Method_X( this, bc, lambda, nu_c, nu_v, u, f &
                                         , bv , i_max, r_red, r_max, ni       )

      class(DG_EllipticOperator_3D),        intent(in)    :: this
      character,                            intent(in)    :: bc(:)
      real(RNP),                            intent(in)    :: lambda
      real(RNP),                  optional, intent(in)    :: nu_c
      real(RNP), contiguous,      optional, intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
      integer,                              intent(in)    :: i_max
      real(RNP),                  optional, intent(in)    :: r_red
      real(RNP),                  optional, intent(in)    :: r_max
      integer,                    optional, intent(out)   :: ni

    end subroutine SchwarzPCG_Method_X

  end interface

contains

  !=============================================================================
  ! Constructor and initialization

  !-----------------------------------------------------------------------------
  !> New diffusion operator
  !>
  !> Skipping `penalty` or passing a negative value yields the default specified
  !> in DG_ElementOptions_1D.
  !>
  !> Skipping `interior_bc` or passing a space results in a direct coupling to
  !> frozen elements. Specicyfing 'D' enforces Dirichlet conditions at interior
  !> boundaries

  function New_DG_EllipticOperator_3D &
        (sem, schwarz_opt, penalty, interior_bc) result(this)

    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_SchwarzOptions_3D),           intent(in) :: schwarz_opt
    real(RNP),                   optional, intent(in) :: penalty
    character,                   optional, intent(in) :: interior_bc

    type(DG_EllipticOperator_3D) :: this

    call Init_DG_EllipticOperator_3D &
            (this, sem, schwarz_opt, penalty, interior_bc)

  end function New_DG_EllipticOperator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of the diffusion operator
  !>
  !> Skipping `penalty` or passing a negative value yields the default specified
  !> in DG_ElementOptions_1D

  subroutine Init_DG_EllipticOperator_3D &
        (this, sem, schwarz_opt, penalty, interior_bc)

    class(DG_EllipticOperator_3D),         intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in)    :: sem
    class(DG_SchwarzOptions_3D),           intent(in)    :: schwarz_opt
    real(RNP),                   optional, intent(in)    :: penalty
    character,                   optional, intent(in)    :: interior_bc

    this % sem => sem
    this % eop =  DG_ElementOperators_1D(sem%std_op, penalty)

    this % schwarz = DG_SchwarzOperator_3D(schwarz_opt, this%eop, sem%mesh)

    if (present(interior_bc)) then
      this % interior_bc = interior_bc
    end if

    if (scan(' D', this%interior_bc) == 0) then
      call Error( 'Init_DG_EllipticOperator_3D'                       &
                , 'interior_bc"'//this%interior_bc//'" not supported' &
                , 'DG__Elliptic_Operator__3D'                         )
    end if

  end subroutine Init_DG_EllipticOperator_3D

  !=============================================================================
  ! Application of the homogeneous elliptic operator, r = Au with bv = 0

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator with constant diffusivity

  subroutine Apply_C(this, bc, lambda, nu, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    character,                     intent(in)  :: bc(:)      !< {P,D,N}
    real(RNP),                     intent(in)  :: lambda     !< λ
    real(RNP),                     intent(in)  :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r(:,:,:,:) !< result
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    if (this % sem % mesh % regular) then
      call Eval_RC(this, bc, lambda, nu, u, r, bv = bv)
    else
      call Eval_DC(this, bc, lambda, nu, u, r, bv = bv)
    end if

  end subroutine Apply_C

  !-----------------------------------------------------------------------------
  !> Application of the elliptic operator with variable diffusivity

  subroutine Apply_V(this, bc, lambda, nu, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    character,                     intent(in)  :: bc(:)       !< {P,D,N}
    real(RNP),                     intent(in)  :: lambda      !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r (:,:,:,:) !< result
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    if (this % sem % mesh % regular) then
      call Eval_RV(this, bc, lambda, nu, u, r, bv = bv)
    else
      call Eval_DV(this, bc, lambda, nu, u, r, bv = bv)
    end if

  end subroutine Apply_V

  !=============================================================================
  ! Residual of the elliptic equation, r = Au - f

  !-----------------------------------------------------------------------------
  !> Residual for constant diffusivity
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine Residual_C(this, bc, lambda, nu, f, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    character,                     intent(in)  :: bc(:)      !< {P,D,N}
    real(RNP),                     intent(in)  :: lambda      !< λ
    real(RNP),                     intent(in)  :: nu          !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: f(:,:,:,:)  !< RHS
    real(RNP), contiguous,         intent(in)  :: u(:,:,:,:)  !< operand
    real(RNP), contiguous,         intent(out) :: r(:,:,:,:)  !< result
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    if (this % sem % mesh % regular) then
      call Eval_RC(this, bc, lambda, nu, u, r, f, bv)
    else
      call Eval_DC(this, bc, lambda, nu, u, r, f, bv)
    end if

  end subroutine Residual_C

  !-----------------------------------------------------------------------------
  !> Residual for variable diffusivity
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine Residual_V(this, bc, lambda, nu, f, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    character,                     intent(in)  :: bc(:)       !< {P,D,N}
    real(RNP),                     intent(in)  :: lambda      !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,         intent(in)  :: f (:,:,:,:) !< RHS
    real(RNP), contiguous,         intent(out) :: r (:,:,:,:) !< result
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    if (this % sem % mesh % regular) then
      call Eval_RV(this, bc, lambda, nu, u, r, f, bv)
    else
      call Eval_DV(this, bc, lambda, nu, u, r, f, bv)
    end if

  end subroutine Residual_V

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with constant ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine CG_Method_C(this, bc, lambda, nu, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)      !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call CG_Method_X( this, bc, lambda, nu_c = nu, u = u, f = f, bv = bv,  &
                      i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine CG_Method_C

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with variable ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine CG_Method_V(this, bc, lambda, nu, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)       !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call CG_Method_X( this, bc, lambda, nu_v = nu, u = u, f = f, bv = bv,  &
                      i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine CG_Method_V

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method with constant ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine Schwarz_Method_C( this, bc, lambda, nu, u, f, bv &
                             , i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)      !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call Schwarz_Method_X( this, bc, lambda, nu_c = nu, u = u, f = f, bv = bv,  &
                           i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine Schwarz_Method_C

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method with variable ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine Schwarz_Method_V( this, bc, lambda, nu, u, f, bv &
                             , i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)       !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call Schwarz_Method_X( this, bc, lambda, nu_v = nu, u = u, f = f, bv = bv,  &
                           i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine Schwarz_Method_V

  !=============================================================================
  ! Schwarz-preconditioned conjugate gradient method

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with constant ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine SchwarzPCG_Method_C( this, bc, lambda, nu, u, f, bv &
                                , i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)       !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP),             intent(in)    :: nu          !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)  !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:)  !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    call SchwarzPCG_Method_X( this, bc, lambda, nu_c = nu &
                            , u = u, f = f, bv = bv       &
                            , i_max = i_max               &
                            , r_red = r_red               &
                            , r_max = r_max               &
                            , ni = ni                     )

  end subroutine SchwarzPCG_Method_C

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with variable ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent

  subroutine SchwarzPCG_Method_V( this, bc, lambda, nu, u, f, bv &
                                , i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D), intent(in) :: this

    character,             intent(in)    :: bc(:)       !< BC {P,D,N}
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    call SchwarzPCG_Method_X( this, bc, lambda, nu_v = nu &
                            , u = u, f = f, bv = bv       &
                            , i_max = i_max               &
                            , r_red = r_red               &
                            , r_max = r_max               &
                            , ni = ni                     )

  end subroutine SchwarzPCG_Method_V

  !=============================================================================
  ! Shared procedures

  !-----------------------------------------------------------------------------
  !> Weak enforcement of boundary conditions
  !>
  !> Input:
  !>
  !>   - `bv`     Dirichlet or Neumann boundary values in component 1
  !>   - `jmp_u`  u⁻    on boundary element faces
  !>   - `avg_q`  n⋅q⁻  on boundary element faces
  !>
  !> Interior and ghost face entries `jmp_u` and `avg_q` will be ignored and
  !> remain unchanged.
  !>
  !> On output, the boundary conditions are applied as follows:
  !>
  !>   - Dirichlet boundaries (bc = 'D'):
  !>
  !>          jmp_u  ←  n⋅[u]  =  2(u⁻ - u_b)
  !>          avg_q  ←  n⋅{q}  =  n⋅q⁻            (unchanged)
  !>
  !>   - Neumann boundaries (bc = 'N'):
  !>
  !>          jmp_u  ←  n⋅[u]  =  0
  !>          avg_q  ←  n⋅{q}  =  q_b

  subroutine EnforceBoundaryConditions(this, bc, bv, jmp_u, avg_q)
    class(DG_EllipticOperator_3D), intent(in) :: this
    character, intent(in) :: bc(:) !< boundary conditions {P,D,N}
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:) !< boundary values
    real(RNP), contiguous, intent(inout) :: jmp_u(:,:,:,:) !< trace of u
    real(RNP), contiguous, intent(inout) :: avg_q(:,:,:,:) !< trace of ν du/dn

    logical :: has_bv
    integer :: b, e, f, l

    if (this % sem % mesh % n_bound < 1) return

    has_bv = present(bv)

    associate(boundary => this % sem % mesh % boundary)

      do b = 1, size(boundary)

        select case(bc(b))

        case('D')  ! avg_q remains unchanged !

          if (has_bv) then
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = 2 * (jmp_u (:,:,f,e) - bv(b) % val(:,:,l,1))
            end do
            !$omp end do nowait
          else
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = 2 * jmp_u (:,:,f,e)
            end do
            !$omp end do nowait
          end if

        case('N')

          if (has_bv) then
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = ZERO
              avg_q(:,:,f,e) = bv(b) % val(:,:,l,1)
            end do
            !$omp end do nowait
          else
            !$omp do
            do l = 1, boundary(b) % n_face
              e = boundary(b) % face(l) % element_id
              f = boundary(b) % face(l) % element_face
              jmp_u(:,:,f,e) = ZERO
              avg_q(:,:,f,e) = ZERO
            end do
            !$omp end do nowait
          end if

        end select

      end do
     !$omp barrier

    end associate

  end subroutine EnforceBoundaryConditions

  !-----------------------------------------------------------------------------
  !> Compose element-boundary fluxes from flux traces -- constant diffusivity

  subroutine GetElementBoundaryFluxes_C( this, element, struct, hom_bc &
                                       , e, f, tr, jmp_u, avg_q        )

    class(DG_EllipticOperator_3D), intent(in) :: this
    class(MeshElement_3D), intent(in) :: element
    logical,   intent(in)  :: struct   !< F/T for un/structured mesh
    logical,   intent(in)  :: hom_bc   !< F/T for in/homogeneous internal BC
    integer,   intent(in)  :: e        !< element ID
    integer,   intent(in)  :: f        !< element face
    real(RNP), contiguous, intent(in)  :: tr(:,:,:,:,:) !< traces of u, q_n
    real(RNP), contiguous, intent(out) :: jmp_u(:,:)    !< normal jump n⋅[u]
    real(RNP), contiguous, intent(out) :: avg_q(:,:)    !< normal flux n⋅{q}

    integer :: i, l, m

    if (element % face(f) % boundary == 0 .and. this % interior_bc == 'D') then
      ! interface to frozen element treated as Dirichlet boundary
      if (hom_bc) then
        jmp_u = 2 * tr(:,:,f,e,1)
      else
        i = element % face(f) % i_neighbor
        l = element % neighbor(i) % id
        m = element % neighbor(i) % component
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,1), jmp_u)
        jmp_u = 2 * (tr(:,:,f,e,1) - jmp_u)
      end if
      avg_q = tr(:,:,f,e,2)

    else if (element % face(f) % i_neighbor > 0) then
      ! active or frozen neighbor
      i = element % face(f) % i_neighbor
      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      if (struct) then
        jmp_u = (tr(:,:,f,e,1) - tr(:,:,m,l,1))
        avg_q = (tr(:,:,f,e,2) - tr(:,:,m,l,2)) * HALF
      else
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,1), jmp_u)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,2), avg_q)
        jmp_u = (tr(:,:,f,e,1) - jmp_u)
        avg_q = (tr(:,:,f,e,2) - avg_q) * HALF
      end if

    else
      ! domain boundary: treated by EnforceBoundaryConditions
      jmp_u = tr(:,:,f,e,1)
      avg_q = tr(:,:,f,e,2)

   end if

  end subroutine GetElementBoundaryFluxes_C


  !-----------------------------------------------------------------------------
  !> Compose element-boundary fluxes from flux traces -- variable diffusivity

  subroutine GetElementBoundaryFluxes_V( this, element, struct, hom_bc  &
                                       , e, f, tr, nu_max, jmp_u, avg_q )

    class(DG_EllipticOperator_3D), intent(in) :: this
    class(MeshElement_3D), intent(in) :: element
    logical,   intent(in)  :: struct   !< F/T for un/structured mesh
    logical,   intent(in)  :: hom_bc   !< F/T for in/homogeneous internal BC
    integer,   intent(in)  :: e        !< element ID
    integer,   intent(in)  :: f        !< element face
    real(RNP), contiguous, intent(in)  :: tr(:,:,:,:,:) !< traces of u, q_n
    real(RNP), contiguous, intent(out) :: nu_max(:,:)   !< max(ν⁻,ν⁺)
    real(RNP), contiguous, intent(out) :: jmp_u(:,:)    !< normal jump n⋅[u]
    real(RNP), contiguous, intent(out) :: avg_q(:,:)    !< normal flux n⋅{q}

    integer :: i, l, m

    if (element % face(f) % boundary == 0 .and. this % interior_bc == 'D') then
      ! interface to frozen element treated as Dirichlet boundary
      nu_max = tr(:,:,f,e,1)
      if (hom_bc) then
        jmp_u = 2 * tr(:,:,f,e,2)
      else
        i = element % face(f) % i_neighbor
        l = element % neighbor(i) % id
        m = element % neighbor(i) % component
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,2), jmp_u)
        jmp_u = 2 * (tr(:,:,f,e,2) - jmp_u)
      end if
      avg_q = tr(:,:,f,e,3)

    else if (element % face(f) % i_neighbor > 0) then
      ! active or frozen neighbor
      i = element % face(f) % i_neighbor
      l = element % neighbor(i) % id
      m = element % neighbor(i) % component
      if (struct) then
        nu_max = max(tr(:,:,f,e,1), tr(:,:,m,l,1))
        jmp_u  = (tr(:,:,f,e,2) - tr(:,:,m,l,2))
        avg_q  = (tr(:,:,f,e,3) - tr(:,:,m,l,3)) * HALF
      else
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,1), nu_max)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,2), jmp_u)
        call element % AlignFromNeighborFace(f, i, tr(:,:,m,l,3), avg_q)
        nu_max = max(tr(:,:,f,e,1), nu_max)
        jmp_u  = (tr(:,:,f,e,2) - jmp_u)
        avg_q  = (tr(:,:,f,e,3) - avg_q) * HALF
      end if

    else
      ! domain boundary: treated by EnforceBoundaryConditions
      nu_max = tr(:,:,f,e,1)
      jmp_u  = tr(:,:,f,e,2)
      avg_q  = tr(:,:,f,e,3)

   end if

  end subroutine GetElementBoundaryFluxes_V

  !=============================================================================

end module DG__Elliptic_Operator__3D
