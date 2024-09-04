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
  use XMPI           , only: XMPI_Bcast
  use Logging_Levels , only: log_level_inner_iteration
  use Array_Assignments
  use Execution_Control
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D

  implicit none
  private

  public :: DG_EllipticOperator_3D


  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM

  type DG_EllipticOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    type(DG_ElementOperators_1D) :: eop
    type(DG_SchwarzOperator_3D)  :: schwarz
    character, allocatable       :: bc(:)  !< boundary conditions {P,D,N}
    real(RNP)                    :: r_nu_s !< ratio νˢ/(νᵖ+νˢ)

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

    procedure :: PhysicalDiffusivity
    procedure :: SpectralDiffusivity

    ! procedures intended for internal use
    procedure :: EnforceBoundaryConditions
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

    module subroutine Eval_RC(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν = νᵖ + νˢ
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_RC

    !---------------------------------------------------------------------------
    !> Evaluation with regular mesh and variable ν

    module subroutine Eval_RV(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda      !< λ
      real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
      real(RNP), contiguous,           intent(in)  :: u (:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r (:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f (:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_RV

    !---------------------------------------------------------------------------
    !> Evaluation with deformed mesh and constant ν

    module subroutine Eval_DC(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν = νᵖ
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_DC

    !---------------------------------------------------------------------------
    !> Evaluation with deformed mesh and variable isotropic ν

    module subroutine Eval_DV(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda      !< λ
      real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:)  !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:)  !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:)  !< RHS
      class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Eval_DV

    !---------------------------------------------------------------------------
    !> Conjugate gradient method with either constant or variable ν

    module subroutine CG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                 , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),   intent(in)    :: this
      real(RNP),                       intent(in)    :: lambda        !< λ
      real(RNP),             optional, intent(in)    :: nu_c          !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:,:,:) !< νᵖ
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine CG_Method_X

    !---------------------------------------------------------------------------
    !> Overlapping Schwarz method with either constant or variable ν

    module subroutine Schwarz_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                      , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),   intent(in)    :: this
      real(RNP),                       intent(in)    :: lambda        !< λ
      real(RNP),             optional, intent(in)    :: nu_c          !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:,:,:) !< νᵖ
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine Schwarz_Method_X

    !---------------------------------------------------------------------------
    !> Schwarz-preconditioned CG method with either constant or variable ν

    module subroutine SchwarzPCG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                         , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),   intent(in)    :: this
      real(RNP),                       intent(in)    :: lambda        !< λ
      real(RNP),             optional, intent(in)    :: nu_c          !< νᵖ+νˢ
      real(RNP), contiguous, optional, intent(in)    :: nu_v(:,:,:,:) !< νᵖ
      real(RNP), contiguous,           intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,           intent(in)    :: f(:,:,:,:)
      class(BoundaryVariable_3D),      intent(in)    :: bv(:)
      integer,                         intent(in)    :: i_max
      real(RNP),             optional, intent(in)    :: r_red
      real(RNP),             optional, intent(in)    :: r_max
      integer,               optional, intent(out)   :: ni

    end subroutine SchwarzPCG_Method_X

  end interface

contains

  !=============================================================================
  ! Constructor and initialization

  !-----------------------------------------------------------------------------
  !> New diffusion operator

  function New_DG_EllipticOperator_3D(sem, dg_opt, schwarz_opt, bc, r_nu_s) &
        result(this)

    class(SpectralElementMesh_3D), target, intent(in) :: sem
    class(DG_ElementOptions_1D),           intent(in) :: dg_opt
    class(DG_SchwarzOptions_3D),           intent(in) :: schwarz_opt
    character,                             intent(in) :: bc(:)
    real(RNP),                   optional, intent(in) :: r_nu_s   !< [0]

    type(DG_EllipticOperator_3D) :: this

    call Init_DG_EllipticOperator_3D(this, sem, dg_opt, schwarz_opt, bc, r_nu_s)

  end function New_DG_EllipticOperator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of the diffusion operator

  subroutine Init_DG_EllipticOperator_3D( this, sem, dg_opt, schwarz_opt, bc, &
                                          r_nu_s                              )

    class(DG_EllipticOperator_3D),         intent(inout) :: this
    class(SpectralElementMesh_3D), target, intent(in)    :: sem
    class(DG_ElementOptions_1D),           intent(in)    :: dg_opt
    class(DG_SchwarzOptions_3D),           intent(in)    :: schwarz_opt
    character,                             intent(in)    :: bc(:)
    real(RNP),                   optional, intent(in)    :: r_nu_s !< [0]

    this % sem => sem
    this % eop =  DG_ElementOperators_1D(dg_opt)
    this % bc  =  bc

    if (present(r_nu_s)) then
      this % r_nu_s = r_nu_s
    else
      this % r_nu_s = 0
    end if

    this % schwarz  =  DG_SchwarzOperator_3D( schwarz_opt, this%eop, sem%mesh &
                                            , bc, this%r_nu_s                 )

  end subroutine Init_DG_EllipticOperator_3D

  !=============================================================================
  ! Application of the homogeneous elliptic operator, r = Au with bv = 0

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator with constant diffusivity

  subroutine Apply_C(this, lambda, nu, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    real(RNP),                     intent(in)  :: lambda     !< λ
    real(RNP),                     intent(in)  :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r(:,:,:,:) !< result

    if (this % sem % mesh % regular) then
      call Eval_RC(this, lambda, nu, u, r)
    else
      call Eval_DC(this, lambda, nu, u, r)
    end if

  end subroutine Apply_C

  !-----------------------------------------------------------------------------
  !> Application of the elliptic operator with variable diffusivity

  subroutine Apply_V(this, lambda, nu, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    real(RNP),                     intent(in)  :: lambda      !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,         intent(out) :: r (:,:,:,:) !< result

    if (this % sem % mesh % regular) then
      call Eval_RV(this, lambda, nu, u, r)
    else
      call Eval_DV(this, lambda, nu, u, r)
    end if

  end subroutine Apply_V

  !=============================================================================
  ! Residual of the elliptic equation, r = Au - f

  !-----------------------------------------------------------------------------
  !> Residual for constant diffusivity

  subroutine Residual_C(this, lambda, nu, f, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    real(RNP),                     intent(in)  :: lambda     !< λ
    real(RNP),                     intent(in)  :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous,         intent(in)  :: f(:,:,:,:) !< RHS
    real(RNP), contiguous,         intent(in)  :: u(:,:,:,:) !< operand
    class(BoundaryVariable_3D),    intent(in)  :: bv(:)      !< boundary values
    real(RNP), contiguous,         intent(out) :: r(:,:,:,:) !< result

    if (this % sem % mesh % regular) then
      call Eval_RC(this, lambda, nu, u, r, f, bv)
    else
      call Eval_DC(this, lambda, nu, u, r, f, bv)
    end if

  end subroutine Residual_C

  !-----------------------------------------------------------------------------
  !> Residual for variable diffusivity

  subroutine Residual_V(this, lambda, nu, f, bv, u, r)
    class(DG_EllipticOperator_3D), intent(in)  :: this
    real(RNP),                     intent(in)  :: lambda      !< λ
    real(RNP), contiguous,         intent(in)  :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous,         intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,         intent(in)  :: f (:,:,:,:) !< RHS
    class(BoundaryVariable_3D),    intent(in)  :: bv(:)       !< boundary values
    real(RNP), contiguous,         intent(out) :: r (:,:,:,:) !< result

    if (this % sem % mesh % regular) then
      call Eval_RV(this, lambda, nu, u, r, f, bv)
    else
      call Eval_DV(this, lambda, nu, u, r, f, bv)
    end if

  end subroutine Residual_V

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with constant ν

  subroutine CG_Method_C(this, lambda, nu, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this

    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(BoundaryVariable_3D), intent(in) :: bv(:)

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call CG_Method_X( this, lambda, nu_c = nu, u = u, f = f, bv = bv,      &
                      i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine CG_Method_C

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with variable ν

  subroutine CG_Method_V(this, lambda, nu, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this

    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(BoundaryVariable_3D), intent(in) :: bv(:)

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call CG_Method_X( this, lambda, nu_v = nu, u = u, f = f, bv = bv,      &
                      i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine CG_Method_V

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method with constant ν

  subroutine Schwarz_Method_C( this, lambda, nu, u, f, bv &
                             , i_max, r_red, r_max, ni    )

    class(DG_EllipticOperator_3D), intent(in) :: this

    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(BoundaryVariable_3D), intent(in) :: bv(:)

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call Schwarz_Method_X( this, lambda, nu_c = nu, u = u, f = f, bv = bv,      &
                           i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine Schwarz_Method_C

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method with variable ν

  subroutine Schwarz_Method_V( this, lambda, nu, u, f, bv &
                             , i_max, r_red, r_max, ni    )

    class(DG_EllipticOperator_3D), intent(in) :: this

    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(BoundaryVariable_3D), intent(in) :: bv(:)

    integer,             intent(in)    :: i_max  !< max num iterations
    real(RNP), optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP), optional, intent(in)    :: r_max  !< max admissible residual
    integer,   optional, intent(out)   :: ni     !< executed num iterations

    call Schwarz_Method_X( this, lambda, nu_v = nu, u = u, f = f, bv = bv,      &
                           i_max = i_max, r_red = r_red, r_max = r_max, ni = ni )

  end subroutine Schwarz_Method_V

  !=============================================================================
  ! Schwarz-preconditioned conjugate gradient method

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with constant ν

  subroutine SchwarzPCG_Method_C( this, lambda, nu, u, f, bv &
                                , i_max, r_red, r_max, ni    )

    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP),             intent(in)    :: nu          !< ν = νᵖ+νˢ
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:)  !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:)  !< right hand side
    class(BoundaryVariable_3D), intent(in) :: bv(:) !< BC
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    call SchwarzPCG_Method_X( this, lambda, nu_c = nu, u = u, f = f, bv = bv &
                            , i_max = i_max, r_red = r_red, r_max = r_max    &
                            , ni = ni                                        )

  end subroutine SchwarzPCG_Method_C

  !-----------------------------------------------------------------------------
  !> Schwarz-preconditioned CG method with variable ν

  subroutine SchwarzPCG_Method_V( this, lambda, nu, u, f, bv &
                                , i_max, r_red, r_max, ni    )

    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP),             intent(in)    :: lambda      !< λ
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν = νᵖ
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side
    class(BoundaryVariable_3D), intent(in) :: bv(:) !< BC
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    call SchwarzPCG_Method_X( this, lambda, nu_v = nu, u = u, f = f, bv = bv &
                            , i_max = i_max, r_red = r_red, r_max = r_max    &
                            , ni = ni                                        )

  end subroutine SchwarzPCG_Method_V

  !-----------------------------------------------------------------------------
  !> Physical diffusivity coefficient νᵖ

  pure real(RNP) function PhysicalDiffusivity(this, nu) result(nu_p)
    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP), intent(in) :: nu !< ν = νᵖ + νˢ

    associate(r => this % r_nu_s)
      nu_p = (1 - r) * nu
    end associate

  end function PhysicalDiffusivity

  !-----------------------------------------------------------------------------
  !> Spectral diffusivity coefficient νˢ

  pure real(RNP) function SpectralDiffusivity(this, nu) result(nu_s)
    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP), intent(in) :: nu !< ν = νᵖ + νˢ

    associate(r => this % r_nu_s)
      nu_s = r * nu
    end associate

  end function SpectralDiffusivity

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

  subroutine EnforceBoundaryConditions(this, bv, jmp_u, avg_q)
    class(DG_EllipticOperator_3D), intent(in) :: this
    class(BoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(inout) :: jmp_u(:,:,:,:) !< trace of u
    real(RNP), contiguous, intent(inout) :: avg_q(:,:,:,:) !< trace of ν du/dn

    logical :: has_bv
    integer :: b, e, f, l

    if (this % sem % mesh % n_bound < 1) return

    has_bv = present(bv)

    associate(boundary => this % sem % mesh % boundary)

      do b = 1, size(boundary)

        select case(this % bc(b))

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

  !=============================================================================

end module DG__Elliptic_Operator__3D
