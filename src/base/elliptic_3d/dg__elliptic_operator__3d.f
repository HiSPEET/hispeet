!> summary:  Abstract 3D elliptic operator
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
  use Execution_Control
  use DG__Element_Operators__1D
  use DG__Schwarz_Operator__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  !-----------------------------------------------------------------------------
  !> Base type for scalar diffusion operators for 3D DG-SEM

  type DG_EllipticOperator_3D

    class(SpectralElementMesh_3D), pointer :: sem => null()
    type(DG_ElementOperators_1D) :: eop
    type(DG_SchwarzOperator_3D)  :: schwarz
    character, allocatable       :: bc(:)  !< boundary conditions {P,D,N}
    real(RNP)                    :: r_nu_s !< ratio νˢ/(ν+νˢ)

  contains

    procedure :: Init_DG_EllipticOperator_3D

    generic :: Apply => Apply_C, Apply_V
    procedure, private :: Apply_C, Apply_V

    generic :: CG_Method => CG_Method_C, CG_Method_V
    procedure, private :: CG_Method_C, CG_Method_V

    generic :: Schwarz_Method => Schwarz_Method_C, Schwarz_Method_V
    procedure, private :: Schwarz_Method_C, Schwarz_Method_V

    generic :: SchwarzPCG_Method => SchwarzPCG_Method_C, SchwarzPCG_Method_V
    procedure, private :: SchwarzPCG_Method_C, SchwarzPCG_Method_V

    procedure :: SpectralDiffusivity
    procedure :: TotalDiffusivity

  end type DG_EllipticOperator_3D

  ! constructors
  interface DG_EllipticOperator_3D
    module procedure New_DG_EllipticOperator_3D
  end interface

  !=============================================================================
  ! Interfaces to separate module procedures

  interface

    !---------------------------------------------------------------------------
    !> Application with regular mesh and constant ν

    module subroutine Apply_RC(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Apply_RC

    !---------------------------------------------------------------------------
    !> Application with regular mesh and variable ν

    module subroutine Apply_RV(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda      !< λ
      real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
      real(RNP), contiguous,           intent(in)  :: u (:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r (:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f (:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Apply_RV

    !---------------------------------------------------------------------------
    !> Application with deformed mesh and constant ν

    module subroutine Apply_DC(this, lambda, nu, u, r, f, bv)
      class(DG_EllipticOperator_3D),   intent(in)  :: this
      real(RNP),                       intent(in)  :: lambda     !< λ
      real(RNP),                       intent(in)  :: nu         !< ν
      real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
      real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
      real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
      class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    end subroutine Apply_DC

    !---------------------------------------------------------------------------
    !> Conjugate gradient method with either constant or variable ν

    module subroutine CG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                 , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),             intent(in)    :: this
      real(RNP),                                 intent(in)    :: lambda
      real(RNP),                       optional, intent(in)    :: nu_c
      real(RNP), contiguous,           optional, intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                     intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv(:)
      integer,                                   intent(in)    :: i_max
      real(RNP),                       optional, intent(in)    :: r_red
      real(RNP),                       optional, intent(in)    :: r_max
      integer,                         optional, intent(out)   :: ni

    end subroutine CG_Method_X

    !---------------------------------------------------------------------------
    !> Overlapping Schwarz method with either constant or variable ν

    module subroutine Schwarz_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                      , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),             intent(in)    :: this
      real(RNP),                                 intent(in)    :: lambda
      real(RNP),                       optional, intent(in)    :: nu_c
      real(RNP), contiguous,           optional, intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                     intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv(:)
      integer,                                   intent(in)    :: i_max
      real(RNP),                       optional, intent(in)    :: r_red
      real(RNP),                       optional, intent(in)    :: r_max
      integer,                         optional, intent(out)   :: ni

    end subroutine Schwarz_Method_X

    !---------------------------------------------------------------------------
    !> Schwarz-preconditioned CG method with either constant or variable ν

    module subroutine SchwarzPCG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                         , i_max, r_red, r_max, ni            )

      class(DG_EllipticOperator_3D),             intent(in)    :: this
      real(RNP),                                 intent(in)    :: lambda
      real(RNP),                       optional, intent(in)    :: nu_c
      real(RNP), contiguous,           optional, intent(in)    :: nu_v(:,:,:,:)
      real(RNP), contiguous,                     intent(inout) :: u(:,:,:,:)
      real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv(:)
      integer,                                   intent(in)    :: i_max
      real(RNP),                       optional, intent(in)    :: r_red
      real(RNP),                       optional, intent(in)    :: r_max
      integer,                         optional, intent(out)   :: ni

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
    real(RNP),                   optional, intent(in) :: r_nu_s !< νˢ/(ν+νˢ) [0]

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

    this % sem      => sem
    this % eop      =  DG_ElementOperators_1D(dg_opt)
    this % bc       =  bc
    this % schwarz  =  DG_SchwarzOperator_3D(schwarz_opt, this%eop, sem%mesh, &
                                             bc, r_nu_s)

  end subroutine Init_DG_EllipticOperator_3D

  !=============================================================================
  ! Application of the diffusion operator, r = Au - f

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator with constant diffusivity

  subroutine Apply_C(this, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    real(RNP),                       intent(in)  :: lambda     !< λ
    real(RNP),                       intent(in)  :: nu         !< ν
    real(RNP), contiguous,           intent(in)  :: u(:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r(:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f(:,:,:,:) !< RHS
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    if (this % sem % mesh % regular) then
      call Apply_RC(this, lambda, nu, u, r, f, bv)
    else
      call Apply_DC(this, lambda, nu, u, r, f, bv)
    end if

  end subroutine Apply_C

  !-----------------------------------------------------------------------------
  !> Application of the diffusion operator with variable diffusivity

  subroutine Apply_V(this, lambda, nu, u, r, f, bv)
    class(DG_EllipticOperator_3D),   intent(in)  :: this
    real(RNP),                       intent(in)  :: lambda      !< λ
    real(RNP), contiguous,           intent(in)  :: nu(:,:,:,:) !< ν
    real(RNP), contiguous,           intent(in)  :: u (:,:,:,:) !< operand
    real(RNP), contiguous,           intent(out) :: r (:,:,:,:) !< result
    real(RNP), contiguous, optional, intent(in)  :: f (:,:,:,:) !< RHS
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    !< boundary values

    if (this % sem % mesh % regular) then
      call Apply_RV(this, lambda, nu, u, r, f, bv)
    else
    ! NOT YET SUPPORTED
    ! call Apply_DV(this, lambda, nu, u, r, f, bv)
    end if

  end subroutine Apply_V

  !-----------------------------------------------------------------------------
  !> Conjugate gradient method with constant ν

  subroutine CG_Method_C(this, lambda, nu, u, f, bv, i_max, r_red, r_max, ni)

    class(DG_EllipticOperator_3D), intent(in) :: this

    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:)

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
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:)

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
    real(RNP),             intent(in)    :: nu         !< ν
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:)

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
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side

    !> boundary values matching the specified conditions
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:)

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
    real(RNP),             intent(in)    :: lambda     !< λ
    real(RNP),             intent(in)    :: nu         !< ν
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f(:,:,:,:) !< right hand side
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:) !< BC
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
    real(RNP), contiguous, intent(in)    :: nu(:,:,:,:) !< ν
    real(RNP), contiguous, intent(inout) :: u (:,:,:,:) !< approximate solution
    real(RNP), contiguous, intent(in)    :: f (:,:,:,:) !< right hand side
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv(:) !< BC
    integer,               intent(in)    :: i_max   !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red   !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max   !< max admissible residual
    integer,     optional, intent(out)   :: ni      !< executed num iterations

    call SchwarzPCG_Method_X( this, lambda, nu_v = nu, u = u, f = f, bv = bv &
                            , i_max = i_max, r_red = r_red, r_max = r_max    &
                            , ni = ni                                        )

  end subroutine SchwarzPCG_Method_V

  !=============================================================================
  ! Helper routines

  !-----------------------------------------------------------------------------
  !> Spectral diffusivity coefficient νˢ

  pure real(RNP) function SpectralDiffusivity(this, nu) result(nu_s)
    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP), intent(in) :: nu

    associate(r_nu_s => this % r_nu_s)
      if (r_nu_s < 1) then
        nu_s = nu * r_nu_s / (1 + r_nu_s)
      else
        nu_s = nu
      end if
    end associate

  end function SpectralDiffusivity

  !-----------------------------------------------------------------------------
  !> Sum of diffusivity coefficients, ν + νˢ

  pure real(RNP) function TotalDiffusivity(this, nu) result(nu_t)
    class(DG_EllipticOperator_3D), intent(in) :: this
    real(RNP), intent(in) :: nu

    associate(r_nu_s => this % r_nu_s)
      if (r_nu_s < 1) then
        nu_t = nu / (1 + r_nu_s)
      else
        nu_t = nu
      end if
    end associate

  end function TotalDiffusivity

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

  subroutine EnforceBoundaryConditions(diffusion_op, bv, jmp_u, avg_q)
    class(DG_EllipticOperator_3D), intent(in) :: diffusion_op
    class(SpectralElementBoundaryVariable_3D), optional, intent(in) :: bv(:)
    real(RNP), contiguous, intent(inout) :: jmp_u(:,:,:,:) !< trace of u
    real(RNP), contiguous, intent(inout) :: avg_q(:,:,:,:) !< trace of ν du/dn

    logical :: has_bv
    integer :: b, e, f, l

    has_bv = present(bv)

    associate(boundary => diffusion_op % sem % mesh % boundary)

      do b = 1, size(boundary)

        select case(diffusion_op % bc(b))

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
