!> summary:  Incompressible Navier-Stokes DG-SEM operator
!> author:   Joerg Stiller
!> date:     2021/12/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - clean treatment of cases: regular/deformed mesh, constant/variable ν
!>   - extension of velocity boundary conditions
!===============================================================================

module INS__Operator__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI

  use TPO__INS_Convection__3D

  use Standard_Operators__1D
  use Embedded_Interpolation__1D
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D

  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  use INS__Problem__3D

  implicit none
  private

  public :: INS_Operator_3D
  public :: INS_Options_3D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for incompressible Navier-Stokes problems

  type INS_Operator_3D

    real(RNP) :: mu_0 !< bulk viscosity,  μ = ζ/ρ
    real(RNP) :: nu_0 !< shear viscosity, ν = η/ρ

    character, allocatable :: bc_v(:) !< velocity BC, copied from problem
    character, allocatable :: bc_p(:) !< pressure BC

    type(DG_ElementOperators_1D) :: eop_v !< DG operators for v
    type(DG_ElementOperators_1D) :: eop_p !< DG operators for p
    type(StandardOperators_1D)   :: sop_q !< quadrature ops for convection

    type(EmbeddedInterpolation_1D) :: iop_vp !< interpolation from v to p points
    type(EmbeddedInterpolation_1D) :: iop_vq !< interpolation from v to q points
    type(EmbeddedInterpolation_1D) :: iop_pv !< interpolation from p to v points

    type(Mesh_3D)                :: mesh  !< local mesh partition
    type(SpectralElementMesh_3D) :: sem_v !< mesh + metrics for v
    type(SpectralElementMesh_3D) :: sem_p !< mesh + metrics for p
    type(SpectralElementMesh_3D) :: sem_q !< mesh + metrics for convection

    type(DG_EllipticOperator_3D) :: laplacian_p  !< negative Laplacian for p
    type(DG_SchwarzOperator_3D)  :: schwarz_v(3) !< Schwarz operators for v

  contains

    procedure :: Init_INS_Operator_3D
    procedure :: SetVelocityBC
    procedure :: PressureSolver
    procedure :: GetConvectionTerm
    generic   :: GetDiffusionTerm       => GetDiffusionTerm_C
    generic   :: ApplyDiffusionOperator => ApplyDiffusionOperator_C
    generic   :: GetDiffusionResidual   => GetDiffusionResidual_C
    procedure :: DiffusionSolver

    procedure, private :: GetDiffusionTerm_C
    procedure, private :: ApplyDiffusionOperator_C
    procedure, private :: GetDiffusionResidual_C

  end type INS_Operator_3D

  ! constructor interface
  interface INS_Operator_3D
    procedure New_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for INS_Operator_3D initialization

  type INS_Options_3D
    real(RNP) :: mu_0 = -2 !< bulk viscosity,  μ = ζ/ρ
    type(DG_ElementOptions_1D) :: eop_v !< DG operator options for v
    type(DG_ElementOptions_1D) :: eop_p !< DG operator options for p
    type(StandardOperatorOptions_1D) :: sop_q !< quadrature opts for convection
    type(DG_SchwarzOptions_3D) :: schwarz_p !< Schwarz options for p-solver
    type(DG_SchwarzOptions_3D) :: schwarz_v !< Schwarz options for v-solver
  contains
    procedure :: Bcast => Bcast_INS_Options_3D
  end type INS_Options_3D

  !=============================================================================
  ! Module procedures

  interface

    !---------------------------------------------------------------------------
    !> Projection-based pressure solver

    module subroutine PressureSolver( this, tau, bv_v, v, f, bv_p, p &
                                    , i_max, r_red, r_max, ni        )

      class(INS_Operator_3D),                    intent(in)    :: this
      real(RNP),                                 intent(in)    :: tau
      class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv_v(:)
      real(RNP), contiguous,                     intent(in)    :: v(:,:,:,:,:)
      real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(inout) :: bv_p(:)
      real(RNP), contiguous,                     intent(inout) :: p(:,:,:,:)
      integer,                                   intent(in)    :: i_max
      real(RNP),                       optional, intent(in)    :: r_red
      real(RNP),                       optional, intent(in)    :: r_max
      integer,                         optional, intent(out)   :: ni

    end subroutine PressureSolver

    !---------------------------------------------------------------------------
    !> Diffusion term with constant viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_DC(this, v, vp, sp, F_d)
      class(INS_Operator_3D), intent(in)    :: this
      real(RNP), contiguous,  intent(in)    :: v(:,:,:,:,:)
      real(RNP), contiguous,  intent(inout) :: vp(:,:,:,:,:)
      real(RNP), contiguous,  intent(inout) :: sp(:,:,:,:,:)
      real(RNP), contiguous,  intent(out)   :: F_d(:,:,:,:,:)
    end subroutine GetDiffusionTerm_DC

    !---------------------------------------------------------------------------
    !> Homogeneous diffusion operator with constant viscosity

    module subroutine ApplyDiffusionOperator_C(this, tau, v, r)
      class(INS_Operator_3D), intent(in)  :: this
      real(RNP),              intent(in)  :: tau
      real(RNP), contiguous,  intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,  intent(out) :: r(:,:,:,:,:)
    end subroutine ApplyDiffusionOperator_C

    !---------------------------------------------------------------------------
    !> Diffusion residual with constant viscosity

    module subroutine GetDiffusionResidual_C(this, tau, f, bv_v, v, r)
      class(INS_Operator_3D), intent(in)  :: this
      real(RNP),              intent(in)  :: tau
      real(RNP), contiguous, intent(in)   :: f(:,:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
      real(RNP), contiguous,  intent(in)  :: v(:,:,:,:,:)
      real(RNP), contiguous,  intent(out) :: r(:,:,:,:,:)
    end subroutine GetDiffusionResidual_C

    !---------------------------------------------------------------------------
    !> Diffusion solver

    module subroutine DiffusionSolver( this, tau, f, bv_v, v   &
                                     , i_max, r_red, r_max, ni )

      class(INS_Operator_3D),                    intent(in)    :: this
      real(RNP),                                 intent(in)    :: tau
      real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:,:)
      class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv_v(:)
      real(RNP), contiguous,                     intent(inout) :: v(:,:,:,:,:)
      integer,                                   intent(in)    :: i_max
      real(RNP),                       optional, intent(in)    :: r_red
      real(RNP),                       optional, intent(in)    :: r_max
      integer,                         optional, intent(out)   :: ni

    end subroutine DiffusionSolver

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of INS_Operator_3D

  function New_INS_Operator_3D(opt, problem, mesh) result(this)
    class(INS_Options_3D), intent(in) :: opt     !< options
    class(INS_Problem_3D), intent(in) :: problem !< INS flow problem
    type(Mesh_3D), intent(in) :: mesh !< local mesh partition, will be copied
    type(INS_Operator_3D) :: this

    call Init_INS_Operator_3D(this, opt, problem, mesh)

  end function New_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of INS_Operator_3D

  subroutine Init_INS_Operator_3D(this, opt, problem, mesh)
    class(INS_Operator_3D),  intent(inout) :: this    !< new INS operator
    class(INS_Options_3D),   intent(in)    :: opt     !< options
    class(INS_Problem_3D),   intent(in)    :: problem !< INS flow problem
    type(Mesh_3D), optional, intent(in)    :: mesh    !< local mesh partition

    integer :: b, d

    this % mu_0 = opt % mu_0
    this % nu_0 = problem % nu_ref
    this % bc_v = problem % bc_v

    this % eop_v = DG_ElementOperators_1D(opt % eop_v)
    this % eop_p = DG_ElementOperators_1D(opt % eop_p)
    this % sop_q = StandardOperators_1D  (opt % sop_q)

    this % iop_vp = EmbeddedInterpolation_1D(this % eop_v, this % eop_p % x)
    this % iop_vq = EmbeddedInterpolation_1D(this % eop_v, this % sop_q % x)
    this % iop_pv = EmbeddedInterpolation_1D(this % eop_p, this % eop_v % x)

    if (present(mesh)) then
      this % mesh = mesh
    else if (this % mesh % n_parts < 1) then
      call Error('Init_INS_Operator_3D', &
                 'this%mesh must be initialized or argument mesh given', &
                 'INS__Operator__3D')
    end if

    this % sem_v = SpectralElementMesh_3D( this % mesh, this % eop_v % po )
    this % sem_p = SpectralElementMesh_3D( this % mesh, this % eop_p % po )
    this % sem_q = SpectralElementMesh_3D( this % mesh          &
                                         , this % sop_q % po    &
                                         , this % sop_q % basis )

    ! pressure BC
    allocate(this % bc_p(this % mesh % n_bound), source = '')
    do b = 1, size(this % bc_p)
      select case(this % bc_v(b))
      case('D')
        this % bc_p(b) = 'N'
      case('P')
        this % bc_p(b) = 'P'
      end select
    end do

    ! pressure operator
    this % laplacian_p = DG_EllipticOperator_3D( sem         = this % sem_p     &
                                               , dg_opt      = opt  % eop_p     &
                                               , schwarz_opt = opt  % schwarz_p &
                                               , bc          = this % bc_p      )

    ! Schwarz operators for the viscous diffusion solver
    do d = 1, 3
      this % schwarz_v(d) = DG_SchwarzOperator_3D( opt  % schwarz_v &
                                                 , this % eop_v     &
                                                 , this % mesh      &
                                                 , this % bc_v      )
    end do

  end subroutine Init_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Inject the velocity boundary conditions into trace variables
  !>
  !> The boundary variable `bv_v` is expected to contain the velocity boundary
  !> values in the first three components. Velocity values will be injected
  !> into `tr_v` and stresses into `tr_s`, if present.
  !>
  !> The traces are stored as element face variables defined as
  !>
  !>     tr_v(np,np,6,nl,3)
  !>     tr_s(np,np,6,nl,3)
  !>
  !> where `np` is the number of velocity element points per direction and
  !> `nl` is the number of elements, possibly including the ghosts.
  !>
  !> @note
  !> So far, only Dirichlet conditions ('D') are supported.

  subroutine SetVelocityBC(this, bv_v, tr_v, tr_s)
    class(INS_Operator_3D) , intent(in) :: this
    !< Navier-Stokes operator
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
    !< boundary values
    real(RNP), optional, intent(inout) :: tr_v(:,:,:,:,:)
    !< velocity on element faces
    real(RNP), optional, intent(inout) :: tr_s(:,:,:,:,:)
    !< stress vector on element faces

    integer :: b

    if (present(tr_v)) then
      do b = 1, this % mesh % n_bound
        if (this % bc_v(b) == 'D') then
          call bv_v(b) % CopyToTraceVariable(tr_v)
        end if
      end do
    end if

    if (present(tr_s)) then
      return ! nothing to do yet
    end if

  end subroutine SetVelocityBC

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity

  subroutine GetConvectionTerm(this, v, vp, F_c)
    class(INS_Operator_3D), intent(in) :: this
    real(RNP), contiguous, intent(in)  :: v(:,:,:,:,:)   !< velocity
    real(RNP), contiguous, intent(in)  :: vp(:,:,:,:,:)  !< outer velocity v⁺
    real(RNP), contiguous, intent(out) :: F_c(:,:,:,:,:) !< convection term

!   if (this % mesh % regular) then
!     not implemented yet
!   else
      call TPO_INS_Convection( nv   = this % eop_v % po + 1       &
                             , nq   = this % sop_q % po + 1       &
                             , ne   = this % mesh % n_elem        &
                             , D_v  = this % eop_v  % D           &
                             , I_vq = this % iop_vq % A           &
                             , w_q  = this % sop_q  % w           &
                             , Jd_q = this % sem_q % metrics % Jd &
                             , Ji_q = this % sem_q % metrics % Ji &
                             , a_q  = this % sem_q % metrics % a  &
                             , n_q  = this % sem_q % metrics % n  &
                             , v    = v                           &
                             , vp   = vp                          &
                             , F_c  = F_c                         )
!   end if

  end subroutine GetConvectionTerm

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity

  subroutine GetDiffusionTerm_C(this, v, vp, sp, F_d)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< velocity (np,np,np,ne,3)

    real(RNP), contiguous, intent(inout) :: vp(:,:,:,:,:)
    !< exterior velocity traces (np,np,6,ne,3)
    !<   - in:  v   at Dirichlet faces, nn⋅v at free-slip faces, undefined else
    !<   - out: v⁺  at all element faces

    real(RNP), contiguous, intent(inout) :: sp(:,:,:,:,:)
    !< viscous flux traces (np,np,6,ne,3)
    !<   - in:  s  = n⋅τ   at free-slip or traction faces, undefined else
    !<   - out: s⁺ = n⋅τ⁺  at all element faces

    real(RNP), contiguous, intent(out) :: F_d(:,:,:,:,:)
    !< diffusion term (np,np,np,ne,3)

!   if (this % mesh % regular) then
!     not implemented yet
!   else
      call GetDiffusionTerm_DC(this, v, vp, sp, F_d)
!   end if

  end subroutine GetDiffusionTerm_C

  !=============================================================================
  ! Type-bound procedures of INS_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type INS_Options_3D

   subroutine Bcast_INS_Options_3D(this, root, comm)
    class(INS_Options_3D), intent(inout) :: this
    integer,               intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),        intent(in)    :: comm !< MPI communicator

    call XMPI_Bcast(this % mu_0, root, comm)

    call this % eop_v     % Bcast(root, comm)
    call this % eop_p     % Bcast(root, comm)
    call this % sop_q     % Bcast(root, comm)
    call this % schwarz_p % Bcast(root, comm)
    call this % schwarz_v % Bcast(root, comm)

  end subroutine Bcast_INS_Options_3D

  !=============================================================================

end module INS__Operator__3D
