!> summary:  Incompressible Navier-Stokes DG-SEM operator
!> author:   Joerg Stiller
!> date:     2021/12/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Operator__3D
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI

  use Standard_Operators__1D
  use Embedded_Interpolation__1D
  use DG__Element_Operators__1D
  use DG__Elliptic_Operator__3D

  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Spectral_Element_Boundary_Variable__3D

  implicit none
  private

  public :: INS_Operator_3D
  public :: INS_Options_3D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for incompressible Navier-Stokes problems

  type INS_Operator_3D

    real(RNP) :: mu !< bulk viscosity, μ = ζ/ρ, may be variable in future

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

!?  type(DG_EllipticnOperator_3D) :: diff_op_p | pressure_op
!< pressure laplacian operator

  contains

    generic :: Init => Init_INS_Operator_3D
    procedure, private :: Init_INS_Operator_3D

    procedure :: SetVelocityBC

    generic :: GetDiffusionTerm => GetDiffusionTerm_C
    procedure, private :: GetDiffusionTerm_C ! diffusion with constant viscosity

  end type INS_Operator_3D

  ! constructor interface
  interface INS_Operator_3D
    procedure New_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for INS_Operator_3D initialization

  type INS_Options_3D
    real(RNP) :: mu = -2 !< bulk viscosity, μ = ζ/ρ
    type(DG_ElementOptions_1D) :: opt_v !< DG operator options for v
    type(DG_ElementOptions_1D) :: opt_p !< DG operator options for p
    type(StandardOperatorOptions_1D) :: opt_q !< quadrature opts for convection
  contains
    procedure :: Bcast => Bcast_INS_Options_3D
  end type INS_Options_3D

  !=============================================================================
  ! Module procedures

  interface

    !---------------------------------------------------------------------------
    !> Diffusion term with constant viscosity on irregular (deformed) mesh

    module subroutine GetDiffusionTerm_DC(this, bc, nu, v, vp, sp, F_d)
      class(INS_Operator_3D), intent(in)    :: this
      character,              intent(in)    :: bc(:)
      real(RNP),              intent(in)    :: nu
      real(RNP), contiguous,  intent(in)    :: v(:,:,:,:,:)
      real(RNP), contiguous,  intent(inout) :: vp(:,:,:,:,:)
      real(RNP), contiguous,  intent(inout) :: sp(:,:,:,:,:)
      real(RNP), contiguous,  intent(out)   :: F_d(:,:,:,:,:)
    end subroutine GetDiffusionTerm_DC

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of INS_Operator_3D

  type(INS_Operator_3D) function New_INS_Operator_3D(opt, mesh) result(this)
    type(INS_Options_3D), intent(in) :: opt !< options
    type(Mesh_3D), intent(in) :: mesh !< local mesh partition, will be copied

    call Init_INS_Operator_3D(this, opt, mesh)

  end function New_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of INS_Operator_3D

  subroutine Init_INS_Operator_3D(this, opt, mesh)
    class(INS_Operator_3D) , intent(inout) :: this !< new Navier-Stokes operator
    type(INS_Options_3D)   , intent(in)    :: opt  !< options
    type(Mesh_3D), optional, intent(in)    :: mesh !< local mesh partition

    this % mu = opt % mu

    this % eop_v = DG_ElementOperators_1D(opt % opt_v)
    this % eop_p = DG_ElementOperators_1D(opt % opt_p)
    this % sop_q = StandardOperators_1D  (opt % opt_q)

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

  subroutine SetVelocityBC(this, bc_v, bv_v, tr_v, tr_s)
    class(INS_Operator_3D) , intent(in) :: this
    !< Navier-Stokes operator
    character, intent(in) :: bc_v(:)
    !< boundary condition types
    class(SpectralElementBoundaryVariable_3D), intent(in) :: bv_v(:)
    !< boundary values
    real(RNP), optional, intent(inout) :: tr_v(:,:,:,:,:)
    !< velocity on element faces
    real(RNP), optional, intent(inout) :: tr_s(:,:,:,:,:)
    !< stress vector on element faces

    integer :: b

    if (present(tr_v)) then
      do b = 1, this % mesh % n_bound
        if (bc_v(b) == 'D') then
          call bv_v(b) % CopyToElementFaceVariable(tr_v)
        end if
      end do
    end if

    if (present(tr_s)) then
      return ! nothing to do yet
    end if

  end subroutine SetVelocityBC

  !-----------------------------------------------------------------------------
  !> Diffusion term with constant viscosity

  subroutine GetDiffusionTerm_C(this, bc, nu, v, vp, sp, F_d)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    character, intent(in) :: bc(:)
    !< boundary conditions ∈ {'D','P'}

    real(RNP), intent(in) :: nu
    !< kinematic shear viscosity

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

    if (this % mesh % regular) then
      ! not implemented yet
      F_d = 0
    else
      call GetDiffusionTerm_DC(this, bc, nu, v, vp, sp, F_d)
    end if

  end subroutine GetDiffusionTerm_C

  !=============================================================================
  ! Type-bound procedures of INS_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type INS_Options_3D

   subroutine Bcast_INS_Options_3D(this, root, comm)
    class(INS_Options_3D), intent(inout) :: this
    integer,               intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),        intent(in)    :: comm !< MPI communicator

    call XMPI_Bcast(this % mu, root, comm)

    call this % opt_v % Bcast(root, comm)
    call this % opt_p % Bcast(root, comm)
    call this % opt_q % Bcast(root, comm)

  end subroutine Bcast_INS_Options_3D

  !=============================================================================

end module INS__Operator__3D
