!> summary:  Incompressible Navier-Stokes DG-SEM operator
!> author:   Joerg Stiller
!> date:     2021/12/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Operator__3D
  use Kind_Parameters
  use Execution_Control
  use XMPI

  use Standard_Operators__1D
  use DG__Element_Operators__1D
  use Embedded_Interpolation__1D

  use Mesh__3D
  use Spectral_Element_Mesh__3D

  implicit none
  private

  public :: INS_Operator_3D
  public :: INS_Options_3D

  !-----------------------------------------------------------------------------
  !> DG-SEM mesh and operators for incompressible Navier-Stokes problems

  type INS_Operator_3D

    real(RNP) :: chi !< non-dimensional first Lamé coefficient

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

  contains

    procedure :: Init_INS_Operator_3D

  end type INS_Operator_3D

  ! constructor interface
  interface INS_Operator_3D
    procedure New_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for INS_Operator_3D initialization

  type INS_Options_3D
    real(RNP) :: chi = -2 !< non-dimensional first Lamé coefficient
    type(DG_ElementOptions_1D) :: opt_v !< DG operator options for v
    type(DG_ElementOptions_1D) :: opt_p !< DG operator options for p
    type(StandardOperatorOptions_1D) :: opt_q !< quadrature opts for convection
  contains
    procedure :: Bcast => Bcast_INS_Options_3D
  end type INS_Options_3D

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of INS_Operator_3D

  type(INS_Operator_3D) function New_INS_Operator_3D(opt, mesh) result(this)
    class(INS_Options_3D), intent(in) :: opt !< options
    type(Mesh_3D), intent(in) :: mesh !< local mesh partition, will be copied

    call Init_INS_Operator_3D(this, opt, mesh)

  end function New_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of INS_Operator_3D

  subroutine Init_INS_Operator_3D(this, opt, mesh)
    class(INS_Operator_3D) , intent(inout) :: this !< new Navier-Stokes operator
    class(INS_Options_3D)  , intent(in)    :: opt  !< options
    type(Mesh_3D), optional, intent(in)    :: mesh !< local mesh partition

    this % chi = opt % chi

    this % eop_v = DG_ElementOperators_1D(opt % opt_v)
    this % eop_p = DG_ElementOperators_1D(opt % opt_p)
    this % sop_q = StandardOperators_1D  (opt % opt_q)

    this % iop_vp = EmbeddedInterpolation_1D(this % eop_v, this % eop_p % x)
    this % iop_vq = EmbeddedInterpolation_1D(this % eop_v, this % sop_q % x)
    this % iop_pv = EmbeddedInterpolation_1D(this % eop_p, this % eop_v % x)

    if (present(mesh)) then
      this % mesh = mesh
    else if (this % mesh % n_part < 1) then
      call Error('Init_INS_Operator_3D', &
                 'this%mesh must be initialized or argument mesh given', &
                 'INS__Operator__3D')
    end if

    this % sem_v = SpectralElementMesh_3D(mesh, this%eop_v%po)
    this % sem_p = SpectralElementMesh_3D(mesh, this%eop_p%po)

    if (this%sop_q%po /= this%eop_v%po .or. this%sop_q%basis /= 'L' ) then
      this%sem_q = SpectralElementMesh_3D(mesh, this%sop_q%po, this%sop_q%basis)
    end if

  end subroutine Init_INS_Operator_3D

  !=============================================================================
  ! Type-bound procedures of INS_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type INS_Options_3D

   subroutine Bcast_INS_Options_3D(this, root, comm)
    class(INS_Options_3D), intent(inout) :: this
    integer,               intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),        intent(in)    :: comm !< MPI communicator

    call XMPI_Bcast(this % chi, root, comm)

    call this % opt_v % Bcast(root, comm)
    call this % opt_p % Bcast(root, comm)
    call this % opt_q % Bcast(root, comm)

  end subroutine Bcast_INS_Options_3D

  !=============================================================================

end module INS__Operator__3D
