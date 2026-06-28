!> summary:  3D multilevel DG operator for incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2025/03/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__INS__Operator__3D
  use Kind_Parameters
  use Execution_Control
  use Logging_Levels
  use XMPI
  use INS__Problem__3D
  use INS__Operator__3D
  use Volume_Integrals__3D
  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__DG__Elliptic_Solver__3D

  implicit none
  private

  public :: ML_INS_Operator_3D
  public :: ML_INS_OperatorOptions_3D

  !=============================================================================
  ! Types

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM operators for incompressible Navier-Stokes problems

  type ML_INS_Operator_3D

    type(ML_Mesh_3D), pointer :: ml_mesh
    type(ML_MeshOperators_3D) :: ml_op_u
    type(ML_MeshOperators_3D), allocatable :: ml_op_p
    type(ML_DG_EllipticSolver_3D), allocatable :: ml_solver_p
    class(INS_Problem_3D), pointer :: problem
    type(INS_Operator_3D), allocatable :: ins_op(:)

    character :: fc_project !< fine-to-coarse projection method
    logical   :: dc_bulk    !< defect correction w/wo bulk diffusion

  contains
    procedure :: Init_ML_INS_Operator_3D
    procedure :: MG_Stokes_Start
    procedure :: MG_Stokes_Cycle
    procedure :: CalibratePressure
  end type ML_INS_Operator_3D

  ! constructor interface
  interface ML_INS_Operator_3D
    procedure New_ML_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM INS operator options

  type ML_INS_OperatorOptions_3D

    logical   :: mixed      = .true.  !< T/F: use mixed/equal order for (v,p)
    integer   :: fc_smooth  =  0      !< fine-to-coarse jump smoothing {0,1,2}
    character :: fc_project = 'I'     !< fine-to-coarse projection method {I,P}
    logical   :: dc_bulk    = .false. !< defect correction w/wo bulk diffusion

    type(ML_DG_EllipticOptions_3D) :: ml_solver_p !< ML pressure solver options
    type(INS_OperatorOptions_3D)   :: ins         !< INS operator options

  contains
    procedure :: Bcast => Bcast_ML_INS_OperatorOptions_3D
  end type ML_INS_OperatorOptions_3D

  !=============================================================================
  ! Interfaces to submodule procedures

  interface

    !---------------------------------------------------------------------------
    !> Computes the pressure and minimizes the divergence of the given velocity

    module subroutine ProjectionStep(this, tau, bv_u, u)
      class(ML_INS_Operator_3D),     intent(in)    :: this
      real(RNP),                     intent(in)    :: tau
      class(ML_BoundaryVariable_3D), intent(in)    :: bv_u
      class(ML_MeshVariable_3D),     intent(inout) :: u
    end subroutine ProjectionStep

    !---------------------------------------------------------------------------
    !> Cascade start procedure for the Stokes part

    module subroutine Stokes_Cascade(this, tau, mu, nu, bv, f, u)
      class(ML_INS_Operator_3D),           intent(in)    :: this
      real(RNP),                           intent(in)    :: tau
      class(ML_MeshVariable_3D), optional, intent(in)    :: mu
      class(ML_MeshVariable_3D), optional, intent(in)    :: nu
      class(ML_BoundaryVariable_3D),       intent(in)    :: bv
      class(ML_MeshVariable_3D),           intent(inout) :: f
      class(ML_MeshVariable_3D),           intent(inout) :: u
    end subroutine Stokes_Cascade

    !---------------------------------------------------------------------------
    !> Cascade and FMG start for the Stokes multigrid solver

    module subroutine MG_Stokes_Start(this, tau, mu, nu, bv, f_d0, f, u, n_cyc)
      class(ML_INS_Operator_3D),     intent(in)    :: this
      real(RNP),                     intent(in)    :: tau
      class(ML_MeshVariable_3D),     intent(in)    :: mu
      class(ML_MeshVariable_3D),     intent(in)    :: nu
      class(ML_BoundaryVariable_3D), intent(in)    :: bv
      class(ML_MeshVariable_3D),     intent(in)    :: f_d0
      class(ML_MeshVariable_3D),     intent(inout) :: f
      class(ML_MeshVariable_3D),     intent(inout) :: u
      integer,             optional, intent(in)    :: n_cyc
    end subroutine MG_Stokes_Start

    !---------------------------------------------------------------------------
    !> Performs one or more FAS-MG V-cycles for the Stokes part

    module subroutine MG_Stokes_Cycle(this, tau, mu, nu, bv, f, u, n_cyc, l_top)
      class(ML_INS_Operator_3D),     intent(in)    :: this
      real(RNP),                     intent(in)    :: tau
      class(ML_MeshVariable_3D),     intent(in)    :: mu
      class(ML_MeshVariable_3D),     intent(in)    :: nu
      class(ML_BoundaryVariable_3D), intent(in)    :: bv
      class(ML_MeshVariable_3D),     intent(inout) :: f
      class(ML_MeshVariable_3D),     intent(inout) :: u
      integer,             optional, intent(in)    :: n_cyc
      integer,             optional, intent(in)    :: l_top
    end subroutine MG_Stokes_Cycle

  end interface

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of ML_INS_Operator_3D

  function New_ML_INS_Operator_3D(ml_mesh, po, problem, opt) result(this)
    class(ML_Mesh_3D), target, intent(in) :: ml_mesh
    integer, intent(in) :: po(size(ml_mesh%mesh))
    class(INS_Problem_3D), target, intent(in) :: problem
    class(ML_INS_OperatorOptions_3D), intent(in) :: opt
    type(ML_INS_Operator_3D) :: this

    call Init_ML_INS_Operator_3D(this, ml_mesh, po, problem, opt)

  end function New_ML_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of ML_INS_Operator_3D

  subroutine Init_ML_INS_Operator_3D(this, ml_mesh, po, problem, opt)
    class(ML_INS_Operator_3D), intent(inout) :: this
    class(ML_Mesh_3D), target, intent(in) :: ml_mesh
    integer, intent(in) :: po(size(ml_mesh%mesh))
    class(INS_Problem_3D), target, intent(in) :: problem
    class(ML_INS_OperatorOptions_3D), intent(in) :: opt

    integer :: l

    ! creating components ......................................................

    this % ml_mesh => ml_mesh
    this % problem => problem

    ! multilevel spectral-element mesh operators and metrics
    this % ml_op_u = ML_MeshOperators_3D( ml_mesh, po, 'L' &
                                        , opt%fc_smooth    )
    if (opt%mixed) then
      this % ml_op_p = ML_MeshOperators_3D( ml_mesh, max(po-1,1), 'L' &
                                          , opt%fc_smooth )
    else
      this % ml_op_p = this % ml_op_u
    end if

    ! fine-to-coarse projection method
    this % fc_project = opt % fc_project

    ! fine-to-coarse defect correction with or without bulk diffusion
    this % dc_bulk = opt % dc_bulk

    ! multilevel pressure solver
    if (opt % ins % pressure_solver == 'MG') then
      if (opt%mixed) then
        this % ml_solver_p = &
            ML_DG_EllipticSolver_3D(this%ml_op_p, opt%ml_solver_p)
      else
        this % ml_solver_p = &
            ML_DG_EllipticSolver_3D(this%ml_op_u, opt%ml_solver_p)
      end if
    else
      if (allocated(this % ml_solver_p)) then
        deallocate(this % ml_solver_p)
      end if
    end if

    ! incompressible Navier-Stokes operator on each level
    allocate(this%ins_op( size(this%ml_mesh%mesh) ))
    do l = 1, size(this%ins_op)
      this % ins_op(l) = INS_Operator_3D( opt%ins                 &
                                        , this % problem          &
                                        , this % ml_op_u % sem(l) &
                                        , this % ml_op_p % sem(l) &
                                        , this % ml_solver_p      &
                                        , level = l               )
    end do

  end subroutine Init_ML_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Remove the mean pressure unless Dirichlet conditions apply
  !>
  !> The multilevel variable `u` must contain the either the flow variables
  !> `[ v(1:3), p, ... ]` or only the pressure `p`.

  subroutine CalibratePressure(this, u)
    class(ML_INS_Operator_3D), intent(in)    :: this !< ML INS operator
    class(ML_MeshVariable_3D), intent(inout) :: u    !< flow variables or p

    real(RNP) :: p_mean, p_mean_l, volume
    integer   :: c, l, l_top

    if (any(this % problem % bc_p == 'D')) return

    ! identify pressure component
    select case(size(u%name))
    case(1)
      c = 1
    case(4:)
      c = 4
    case default
      call Error( 'CalibratePressure'                                   &
                , 'u must contain either the flow variables or p alone' &
                , 'ML__INS__Operator__3D'                               )
    end select

    call this % ml_op_u % Get_Volume(volume)

    l_top = size(u % level)

    p_mean = 0
    do l = 1, l_top
      associate(p_l => u % level(l) % val(:,:,:,:,c))
        call GetVolumeIntegral(this%ins_op(l)%sem_u, p_l, p_mean_l, leaf=.true.)
        p_mean = p_mean + p_mean_l
      end associate
    end do
    p_mean = p_mean / volume

    do l = 1, l_top
      associate(p_l => u % level(l) % val(:,:,:,:,c))
        p_l = p_l - p_mean
      end associate
    end do

  end subroutine CalibratePressure

  !=============================================================================
  ! Type-bound procedures of ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type ML_INS_OperatorOptions_3D

   subroutine Bcast_ML_INS_OperatorOptions_3D(this, root, comm)
    class(ML_INS_OperatorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % mixed     , root, comm)
    call XMPI_Bcast(this % fc_smooth , root, comm)
    call XMPI_Bcast(this % fc_project, root, comm)
    call XMPI_Bcast(this % dc_bulk   , root, comm)

    call this % ml_solver_p % Bcast(root, comm)
    call this % ins         % Bcast(root, comm)

  end subroutine Bcast_ML_INS_OperatorOptions_3D

  !=============================================================================

end module ML__INS__Operator__3D
