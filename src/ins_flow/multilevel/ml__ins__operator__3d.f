!> summary:  3D multilevel DG operator for incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2025/03/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__INS__Operator__3D
  use Execution_Control
  use XMPI
  use INS__Problem__3D
  use INS__Operator__3D
  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__DG__Elliptic_Solver__3D

  implicit none
  private

  public :: ML_INS_Operator_3D
  public :: ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM operators for incompressible Navier-Stokes problems

  type ML_INS_Operator_3D
    type(ML_Mesh_3D), pointer :: ml_mesh
    type(ML_MeshOperators_3D) :: ml_op_u
    type(ML_MeshOperators_3D), allocatable :: ml_op_p
    type(ML_DG_EllipticSolver_3D), allocatable :: ml_solver_p
    class(INS_Problem_3D), pointer :: problem
    type(INS_Operator_3D), allocatable :: ins_op(:)
  contains
    procedure :: Init_ML_INS_Operator_3D
  end type ML_INS_Operator_3D

  ! constructor interface
  interface ML_INS_Operator_3D
    procedure New_ML_INS_Operator_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM INS operator options

  type ML_INS_OperatorOptions_3D
    integer, allocatable :: po(:) !< sequence of polynomial orders
    logical :: mixed  = .true.    !< T/F: use mixed/equal order for (v,p)
    integer :: smooth = 0         !< FC discontinuity smoothing {0,1,2}
    type(ML_DG_EllipticOptions_3D) :: ml_solver_p !< ML pressure solver options
    type(INS_OperatorOptions_3D)   :: ins         !< INS operator options
  contains
    procedure :: Bcast => Bcast_ML_INS_OperatorOptions_3D
  end type ML_INS_OperatorOptions_3D

contains

  !=============================================================================
  ! Type-bound procedures of INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Constructor of ML_INS_Operator_3D

  function New_ML_INS_Operator_3D(opt, ml_mesh, problem) result(this)
    class(ML_INS_OperatorOptions_3D), intent(in) :: opt
    class(ML_Mesh_3D), target, intent(in) :: ml_mesh
    class(INS_Problem_3D), target, intent(in) :: problem
    type(ML_INS_Operator_3D) :: this

    call Init_ML_INS_Operator_3D(this, opt, ml_mesh, problem)

  end function New_ML_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of ML_INS_Operator_3D

  subroutine Init_ML_INS_Operator_3D(this, opt, ml_mesh, problem)
    class(ML_INS_Operator_3D), intent(inout) :: this
    class(ML_INS_OperatorOptions_3D), intent(in) :: opt
    class(ML_Mesh_3D), target, intent(in) :: ml_mesh
    class(INS_Problem_3D), target, intent(in) :: problem

    integer :: l

    ! sanity check .............................................................

    if (.not. allocated(opt%po)) then
      call Error('Init_ML_INS_Operator_3D', 'opt%po not allocated')
    end if

    if (.not. allocated(ml_mesh%mesh)) then
      call Error('Init_ML_INS_Operator_3D', 'ml_mesh%mesh not allocated')
    end if

    if (size(opt%po) < size(ml_mesh%mesh)) then
      call Error('Init_ML_INS_Operator_3D', 'opt%po does not match ml_mesh')
    end if

    ! creating components ......................................................

    associate(po => opt%po, smooth => opt%smooth)

      this % ml_mesh => ml_mesh
      this % problem => problem

      ! multilevel spectral-element mesh operators and metrics
      this % ml_op_u = ML_MeshOperators_3D(ml_mesh, po, 'L', smooth)
      if (opt%mixed) then
        this % ml_op_p = ML_MeshOperators_3D(ml_mesh, max(po-1,1), 'L', smooth)
      end if

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

    end associate

  end subroutine Init_ML_INS_Operator_3D

  !=============================================================================
  ! Type-bound procedures of ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type ML_INS_OperatorOptions_3D

   subroutine Bcast_ML_INS_OperatorOptions_3D(this, root, comm)
    class(ML_INS_OperatorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    integer :: n, rank

    call MPI_Comm_rank(comm, rank)

    if (rank == root) then
      n = size(this%po)
    end if

    call XMPI_Bcast(n, root, comm)

    if (rank /= root) then
      if (allocated(this%po)) then
        deallocate(this%po)
      end if
      allocate(this%po(n))
    end if

    call XMPI_Bcast(this % po    , root, comm)
    call XMPI_Bcast(this % mixed , root, comm)
    call XMPI_Bcast(this % smooth, root, comm)

    call this % ml_solver_p % Bcast(root, comm)
    call this % ins         % Bcast(root, comm)

  end subroutine Bcast_ML_INS_OperatorOptions_3D

  !=============================================================================

end module ML__INS__Operator__3D
