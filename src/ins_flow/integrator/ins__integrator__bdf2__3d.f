!> summary:  BDF2 method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/11/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Integrator__BDF2__3D
  use Kind_Parameters
  use Constants
  use XMPI
  use Boundary_Variable__3D
  use INS__Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D
  use INS__Integrator__BDF2__PrepStep__3D

  implicit none
  private

  public :: INS_Integrator_BDF2_3D
  public :: INS_Integrator_BDF2_Options_3D

  !-----------------------------------------------------------------------------
  !> BDF2 method for incompressible flows

  type, extends(INS_Integrator_3D) :: INS_Integrator_BDF2_3D
  contains
    procedure, non_overridable :: Init_INS_Integrator_BDF2_3D
    procedure :: TimeStep
  end type INS_Integrator_BDF2_3D

  ! constructor
  interface INS_Integrator_BDF2_3D
    module procedure New_INS_Integrator_BDF2_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing BDF2 time-integrator options (none, so far)

  type, extends(INS_IntegratorOptions_3D) :: &
    INS_Integrator_BDF2_Options_3D
  contains
    procedure :: Bcast => Bcast_INS_Integrator_BDF2_Options
  end type INS_Integrator_BDF2_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_Integrator_BDF2_3D

  function New_INS_Integrator_BDF2_3D(problem, ins_op, opt) result(this)
    class(INS_Problem_3D),                 intent(in) :: problem
    class(INS_Operator_3D),                intent(in) :: ins_op
    class(INS_Integrator_BDF2_Options_3D), intent(in) :: opt
    type(INS_Integrator_BDF2_3D) :: this

    call Init_INS_Integrator_BDF2_3D(this, problem, ins_op, opt)

  end function New_INS_Integrator_BDF2_3D

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_Integrator_BDF2_3D object

  subroutine Init_INS_Integrator_BDF2_3D(this, problem, ins_op, opt)
    class(INS_Integrator_BDF2_3D),         intent(inout) :: this
    class(INS_Problem_3D),                 intent(in)    :: problem
    class(INS_Operator_3D),                intent(in)    :: ins_op
    class(INS_Integrator_BDF2_Options_3D), intent(in)    :: opt

    ! intialize parent type
    call this % Init_INS_Integrator_3D(problem, ins_op, opt)

    this % name = 'BDF2 method'

  end subroutine Init_INS_Integrator_BDF2_3D

  !-----------------------------------------------------------------------------
  !> Execution of a BDF2 time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_Integrator_BDF2_3D), intent(inout) :: this
    real(RNP),             intent(inout) :: t            !< time t₀ → t
    real(RNP),             intent(in)    :: dt           !< step size ∆t = t-t₀
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! internal variables .......................................................

    real(RNP), allocatable, save :: f  (:,:,:,:,:) ! unweighted RHS
    real(RNP), allocatable, save :: f_d(:,:,:,:,:) ! unweighted diffusion term
    real(RNP), allocatable, save :: mu (:,:,:,:)   ! bulk diffusivity μ
    real(RNP), allocatable, save :: nu (:,:,:,:)   ! shear diffusivity ν

    ! saved at time t₁ = t₀-∆t
    real(RNP), allocatable, save :: u1  (:,:,:,:,:) ! flow variables
    real(RNP), allocatable, save :: f_c1(:,:,:,:,:) ! unweighted convection term
    real(RNP), allocatable, save :: f_d1(:,:,:,:,:) ! unweighted diffusion term

    type(BoundaryVariable_3D), allocatable, save :: bv_u(:) ! boundary values

    logical, save :: first
    real(RNP) :: tau
    integer   :: b, np

    associate( ins_op => this % ins_op              &
             , mesh   => this % ins_op % mesh       &
             , po     => this % ins_op % eop_u % po &
             , nc     => this % problem % nc        )

      ! initialization .........................................................

      np = po + 1

      !$omp master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (allocated(f)) then
        first = any(shape(f) /= shape(u))
        if (first) then
          deallocate(f, bv_u, f_d, u1, f_c1, f_d1)
          if (allocated(mu)) deallocate(mu)
          if (allocated(nu)) deallocate(nu)
        end if
      else
        first = .true.
      end if

      if (first) then

        allocate(f(np, np, np, mesh%n_elem, nc), source = ZERO)
        allocate(f_d, u1, f_c1, f_d1, source = f)

        if (ins_op % HasVariableViscosity()) then
          allocate(mu(np, np, np, mesh%n_elem))
          allocate(nu(np, np, np, mesh%n_elem))
        end if

        allocate(bv_u(mesh % n_bound))
        do b = 1, mesh % n_bound
          call bv_u(b) % Init(mesh % boundary(b), po, nc = nc)
        end do

      end if

      !$omp end master !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!omp barrier    !! not required because of barrier in BDF2_PrepStep_3D

      ! time step ..............................................................

      call INS_Integrator_BDF2_PrepStep_3D( ins_op, t, dt, u, u1, f_c1, f_d1 &
                                          , tau, f, f_d, bv_u, mu, nu, first )

      call ins_op % StokesSolver(tau, f, bv_u, mu, nu, u, f_d0 = f_d)

      ! cleanup ................................................................

      !$omp master
      t = t + dt
      !$omp end master

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(f, bv_u, f_d, u1, f_c1, f_d1)
      if (allocated(mu)) deallocate(mu)
      if (allocated(nu)) deallocate(nu)
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================
  ! TBP of INS_Integrator_BDF2_Options_3D

  !-----------------------------------------------------------------------------
  !> MPI broadcasting of BDF2 time-integrator options

  subroutine Bcast_INS_Integrator_BDF2_Options(this, root, comm)
    class(INS_Integrator_BDF2_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call this % INS_IntegratorOptions_3D % Bcast(root, comm)

  end subroutine Bcast_INS_Integrator_BDF2_Options

  !=============================================================================

end module INS__Integrator__BDF2__3D
