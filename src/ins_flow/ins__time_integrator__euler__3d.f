!> summary:  Euler method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @note
!>   - The present version implements the IMEX Euler method
!===============================================================================

module INS__Time_Integrator__Euler__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use XMPI

  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  use Boundary_Variable__3D

  use INS__Time_Integrator__3D
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_TimeIntegrator_Euler_3D
  public :: INS_TimeIntegrator_Euler_Options_3D

  !-----------------------------------------------------------------------------
  !> Euler method for incompressible flows

  type, extends(INS_TimeIntegrator_3D) :: INS_TimeIntegrator_Euler_3D
    integer   :: i_max_p !< max num p-iterations   in projection step
    integer   :: i_max_v !< max num v-iterations   in projection step
    real(RNP) :: r_red   !< min residual reduction in projection step, if > 0
    real(RNP) :: r_max   !< max residual to reach  in projection step, if > 0
  contains
    procedure, non_overridable :: Init_INS_TimeIntegrator_Euler_3D
    procedure :: TimeStep
  end type INS_TimeIntegrator_Euler_3D

  ! constructor
  interface INS_TimeIntegrator_Euler_3D
    module procedure New_INS_TimeIntegrator_Euler_3D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(INS_TimeIntegratorOptions_3D) :: &
    INS_TimeIntegrator_Euler_Options_3D
    integer   :: i_max_p = 5 !< max num p-iterations   in projection step
    integer   :: i_max_v = 2 !< max num v-iterations   in projection step
    real(RNP) :: r_red   = 0 !< min residual reduction in projection step, if > 0
    real(RNP) :: r_max   = 0 !< max residual to reach  in projection step, if > 0
  contains
    procedure :: Bcast => Bcast_INS_TimeIntegrator_Euler_3D
  end type INS_TimeIntegrator_Euler_Options_3D

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type INS_TimeIntegrator_Euler_3D

  function New_INS_TimeIntegrator_Euler_3D(problem, ins_op, opt) result(this)
    class(INS_Problem_3D),  intent(in) :: problem
    class(INS_Operator_3D), intent(in) :: ins_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt
    type(INS_TimeIntegrator_Euler_3D) :: this

    call Init_INS_TimeIntegrator_Euler_3D(this, problem, ins_op, opt)

  end function New_INS_TimeIntegrator_Euler_3D

  !-----------------------------------------------------------------------------
  !> Initialization of a INS_TimeIntegrator_Euler_3D object

  subroutine Init_INS_TimeIntegrator_Euler_3D(this, problem, ins_op, opt)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    class(INS_Problem_3D),              intent(in)    :: problem
    class(INS_Operator_3D),             intent(in)    :: ins_op
    class(INS_TimeIntegrator_Euler_Options_3D), optional, intent(in) :: opt

    ! intialize parent type
    call this % Init_INS_TimeIntegrator_3D(problem, ins_op, opt)

    this % name    = 'Euler method'
    this % i_max_p = opt % i_max_p
    this % i_max_v = opt % i_max_v
    this % r_red   = opt % r_red
    this % r_max   = opt % r_max

  end subroutine Init_INS_TimeIntegrator_Euler_3D

  !-----------------------------------------------------------------------------
  !> Execution of an Euler time step

  subroutine TimeStep(this, t, dt, u, standby)
    class(INS_TimeIntegrator_Euler_3D), intent(inout) :: this
    real(RNP),             intent(inout) :: t            !< time t₀ → t
    real(RNP),             intent(in)    :: dt           !< step size ∆t = t-t₀
    real(RNP), contiguous, intent(inout) :: u(:,:,:,:,:) !< u(x,t₀) → u(x,t)
    logical,     optional, intent(in)    :: standby      !< reuse workspace [F]

    ! internal variables .......................................................

    real(RNP), allocatable, save :: inv_mm(:,:,:,:) ! inverse diagonal mass matrix
    real(RNP), allocatable, save :: v_0(:,:,:,:,:)  ! initial velocity
    real(RNP), allocatable, save :: F_c(:,:,:,:,:)  ! convection term
    real(RNP), allocatable, save :: F_d(:,:,:,:,:)  ! viscous diffusion term
    real(RNP), allocatable, save :: Q  (:,:,:,:,:)  ! source term
    real(RNP), allocatable, save :: vp (:,:,:,:,:)  ! outer velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:)  ! outer viscous flux traces s⁺

    ! boundary points and values
    type(BoundaryVariable_3D), allocatable, save :: bv_x(:), bv_u(:), bv_v(:)

    ! control
    real(RNP), save :: t_0 = -huge(ONE)
    logical,   save :: fresh

    ! auxiliary
    integer :: b, e, d, np, po

    associate( problem => this % problem        &
             , ins_op  => this % ins_op         &
             , mesh    => this % ins_op % mesh  &
             , sem_v   => this % ins_op % sem_v &
             , v       => u(:,:,:,:,1:3)        )

      ! initialization .........................................................

      po = ins_op % eop_v % po
      np = po + 1

      !$omp master

      if (allocated(v_0)) then
        if (any(shape(v_0) /= shape(v))) then
          deallocate(inv_mm, v_0, F_c, F_d, Q, vp, sp)
          deallocate(bv_x, bv_u, bv_v)
        end if
      end if

      if (allocated(v_0)) then

        fresh = abs(t_0 - t) > epsilon(ONE)

      else

        fresh = .true.

        allocate( inv_mm (np, np, np, mesh % n_elem   ), source = ZERO )
        allocate( v_0    (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( F_c    (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( F_d    (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( Q      (np, np, np, mesh % n_elem, 3), source = ZERO )
        allocate( vp     (np, np,  6, mesh % n_elem, 3), source = ZERO )
        allocate( sp     (np, np,  6, mesh % n_elem, 3), source = ZERO )

        allocate(bv_x (mesh % n_bound) )
        allocate(bv_u (mesh % n_bound) )
        allocate(bv_v (mesh % n_bound) )
        do b = 1, mesh % n_bound
          bv_u(b) = BoundaryVariable_3D(mesh % boundary(b), po, nc = 4)
          bv_v(b) = BoundaryVariable_3D(bv_u(b), first=1, last=3)
          bv_x(b) = BoundaryVariable_3D(mesh % boundary(b), po, nc = 3)
          call bv_x(b) % Extract(sem_v % metrics % x)
        end do

      end if

      t_0 = t
      t   = t + dt

      !$omp end master
      !$omp barrier

      ! inverse diagonal mass matrix
      call sem_v % Get_DG_DiagonalMassMatrix(inv_mm)
      !$omp workshare
      inv_mm = 1 / inv_mm
      !$omp end workshare nowait

      ! convective terms at time t₀ ............................................

      if (fresh) then
        do b = 1, mesh % n_bound
          call problem % GetBoundaryValues(b, bv_x(b) % val, t_0, bv_u(b) % val)
        end do
        call GetOuterTraces_3D(mesh, v, vp)
        call ins_op % SetVelocityBC(bv_v, vp, sp)
      end if

      call ins_op % GetConvectionTerm(v, vp, F_c)

      ! viscous diffusion term at time t .......................................
      ! so far ν is constant and boundaries are periodic or have Dirichlet BC

      do b = 1, mesh % n_bound
        call problem % GetBoundaryValues(b, bv_x(b) % val, t, bv_u(b) % val)
      end do
      call ins_op % SetVelocityBC(bv_v, vp, sp)
      call ins_op % GetDiffusionTerm(v, vp, sp, F_d)

      ! viscous and convective RHS .............................................

      !$omp do
      do e = 1, mesh % n_elem
        do d = 1, 3
          F_c(:,:,:,e,d) = inv_mm(:,:,:,e) * F_c(:,:,:,e,d)
          F_d(:,:,:,e,d) = inv_mm(:,:,:,e) * F_d(:,:,:,e,d)
        end do
      end do

      ! sources at time t ......................................................

      call problem % GetExternalSources(sem_v % metrics % x, t, Q)

      ! initial velocity .......................................................

      call SetArray(v_0, v, multi = .true.)

      ! extrapolation-projection-diffusion step ................................

      call this % ProjectionStep( dt, v_0, F_c, F_d, Q, bv_u, u  &
                                , this % i_max_p, this % i_max_v &
                                , this % r_red  , this % r_max   )

      ! cleanup ................................................................

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(inv_mm, v_0, F_c, F_d, Q, vp, sp)
      deallocate(bv_x, bv_u, bv_v)
      !$omp end master

    end associate

  end subroutine TimeStep

  !=============================================================================
  ! TBP of INS_TimeIntegrator_Euler_Options_3D


  !-----------------------------------------------------------------------------
  !> MPI broadcasting of Euler time-integrator options

  subroutine Bcast_INS_TimeIntegrator_Euler_3D(this, root, comm)
    class(INS_TimeIntegrator_Euler_Options_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    type(MPI_Request) :: request(4)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    call this % INS_TimeIntegratorOptions_3D % Bcast(root, comm)

    n = 1
    call XMPI_Ibcast( this % i_max_p, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % i_max_v, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % r_red  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % r_max  , root, comm, request(n) )

    call MPI_Waitall( n, request, stat )

  end subroutine Bcast_INS_TimeIntegrator_Euler_3D

  !=============================================================================

end module INS__Time_Integrator__Euler__3D
