!> summary:  Euler method for incompressible flows
!> author:   Joerg Stiller
!> date:     2022/09/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Integrator__Euler__3D
  use Kind_Parameters
  use Constants
  use Array_Assignments
  use XMPI

  use Trace_Operators__3D
  use Element_Face_Transfer_Buffer__3D
  use SEM__Boundary_Variable__3D

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

    real(RNP), allocatable, save :: inv_M(:,:,:,:) ! inverse diagonal mass matrix
    real(RNP), allocatable, save :: v_0(:,:,:,:,:) ! initial velocity
    real(RNP), allocatable, save :: F_c(:,:,:,:,:) ! convection term
    real(RNP), allocatable, save :: F_d(:,:,:,:,:) ! viscous diffusion term
    real(RNP), allocatable, save :: Q  (:,:,:,:,:) ! source term
    real(RNP), allocatable, save :: vp (:,:,:,:,:) ! outer velocity traces v⁺
    real(RNP), allocatable, save :: sp (:,:,:,:,:) ! outer viscous flux traces s⁺

    ! boundary points and values
    type(SEM_BoundaryVariable_3D), allocatable, save :: bv_x(:), bv_u(:), bv_v(:)

    integer :: b, e, d, np
!### CHECK
print *, '# 00'
!### CHECK END

    associate( problem => this % problem        &
             , ins_op  => this % ins_op         &
             , mesh    => this % ins_op % mesh  &
             , sem_v   => this % ins_op % sem_v &
             , v       => u(:,:,:,:,1:3)        )

      ! initialization .........................................................

      np = ins_op % eop_v % po + 1  ! = size(u,1)

      !$omp master
      allocate( inv_M (np, np, np, mesh % n_elem   ), source = ZERO )
      allocate( v_0   (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( F_c   (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( F_d   (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( Q     (np, np, np, mesh % n_elem, 3), source = ZERO )
      allocate( vp    (np, np,  6, mesh % n_elem, 3), source = ZERO )
      allocate( sp    (np, np,  6, mesh % n_elem, 3), source = ZERO )
      bv_x = SEM_BoundaryVariable_3D(sem_v, mesh % boundary, nc = 3)
      bv_u = SEM_BoundaryVariable_3D(sem_v, mesh % boundary, nc = 4)
      bv_v = SEM_BoundaryVariable_3D(bv_u, first=1, last=3)
      !$omp end master
      !$omp barrier
!### CHECK
print *, '# 01'
print *, '# 01, shape(this % ins_op % mesh % boundary)   =', shape(this % ins_op % mesh % boundary)
print *, '# 01, size(bv_x)                               =',  size(bv_x)
print *, '# 01, size(bv_u)                               =',  size(bv_u)
print *, '# 01, size(bv_v)                               =',  size(bv_v)
do b = 1, mesh % n_bound
print *, '# 01, b ========================================', b
print *, '# 01, associated(bv_x(b) % val)                =', associated(bv_x(b) % val)
print *, '# 01, associated(bv_x(b) % sem)                =', associated(bv_x(b) % sem)
print *, '# 01, associated(bv_x(b) % sem % mesh)         =', associated(bv_x(b) % sem % mesh)
print *, '# 01, shape(bv_x(b) % sem % mesh % boundary)   =', shape(bv_x(b) % sem % mesh % boundary)
print *, '# 01, associated(bv_u(b) % val)                =', associated(bv_u(b) % val)
print *, '# 01, associated(bv_u(b) % sem)                =', associated(bv_u(b) % sem)
print *, '# 01, associated(bv_u(b) % sem % mesh)         =', associated(bv_u(b) % sem % mesh)
print *, '# 01, shape(bv_u(b) % sem % mesh % boundary)   =', shape(bv_u(b) % sem % mesh % boundary)
print *, '# 01, associated(bv_v(b) % val)                =', associated(bv_v(b) % val)
print *, '# 01, associated(bv_v(b) % sem)                =', associated(bv_v(b) % sem)
print *, '# 01, associated(bv_v(b) % sem % mesh)         =', associated(bv_v(b) % sem % mesh)
print *, '# 01, shape(bv_v(b) % sem % mesh % boundary)   =', shape(bv_v(b) % sem % mesh % boundary)
end do
!### CHECK END

      ! inverse diagonal mass matrix
      call sem_v % Get_DG_DiagonalMassMatrix(inv_M)
      !$omp workshare
      inv_M = 1 / inv_M
      !$omp end workshare nowait
!### CHECK
print *, '# 02'
!### CHECK END

      ! initial velocity
      call SetArray(v_0, v, multi = .true.)
!### CHECK
print *, '# 03'
!### CHECK END

      ! final time
      t = t + dt

      ! sources ................................................................
      ! think about reusing

      call problem % GetExternalSources(sem_v % metrics % x, t, Q)
!### CHECK
print *, '# 04'
do b = 1, mesh % n_bound
print *, '# 04, b ========================================', b
print *, '# 04, associated(bv_x(b) % val)                =', associated(bv_x(b) % val)
print *, '# 04, associated(bv_x(b) % sem)                =', associated(bv_x(b) % sem)
print *, '# 04, associated(bv_x(b) % sem % mesh)         =', associated(bv_x(b) % sem % mesh)
print *, '# 04, shape(bv_x(b) % sem % mesh % boundary)   =', shape(bv_x(b) % sem % mesh % boundary)
print *, '# 04, associated(bv_u(b) % val)                =', associated(bv_u(b) % val)
print *, '# 04, associated(bv_u(b) % sem)                =', associated(bv_u(b) % sem)
print *, '# 04, associated(bv_u(b) % sem % mesh)         =', associated(bv_u(b) % sem % mesh)
print *, '# 04, shape(bv_u(b) % sem % mesh % boundary)   =', shape(bv_u(b) % sem % mesh % boundary)
print *, '# 04, associated(bv_v(b) % val)                =', associated(bv_v(b) % val)
print *, '# 04, associated(bv_v(b) % sem)                =', associated(bv_v(b) % sem)
print *, '# 04, associated(bv_v(b) % sem % mesh)         =', associated(bv_v(b) % sem % mesh)
print *, '# 04, shape(bv_v(b) % sem % mesh % boundary)   =', shape(bv_v(b) % sem % mesh % boundary)
end do
!### CHECK END

      ! boundary values ........................................................
      ! also think about reusing

      do b = 1, mesh % n_bound
        call bv_x(b) % Extract(sem_v % metrics % x, mesh % boundary(b))
!### CHECK
print *, '# 04+, b ========================================', b
print *, '# 04+, associated(bv_x(b) % val)                =', associated(bv_x(b) % val)
print *, '# 04+, associated(bv_x(b) % sem)                =', associated(bv_x(b) % sem)
print *, '# 04+, associated(bv_x(b) % sem % mesh)         =', associated(bv_x(b) % sem % mesh)
print *, '# 04+, shape(bv_x(b) % sem % mesh % boundary)   =', shape(bv_x(b) % sem % mesh % boundary)
print *, '# 04+, associated(bv_u(b) % val)                =', associated(bv_u(b) % val)
print *, '# 04+, associated(bv_u(b) % sem)                =', associated(bv_u(b) % sem)
print *, '# 04+, associated(bv_u(b) % sem % mesh)         =', associated(bv_u(b) % sem % mesh)
print *, '# 04+, shape(bv_u(b) % sem % mesh % boundary)   =', shape(bv_u(b) % sem % mesh % boundary)
!### CHECK END
        call problem % GetBoundaryValues(b, bv_x(b) % val, t, bv_u(b) % val)
      end do
!### CHECK
print *, '# 05'
print *, '# 05, shape(this % ins_op % mesh % boundary)   =', shape(this % ins_op % mesh % boundary)
do b = 1, size(bv_u)
print *, '# 05, b ==========================================', b
print *, '# 05, associated(bv_u(b) % sem)                  =', associated(bv_u(b) % sem)
print *, '# 05, associated(bv_u(b) % sem % mesh)           =', associated(bv_u(b) % sem % mesh)
print *, '# 05, shape(bv_u(b) % sem % mesh % boundary)     =', shape(bv_u(b) % sem % mesh % boundary)
print *, '# 05, allocated(bv_u(b) % sem % mesh % boundary) =', allocated(bv_u(b) % sem % mesh % boundary)
print *, '# 05, associated(bv_v(b) % sem)                  =', associated(bv_v(b) % sem)
print *, '# 05, associated(bv_v(b) % sem % mesh)           =', associated(bv_v(b) % sem % mesh)
print *, '# 05, shape(bv_v(b) % sem % mesh % boundary)     =', shape(bv_v(b) % sem % mesh % boundary)
print *, '# 05, allocated(bv_v(b) % sem % mesh % boundary) =', allocated(bv_v(b) % sem % mesh % boundary)
end do
!### CHECK END


      ! viscous and convective RHS .............................................
      ! so far ν is constant and boundaries are periodic or have Dirichlet BC

      call ins_op % SetVelocityBC(bv_v, vp, sp)
!### CHECK
print *, '# 06'
b = 1
print *, '# 06, shape(this % ins_op % mesh % boundary)   =', shape(this % ins_op % mesh % boundary)
print *, '# 06, associated(bv_u(b) % sem)                =', associated(bv_u(b) % sem)
print *, '# 06, associated(bv_u(b) % sem % mesh)         =', associated(bv_u(b) % sem % mesh)
print *, '# 06, shape(bv_u(b) % sem % mesh % boundary)   =', shape(bv_u(b) % sem % mesh % boundary)
print *, '# 06, associated(bv_v(b) % sem)                =', associated(bv_v(b) % sem)
print *, '# 06, associated(bv_v(b) % sem % mesh)         =', associated(bv_v(b) % sem % mesh)
print *, '# 06, shape(bv_v(b) % sem % mesh % boundary)   =', shape(bv_v(b) % sem % mesh % boundary)
!### CHECK END
!###!      call ins_op % GetDiffusionTerm(v, vp, sp, F_d)
!### CHECK
print *, '# 07'
!### CHECK END
!###!      call ins_op % GetConvectionTerm(v, vp, F_c)
!### CHECK
print *, '# 08'
!### CHECK END

      !$omp do
      do e = 1, mesh % n_elem
        do d = 1, 3
          F_c(:,:,:,e,d) = inv_M(:,:,:,e) * F_c(:,:,:,e,d)
          F_d(:,:,:,e,d) = inv_M(:,:,:,e) * F_d(:,:,:,e,d)
        end do
      end do
!### CHECK
print *, '# 09'
b = 1
print *, '# 09, shape(this % ins_op % mesh % boundary)   =', shape(this % ins_op % mesh % boundary)
print *, '# 09, associated(bv_u(b) % sem)                =', associated(bv_u(b) % sem)
print *, '# 09, associated(bv_u(b) % sem % mesh)         =', associated(bv_u(b) % sem % mesh)
print *, '# 09, shape(bv_u(b) % sem % mesh % boundary)   =', shape(bv_u(b) % sem % mesh % boundary)
print *, '# 09, associated(bv_v(b) % sem)                =', associated(bv_v(b) % sem)
print *, '# 09, associated(bv_v(b) % sem % mesh)         =', associated(bv_v(b) % sem % mesh)
print *, '# 09, shape(bv_v(b) % sem % mesh % boundary)   =', shape(bv_v(b) % sem % mesh % boundary)
!### CHECK END

      ! extrapolation-projection-diffusion step ................................

      call this % ProjectionStep( dt, v_0, F_c, F_d, Q, bv_u, u  &
                                , this % i_max_p, this % i_max_v &
                                , this % r_red  , this % r_max   )
!### CHECK
print *, '# 10'
!### CHECK END

      ! cleanup ................................................................

      if (present(standby)) then
        if (standby) return
      end if

      !$omp master
      deallocate(inv_M, v_0, F_c, F_d, Q, vp, sp)
      !$omp end master

    end associate
!### CHECK
print *, '# XX'
!### CHECK END

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
