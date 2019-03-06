!> summary:  Polynomial multigrid for use with elliptic solvers
!> author:   Joerg Stiller
!> date:     2019/02/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid for use with elliptic solvers
!===============================================================================

module CART__Elliptic_PMG
  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Execution_Control, only: Error
  use XMPI

  use IP_Element_Operators_1D
  use Array_Assignments
  use Array_Reductions

  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Schwarz_Operator
  use CART__Elliptic_PMG_Level

  implicit none
  private

  public :: PMG_Method3D
  public :: PMG_Options3D

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid operators and related procedures

  type PMG_Method3D

    ! mesh, levels and operators
    type(MeshPartition), pointer     :: mesh      !< mesh partition
    type(PMG_Level),     allocatable :: level(:)  !< levels

    ! p-MG/CG solver settings
    integer   :: i_max   !< max number iterations (cycles)
    real(RNP) :: r_red   !< min residual reduction
    real(RNP) :: r_max   !< max admissible residual
    real(RNP) :: dr_min  !< termination threshold for Δr

    ! coarse grid solver settings
    character :: solver  !< coarse grid solver, 'C': CG, 'S': Schwarz
    integer   :: i0_max  !< max num of iterations on coarse grid
    real(RNP) :: r0_red  !< min residual reduction on coarse grid

    ! control
    logical   :: monitor !< switch for monitoring

  contains

    generic :: Init_PMG_Method3D => Init_IP
    procedure, private :: Init_IP

    generic :: SetProblem => SetProblem_CI, SetProblem_VI
    procedure, private :: SetProblem_CI
    procedure, private :: SetProblem_VI

    procedure :: GetWorkspace
    procedure :: FreeWorkspace
    procedure :: BCtoRHS

    procedure :: MG_Solver
    procedure :: MG_CG_Solver

  end type PMG_Method3D

  ! constructor interface
  interface PMG_Method3D
    module procedure New_IP
  end interface

  !-----------------------------------------------------------------------------
  !> Type bundling polynomial multigrid options

  type PMG_Options3D

    ! levels
    integer   :: po_top    = -1  !< polynomial order at top level
    integer   :: po_bot    =  1  !< polynomial order at bottom level
    real(RNP) :: cr        =  2  !< target coarsening ratio
    real(RNP) :: cr_max    =  2  !< max coarsening ratio

    ! p-MG/CG solver settings
    integer   :: i_max     =  1  !< max number iterations (cycles)
    real(RNP) :: r_red     = -1  !< min residual reduction
    real(RNP) :: r_max     = -1  !< max admissible residual
    real(RNP) :: dr_min    = -1  !< termination threshold for Δr

    ! smoothing settings
    integer   :: ns1       =  1  !< num pre-smoothing  steps on top level
    integer   :: ns2       =  1  !< num post-smoothing steps on top level
    integer   :: mvs       =  1  !< multiplier for variable smoothing

    type(SchwarzOptions3D) :: schwarz !< Schwarz method

    ! coarse grid solver settings
    character :: solver    = 'C' !< coarse grid solver, 'C': CG, 'S': Schwarz
    integer   :: i0_max    =  1  !< max number coarse grid iterations
    real(RNP) :: r0_red    = -1  !< min coarse grid residual reduction

    ! control
    logical   :: monitor   = .false. !< switch for monitoring

  contains

    procedure :: Bcast => PMG_Options3D_Bcast

  end type PMG_Options3D

  !=============================================================================
  ! Separate procedures

  interface

    !---------------------------------------------------------------------------
    !> PMG_Method3D initialization: IP/DG-SEM, no problem data

    module subroutine Init_IP( this, mesh, ip_opt, pmg_opt )
      class(PMG_Method3D),          intent(inout) :: this
      class(MeshPartition), target, intent(in)    :: mesh    !< mesh partition
      class(IP_ElementOptions1D),   intent(in)    :: ip_opt  !< IP/DG options
      class(PMG_Options3D),         intent(in)    :: pmg_opt !< PMG options
    end subroutine Init_IP

    !---------------------------------------------------------------------------
    !> Initialize problem with constant isotropic diffusivity

    module subroutine SetProblem_CI(this, lambda, nu, bc)
      class(PMG_Method3D), intent(inout) :: this
      real(RNP),           intent(in)    :: lambda  !< Helmholtz parameter
      real(RNP),           intent(in)    :: nu      !< diffusivity
      character,           intent(in)    :: bc(:)   !< boundary conditions
    end subroutine SetProblem_CI

    !---------------------------------------------------------------------------
    !> Initialize problem with variable isotropic diffusivity

    module subroutine SetProblem_VI(this, lambda, nu, bc)
      class(PMG_Method3D), intent(inout) :: this
      real(RNP),           intent(in)    :: lambda         !< Helmholtz parameter
      real(RNP),           intent(in)    :: nu(0:,0:,0:,:) !< diffusivity
      character,           intent(in)    :: bc(:)          !< boundary conditions
    end subroutine SetProblem_VI

  end interface

contains

!===============================================================================
! PMG_Method3D: constructors

!-------------------------------------------------------------------------------
!> PMG_Method3D constructor: IP/DG-SEM

type(PMG_Method3D) function New_IP(mesh, ip_opt, pmg_opt) result(this)
  class(MeshPartition), target, intent(in) :: mesh    !< mesh partition
  class(IP_ElementOptions1D),   intent(in) :: ip_opt  !< IP/DG options
  class(PMG_Options3D),         intent(in) :: pmg_opt !< PMG options

  call Init_IP(this, mesh, ip_opt, pmg_opt)

end function New_IP

!===============================================================================
! PMG_Method3D: type-bound procedures

!-------------------------------------------------------------------------------
!> Allocate workspace

subroutine GetWorkspace(this)
  class(PMG_Method3D), intent(inout) :: this
  integer :: l

  do l = 0, ubound(this % level,1)
    call this % level(l) % GetWorkspace()
  end do

end subroutine GetWorkspace

!-------------------------------------------------------------------------------
!> Delete workspace

subroutine FreeWorkspace(this)
  class(PMG_Method3D), intent(inout) :: this
  integer :: l

  do l = 0, ubound(this % level,1)
    call this % level(l) % FreeWorkspace()
  end do

end subroutine FreeWorkspace

!-------------------------------------------------------------------------------
!>  Adds the boundary contributions of the right hand side

subroutine BCtoRHS(this, bv, f)
  class(PMG_Method3D),    intent(in)    :: this
  type(BoundaryVariable), intent(in)    :: bv(:)      !< boundary values
  real(RNP),              intent(inout) :: f(:,:,:,:) !< RHS

  integer :: l_top

  l_top = ubound(this % level, 1)
  call this % level(l_top) % elliptic_op % BcToRHS(bv, f)

end subroutine BCtoRHS

!-------------------------------------------------------------------------------
!> p-MG solver

subroutine MG_Solver(this, u, f, ni, r_2)
  class(PMG_Method3D), intent(inout) :: this
  real(RNP),           intent(inout) :: u(:,:,:,:) !< approx/final solution
  real(RNP),           intent(in)    :: f(:,:,:,:) !< RHS
  integer,   optional, intent(out)   :: ni         !< number of executed cycles
  real(RNP), optional, intent(out)   :: r_2        !< L2 norm of residual

  ! local data .................................................................

  logical, save :: converged

  logical   :: check_convergence
  real(RNP) :: dr_min, rr, r_max, r_new, r_old
  integer   :: i, l_top

  ! prerequisites ..............................................................

  l_top = ubound(this%level, 1)

  check_convergence = max(this%r_red, this%r_max, this%dr_min) > 0

  associate( elliptic_op => this % level( l_top ) % elliptic_op, &
             u_top       => this % level( l_top ) % u,           &
             f_top       => this % level( l_top ) % f,           &
             r_top       => this % level( l_top ) % v            )

    ! initialization ...........................................................

    call SetArray(u_top, u)
    call SetArray(f_top, f)

    ! termination conditions
    if (check_convergence) then
      call elliptic_op % Residual(u_top, f_top, r_top)
      rr = ScalarProduct(r_top, r_top, this%mesh%comm)
      r_old  = sqrt(rr)
      r_max  = max(r_old * this%r_red, this%r_max)
      dr_min = this%dr_min
    end if

    !$omp single
    converged = .false.
    !$omp end single

    ! MG cycles ................................................................

    do i = 1, this % i_max

      call V_Cycle(this)

      if (check_convergence) then

        call elliptic_op % Residual(u_top, f_top, r_top)
        rr = ScalarProduct(r_top, r_top, this%mesh%comm)
        r_new = sqrt(rr)

        converged = r_new <= r_max .or. abs(r_new - r_old) <= dr_min

        !$omp master
        call XMPI_Bcast(converged, root=0, comm=this%mesh%comm)
        !$omp end master
        !$omp barrier

        if (converged) exit
        r_old = r_new

      end if

    end do

    ! finalization .............................................................

    call SetArray(u, u_top)

    if (present(ni)) then
      !$omp master
      ni = min(i, this%i_max)
      !$omp end master
    end if

    if (present(r_2)) then
      if (r_max <= 0) then
        call elliptic_op % Residual(u_top, f_top, r_top)
        rr = ScalarProduct(r_top, r_top, this%mesh%comm)
      end if
      !$omp master
      r_2 = sqrt(rr)
      !$omp end master
    end if

  end associate

end subroutine MG_Solver

!-------------------------------------------------------------------------------
!> p-MG/CG solver

subroutine MG_CG_Solver(this, u, f, ni, r_2, i_max)
  class(PMG_Method3D), intent(inout) :: this
  real(RNP),           intent(inout) :: u(:,:,:,:) !< approx/final solution
  real(RNP),           intent(in)    :: f(:,:,:,:) !< RHS
  integer,   optional, intent(out)   :: ni         !< number of executed cycles
  real(RNP), optional, intent(out)   :: r_2        !< L2 norm of residual
  integer,   optional, intent(in)    :: i_max      !< overrides preset num cycles

  ! local data .................................................................

  real(RNP), dimension(:,:,:,:), allocatable, save :: p, q, r, s, z
  logical,   save :: converged

  logical   :: check_convergence, singular
  real(RNP) :: dr_min, r_max, r_new, r_old, rr
  real(RNP) :: alpha, beta, delta
  integer   :: i, i_max_, l_top

  ! prerequisites ..............................................................

  l_top = ubound(this%level,1)

  check_convergence = max(this%r_red, this%r_max, this%dr_min) > 0

  ! workspace
  !$omp single
  allocate(p, mold=u)
  allocate(q, mold=u)
  allocate(r, mold=u)
  allocate(s, mold=u)
  allocate(z, mold=u)
  !$omp end single

  associate( elliptic_op => this % level( l_top ) % elliptic_op, &
             u_top       => this % level( l_top ) % u,           &
             f_top       => this % level( l_top ) % f,           &
             mesh        => this % mesh                          )

    ! initialization ...........................................................

    singular = elliptic_op % lambda == 0 .and. all(elliptic_op % bc == 'D')

    ! initial residual
    call elliptic_op % Residual(u, f, r)

    ! termination conditions
    if (check_convergence) then
      rr = ScalarProduct(r, r, mesh%comm)
      r_old  = sqrt(rr)
      r_max  = max(r_old * this%r_red, this%r_max)
      dr_min = this%dr_min
      !$omp master
      converged = r_old < r_max
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier
    else
      !$omp single
      converged = .false.
      !$omp end single
    end if

    if (converged) then
      i_max_ = 0
      i      = 0
      r_new  = r_old
    else if (present(i_max)) then
      i_max_ = i_max
    else
      i_max_ = this % i_max
    end if

    ! iteration ................................................................

    do i = 1, i_max_

      ! MG preconditioner: z = MG(r, 0)
      call SetArray(u_top, ZERO)                          ! u_L = 0
      call SetArray(f_top, r)                             ! f_L = r
      call V_Cycle(this)                                  ! u_L = MG(r, 0)
      call SetArray(z, u_top)                             ! z = u_L

      ! set/update search vector
      if (i == 1) then
        if (singular) then
          call CalibrateArray(z, mesh%comm)
        end if
        call SetArray(p, z)                               ! p = z
      else
        call SetArray(q, r)                               ! q = r
        call MergeArrays(ONE, q, -ONE, s)                 ! q = r - s
        beta = ScalarProduct(q, z, mesh%comm) / delta
        call MergeArrays(beta, p, ONE, z)                 ! p = beta p + z
      end if

      ! save old residual
      call SetArray(s, r)

      ! correction
      call elliptic_op % Apply(p, q)
      delta = ScalarProduct(r, z, mesh%comm)
      alpha = delta / ScalarProduct(p, q, mesh%comm)
      call MergeArrays(ONE, u,  alpha, p)                 ! u = u + alpha p
      call MergeArrays(ONE, r, -alpha, q)                 ! r = r - alpha q

      if (check_convergence) then

        r_new = sqrt( ScalarProduct(r, r, mesh%comm) )

        converged = r_new <= r_max .or. abs(r_new - r_old) <= dr_min

        !$omp master
        call XMPI_Bcast(converged, root=0, comm=mesh%comm)
        !$omp end master
        !$omp barrier

        r_old = r_new
      end if

      if (converged .or. i == this%i_max) exit

    end do

    ! finalization .............................................................

    if (present(ni)) then
      !$omp master
      ni = i
      !$omp end master
    end if

    if (present(r_2)) then
      if (.not. check_convergence) then
        r_new = sqrt( ScalarProduct(r, r, mesh%comm) )
      end if
      !$omp master
      r_2 = r_new
      !$omp end master
    end if

  end associate

  !$omp barrier
  !$omp master
  deallocate(p, q, r, s, z)
  !$omp end master

end subroutine MG_CG_Solver

!===============================================================================
! PMG_Method3D: helpers

!-------------------------------------------------------------------------------
!> PMG correction scheme V-cycle

subroutine V_Cycle(this)
  class(PMG_Method3D), intent(inout) :: this

  integer :: l, l_top

  associate(level => this%level)

    ! prerequisites ............................................................

    l_top = ubound(level,1)

    ! downward leg .............................................................

    if (this%monitor) call Monitoring(this, l_top, '0')

    do l = l_top, 1, -1

      ! initialize correction: u_0 = 0
      if (l < l_top) then
        call SetArray(level(l)%u, ZERO)
      end if

      ! pre-smoothing
      call level(l) % elliptic_op % SchwarzMethod( level(l)%u,  &
                                                   level(l)%f,  &
                                                   level(l)%ns1 )
      if (this%monitor) call Monitoring(this, l, '1')

      ! residual evaluation: v_l = f_l - A_l u_l
      call level(l) % elliptic_op % Residual(level(l)%u, level(l)%f, level(l)%v)

      ! restriction: f_l-1 = R v_l
      call level(l) % Restrict(level(l)%v, level(l-1)%f)

    end do

    ! coarse mesh solution .....................................................

    associate( u_0 => level(0)%u, i_max => this%i0_max, &
               f_0 => level(0)%f, r_red => this%r0_red  )

      call SetArray(u_0, ZERO)
      if (this%monitor) call Monitoring(this, 0, '0')
      select case(this % solver)
      case('C')
        call level(0) % elliptic_op % ConjugateGradients(u_0, f_0, i_max, r_red)
      case('S')
        call level(0) % elliptic_op % SchwarzMethod(u_0, f_0, i_max, r_red)
      end select
      if (this%monitor) call Monitoring(this, 0, 's')

    end associate

    ! upward leg ...............................................................

    do l = 1, l_top

      ! prolongation: v_l = I u_l-1
      call level(l-1) % Prolongate(level(l-1)%u, level(l)%v)

      ! correction: u_l = u_l + v_l
      call MergeArrays(ONE, level(l)%u, ONE, level(l)%v)
      if (this%monitor) call Monitoring(this, l, 'c')

      ! post-smoothing
      call level(l) % elliptic_op % SchwarzMethod( level(l)%u,  &
                                                   level(l)%f,  &
                                                   level(l)%ns2 )
      if (this%monitor) call Monitoring(this, l, '2')

    end do

  end associate

end subroutine V_Cycle

!-------------------------------------------------------------------------------
!> Monitoring of residual

subroutine Monitoring(this, l, step)
  type(PMG_Method3D), intent(inout) :: this
  integer,            intent(in)    :: l    !< level
  character(len=*),   intent(in)    :: step !< current step

  real(RNP) :: rr

  associate(u => this%level(l)%u, f => this%level(l)%f, r => this%level(l)%v)
    call this % level(l) % elliptic_op % Residual(u, f, r)
    rr = ScalarProduct(r, r, this%mesh%comm)
  end associate

  !$omp barrier
  !$omp master
  if (this%mesh%part == 0) then
    print '(2X,A,I2,3A,ES10.3)', 'level ', l, ': r[', step, '] =', sqrt(rr)
  end if
  !$omp end master

end subroutine Monitoring

!===============================================================================
! PMG_Options3D: type-bound procedures

subroutine PMG_Options3D_Bcast(this, root, comm)
  class(PMG_Options3D), intent(inout) :: this
  integer,        intent(in) :: root !< rank of broadcast root
  type(MPI_Comm), intent(in) :: comm !< MPI communicator

  type(MPI_Request) :: request(15)
  type(MPI_Status)  :: stat(size(request))
  integer :: n

  n = 1
  call XMPI_Ibcast( this % po_top    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % po_bot    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % cr        , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % cr_max    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % i_max     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r_red     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r_max     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % dr_min    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % ns1       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % ns2       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % mvs       , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % solver    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % i0_max    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % r0_red    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % monitor   , root, comm, request(n) )

  call MPI_Waitall(n, request, stat)

  call this % Schwarz % Bcast(root, comm)

end subroutine PMG_Options3D_Bcast

!===============================================================================

end module CART__Elliptic_PMG
