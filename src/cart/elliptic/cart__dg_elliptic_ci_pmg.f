!> summary:  Polynomial multigrid for use with DG elliptic solvers
!> author:   Joerg Stiller
!> date:     2018/01/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid for use with DG elliptic solvers
!===============================================================================

module CART__DG_Elliptic_CI_PMG

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE
  use Array_Assignments
  use Array_Reductions

  use XMPI

  use CART__Mesh_Partition
  use CART__DG_PMG_Transfer_Operators
  use CART__DG_Elliptic_CI_PMG_Level
  use CART__DG_Elliptic_CI_Operator
  use CART__DG_Elliptic_CI_Residual
  use CART__DG_Elliptic_CI_Conj_Grad

  implicit none
  private

  public :: PolynomialMultigrid_Options
  public :: PolynomialMultigrid

  !-----------------------------------------------------------------------------
  !> Type bundling polynomial multigrid options

  type PolynomialMultigrid_Options

    ! levels
    integer   :: po_top    = -1        !< polynomial order at top level
    integer   :: po_bot    =  1        !< polynomial order at bottom4 level
    real(RNP) :: cr        =  2        !< target coarsening ratio
    real(RNP) :: cr_max    =  2        !< max coarsening ratio

    ! p-MG/CG solver settings
    integer   :: i_max     =  1        !< max number iterations (cycles)
    real(RNP) :: r_red     = -1        !< min residual reduction
    real(RNP) :: r_max     = -1        !< max admissible residual
    real(RNP) :: dr_min    = -1        !< termination threshold for Δr

    ! smoothing settings
    integer   :: ns1       =  1        !< num pre-smoothing  steps on top level
    integer   :: ns2       =  1        !< num post-smoothing steps on top level
    integer   :: mvs       =  1        !< multiplier for variable smoothing

    ! coarse grid solver settings
    character(len=20) :: &
                 solver    = 'CG'      !< coarse grid solver {'CG','Schwarz}
    integer   :: i0_max    =  1        !< max number coarse grid iterations
    real(RNP) :: r0_red    = -1        !< min coarse grid residual reduction

    ! Schwarz smoother/solver settings
    real(RNP) :: delta(3)  =  0.08     !< relative element overlap
    integer   :: no_min    =  1        !< min overlap in points >= 0
    integer   :: weighting =  5        !< weighting method {0,1,3,5,7,9}

    ! control
    logical   :: monitor   = .false.   !< switch for monitoring

  contains

    procedure :: Bcast => PMG_Options_Bcast

  end type PolynomialMultigrid_Options

  !-----------------------------------------------------------------------------
  !> Polynomial multigrid operators and related procedures

  type PolynomialMultigrid

    ! mesh
    type(MeshPartition),         pointer     :: mesh        !< mesh partition
    type(PMG_Level),             allocatable :: level(:)    !< levels
    type(PMG_TransferOperators), allocatable :: transfer(:) !< transfer ops

    ! p-MG/CG solver settings
    integer   :: i_max            !< max number iterations (cycles)
    real(RNP) :: r_red            !< min residual reduction
    real(RNP) :: r_max            !< max admissible residual
    real(RNP) :: dr_min           !< termination threshold for Δr

    ! coarse grid solver settings
    character(len=20) :: solver   !< coarse grid solver {'CG','Schwarz}
    integer   :: i0_max           !< max num of iterations on coarse grid
    real(RNP) :: r0_red           !< min residual reduction on coarse grid

    ! control
    logical   :: monitor          !< switch for monitoring

  contains

    procedure :: New => New_PolynomialMultigrid

    generic   :: SetProblem => SetProblem_C
    procedure, private :: SetProblem_C

    procedure :: MG_Solver
    procedure :: MG_CG_Solver
    final     :: Delete_PolynomialMultigrid

  end type PolynomialMultigrid

contains

!===============================================================================
! PolynomialMultigrid_Options: type-bound procedures

subroutine PMG_Options_Bcast(this, root, comm)
  class(PolynomialMultigrid_Options), intent(inout) :: this
  integer,        intent(in) :: root !< rank of broadcast root
  type(MPI_Comm), intent(in) :: comm !< MPI communicator

  type(MPI_Request) :: request(18)
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
  call XMPI_Ibcast( this % delta     , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % no_min    , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % weighting , root, comm, request(n) );  n = n + 1
  call XMPI_Ibcast( this % monitor   , root, comm, request(n) )

  call MPI_Waitall(n, request, stat)

end subroutine PMG_Options_Bcast

!===============================================================================
! PolynomialMultigrid: type-bound procedures

!-------------------------------------------------------------------------------
!> Initialization of a PolynomialMultigrid object with all levels given

subroutine New_PolynomialMultigrid(this, opt, mesh, penalty)

  ! arguments ..................................................................

  class(PolynomialMultigrid),         intent(inout) :: this
  class(PolynomialMultigrid_Options), intent(in)    :: opt     !< PMG options
  class(MeshPartition), target,       intent(in)    :: mesh    !< mesh partition
  real(RNP),                          intent(in)    :: penalty !< IP penalty

  ! local variables ............................................................

  integer, allocatable :: po(:)
  integer :: l, l_top, ns1, ns2

  ! prerequisites...............................................................

  call Delete_PolynomialMultigrid(this)

  call CreatePolynomialLevels(opt%po_top, opt%po_bot, opt%cr, opt%cr_max, po)
  l_top = ubound(po,1)

  ! initialization of components ...............................................

  this % mesh    => mesh

  this % i_max   = opt % i_max
  this % r_red   = opt % r_red
  this % r_max   = opt % r_max
  this % dr_min  = opt % dr_min

  this % solver  = opt % solver
  this % i0_max  = opt % i0_max
  this % r0_red  = opt % r0_red
  this % monitor = opt % monitor

  allocate(this % level(0:l_top))
  ns1 = opt % ns1
  ns2 = opt % ns2
  do l = l_top, 0, -1
    call this % level(l) % New( po(l), l == 0, l == l_top  &
                              , mesh % ne                  &
                              , mesh % dx                  &
                              , penalty                    &
                              , ns1, ns2                   &
                              , opt % delta                &
                              , opt % no_min               &
                              , opt % weighting            &
                              )
    ns1 = ns1 * opt % mvs
    ns2 = ns2 * opt % mvs
  end do

  allocate(this % transfer(1:l_top))
  do l = 1, l_top
    call this % transfer(l) % New(po(l-1), po(l))
  end do

end subroutine New_PolynomialMultigrid

!-------------------------------------------------------------------------------
!> Creates a set of ascending polynomial orders

subroutine CreatePolynomialLevels(po_top, po_bot, cr, cr_max, po)
  integer,   intent(in) :: po_top  !< polynomial order at top level
  integer,   intent(in) :: po_bot  !< polynomial order at bottom level
  real(RNP), intent(in) :: cr      !< target coarsening ratio
  real(RNP), intent(in) :: cr_max  !< max coarsening ratio
  integer, allocatable, intent(out) :: po(:)  !< polynomial levels

  integer, allocatable :: q(:)
  integer   :: lb, lt

  if (po_top > po_bot) then

    allocate(q(po_bot:po_top), source = po_top)
    lt = po_top
    lb = lt
    do
      lb = lb -1
      q(lb) = max( nint(q(lb+1)/cr), ceiling(q(lb+1)/cr_max), 1 )
      q(lb) = min( q(lb), q(lb+1) - 1 )
      if (q(lb) <= max(po_bot, 1)) then
        q(lb) = max(q(lb), po_bot)
        exit
      end if
    end do
    allocate(po(0:lt-lb), source=q(lb:lt))

  else if (po_top == po_bot) then

    allocate(po(0:0), source=po_top)

  else

    allocate(po(0:-1))

  end if

end subroutine CreatePolynomialLevels

!-------------------------------------------------------------------------------
!> (Re)initialize problem with constant diffusivity

subroutine SetProblem_C(this, lambda, nu, bc)
  class(PolynomialMultigrid), intent(inout) :: this
  real(RNP), intent(in) :: lambda    !< Helmholtz parameter
  real(RNP), intent(in) :: nu        !< constant diffusivity
  character, intent(in) :: bc(:)     !< boundary conditions {P,D,N}

  integer :: l

  do l = ubound(this%level,1), lbound(this%level,1), -1
    call this % level(l) % SetProblem(lambda, nu, bc)
  end do

end subroutine SetProblem_C

!-------------------------------------------------------------------------------
!> PMG correction scheme V-cycle

subroutine V_Cycle(this)
  class(PolynomialMultigrid), intent(inout) :: this

  integer :: l, l_top

  associate(level => this%level, transfer => this%transfer)

    ! prerequisites ............................................................

    l_top = ubound(level,1)

    ! downward leg .............................................................

    if (this%monitor) call Monitoring(this, l_top, '0')

    do l = l_top, 1, -1

      ! initialize correction: u_0 = 0
      if (l < l_top) then
        call AssignScalar(level(l)%u, ZERO)
      end if

      ! pre-smoothing
      call Smoother(this, l, level(l)%ns1)
      if (this%monitor) call Monitoring(this, l, '1')

      ! residual evaluation: v_l = f_l - A_l u_l
      call Residual(this, l)

      ! restriction: f_l-1 = R v_l
      call transfer(l) % Restrict(level(l)%v, level(l-1)%f)

    end do

    ! coarse mesh solution .....................................................

    call AssignScalar(level(0)%u, ZERO)
    if (this%monitor) call Monitoring(this, 0, '0')
    call CoarseGridSolver(this)
    if (this%monitor) call Monitoring(this, 0, 's')

    ! upward leg ...............................................................

    do l = 1, l_top

      ! prolongation: v_l = I u_l-1
      call transfer(l) % Prolongate(level(l-1)%u, level(l)%v)

      ! correction: u_l = u_l + v_l
      call MergeArrays(ONE, level(l)%u, ONE, level(l)%v)
      if (this%monitor) call Monitoring(this, l, 'c')

      ! post-smoothing
      call Smoother(this, l, level(l)%ns2)
      if (this%monitor) call Monitoring(this, l, '2')

    end do

  end associate

end subroutine V_Cycle

!-------------------------------------------------------------------------------
!> p-MG solver

subroutine MG_Solver( this, lambda, nu, bc, f, u, ni, r_2 )

  ! arguments ..................................................................

  class(PolynomialMultigrid), intent(inout) :: this

  real(RNP),           intent(in)    :: lambda     !< Helmholtz parameter
  real(RNP),           intent(in)    :: nu         !< diffusivity
  character,           intent(in)    :: bc(:)      !< boundary conditions {P,D,N}
  real(RNP),           intent(in)    :: f(:,:,:,:) !< RHS
  real(RNP),           intent(inout) :: u(:,:,:,:) !< approx/final solution
  integer,   optional, intent(out)   :: ni         !< number of executed cycles
  real(RNP), optional, intent(out)   :: r_2        !< L2 norm of residual

  ! local data .................................................................

  logical, save :: converged

  logical   :: check_convergence
  real(RNP) :: dr_min, rr, r_max, r_new, r_old
  integer   :: i, l_top

  ! prerequisites ..............................................................

  l_top = ubound(this%level,1)

  check_convergence = max(this%r_red, this%r_max, this%dr_min) > 0

  call this % SetProblem(lambda, nu, bc)

  associate( u_top => this % level( l_top ) % u  &
           , f_top => this % level( l_top ) % f  &
           , r_top => this % level( l_top ) % v  )

    ! initialization ...........................................................

    call AssignArray(u_top, u)
    call AssignArray(f_top, f)

    ! termination conditions
    if (check_convergence) then
      call Residual(this, l_top)
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

        call Residual(this, l_top)
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

    call AssignArray(u, u_top)

    if (present(ni)) then
      !$omp master
      ni = min(i, this%i_max)
      !$omp end master
    end if

    if (present(r_2)) then
      if (r_max <= 0) then
        call Residual(this, l_top)
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

subroutine MG_CG_Solver( this, lambda, nu, bc, f, u, ni, r_2, i_max )

  ! arguments ..................................................................

  class(PolynomialMultigrid), intent(inout) :: this

  real(RNP),           intent(in)    :: lambda     !< Helmholtz parameter
  real(RNP),           intent(in)    :: nu         !< diffusivity
  character,           intent(in)    :: bc(:)      !< boundary conditions {P,D,N}
  real(RNP),           intent(in)    :: f(:,:,:,:) !< RHS
  real(RNP),           intent(inout) :: u(:,:,:,:) !< approx/final solution
  integer,   optional, intent(out)   :: ni         !< number of executed cycles
  real(RNP), optional, intent(out)   :: r_2        !< L2 norm of residual
  integer,   optional, intent(in)    :: i_max      !< overrides preset num cycles

  ! local data .................................................................

  real(RNP), dimension(:,:,:,:), allocatable, save :: g, p, q, r, s, z
  logical,   save :: converged

  logical   :: check_convergence
  real(RNP) :: dr_min, r_max, r_new, r_old
  real(RNP) :: alpha, beta, delta
  integer   :: i, i_max_, l_top

  ! prerequisites ..............................................................

  l_top = ubound(this%level,1)

  check_convergence = max(this%r_red, this%r_max, this%dr_min) > 0

  call this % SetProblem(lambda, nu, bc)

  ! workspace
  !$omp single
  allocate(g, mold=u)
  allocate(p, mold=u)
  allocate(q, mold=u)
  allocate(r, mold=u)
  allocate(s, mold=u)
  allocate(z, mold=u)
  !$omp end single

  associate( mesh   => this % mesh                  &
           , eop    => this % level( l_top ) % eop  &
           , u_top  => this % level( l_top ) % u    &
           , f_top  => this % level( l_top ) % f    &
           )

    ! initialization ...........................................................

    ! calibrate RHS of singular problem
    call AssignArray(g, f)
    if (abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')) then
      call CalibrateArray(g, mesh%comm)
    end if

    ! initial residual
    call EllipticResidual(mesh, eop, lambda, nu, bc, u, g, r)

    ! termination conditions
    if (check_convergence) then
      call Residual(this, l_top)
      r_old  = sqrt( ScalarProduct(r, r, mesh%comm) )
      r_max  = max(r_old * this%r_red, this%r_max)
      dr_min = this%dr_min
      !$omp master
      converged = r_old < r_max
      call XMPI_Bcast(converged, root=0, comm=this%mesh%comm)
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
      call AssignScalar(u_top, ZERO)                      ! u_L = 0
      call AssignArray(f_top, r)                          ! f_L = r
      call V_Cycle(this)                                  ! u_L = MG(r, 0)
      call AssignArray(z, u_top)                          ! z = u_L

      ! set/update search vector
      if (i == 1) then
        if (abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')) then
          call CalibrateArray(z, mesh%comm)
        end if
        call AssignArray(p, z)                            ! p = z
      else
        call AssignArray(q, r)                            ! q = r
        call MergeArrays(ONE, q, -ONE, s)                 ! q = r - s
        beta = ScalarProduct(q, z, mesh%comm) / delta
        call MergeArrays(beta, p, ONE, z)                 ! p = beta p + z
      end if

      ! save old residual
      call AssignArray(s, r)

      ! correction
      call EllipticOperator(mesh, eop, lambda, nu, bc, p, q)
      delta = ScalarProduct(r, z, mesh%comm)
      alpha = delta / ScalarProduct(p, q, mesh%comm)
      call MergeArrays(ONE, u,  alpha, p)                 ! u = u + alpha p
      call MergeArrays(ONE, r, -alpha, q)                 ! r = r - alpha q

      if (check_convergence) then

        r_new = sqrt( ScalarProduct(r, r, mesh%comm) )

        converged = r_new <= r_max .or. abs(r_new - r_old) <= dr_min

        !$omp master
        call XMPI_Bcast(converged, root=0, comm=this%mesh%comm)
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
  deallocate(g, p, q, r, s, z)
  !$omp end master

end subroutine MG_CG_Solver

!-------------------------------------------------------------------------------
!> Finalization of a PolynomialMultigrid object

subroutine Delete_PolynomialMultigrid(this)
  type(PolynomialMultigrid), intent(inout) :: this

  nullify(this % mesh)

  if (allocated(this % level   )) deallocate(this % level   )
  if (allocated(this % transfer)) deallocate(this % transfer)

end subroutine Delete_PolynomialMultigrid

!===============================================================================
! Helper routines

!-------------------------------------------------------------------------------
!> Residual

subroutine Residual(this, l)
  type(PolynomialMultigrid), intent(inout) :: this
  integer,                   intent(in)    :: l !< level

  if (this % level(l) % has_var_nu) then

  else

    call EllipticResidual( this % mesh               &
                          , this % level(l) % eop     &
                          , this % level(l) % lambda  &
                          , this % level(l) % nu_c    &
                          , this % level(l) % bc      &
                          , this % level(l) % u       &
                          , this % level(l) % f       &
                          , this % level(l) % v       &
                          )
  end if

end subroutine Residual

!-------------------------------------------------------------------------------
!> Smoother

subroutine Smoother(this, l, ns)
  type(PolynomialMultigrid), intent(inout) :: this
  integer,                   intent(in)    :: l  !< level
  integer,                   intent(in)    :: ns !< number of smoothing steps

  if (this % level(l) % has_var_nu) then

  else

    call this % level(l) % schwarz % Iteration( this % mesh               &
                                              , this % level(l) % lambda  &
                                              , this % level(l) % nu_c    &
                                              , this % level(l) % bc      &
                                              , this % level(l) % u       &
                                              , this % level(l) % f       &
                                              , ns                        &
                                              )
  end if

end subroutine Smoother

!-------------------------------------------------------------------------------
!> Coarse grid solver

subroutine CoarseGridSolver(this)
  type(PolynomialMultigrid), intent(inout) :: this

  if (this % level(0) % has_var_nu) then
  else

    if (trim(this%solver) == 'Schwarz') then
      call this % level(0) % schwarz % Iteration( this % mesh               &
                                                , this % level(0) % lambda  &
                                                , this % level(0) % nu_c    &
                                                , this % level(0) % bc      &
                                                , this % level(0) % u       &
                                                , this % level(0) % f       &
                                                , this % i0_max             &
                                                , this % r0_red             &
                                                )
    else
      call ConjugateGradients( this % mesh               &
                             , this % level(0) % eop     &
                             , this % level(0) % lambda  &
                             , this % level(0) % nu_c    &
                             , this % level(0) % bc      &
                             , this % level(0) % u       &
                             , this % level(0) % f       &
                             , this % i0_max             &
                             , this % r0_red             &
                             )
    end if
  end if

end subroutine CoarseGridSolver

!-------------------------------------------------------------------------------
!> Monitoring of residual

subroutine Monitoring(this, l, step)
  type(PolynomialMultigrid), intent(inout) :: this
  integer,                   intent(in)    :: l    !< level
  character(len=*),          intent(in)    :: step !< current step

  real(RNP) :: rr

  call Residual(this, l)

  associate(r => this % level(l) % v)
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

end module CART__DG_Elliptic_CI_PMG
