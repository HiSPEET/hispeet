!> summary:  Validation of the tensor-product Schwarz operator
!> author:   Joerg Stiller
!> date:     2018/11/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__CART__TPO_Schwarz_Cubic
  use Kind_Parameters,   only: IXL, RNP
  use Constants,         only: ONE, THIRD, ZERO
  use Array_Assignments, only: AssignScalar
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use Execution_Control
  use IP_Element_Operators_1D

  use XMPI

  use CART__Mesh_Partition
  use CART__Generate_Structured_Mesh
  use CART__Schwarz_Operator
  use CART__TPO_Schwarz_Cubic
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer   :: po          =  8       ! polynomial order
  integer   :: np(3)       =  1       ! number of partitions in directions 1:3
  integer   :: ep(3)       =  10      ! elements per partition and direction
  real(RNP) :: delta(3)    = -1       ! relative overlap in directions 1:3
  integer   :: nt          =  1       ! number of test runs

  namelist /input/ po, ep, delta, nt

  ! MPI ........................................................................

  type(MPI_Comm) :: comm              ! communicator
  integer        :: comm_size         ! number of processes
  integer        :: rank              ! local rank

  ! mesh .......................................................................

  integer   :: ne
  logical   :: periodic(3) = .false.
  real(RNP) :: xo(3) = 0
  real(RNP) :: lx(3) = 1
  real(RNP) :: dx(3)

  type(MeshPartition) :: mesh

  ! operators and variables ....................................................

  real(RNP) :: lambda = 1                         ! Helmholtz parameter
  real(RNP) :: nu     = 2                         ! diffusivity
  character :: bc(6)  = ['D','D','D','N','N','N'] ! boundary conditions

  type(IP_ElementOperators1D) :: element_op ! IP/DG element oprators
  type(SchwarzOperator3D)     :: schwarz_op ! Schwarz operator

  procedure(TPO_Schwarz_Cubic_Proc), pointer :: SchwarzOp_Gen
  procedure(TPO_Schwarz_Cubic_Proc), pointer :: SchwarzOp_Par

  real(RNP), allocatable :: f(:,:,:,:), u(:,:,:,:), r(:,:,:,:)

  ! auxiliary ..................................................................

  character(len=80) :: input_file = 'validate__cart__tpo_schwarz_cubic.prm'

  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_par, mflops_par, mlups_par

  logical :: exists, parametrized
  integer :: nflop, nop, prm
  integer :: n1, n2, n3, nc, nd
  integer :: i

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! MPI ........................................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_size(comm, comm_size)
  call MPI_Comm_rank(comm, rank)

  if (comm_size > 1) then
    call Error('Validate__CART__TPO_Schwarz', &
               'Please run this program with one MPI process only')
  end if

  ! read test parameters .......................................................

  inquire(file=input_file, exist=exists)
  if (exists) then
    open(newunit=prm, file=input_file)
    read(prm, nml=input)
    close(prm)
  end if

  ! enforce uniform subdomains
  delta = delta(1)

  ! mesh .......................................................................

  ne = product(ep)
  dx = lx / (np * ep)

  call GenerateStructuredMesh(mesh, np, ep, xo, dx, periodic, comm)

  ! operators ..................................................................

  call element_op % New(po)
  call schwarz_op % New(element_op, mesh, lambda, nu, bc, delta)

  n1 = schwarz_op % n1
  n2 = schwarz_op % n2
  n3 = schwarz_op % n3
  nc = schwarz_op % nc
  nd = ne

  call TPO_Schwarz_Cubic_Assign(-1, SchwarzOp_Gen)  ! generic, for reference
  call TPO_Schwarz_Cubic_Assign(n1, SchwarzOp_Par)  ! parametrized

  parametrized = .not. associated( SchwarzOp_Par, &
                                   SchwarzOp_Gen  )

  ! workspace ..................................................................

  nop   = n1 * n2 * n3
  nflop = nop * (4 * (n1 + n2 + n3) + 1)

  allocate(f(n1,n2,n3,nd), u(n1,n2,n3,nd), r(n1,n2,n3,nd))
  call AssignScalar(u, ZERO)
  call AssignScalar(r, ZERO)
  call AssignScalar(f, ZERO)
  call random_number(f)

  !-----------------------------------------------------------------------------

  associate( S1  => schwarz_op % S1,  W1 => schwarz_op % W1,      &
             cfg => schwarz_op % cfg, D_inv => schwarz_op % D_inv )

    ! generic implementation ...................................................

    !$omp parallel
    !$acc data copyin(S1, W1, cfg, D_inv, f) copyout(r, u)

    ! r = reference result
    call SchwarzOp_Gen(n1, nc, nd, S1, W1, cfg, D_inv, f, r)

    call system_clock(count0, rate)
    do i = 1, nt
      call SchwarzOp_Gen(n1, nc, nd, S1, W1, cfg, D_inv, f, u)
      !$acc wait
    end do
    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_gen  = maxval(abs(u - r))
    mflops_gen = 1E-6 / time * nd * nflop
    mlups_gen  = 1E-6 / time * nd * nop

    ! parametrized implementation ..............................................

    !$omp parallel
    !$acc data copyin(S1, W1, cfg, D_inv, f) copyout(u)

    call SchwarzOp_Par(n1, nc, nd, S1, W1, cfg, D_inv, f, u)
    !$acc wait

    call system_clock(count0, rate)
    do i = 1, nt
      call SchwarzOp_Par(n1, nc, nd, S1, W1, cfg, D_inv, f, u)
      !$acc wait
    end do
    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_par  = maxval(abs(u - r))
    mflops_par = 1E-6 / time * ne * nflop
    mlups_par  = 1E-6 / time * ne * nop

  end associate

  !-----------------------------------------------------------------------------
  ! print results

  write(*,'(/,A,/)') 'Uniform (cubic) Schwarz operator'

  write(*,'(3A)') '#                                  ',   &
                  '   ------------ generic ------------',  &
                  '   --------- parametrized  ---------'
  write(*,'(3A)') '#  n1   n2   n3        nd        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',    &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(3I5,2(2X,I8))', advance='NO') n1, n2, n3, nd, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  if (parametrized) then
    write(*,'(3(2X,ES10.3))') error_par, mflops_par, mlups_par
  else
    write(*,'(3(8X,A))') 'None', 'None', 'None'
  end if
  write(*,*)

!===============================================================================

end program Validate__CART__TPO_Schwarz_Cubic
