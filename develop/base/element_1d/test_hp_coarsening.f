!> summary:  Test of the hp coarsening operator
!> author:   Joerg Stiller
!> date:     2024/12/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Test_HP_Coarsening
  use Kind_Parameters
  use Constants
  use Standard_Element_Operators__1D
  use Embedded_Interpolation_Operator__1D
  use HP__Coarsening_Operator__1D
  implicit none

  character(len=80) :: test_case = ''

  integer      :: mode    =  2    ! p coarsening (1) or hp coarsening (2)
  integer      :: po_f    =  7    ! fine element order
  integer      :: po_c    =  7    ! coarse element order
  character(2) :: nodes   = 'L'   ! node type
  character    :: method  = 'I'   ! use interpolation 'I' or L²-projection 'P'
  integer      :: smooth  =  1    ! jump smoothing 0/1/2: none/linear/average
  real(RNP)    :: jmp_rot =  0.1  ! jump caused by rotation about endpoints
  real(RNP)    :: jmp_lft =  0.1  ! jump caused by lifting
  integer      :: n_s     =  80   ! number of sampling intervals per element

  namelist/input/ mode, po_f, po_c, nodes, method, smooth, jmp_rot, jmp_lft, n_s

  type(StandardElementOperators_1D)      :: eop_c
  type(StandardElementOperators_1D)      :: eop_f
  type(HP_CoarseningOperator_1D)         :: pop_fc
  type(HP_CoarseningOptions_1D)          :: opt_fc
  type(EmbeddedInterpolationOperator_1D) :: iop_fs
  type(EmbeddedInterpolationOperator_1D) :: iop_cs

  real(RNP), allocatable :: x_f (:,:) ! fine element nodes
  real(RNP), allocatable :: u_f (:,:) ! fine element coefficients
  real(RNP), allocatable :: u_c (:)   ! coarse element coefficients
  real(RNP), allocatable :: xi_s(:)   ! sampling point standard coordinates
  real(RNP), allocatable :: x_s (:)   ! sampling point coordinates
  real(RNP), allocatable :: u_fs(:)   ! fine solution
  real(RNP), allocatable :: u_cs(:)   ! coarse solution

  integer :: i, io, n_s2, stat

  ! initialization .............................................................

  ! read parameters
  call get_command_argument(1, test_case, status=stat)
  if (stat == 0) then
    open(newunit=io, file=trim(test_case)//'.prm')
    read(io,input)
    close(io)
  else
    test_case = ''
  end if

  eop_c = StandardElementOperators_1d(po_c, nodes)
  eop_f = StandardElementOperators_1d(po_f, nodes)

  opt_fc % nodes  = nodes
  opt_fc % po_f   = po_f
  opt_fc % po_c   = po_c
  opt_fc % mode   = mode
  opt_fc % method = method
  opt_fc % smooth = smooth

  pop_fc = HP_CoarseningOperator_1D(opt_fc)

  ! fine element nodes and coefficients ........................................

  select case(mode)
  case(1)
    allocate(x_f(0:po_f,1), source = reshape(eop_f % x, [po_f+1,1]))
    allocate(u_f(0:po_f,1), source = real(sin(PI * x_f), RNP))

  case(2)
    allocate(x_f(0:po_f,2))
    x_f(:,1) = HALF * (eop_f % x + 1) - 1
    x_f(:,2) = HALF * (eop_f % x + 1)
    allocate(u_f(0:po_f,2), source = real(sin(PI * x_f), RNP))
    ! jump by rotation
    u_f(:,1) = u_f(:,1) + HALF * jmp_rot * (x_f(:,1) + 1)
    u_f(:,2) = u_f(:,2) - HALF * jmp_rot *  x_f(:,2)
    ! jump by lifting
    u_f(:,1) = u_f(:,1) + HALF * jmp_lft
    u_f(:,2) = u_f(:,2) - HALF * jmp_lft

  end select

  ! fine-to-coarse projection ..................................................

  allocate(u_c(0:po_c), source = ZERO)
  do i = 1, mode
    u_c = u_c + matmul(pop_fc % A(:,:,i), u_f(:,i))
  end do

  ! interpolation to sampling points ...........................................

  allocate(xi_s(0:n_s), source = [(real(TWO*i/n_s - ONE, RNP), i = 0, n_s)])

  select case(mode)
  case(1)
    allocate(x_s(0:n_s), source = xi_s)
    iop_fs = EmbeddedInterpolationOperator_1D(eop_f, xi_s)
    iop_cs = EmbeddedInterpolationOperator_1D(eop_c, xi_s)
    allocate(u_fs, u_cs, mold = x_s)
    u_fs = matmul(iop_fs % A, u_f(:,1))
    u_cs = matmul(iop_cs % A, u_c)

  case(2)
    n_s2 = 2 * n_s + 1
    allocate(x_s(0:n_s2))
    x_s( 0    :n_s  ) = HALF * (xi_s + 1) - 1
    x_s( n_s+1:n_s2 ) = HALF * (xi_s + 1)
    iop_fs = EmbeddedInterpolationOperator_1D(eop_f, xi_s)
    iop_cs = EmbeddedInterpolationOperator_1D(eop_c, x_s)
    allocate(u_fs, u_cs, mold = x_s)
    u_fs( 0    :n_s  ) = matmul(iop_fs % A, u_f(:,1))
    u_fs( n_s+1:n_s2 ) = matmul(iop_fs % A, u_f(:,2))
    u_cs = matmul(iop_cs % A, u_c)

  end select

  ! output .....................................................................

  open(newunit = io, file = trim(test_case)//'.dat')
  write(io,'(A)') '#    x             u_f           u_c'
  do i = 0, ubound(x_s,1)
    write(io,'(3(ES12.5,2X))') x_s(i), u_fs(i), u_cs(i)
  end do
  close(io)

  !=============================================================================

end program Test_HP_Coarsening
