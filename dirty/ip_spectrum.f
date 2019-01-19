program IP_Spectrum
  use Kind_Parameters, only: RNP, RDP
  use Eigenproblems,   only: SolveSymmetricEigenproblem
  use IP_Element_Operators_1D
  implicit none

  character :: bc(2) = 'D'
  integer   :: ne = 2
  integer   :: po = 2
  real(RNP) :: lx = 1
  real(RNP) :: penalty = 2
  logical   :: hybrid  = .false.
  logical   :: write_A = .false.
  logical   :: write_V = .false.

  namelist /input/ bc, ne, po, lx, penalty, hybrid, write_a, write_v

  type(IP_ElementOperators1D) :: eop

  real(RNP), allocatable :: Le(:,:,:) ! stiffness matrix, L(0:po,0:po,-1:1)
  real(RDP), allocatable :: A(:,:)    ! system matrix, A(na,na), na = (po+1)*ne
  real(RDP), allocatable :: V(:,:)    ! eigenvectors, V(na,na)
  real(RDP), allocatable :: lambda(:) ! eigenvalues, lambda(na)

  real(RNP) :: dx(-1:1), delta
  logical   :: exists
  integer   :: na, np, io
  integer   :: i, i1, i2

  ! read input -----------------------------------------------------------------

  inquire(file='ip-input.prm', exist=exists)
  if (exists) then
    open(newunit=io, file='ip-input.prm')
    read(io, nml=input)
    close(io)
  end if

  ! initialization -------------------------------------------------------------

  np = po + 1
  na = np * ne

  call eop % New(po, penalty, hybrid)

  allocate(Le(0:po,0:po,-1:1))
  allocate(A(na,na), V(na,na), lambda(na))
  A = 0

  dx = lx / ne

  ! system matrix --------------------------------------------------------------

  select case(ne)

  case(1) ! single element .....................................................

    call eop % GetStiffnessMatrix(dx, bc, Le)
    if (all(bc == 'P')) then
      A = Le(:,:,-1) + Le(:,:,0) + Le(:,:,1)
    else
      A = Le(:,:,0)
    end if

  case(2) ! two elements .......................................................

    if (all(bc == 'P')) then

      call eop % GetStiffnessMatrix(dx, bc, Le)
      i1 = 1
      i2 = np
      A( i1:i2,  i1   :i2    ) = Le(:,:, 0)
      A( i1:i2,  i1+np:i2+np ) = Le(:,:,-1) + Le(:,:, 1)

      i1 = i1 + np
      i2 = i2 + np
      A( i1:i2,  i1-np:i2-np ) = Le(:,:,-1) + Le(:,:, 1)
      A( i1:i2,  i1   :i2    ) = Le(:,:, 0)

    else

      call eop % GetStiffnessMatrix(dx, [bc(1), ' '], Le)
      i1 = 1
      i2 = np
      A( i1:i2,  i1   :i2    ) = Le(:,:, 0)
      A( i1:i2,  i1+np:i2+np ) = Le(:,:, 1)

      call eop % GetStiffnessMatrix(dx, [' ', bc(2)], Le)
      i1 = i1 + np
      i2 = i2 + np
      A( i1:i2,  i1-np:i2-np ) = Le(:,:,-1)
      A( i1:i2,  i1   :i2    ) = Le(:,:, 0)

    end if

  case(3:) ! 3 or more elements ................................................

    ! first
    call eop % GetStiffnessMatrix(dx, [bc(1), ' '], Le)
    i1 = 1
    i2 = np
    A( i1:i2,  i1   :i2    ) = Le(:,:, 0)
    A( i1:i2,  i1+np:i2+np ) = Le(:,:, 1)
    if (bc(1) == 'P') then
      A( i1:i2, na-po:na ) = Le(:,:,-1)
    end if

    ! interior
    call eop % GetStiffnessMatrix(dx, [' ', ' '], Le)
    do i = 2, ne-1
      i1 = i1 + np
      i2 = i2 + np
      A( i1:i2,  i1-np:i2-np ) = Le(:,:,-1)
      A( i1:i2,  i1   :i2    ) = Le(:,:, 0)
      A( i1:i2,  i1+np:i2+np ) = Le(:,:, 1)
    end do

    ! last
    call eop % GetStiffnessMatrix(dx, [' ', bc(2)], Le)
    i1 = i1 + np
    i2 = i2 + np
    A( i1:i2,  i1-np:i2-np ) = Le(:,:,-1)
    A( i1:i2,  i1   :i2    ) = Le(:,:, 0)
    if (bc(2) == 'P') then
      A( i1:i2, 1:np ) = Le(:,:,1)
    end if

  end select

  ! check symmetry .............................................................

  delta = maxval(abs(A - transpose(A)))
  print '(A,ES10.3)', '|A - Aᵀ| ≤', delta

  ! solve eigenproblem ---------------------------------------------------------

  call SolveSymmetricEigenproblem(A, lambda, V)

  if (any(bc == 'D')) then
    i = 1
  else
    i = 2
  end if
  print '(A,ES10.3)', 'cond(A)  =', lambda(na) / lambda(i)

  ! save results ---------------------------------------------------------------

  open(newunit=io, file='ip-lambda.dat')
  do i = 1, na
    write(io,'(ES15.8)') lambda(i)
  end do
  close(io)

end program IP_Spectrum
