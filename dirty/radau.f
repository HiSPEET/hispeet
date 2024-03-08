program Radau
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  implicit none

  real(RNP), allocatable :: xr(:), wr(:), Dr(:,:), fr(:), Dfr(:)
  real(RNP), allocatable :: x(:), l(:,:), Dl(:,:), f(:), Df(:)

  logical :: right
  integer :: po, ns, io
  integer :: i, j, k

  ! initialization .............................................................

  write(*,'(A,/)') 'Testing Radau quadrature and polynomials'
  write(*,'(2X,A)',advance='NO') 'degree:          po = '; read *, po
  write(*,'(2X,A)',advance='NO') 'sampling:        ns = '; read *, ns
  write(*,'(2X,A)',advance='NO') 'right (T/F):  right = '; read *, right

  allocate(xr(0:po), wr(0:po), Dr(0:po,0:po), fr(0:po), Dfr(0:po),  source=ZERO)
  allocate(x(0:ns), l(0:ns,0:po), Dl(0:ns,0:po), f(0:ns), Df(0:ns), source=ZERO)

  xr = RadauPoints(po, right)
  wr = RadauWeights(xr)
  Dr = RadauDiffMatrix(xr)

  ! basis test .................................................................

  write(*,'(/,A)') 'Quadrature points and weights:'
  do k = 0, po
    write(*,'(2(2X,ES16.9))') xr(k), wr(k)
  end do

  do i = 0, ns
    x(i) = TWO*i/ns - ONE
    do k = 0, po
      l(i,k) = RadauPolynomial(k, xr, x(i))
      do j = 0, po
        Dl(i,j) = Dl(i,j) + Dr(k,j) * l(i,k)
      end do
    end do
  end do

  open(newunit = io, file = 'radau_basis.dat')

  write(io,'(A)', advance='NO') '# x'
  do k = 0, po
    write(io,'(A,I0)', advance='NO') ', L_', k
  end do
  do k = 0, po
    write(io,'(A,I0)', advance='NO') ', dL_', k
  end do
  write(io, '(A)', advance='YES')

  do i = 0, ns
    write(io,'(99(ES16.9,2X))') x(i), l(i,:), Dl(i,:)
  end do

  close(io)

  ! integration test ...........................................................

  write(*,'(/,A)') 'Quadrature:'

  do k = 0, po
    fr(k) = sin(PI * xr(k))
  end do

  write(*,'(2X,A,ES0.3)') 'int(sin(x)) = ', sum(wr * fr)

  ! interpolation and differentiation test .....................................

  write(*,'(/,A)') 'Interpolation and differentiation:'

  Dfr = matmul(Dr, fr)

  f(0:ns)  = matmul(l, fr)
  Df(0:ns) = matmul(l, Dfr)

  open(newunit = io, file = 'radau_test.dat')
  write(io,'(A)') '# x, f, Df'
  do i = 0, ns
    write(io,'(99(ES16.9,2X))') x(i), f(i), Df(i)
  end do
  close(io)

  write(*,'(2X,A,ES0.3)') "err(f)  = ", maxval(abs(f  - sin(PI * x)))
  write(*,'(2X,A,ES0.3)') "err(f') = ", maxval(abs(Df - PI * cos(PI * x)))

end program Radau
