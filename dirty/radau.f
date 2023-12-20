program Radau
  use Kind_Parameters
  use Constants
  use Gauss_Jacobi
  implicit none

  real(RNP), allocatable :: xr(:), wr(:), Dr(:,:)
  real(RNP), allocatable :: x(:), lr(:,:), d_lr(:,:)

  logical :: right
  integer :: po, ns, io
  integer :: i, j, k

  write(*,'(A,/)') 'Testing Radau quadrature and polynomials'
  write(*,'(2X,A)',advance='NO') 'degree:          po = '; read *, po
  write(*,'(2X,A)',advance='NO') 'sampling:        ns = '; read *, ns
  write(*,'(2X,A)',advance='NO') 'right (T/F):  right = '; read *, right

  allocate(xr(0:po), wr(0:po), Dr(0:po,0:po), source = ZERO)
  allocate(x(0:ns), lr(0:ns,0:po), d_lr(0:ns,0:po), source = ZERO)

  xr = RadauPoints(po, right)
  wr = RadauWeights(xr)
  Dr = RadauDiffMatrix(xr)

  write(*,'(/,A)') 'Quadrature points and weights:'
  do k = 0, po
    write(*,'(2(2X,ES16.9))') xr(k), wr(k)
  end do

  do i = 0, ns
    x(i) = TWO*i/ns - ONE
    do k = 0, po
      lr(i,k) = RadauPolynomial(k, xr, x(i))
      do j = 0, po
        d_lr(i,j) = d_lr(i,j) + Dr(k,j) * lr(i,k)
      end do
    end do
  end do

  open(newunit = io, file = 'radau.dat')

  write(io,'(A)', advance='NO') '# x'
  do k = 0, po
    write(io,'(A,I0)', advance='NO') ', L_', k
  end do
  do k = 0, po
    write(io,'(A,I0)', advance='NO') ', dL_', k
  end do
  write(io, '(A)', advance='YES')

  do i = 0, ns
    write(io,'(99(ES16.9,2X))') x(i), lr(i,:), d_lr(i,:)
  end do

  close(io)

end program Radau
