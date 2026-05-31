program Schwarz_Overlap
  use Kind_Parameters
  use Gauss_Jacobi
  implicit none

  real(RNP) :: delta
  integer   :: p

  write(*,'(/,A,2/,2X,A)',advance='NO') 'Schwarz overlap', 'delta = '
  read *, delta

  write(*,'(/,2A5,/)') 'po', 'no'

  do p = 1, 32
    write(*,'(2I5)') p, count(LobattoPoints(p) <= 2 * delta - 1)
  end do

end program Schwarz_Overlap
