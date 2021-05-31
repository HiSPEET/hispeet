!> summary:  Validation of 3d constant anisotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__Schwarz__3D_CA
  use Kind_Parameters,   only: IXL, RNP
  use Constants,         only: ZERO
  use Array_Assignments, only: SetArray

  use TPO__Schwarz__3D_CA
  use TPO__Schwarz__3D_CA__Gen

  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: n1 =  6      ! subdomain points in direction 1
  integer :: n2 =  7      ! subdomain points in direction 2
  integer :: n3 =  8      ! subdomain points in direction 3
  integer :: nc =  9      ! number of subdomain configurations
  integer :: nd =  1      ! number of subdomains
  integer :: nt =  1      ! number of test runs

  namelist /input/ n1, n2, n3, nc, nd, nt

  ! operators and operands .....................................................

  real(RNP), allocatable :: S1(:,:,:), W1(:,:)
  real(RNP), allocatable :: S2(:,:,:), W2(:,:)
  real(RNP), allocatable :: S3(:,:,:), W3(:,:)
  real(RNP), allocatable :: D_inv(:,:,:,:)
  integer,   allocatable :: cfg(:,:)

  real(RNP), allocatable :: f(:,:,:,:), u(:,:,:,:), r(:,:,:,:)

  ! auxiliary ..................................................................

  character(len=80) :: input_file = 'validate__tpo__schwarz_3d_ca.prm'

  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_opt, mflops_opt, mlups_opt
  real(RNP) :: c(3)

  logical :: exists
  integer :: nflop, nop, prm
  integer :: i

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file=input_file, exist=exists)
  if (exists) then
    open(newunit=prm, file=input_file)
    read(prm, nml=input)
    close(prm)
  end if

  ! operators and operands .....................................................

  allocate(S1(n1,n1,nc), W1(n1,nc))
  allocate(S2(n2,n2,nc), W2(n2,nc))
  allocate(S3(n3,n3,nc), W3(n3,nc))
  allocate(cfg(3,nd))
  allocate(D_inv(n1,n2,n3,nd))
  allocate(f, u, r, mold = D_inv)

  !$omp parallel
  !$omp do
  do i = 1, nd
    call random_number(c)
    cfg(:,i) = max(1, min(nc, int(nc*c + 1)))
  end do
  call SetArray(D_inv, ZERO)
  call SetArray(u    , ZERO)
  call SetArray(r    , ZERO)
  call SetArray(f    , ZERO)
  !$omp end parallel


  call random_number(S1)
  call random_number(W1)

  call random_number(S2)
  call random_number(W2)

  call random_number(S3)
  call random_number(W3)

  call random_number(D_inv)
  call random_number(f)

  ! problem dimensions .........................................................

  nop   = n1 * n2 * n3
  nflop = nop * (4 * (n1 + n2 + n3) + 1)

  !-----------------------------------------------------------------------------
  ! test of generic implementation

  !$omp parallel
  !$acc data copyin(S1, S2, S3, W1, W2, W3, cfg, D_inv, f) copyout(u)
  !$acc data create(r)

  ! r = reference result
  call TPO_Schwarz_CA_Gen(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, r)
  !$acc wait

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz_CA_Gen(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
    !$acc wait
  end do
  call system_clock(count)

  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(u - r))
  mflops_gen = 1E-6 / time * nd * nflop
  mlups_gen  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! test of optimized implementation

  call random_number(u)

  !$omp parallel
  !$acc data copyin(S1, S2, S3, W1, W2, W3, cfg, D_inv, f) copyout(u) 
  !$acc data create(r)

  call TPO_Schwarz_CA(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
  !$acc wait

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz_CA(S1, S2, S3, W1, W2, W3, cfg, D_inv, f, u)
    !$acc wait
  end do
  call system_clock(count)

  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_opt  = maxval(abs(u - r))
  mflops_opt = 1E-6 / time * nd * nflop
  mlups_opt  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! print results

  write(*,'(/,A,/)') 'Constant anisotropic Schwarz operator'

  write(*,'(3A)') '#                                  ',   &
                  '   ------------ generic ------------',  &
                  '   ----------- optimized -----------'
  write(*,'(3A)') '#  n1   n2   n3        nd        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',    &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(3I5,2(2X,I8))', advance='NO') n1, n2, n3, nd, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

!===============================================================================

end program Validate__TPO__Schwarz__3D_CA
