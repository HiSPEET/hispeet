!> summary:  Validation of the 3D anisotropic spectral TPO
!> author:   Joerg Stiller, Erik Pfister
!> date:     2020/05/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__Spectral__3D_CA
  use Kind_Parameters, only: IXL, RNP
  use TPO__Spectral__3D_CA
  use TPO__Spectral__3D_CA__Gen
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: n1 = 4   ! number of element-domain points in direction 1
  integer :: n2 = 4   ! number of element-domain points in direction 2
  integer :: n3 = 4   ! number of element-domain points in direction 3
  integer :: ne = 1   ! number of elements
  integer :: nt = 1   ! number of test runs

  namelist /input/ n1, n2, n3, ne, nt

  ! operators and variables ....................................................

  real(RNP), dimension(:,:,:,:), allocatable :: u, v, w
  real(RNP), dimension(:,:,:),   allocatable :: Lambda
  real(RNP), dimension(:,:),     allocatable :: S1, S2, S3

  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_opt, mflops_opt, mlups_opt

  logical :: exists
  integer :: nflop, npop, prm
  integer :: i

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file='validate__tpo__spectral_3d_ca.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__tpo__spectral_3d_ca.prm')
    read(prm, nml=input)
    close(prm)
  end if

  ! problem dimensions
  nflop = n1*n2*n3 * (4*(n1 + n2 + n3) + 1)
  npop  = n1*n2*n3

  ! workspace ..................................................................

  allocate(u(n1,n2,n3,ne), v(n1,n2,n3,ne), w(n1,n2,n3,ne))
  allocate(Lambda(n1,n2,n3), S1(n1,n1), S2(n2,n2), S3(n3,n3))

  !$omp do
  do i = 1, ne
    call random_number(u(:,:,:,i))
    call random_number(v(:,:,:,i))
    call random_number(w(:,:,:,i))
  end do

  call random_number(Lambda)
  call random_number(S1)
  call random_number(S2)
  call random_number(S3)

  !-----------------------------------------------------------------------------
  ! test generic operator

  !$acc data copyin(u,S1,S2,S3,Lambda) copyout(v) create(w)
  !$omp parallel

  call TPO_Spectral_CA_Gen(S1, S2, S3, Lambda, u, w)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_Spectral_CA_Gen(S1, S2, S3, Lambda, u, v)
    !$acc wait
  end do

  call system_clock(count)
  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(v - w))
  mflops_gen = 1E-6 / time * ne * nflop
  mlups_gen  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! test optimized operator

  call random_number(v)

  !$acc data copyin(u,S1,S2,S3,Lambda) copyout(v) create(w)
  !$omp parallel

  call TPO_Spectral(S1, S2, S3, Lambda, u, v)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_Spectral(S1, S2, S3, Lambda, u, v)
    !$acc wait
  end do

  call system_clock(count)

  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_opt  = maxval(abs(v - w))
  mflops_opt = 1E-6 / time * ne * nflop
  mlups_opt  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                                    ',   &
                  ' ------- generic operator --------  ', &
                  ' ------- optimized operator ------'

  write(*,'(3A)') '#  n1   n2   n3        ne        nt    ',  &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s    '

  write(*,'(3I5,2(2X,I8))', advance='NO') n1, n2, n3, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

!===============================================================================

end program Validate__TPO__Spectral__3D_CA
