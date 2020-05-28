!> summary:  Validation of the tensor-product operator for 3d interpolation
!> author:   Joerg Stiller, Erik Pfister
!> date:     2020/05/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__AAA_3d
  use Kind_Parameters, only: IXL, RNP
  use TPO__AAA_3d
  use TPO__AAA_3d__Gen
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: np = 4   ! number of element points per direction in u
  integer :: nq = 4   ! number of element points per direction in v and w
  integer :: ne = 1   ! number of elements
  integer :: nt = 1   ! number of test runs

  namelist /input/ np, nq, ne, nt

  ! operators and variables ....................................................

  real(RNP), dimension(:,:,:,:), allocatable :: u, v, w
  real(RNP), dimension(:,:),     allocatable :: A

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

  inquire(file='validate__tpo__aaa_3d.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__tpo__aaa_3d.prm')
    read(prm, nml=input)
    close(prm)
  end if

  ! problem dimensions
  nflop = np**3 * 2 * np * nq * (np*np + np*nq + nq*nq)
  npop  = np**3

  ! workspace ..................................................................

  allocate(u(np,np,np,ne), v(nq,nq,nq,ne), w(nq,nq,nq,ne), A(nq,np))

  !$omp do
  do i = 1, ne
    call random_number(u(:,:,:,i))
    call random_number(v(:,:,:,i))
    call random_number(w(:,:,:,i))
  end do

  call random_number(A)

  !-----------------------------------------------------------------------------
  ! test generic operator

  !$omp parallel
  !$acc data copyin(u) copyout(v) create(w)

  call TPO_AAA_Gen(A, u, w)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_AAA_Gen(A, u, v)
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

  call random_number(v(:,:,:,i))

  !$omp parallel
  !$acc data copyin(u) copyout(v)

  call TPO_AAA(A, u, v)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_AAA(A, u, v)
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
  write(*,'(3A)') '#                                 ',   &
                  ' ------- generic operator --------  ', &
                  ' ------- optimized operator ------'

  write(*,'(3A)') '#  np     nq        ne        nt    ',  &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s    '

  write(*,'(I5,2X,I5,2(2X,I8))',  advance='NO') np, nq, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

!===============================================================================

end program Validate__TPO__AAA_3d
