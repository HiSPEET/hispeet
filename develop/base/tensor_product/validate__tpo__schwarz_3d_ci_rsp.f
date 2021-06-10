!> summary:  Validation of 3d constant isotropic Schwarz TPO for real(RSP)
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__Schwarz__3D_CI_RSP
  use Kind_Parameters,   only: IXL, RDP, RSP
  use Constants,         only: ZERO
  use Array_Assignments, only: SetArray

  use TPO__Schwarz__3D_CI
  use TPO__Schwarz__3D_CI__Gen

  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: np =  7      ! subdomain points per direction
  integer :: nc =  9      ! number of subdomain configurations
  integer :: nd =  1      ! number of subdomains
  integer :: nt =  1      ! number of test runs

  namelist /input/ np, nc, nd, nt

  ! operators and operands .....................................................

  real(RSP), allocatable :: S(:,:,:), W(:,:)
  real(RSP), allocatable :: D_inv(:,:,:,:)
  integer,   allocatable :: cfg(:,:)

  real(RSP), allocatable :: f(:,:,:,:), u(:,:,:,:), r(:,:,:,:)

  ! auxiliary ..................................................................

  character(len=80) :: input_file = 'validate__tpo__schwarz_3d_ci.prm'

  real(RDP) :: time
  real(RDP) :: error_gen, mflops_gen, mlups_gen
  real(RDP) :: error_opt, mflops_opt, mlups_opt
  real(RDP) :: c(3)

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

  allocate(S(np,np,nc), W(np,nc))
  allocate(cfg(3,nd))
  allocate(D_inv(np,np,np,nd))
  allocate(f, u, r, mold = D_inv)

  !$omp parallel
  !$omp do
  do i = 1, nd
    call random_number(c)
    cfg(:,i) = max(1, min(nc, int(nc*c + 1)))
    D_inv (:,:,:,i) = ZERO
    u     (:,:,:,i) = ZERO
    r     (:,:,:,i) = ZERO
    f     (:,:,:,i) = ZERO
  end do
  !$omp end parallel

  call random_number(S)
  call random_number(W)
  call random_number(D_inv)
  call random_number(f)

  ! problem dimensions .........................................................

  nop   = np **3
  nflop = nop * (12*np + 1)

  !-----------------------------------------------------------------------------
  ! test of generic implementation

  !$omp parallel
  !$acc data copyin(S, W, cfg, D_inv, f) copyout(u) create(r)

  ! r = reference result
  call TPO_Schwarz_CI_Gen(S, W, cfg, D_inv, f, r)
  !$acc wait

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz_CI_Gen(S, W, cfg, D_inv, f, u)
    !$acc wait
  end do
  call system_clock(count)

  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RDP) / nt

  error_gen  = maxval(abs(u - r))
  mflops_gen = 1E-6 / time * nd * nflop
  mlups_gen  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! test of optimized implementation

  call random_number(u)

  !$omp parallel
  !$acc data copyin(S, W, cfg, D_inv, f) copyout(u) create(r)

  call TPO_Schwarz(S, W, cfg, D_inv, f, u)
  !$acc wait

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz(S, W, cfg, D_inv, f, u)
    !$acc wait
  end do
  call system_clock(count)

  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RDP) / nt

  error_opt  = maxval(abs(u - r))
  mflops_opt = 1E-6 / time * nd * nflop
  mlups_opt  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! print results

  write(*,'(/,A,/)') 'Constant isotropic Schwarz operator, real(RSP)'

  write(*,'(3A)') '#                        ',   &
                  '   ------------ generic ------------',  &
                  '   ----------- optimized -----------'
  write(*,'(3A)') '#  np        nd        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',    &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I8))', advance='NO') np, nd, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

!===============================================================================

end program Validate__TPO__Schwarz__3D_CI_RSP
