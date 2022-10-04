!> summary:  Validation of 3d constant isotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__Schwarz__3D_I_RDP
  use Kind_Parameters,   only: IXL, RDP
  use Constants,         only: ZERO
  use Array_Assignments, only: SetArray

  use TPO__Schwarz__3D_I
  use TPO__Schwarz__3D_I__Gen

  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: np =  7      ! subdomain points per direction
  integer :: nc =  9      ! number of subdomain configurations
  integer :: nd =  1      ! number of subdomains
  integer :: nt =  1      ! number of test runs

  namelist /input/ np, nc, nd, nt

  ! operators, parameters and operands .........................................

  real(RDP) :: lambda = 1
  real(RDP), allocatable :: f(:,:,:,:), r(:,:,:,:), u(:,:,:,:), nu(:)
  real(RDP), allocatable :: S(:,:,:), V(:,:), W(:,:), g(:,:)
  integer,   allocatable :: cfg(:,:)

  ! auxiliary ..................................................................

  character(len=80) :: input_file = 'validate__tpo__schwarz_3d_i.prm'

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

  ! operators, parameters and operands .........................................

  allocate(f(np,np,np,nd))
  allocate(u, r, mold = f)
  allocate(S(np,np,nc), V(np,nc), W(np,nc), g(4,nd), nu(nd))
  allocate(cfg(3,nd))

  ! initialization performed in parallel to exploit first touch principle

  !$omp parallel
  !$omp do
  do i = 1, nd
    call random_number(c)
    cfg(:,i) = max(1, min(nc, int(nc*c + 1)))
  end do
  call SetArray(u, ZERO)
  call SetArray(r, ZERO)
  call SetArray(f, ZERO)
  !$omp end parallel

  call random_number(S)
  call random_number(V)
  call random_number(W)
  call random_number(g)
  call random_number(f)
  call random_number(nu)

  ! enforce positivity
  V  = V  + 1
  g  = g  + 1
  nu = nu + 1

  ! problem dimensions .........................................................

  nop   = np ** 3
  nflop = nop * (12*np + 7) + 4

  !-----------------------------------------------------------------------------
  ! test of generic implementation

  !$omp parallel

  ! r = reference result
  call TPO_Schwarz_I_Gen(S, V, W, g, cfg, lambda, nu, f, r)

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz_I_Gen(S, V, W, g, cfg, lambda, nu, f, u)
  end do
  call system_clock(count)

  !$omp end parallel

  time = (count - count0) / real(rate, RDP) / nt

  error_gen  = maxval(abs(u - r))
  mflops_gen = 1E-6 / time * nd * nflop
  mlups_gen  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! test of optimized implementation

  call random_number(u)

  !$omp parallel

  call TPO_Schwarz(S, V, W, g, cfg, lambda, nu, f, u)

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz(S, V, W, g, cfg, lambda, nu, f, u)
  end do
  call system_clock(count)

  !$omp end parallel

  time = (count - count0) / real(rate, RDP) / nt

  error_opt  = maxval(abs(u - r))
  mflops_opt = 1E-6 / time * nd * nflop
  mlups_opt  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! print results

  write(*,'(/,A,/)') 'Constant isotropic Schwarz operator, real(RDP)'

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

end program Validate__TPO__Schwarz__3D_I_RDP
