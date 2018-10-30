!> \file       validate__cart__tpo_schwarz_iso.f
!> \brief      Validation of the isotropic tensor-product spectral operator
!> \author     Joerg Stiller
!> \date       2017/12/13
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__CART__TPO_Schwarz_Iso
  use Kind_Parameters,   only: IXL, RNP
  use Constants,         only: ONE, THIRD, ZERO
  use Array_Assignments, only: AssignScalar
  use Eigenproblems,     only: SolveGeneralizedEigenproblem
  use CART__TPO_Schwarz_Iso
  use CART__DG_Element_Operators
  use CART__DG_Diffusion_CI_Schwarz
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: po = -1       ! polynomial order
  integer :: np =  8       ! operator dimension = po + 1
  integer :: ne =  1       ! number of elements
  integer :: nt =  1       ! number of test runs
  integer :: scenario = 1  ! 1/2/3 = interior / cube with mixed BC / random

  namelist /input/ po, np, ne, nt
  namelist /control/ scenario

  ! operators and variables ....................................................

  integer   :: nc
  real(RNP) :: lambda   =  1  ! Helmholtz parameter
  real(RNP) :: nu       =  1  ! diffusivity
  real(RNP) :: dx(3)    =  1  ! element extensions
  real(RNP) :: penalty  =  2  ! penalty coefficient
  real(RNP) :: delta(3) = -1  ! enforce zero overlap

  type(ElementOperators) :: eop     ! DG element oprators
  type(SchwarzOperator)  :: schwarz ! Schwarz operator and procedures

  procedure(TPO_Schwarz_Iso_Proc), pointer :: SchwarzOp_Gen
  procedure(TPO_Schwarz_Iso_Proc), pointer :: SchwarzOp_Par

  integer,   allocatable :: ec(:,:)
  real(RNP), allocatable :: nu_e(:)
  real(RNP), allocatable :: f(:,:,:,:), u(:,:,:,:), r(:,:,:,:)

  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_par, mflops_par, mlups_par

  logical :: exists, parametrized
  integer :: nflop, nop, prm
  integer :: i

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file='validate__cart__tpo_schwarz_iso.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__cart__tpo_schwarz_iso.prm')
    read(prm, nml=input)
    read(prm, nml=control)
    close(prm)
  end if

  if (po > 0) then
    np = po + 1
  else
    po = np - 1
  end if

  ! number of operands and floating point operations per element
  nop   = np * np * np
  nflop = nop * (12*np + 7) + 4

  ! operators, procedures, data ................................................

  ! element-configurations
  allocate(ec(3,ne))
  select case(scenario)
  case(2)
    call CubicDomainConfiguration(ec)
  case(3)
    call RandomConfiguration(ec)
  case default
    ec = 1
  end select

  ! operators
  call eop % New(po, dx, penalty)
  call schwarz % New(eop, delta)
  nc = size(schwarz % S1, 3)

  ! procedures
  call TPO_Schwarz_Iso_Assign(-1, SchwarzOp_Gen)  ! generic, for reference
  call TPO_Schwarz_Iso_Assign(np, SchwarzOp_Par)  ! parametrized

  parametrized = .not. associated( SchwarzOp_Par, &
                                   SchwarzOp_Gen  )

  ! element-mean diffusivity
  allocate(nu_e(ne))
  call AssignScalar(ne, nu_e, nu)

  ! workspace
  allocate(f(np,np,np,ne), u(np,np,np,ne), r(np,np,np,ne))
  call AssignScalar(size(f), f, ZERO)
  call AssignScalar(size(u), u, ZERO)
  call AssignScalar(size(r), r, ZERO)

  associate( S  => schwarz % S1, &
             V1 => schwarz % V1, &
             V2 => schwarz % V2, &
             V3 => schwarz % V3, &
             W  => schwarz % W1  )

    ! operand and result
    call random_number(f)
    call SchwarzOp_Gen(np, nc, ne, S, V1, V2, V3, W, ec, dx, lambda, nu_e, f, r)

    ! generic implementation ...................................................

    !$omp parallel
    !$acc data copyin(S,V1,V2,V3,W,ec,nu_e,f) copyout(u)

    call SchwarzOp_Gen(np, nc, ne, S, V1, V2, V3, W, ec, dx, lambda, nu_e, f, u)
    !$acc wait

    call system_clock(count0, rate)
    do i = 1, nt
      call SchwarzOp_Gen(np, nc, ne, S, V1, V2, V3, W, ec, dx, lambda, nu_e, f, u)
      !$acc wait
    end do
    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_gen  = maxval(abs(u - r))
    mflops_gen = 1E-6 / time * ne * nflop
    mlups_gen  = 1E-6 / time * ne * nop

    ! parametrized implementation ..............................................

    !$omp parallel
    !$acc data copyin(S,V1,V2,V3,W,ec,nu_e,f) copyout(u)

    call SchwarzOp_par(np, nc, ne, S, V1, V2, V3, W, ec, dx, lambda, nu_e, f, u)
    !$acc wait

    call system_clock(count0, rate)
    do i = 1, nt
      call SchwarzOp_par(np, nc, ne, S, V1, V2, V3, W, ec, dx, lambda, nu_e, f, u)
      !$acc wait
    end do
    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_par  = maxval(abs(u - r))
    mflops_par = 1E-6 / time * ne * nflop
    mlups_par  = 1E-6 / time * ne * nop

  end associate

  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                        ',             &
                  '   ------------ generic ------------',  &
                  '   --------- parametrized  ---------'
  write(*,'(3A)') '#  np        ne        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',        &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I8))',  advance='NO') np, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  if (parametrized) then
    write(*,'(3(2X,ES10.3))') error_par, mflops_par, mlups_par
  else
    write(*,'(3(8X,A))') 'None', 'None', 'None'
  end if
  write(*,*)

contains

  !-----------------------------------------------------------------------------
  !> Element configuration mimicking cubic domain

  subroutine CubicDomainConfiguration(ec)
    integer, intent(out) :: ec(:,:)

    integer :: i, j, k, l, m, n

    n = size(ec,2)
    m = int( real(n, RNP)**THIRD )

    ec = 1

    ! west
    i = 1
    do k = 1, m
    do j = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(1,l) = 2   ! Dirichlet - Interior
      end if
    end do
    end do

    ! east
    i = m
    do k = 1, m
    do j = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(1,l) = 3   ! Interior  - Dirichlet
      end if
    end do
    end do

    ! south
    j = 1
    do k = 1, m
    do i = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(2,l) = 4   ! Neumann   - Interior
      end if
    end do
    end do

    ! north
    j = m
    do k = 1, m
    do i = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(2,l) = 5   ! Interior  - Neumann
      end if
    end do
    end do

    ! bottom
    k = 1
    do j = 1, m
    do i = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(3,l) = 2   ! Dirichlet - Interior
      end if
    end do
    end do

    ! top
    k = m
    do j = 1, m
    do i = 1, m
      l = i + m * (j-1 + m * (k-1))
      if (l <= n) then
        ec(3,l) = 5   ! Interior  - Neumann
      end if
    end do
    end do

  end subroutine CubicDomainConfiguration

  !-----------------------------------------------------------------------------
  !> Random element configuration

  subroutine RandomConfiguration(ec)
    integer, intent(out) :: ec(:,:)

    real, allocatable :: rc(:,:)

    allocate(rc( size(ec,1), size(ec,2) ))
    call random_number(rc)
    ec = nint(9 * rc - 0.5)
    ec = max(1, min(9, ec))

  end subroutine RandomConfiguration

!===============================================================================

end program Validate__CART__TPO_Schwarz_Iso
