!> \file       validate__cart__tpo_spectral.f
!> \brief      Validation of the tensor-product spectral operator
!> \author     Joerg Stiller
!> \date       2017/01/27
!> \copyright  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__CART__TPO_Spectral
  use Kind_Parameters, only: IXL, RNP
  use Constants,       only: ONE
  use Eigenproblems,   only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  use CART__TPO_Spectral
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: na = 7   ! operator dimension for direction 1
  integer :: nb = 8   ! operator dimension for direction 2
  integer :: nc = 9   ! operator dimension for direction 3
  integer :: ne = 1   ! number of elements
  integer :: nt = 1   ! number of test runs

  namelist /input/ na, nb, nc, ne, nt

  ! operators and variables ....................................................

  type(StandardOperators1D) :: sop(3)

  procedure(TPO_Spectral_Proc), pointer :: SpectralOperator_Gen ! generic
  procedure(TPO_Spectral_Proc), pointer :: SpectralOperator_Par ! parametrized

  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP), allocatable :: A(:,:), B(:,:), C(:,:), Lambda(:,:,:)
  real(RNP), allocatable :: lambda_A(:), lambda_B(:), lambda_C(:)
  real(RNP) :: tmp
  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_par, mflops_par, mlups_par

  logical :: exists, parametrized
  integer :: pa, pb, pc
  integer :: nflop, npop, prm
  integer :: e, i, j, k, l

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file='validate__cart__tpo_spectral.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__cart__tpo_spectral.prm')
    read(prm, nml=input)
    close(prm)
  end if

  ! dimensions
  npop  = na * nb * nc
  nflop = npop * (4*(na + nb + nc) + 1)

  ! corresponding polynomial orders
  pa = na - 1
  pb = nb - 1
  pc = nc - 1

  ! operators ..................................................................

  ! generic operator procedure
  call TPO_Spectral_Assign(-1, -1, -1, SpectralOperator_Gen)

  ! try parametrized operator procedure
  call TPO_Spectral_Assign(na, nb, nc, SpectralOperator_Par)

  parametrized = .not. associated( SpectralOperator_Par, &
                                   SpectralOperator_Gen  )

  ! workspace ..................................................................

  allocate( u(0:pa,0:pb,0:pc,ne), v(0:pa,0:pb,0:pc,ne), w(0:pa,0:pb,0:pc,ne) )
  allocate( A(na,na), B(nb,nb), C(nc,nc), Lambda(na,nb,nc) )
  allocate( lambda_A(na), lambda_B(nb), lambda_C(nc) )

  ! 1D operators ...............................................................

  ! standard operators
  call sop(1) % New(pa)
  call sop(2) % New(pb)
  call sop(3) % New(pc)

  ! 1D eigenvalues and eigenvectors
  call SolveGeneralizedEigenproblem(sop(1)%L, sop(1)%w, lambda_A, A)
  call SolveGeneralizedEigenproblem(sop(2)%L, sop(2)%w, lambda_B, B)
  call SolveGeneralizedEigenproblem(sop(3)%L, sop(3)%w, lambda_C, C)

  ! diagonal operator
  do k = 1, nc
  do j = 1, nb
  do i = 1, na
    Lambda(i,j,k) = ONE / (1 + lambda_A(i) + lambda_B(j) + lambda_C(k))
  end do
  end do
  end do

  ! result and operand .........................................................

  ! result
  call random_number(w)

  ! operand: u = H^e w
  associate( M1 => sop(1)%w, L1 => sop(1)%L, &
             M2 => sop(2)%w, L2 => sop(2)%L, &
             M3 => sop(3)%w, L3 => sop(3)%L  )

    do e = 1, ne

      do k = 0, pc
      do j = 0, pb
      do i = 0, pa
        tmp = M1(i) * M2(j) * M3(k) * w(i,j,k,e)
        do l = 0, pa
          tmp = tmp + L1(i,l) * M2(j) * M3(k) * w(l,j,k,e)
        end do
        do l = 0, pb
          tmp = tmp + M1(i) * L2(j,l) * M3(k) * w(i,l,k,e)
        end do
        do l = 0, pc
          tmp = tmp + M1(i) * M2(j) * L3(k,l) * w(i,j,l,e)
        end do
        u(i,j,k,e) = tmp
      end do
      end do
      end do

    end do

  end associate

  !-----------------------------------------------------------------------------
  ! test generic procedure

  !$omp parallel
  !$acc data copyin(u) copyout(v)

  call SpectralOperator_Gen(na, nb, nc, ne, A, B, C, Lambda, u, v)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call SpectralOperator_Gen(na, nb, nc, ne, A, B, C, Lambda, u, v)
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
  ! test parametrized procedure

  if (parametrized) then

    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call SpectralOperator_Par(na, nb, nc, ne, A, B, C, Lambda, u, v)
    !$acc wait

    call system_clock(count0, rate)

    do i = 1, nt
      call SpectralOperator_Par(na, nb, nc, ne, A, B, C, Lambda, u, v)
      !$acc wait
    end do

    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_par  = maxval(abs(v - w))
    mflops_par = 1E-6 / time * ne * nflop
    mlups_par  = 1E-6 / time * ne * npop

  end if

  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                                      ', &
                  '   ------------ generic ------------',    &
                  '   --------- parametrized  ---------'
  write(*,'(3A)') '#  na     nb     nc        ne        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',        &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I5),2(2X,I8))',  advance='NO') na, nb, nc, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  if (parametrized) then
    write(*,'(3(2X,ES10.3))') error_par, mflops_par, mlups_par
  else
    write(*,'(3(8X,A))') 'None', 'None', 'None'
  end if
  write(*,*)

!===============================================================================

end program Validate__CART__TPO_Spectral
