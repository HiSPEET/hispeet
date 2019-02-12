!> summary:  Validation of the isotropic tensor-product spectral operator
!> author:   Joerg Stiller
!> date:     2017/07/21
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__CART__TPO_Spectral_Iso
  use Kind_Parameters, only: IXL, RNP
  use Constants,       only: ONE
  use Eigenproblems,   only: SolveGeneralizedEigenproblem
  use Standard_Operators_1D
  use CART__TPO_Spectral_Iso
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: po = -1   ! polynomial order
  integer :: na =  8   ! operator dimension = po + 1
  integer :: ne =  1   ! number of elements
  integer :: nt =  1   ! number of test runs

  namelist /input/ po, na, ne, nt

  ! operators and variables ....................................................

  type(StandardOperators1D) :: standard_op

  procedure(TPO_Spectral_Iso_Proc), pointer :: SpectralOperator_Gen ! generic
  procedure(TPO_Spectral_Iso_Proc), pointer :: SpectralOperator_Par ! param.

  real(RNP), allocatable :: u(:,:,:,:), v(:,:,:,:), w(:,:,:,:)
  real(RNP), allocatable :: A(:,:), Lambda(:,:,:)
  real(RNP), allocatable :: lambda_A(:)
  real(RNP) :: tmp
  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_par, mflops_par, mlups_par

  logical :: exists, parametrized
  integer :: nflop, nop, prm
  integer :: e, i, j, k, p

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file='validate__cart__tpo_spectral_iso.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__cart__tpo_spectral_iso.prm')
    read(prm, nml=input)
    close(prm)
  end if

  if (po > 0) then
    na = po + 1
  else
    po = na - 1
  end if

  ! dimensions
  nop   = na * na * na
  nflop = nop * (12*na + 1)

  ! operators ..................................................................

  ! generic operator procedure
  call TPO_Spectral_Iso_Assign(-1, SpectralOperator_Gen)

  ! try parametrized operator procedure
  call TPO_Spectral_Iso_Assign(na, SpectralOperator_Par)

  parametrized = .not. associated( SpectralOperator_Par, &
                                   SpectralOperator_Gen  )

  ! workspace ..................................................................

  allocate( u(0:po,0:po,0:po,ne), v(0:po,0:po,0:po,ne), w(0:po,0:po,0:po,ne) )
  allocate( A(na,na), lambda_A(na), Lambda(na,na,na) )

  ! 1D operators ...............................................................

  ! standard operators
  standard_op = StandardOperators1D(po)

  ! 1D eigenvalues and eigenvectors
  call SolveGeneralizedEigenproblem(standard_op%L, standard_op%w, lambda_A, A)

  ! diagonal operator
  do k = 1, na
  do j = 1, na
  do i = 1, na
    Lambda(i,j,k) = ONE / (1 + lambda_A(i) + lambda_A(j) + lambda_A(k))
  end do
  end do
  end do

  ! result and operand .........................................................

  ! result
  call random_number(w)

  ! operand: u = H^e w
  associate( M => standard_op%w, L => standard_op%L )

    do e = 1, ne

      do k = 0, po
      do j = 0, po
      do i = 0, po
        tmp = M(i) * M(j) * M(k) * w(i,j,k,e)
        do p = 0, po
          tmp = tmp + L(i,p) * M(j) * M(k) * w(p,j,k,e)
        end do
        do p = 0, po
          tmp = tmp + M(i) * L(j,p) * M(k) * w(i,p,k,e)
        end do
        do p = 0, po
          tmp = tmp + M(i) * M(j) * L(k,p) * w(i,j,p,e)
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

  call SpectralOperator_Gen(na, ne, A, Lambda, u, v)
  !$acc wait

  call system_clock(count0, rate)

  do i = 1, nt
    call SpectralOperator_Gen(na, ne, A, Lambda, u, v)
    !$acc wait
  end do

  call system_clock(count)
  !$acc end data
  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(v - w))
  mflops_gen = 1E-6 / time * ne * nflop
  mlups_gen  = 1E-6 / time * ne * nop

  !-----------------------------------------------------------------------------
  ! test parametrized procedure

  if (parametrized) then

    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call SpectralOperator_Par(na, ne, A, Lambda, u, v)
    !$acc wait

    call system_clock(count0, rate)

    do i = 1, nt
      call SpectralOperator_Par(na, ne, A, Lambda, u, v)
      !$acc wait
    end do

    call system_clock(count)

    !$acc end data
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_par  = maxval(abs(v - w))
    mflops_par = 1E-6 / time * ne * nflop
    mlups_par  = 1E-6 / time * ne * nop

  end if

  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                        ',             &
                  '   ------------ generic ------------',  &
                  '   --------- parametrized  ---------'
  write(*,'(3A)') '#  na        ne        nt    ', &
                  '   error     MFLOP/s      MLUP/s    ',        &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I8))',  advance='NO') na, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen

  if (parametrized) then
    write(*,'(3(2X,ES10.3))') error_par, mflops_par, mlups_par
  else
    write(*,'(3(8X,A))') 'None', 'None', 'None'
  end if
  write(*,*)

!===============================================================================

end program Validate__CART__TPO_Spectral_Iso
