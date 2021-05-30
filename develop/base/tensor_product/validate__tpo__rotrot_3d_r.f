!> summary:  Validation of the tensor-product rotrotrotrot operator
!> author:   Jörg Stiller, Erik Pfister
!> date:     2020/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__RotRot_3d_R
  use Kind_Parameters, only: IXL, RNP
  use Standard_Operators__1D
  use TPO__RotRot_3d_R
  use TPO__RotRot_3d_R__Gen
  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! test parameters ............................................................

  integer :: po = 4   ! polynomial order of elements
  integer :: ne = 1   ! number of elements
  integer :: nt = 1   ! number of test runs

  namelist /input/ po, ne, nt

  ! operators and variables ....................................................

  type(StandardOperators_1D) :: standard_op

  real(RNP), dimension(:,:,:,:,:), allocatable :: u
  real(RNP), dimension(:,:,:,:,:), allocatable :: v, w
  real(RNP), dimension(:),         allocatable :: x, y, z

  real(RNP) :: x0, y0, z0
  real(RNP) :: dx(3) = 2
  real(RNP) :: time
  real(RNP) :: error_gen, mflops_gen, mlups_gen
  real(RNP) :: error_opt, mflops_opt, mlups_opt

  logical :: exists
  integer :: np, nflop, npop, prm
  integer :: p, pm1, pm2, pm3
  integer :: i, j, k, e

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  inquire(file='validate__tpo__rotrot_3d_r.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='validate__tpo__rotrot_3d_r.prm')
    read(prm, nml=input)
    close(prm)
  end if

  ! operator dimension
  np = po + 1

  ! problem dimensions
  nflop = np**3 * (30*np + 35)
  npop  = np**3

  ! operators ..................................................................

  standard_op = StandardOperators_1D(po)

  ! workspace ..................................................................

  allocate( u(0:po,0:po,0:po,ne,3), &
            v(0:po,0:po,0:po,ne,3), &
            w(0:po,0:po,0:po,ne,3)  )

  allocate( x(0:po), y(0:po), z(0:po) )

  ! order of test function .....................................................

  p = min(po, 5)

  pm1 = max(p - 1, 0)
  pm2 = max(p - 2, 0)
  pm3 = max(p - 3, 0)

  !-----------------------------------------------------------------------------
  ! operand und exact result

  associate( xs => standard_op % x )

    do e = 1, ne

      ! element points .........................................................

      call random_number(x0)
      call random_number(y0)
      call random_number(z0)

      x = x0 + xs
      y = y0 + xs
      z = z0 + xs

      ! operand ................................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po

        u(i,j,k,e,1) =  x(i) ** p  *  y(j) ** pm1
        u(i,j,k,e,2) =  y(j) ** p  *  z(k) ** pm1
        u(i,j,k,e,3) =  z(k) ** p  *  x(i) ** pm1

      end do
      end do
      end do

      ! exact result ...........................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po

        w(i,j,k,e,1) =  pm1 * ( p    *  z(k) ** pm1  *  x(i) ** pm2  &
                              - pm2  *  x(i) ** p    *  y(j) ** pm3  )
        w(i,j,k,e,2) =  pm1 * ( p    *  x(i) ** pm1  *  y(j) ** pm2  &
                              - pm2  *  y(j) ** p    *  z(k) ** pm3  )
        w(i,j,k,e,3) =  pm1 * ( p    *  y(j) ** pm1  *  z(k) ** pm2  &
                              - pm2  *  z(k) ** p    *  x(i) ** pm3  )

      end do
      end do
      end do

    end do

  end associate

  !-----------------------------------------------------------------------------
  ! test generic procedure

  associate( Ds => standard_op % D )

    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call TPO_RotRot_R_Gen(np, ne, Ds, dx, u, v)
    !$acc wait

    call system_clock(count0, rate)

    do i = 1, nt
    call TPO_RotRot_R_Gen(np, ne, Ds, dx, u, v)
    !$acc wait
    end do

    call system_clock(count)
    !$acc end data
    !$omp end parallel

  end associate

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(v - w))
  mflops_gen = 1E-6 / time * ne * nflop
  mlups_gen  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! test optimized procedure

  associate( Ds => standard_op % D )

    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call TPO_RotRot_R(Ds, dx, u, v)
    !$acc wait

    call system_clock(count0, rate)
    do i = 1, nt
      call TPO_RotRot_R(Ds, dx, u, v)
      !$acc wait
    end do

    call system_clock(count)
    !$acc end data
    !$omp end parallel

  end associate

  time = (count - count0) / real(rate, RNP) / nt

  error_opt  = maxval(abs(v - w))
  mflops_opt = 1E-6 / time * ne * nflop
  mlups_opt  = 1E-6 / time * ne * npop


  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                        ',             &
                  '   ------------ generic ------------',  &
                  '   ----------- optimized  ----------'
  write(*,'(3A)') '#  np        ne        nt    ',         &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I8))',  advance='NO') np, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

!===============================================================================

end program Validate__TPO__RotRot_3d_R
