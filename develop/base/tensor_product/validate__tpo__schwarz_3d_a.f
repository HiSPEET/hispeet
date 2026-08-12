!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Validation of 3d constant anisotropic Schwarz TPO
!> author:   Joerg Stiller
!> date:     2020/06/03
!===============================================================================

program Validate__TPO__Schwarz__3D_A
  use Kind_Parameters,   only: IXL, RNP
  use Constants,         only: ZERO
  use Array_Assignments, only: SetArray

  use TPO__Schwarz__3D_A
  use TPO__Schwarz__3D_A__Gen

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

  ! operators, parameters and operands .........................................

  real(RNP) :: lambda = 1
  real(RNP), allocatable :: S1(:,:,:), V1(:,:), W1(:,:)
  real(RNP), allocatable :: S2(:,:,:), V2(:,:), W2(:,:)
  real(RNP), allocatable :: S3(:,:,:), V3(:,:), W3(:,:)
  real(RNP), allocatable :: g(:,:)
  integer,   allocatable :: cfg(:,:)

  real(RNP), allocatable :: f(:,:,:,:), u(:,:,:,:), r(:,:,:,:), nu(:)

  ! auxiliary ..................................................................

  character(len=80) :: input_file = 'validate__tpo__schwarz_3d_a.prm'

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

  allocate(S1(n1,n1,nc), V1(n1,nc), W1(n1,nc))
  allocate(S2(n2,n2,nc), V2(n2,nc), W2(n2,nc))
  allocate(S3(n3,n3,nc), V3(n3,nc), W3(n3,nc))
  allocate(cfg(3,nd), g(4,nd), nu(nd))
  allocate(f(n1,n2,n3,nd))
  allocate(u, r, mold = f)

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


  call random_number(S1)
  call random_number(V1)
  call random_number(W1)

  call random_number(S2)
  call random_number(V2)
  call random_number(W2)

  call random_number(S3)
  call random_number(V3)
  call random_number(W3)

  call random_number(g)
  call random_number(f)
  call random_number(nu)

  ! enforce positivity
  V1 = V1 + 1
  V2 = V2 + 1
  V3 = V3 + 1
  g  = g  + 1
  nu = nu + 1

  ! problem dimensions .........................................................

  nop   = n1 * n2 * n3
  nflop = nop * (4 * (n1 + n2 + n3) + 7) + 4

  !-----------------------------------------------------------------------------
  ! test of generic implementation

  !$omp parallel

  ! r = reference result
  call TPO_Schwarz_A_Gen(S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                         lambda, nu, f, r)

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz_A_Gen(S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, &
                           lambda, nu, f, u)
  end do
  call system_clock(count)

  !$omp end parallel

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(u - r))
  mflops_gen = 1E-6 / time * nd * nflop
  mlups_gen  = 1E-6 / time * nd * nop

  !-----------------------------------------------------------------------------
  ! test of optimized implementation

  call random_number(u)

  !$omp parallel

  call TPO_Schwarz(S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, lambda, nu, f, u)

  call system_clock(count0, rate)
  do i = 1, nt
    call TPO_Schwarz(S1, S2, S3, V1, V2, V3, W1, W2, W3, g, cfg, lambda, nu, f, u)
  end do
  call system_clock(count)

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

end program Validate__TPO__Schwarz__3D_A
