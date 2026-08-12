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

!> summary:  Validation of the curvilinear TPO for 3D convection
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/01/13
!===============================================================================

program Validate__TPO_INS_Convection_D
  use Kind_Parameters
  use TPO__INS_Convection__3D_D__Gen
  use TPO__INS_Convection__3D_D__XSMM
  implicit none

  ! input parameters ...........................................................

  character(len=*), parameter :: &
         input_file = 'validate__tpo__ins_convection_3d_d.prm'

  integer :: po_v = 7       ! polynomial order of velocity
  integer :: po_q = 12      ! polynomial order used for quadrature
  integer :: ne   = 1024    ! number of elements
  integer :: nt   = 1       ! number of test loops ≥ 1

  namelist/input/ po_v, po_q, ne, nt

  ! operators and variables ....................................................

  character, allocatable, save :: bc    (:,:)         !< face BC
  real(RNP), allocatable, save :: D_v   (:,:)         !< 1D diff operator
  real(RNP), allocatable, save :: I_vq  (:,:)         !< 1D interpolation
  real(RNP), allocatable, save :: w_q   (:)           !< 1D quadrature weights
  real(RNP), allocatable, save :: Jd_q  (:,:,:,:)     !< Jacobian determinant
  real(RNP), allocatable, save :: Ji_q  (:,:,:,:,:,:) !< inverse Jacobian matrix
  real(RNP), allocatable, save :: a_q   (:,:,:,:)     !< face area coefficients
  real(RNP), allocatable, save :: n_q   (:,:,:,:,:)   !< face unit normal vectors
  real(RNP), allocatable, save :: v     (:,:,:,:,:)   !< velocity
  real(RNP), allocatable, save :: vp    (:,:,:,:,:)   !< velocity v⁺ at faces
  real(RNP), allocatable, save :: F_ref (:,:,:,:,:)   !< result: reference
  real(RNP), allocatable, save :: F_gen (:,:,:,:,:)   !< result: generic
  real(RNP), allocatable, save :: F_opt (:,:,:,:,:)   !< result: libxsmm

  real(RNP) :: alpha = 0.5
  real(RNP) :: time  = 0
  real(RNP) :: error_gen = 0, mflops_gen = -1, mlups_gen = -1
  real(RNP) :: error_opt = 0, mflops_opt = -1, mlups_opt = -1

  integer(IXL) :: count, count0, rate
  integer :: nv, nq, nvv, nvq, nqq, nvvv, nvvq, nvqq, nqqq
  integer :: nope, nflop, nflop_f
  integer :: i, io

  !-----------------------------------------------------------------------------
  ! initialization

  ! read test parameters .......................................................

  open(newunit = io, file = input_file)
  read(io, nml = input)
  close(io)

  ! workspace ..................................................................

  nv = po_v + 1
  nq = po_q + 1

  allocate( bc(6,ne), source = ' ')

  allocate( D_v   (nv,nv)           )
  allocate( I_vq  (nq,nv)           )
  allocate( w_q   (nq)              )
  allocate( Jd_q  (nq,nq,nq,ne)     )
  allocate( Ji_q  (nq,nq,nq,ne,3,3) )
  allocate( a_q   (nq,nq,6,ne)      )
  allocate( n_q   (nq,nq,6,ne,3)    )
  allocate( v     (nv,nv,nv,ne,3)   )
  allocate( vp    (nv,nv,6,ne,3)    )
  allocate( F_ref (nv,nv,nv,ne,3)   )
  allocate( F_gen (nv,nv,nv,ne,3)   )
  allocate( F_opt (nv,nv,nv,ne,3)   )

  call random_number( D_v   )
  call random_number( I_vq  )
  call random_number( w_q   )
  call random_number( Jd_q  )
  call random_number( Ji_q  )
  call random_number( a_q   )
  call random_number( n_q   )
  call random_number( v     )
  call random_number( vp    )

  if (nv == nq) then
    I_vq( 1, 1) = 1
    I_vq(nq,nv) = 1
  else
    I_vq( 1, 1) = 0.1
    I_vq(nq,nv) = 0.2
  end if

  F_ref = 0
  F_gen = 0
  F_opt = 0

  ! problem dimensions .........................................................

  ! number of operands per element
  nope = nv**3

  nvv  = nv**2
  nvq  = nv * nq
  nqq  = nq**2

  nvvv = nv**3
  nvvq = nv**2 * nq
  nvqq = nv * nq**2
  nqqq = nq**3

  ! number of FLOPs ............................................................

  ! v-interpolation
  if (nv == nq) then
    nflop = 0
  else
    nflop = 3 * (nvvq * 2 * nv + nvqq * 2 * nv + nqqq * 2 * nv)
  end if

  ! mass matrix
  nflop = nflop + 3 * nqqq

  ! MJ⁻¹⋅vv
  nflop = nflop + nqqq * (6 + 9 * 6)

  ! volume integral
  if (nv == nq) then
    nflop = nflop + 3 * (nvvv * 2 * nv) + 6 * (nvvv * (2 * nv + 1)) &
                  + 3 * (nvvv * 2 * nv)                             &
                  + 3 * (nvvv * 3) + nvvv * 18
  else
    nflop = nflop + 3 * (nvqq * 2 * nq + nvvq * 2 * nq + nvvv * 2 * nq)       &
                  + 6 * (nvqq * 2 * nq + nvvq * 2 * nq + nvvv * (2 * nq + 1)) &
                  + 3 * (nvvv * 4 * nv + nvqq * 6 * nv + nqqq * 6)            &
                  + 3 * (nvqq * 3 + nvvq * 3 + nvvv * 3) * nq + nqqq * 18
  end if

  ! face integrals
  nflop_f = 8 + nqq * (5 + 5 + 2 + 3 * 8)
  if (nv == nq) then
    nflop_f = nflop_f + 3 * nvv
  else
    nflop_f = nflop_f + 6 * (nvq * 2 * nv + nqq * 2 * nv)       &
                      + 3 * (nvq * 2 * nq + nvv * (2 * nq + 1))
  end if

  ! number of FLOPs per element
  nflop = nflop + 6 * nflop_f

  !-----------------------------------------------------------------------------
  ! Test generic procedure

  !$omp parallel

  call TPO_INS_Convection_D_Gen(nv, nq, ne, bc, alpha, &
           D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_ref)

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_INS_Convection_D_Gen(nv, nq, ne, bc, alpha, &
             D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_gen)
  end do

  call system_clock(count)

  !$omp end parallel

  if (nt > 0) then
    time = (count - count0) / real(rate, RNP) / max(nt, 1)
    error_gen  = maxval(abs(F_gen - F_ref))
    mflops_gen = 1E-6 / time * ne * nflop
    mlups_gen  = 1E-6 / time * ne * nope
  end if

  !-----------------------------------------------------------------------------
  ! Test optimized procedure

#ifdef __LIBXSMM__
  !$omp parallel

  call TPO_INS_Convection_D_XSMM(nv, nq, ne, bc, alpha, &
           D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_opt)

  call system_clock(count0, rate)

  do i = 1, nt
    call TPO_INS_Convection_D_XSMM(nv, nq, ne, bc, alpha, &
             D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_opt)
  end do

  call system_clock(count)

  !$omp end parallel

  if (nt > 0) then
    time = (count - count0) / real(rate, RNP) / nt
    error_opt  = maxval(abs(F_opt - F_ref))
    mflops_opt = 1E-6 / time * ne * nflop
    mlups_opt  = 1E-6 / time * ne * nope
  end if
#endif

  !-----------------------------------------------------------------------------
  ! Print results

  write(*,*)
  write(*,'(3A)') '                                ',      &
                  '   ------ generic operator ---------',  &
                  '   ------ optimized operator -------'
  write(*,'(3A)') '   np_v   np_q      ne        nt    ',  &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s'
  write(*,'(2(2X,I5),I8,2X,I8)',  advance='NO') nv, nq, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))', advance='NO') error_opt, mflops_opt, mlups_opt
  write(*,'(2/)')

  !=============================================================================

end program Validate__TPO_INS_Convection_D
