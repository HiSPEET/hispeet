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

!> summary:  Evaluation of INS boundary fluxes
!> author:   Joerg Stiller
!> date:     2026/08/20
!===============================================================================

module INS__Boundary_Fluxes__3D
  use Kind_Parameters
  use Constants

  use Boundary_Variable__3D
  use Surface_Integrals__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_BoundaryFluxes_3D

  !-----------------------------------------------------------------------------
  !> Incompressible Navier-Stokes boundary fluxes

  type INS_BoundaryFluxes_3D

    real(RNP)              :: t         !< time
    real(RNP), allocatable :: m(:)      !< mass               (n_nound)
    real(RNP), allocatable :: f_p(:,:)  !< pressure           (n_bound,3)
    real(RNP), allocatable :: f_d(:,:)  !< viscous diffusion  (n_bound,3)

    integer, private :: part = -1

  contains
    procedure :: Evaluate
    procedure :: PrintValues
  end type INS_BoundaryFluxes_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of INS boundary fluxes

  subroutine Evaluate(this, ins_op, t, mu, nu, u, leaf)
    class(INS_BoundaryFluxes_3D), intent(inout) :: this
    class(INS_Operator_3D), intent(in) :: ins_op
    real(RNP),              intent(in) :: t
    real(RNP), contiguous,  intent(in) :: mu(:,:,:,:)
    real(RNP), contiguous,  intent(in) :: nu(:,:,:,:)
    real(RNP), contiguous,  intent(in) :: u(:,:,:,:,:)
    logical,   optional,    intent(in) :: leaf   !< only leaf elements [F]

    !  internal variables ......................................................

    type(BoundaryVariable_3D), allocatable, save :: bv(:)
    ! boundary variables (1:n_bound)
    !   - var 1   : normal velocity, v⋅n
    !   - var 2:4 : pressure times normal, pn
    !   - var 5:7 : viscous stress vector, s = n⋅τ

    type(BoundaryVariable_3D), allocatable, save :: bv_vn(:) ! handle for v⋅n
    type(BoundaryVariable_3D), allocatable, save :: bv_pn(:) ! handle for pn
    type(BoundaryVariable_3D), allocatable, save :: bv_s (:) ! handle for s

    real(RNP), allocatable, save :: int_bv(:,:) ! boundary integrals

    integer :: n_bound, n_flux, po
    integer :: b, e, f, i, k

    ! initialization ...........................................................

    po      = ins_op % eop_u % po
    n_bound = ins_op % mesh % n_bound
    n_flux  = 7

    !$omp master
    allocate(bv(n_bound), bv_vn(n_bound), bv_pn(n_bound), bv_s(n_bound))
    do b = 1, n_bound
      call bv(b) % Init(ins_op % mesh % boundary(b), po, nc = n_flux)
      call bv(b) % GetSlice(first=1, last=1, slice=bv_vn(b))
      call bv(b) % GetSlice(first=2, last=4, slice=bv_pn(b))
      call bv(b) % GetSlice(first=5, last=7, slice=bv_s (b))
    end do
    allocate(int_bv(n_flux,n_bound), source = ZERO)
    !$omp end master
    !$omp barrier

    ! normal velocity ..........................................................

    associate(v => u(:,:,:,:,1:3))
      do b = 1, n_bound
        call bv_vn(b) % ExtractNormalComponent(ins_op % sem_u, v)
      end do
    end associate

    ! pressure x normal ........................................................

    associate( p => u(:,:,:,:,4)                    &
             , n => ins_op % sem_u % metrics % n    &
             , boundary => ins_op % mesh % boundary )

      k = size(p,1)

      !$omp do collapse(2)
      do b = 1, n_bound
      do f = 1, boundary(b) % n_face

        associate(pn => bv_pn(b) % val(:,:,:,1))

          e = boundary(b) % face(f) % element_id
          i = boundary(b) % face(f) % element_face
          select case(i)
          case(1)
            pn(:,:,1) = p(1,:,:,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(1,:,:,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(1,:,:,e) * n(:,:,i,e,3)
          case(2)
            pn(:,:,1) = p(k,:,:,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(k,:,:,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(k,:,:,e) * n(:,:,i,e,3)
          case(3)
            pn(:,:,1) = p(:,1,:,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(:,1,:,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(:,1,:,e) * n(:,:,i,e,3)
          case(4)
            pn(:,:,1) = p(:,k,:,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(:,k,:,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(:,k,:,e) * n(:,:,i,e,3)
          case(5)
            pn(:,:,1) = p(:,:,1,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(:,:,1,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(:,:,1,e) * n(:,:,i,e,3)
          case(6)
            pn(:,:,1) = p(:,:,k,e) * n(:,:,i,e,1)
            pn(:,:,2) = p(:,:,k,e) * n(:,:,i,e,2)
            pn(:,:,3) = p(:,:,k,e) * n(:,:,i,e,3)
          end select

        end associate
      end do
      end do
    end associate

    ! viscous stress vector ....................................................

    do b = 1, n_bound
      call ins_op % GetViscousBoundaryStress( b, mu, nu, u       &
                                            , sb = bv_s(b) % val &
                                            , xout = .true.      )
    end do

    ! integrals ................................................................

    call GetSurfaceIntegrals(ins_op % sem_u, bv, int_bv, leaf)

    ! assignment and finalization ..............................................

    !$omp master
    this % t    = t
    this % part = ins_op % mesh % part
    this % m    = int_bv( 1 ,:)
    this % f_p  = int_bv(2:4,:)
    this % f_d  = int_bv(5:7,:)
    deallocate(bv, bv_vn, bv_pn, bv_s, int_bv)
    !$omp end master

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print INS boundary fluxes

  subroutine PrintValues(this)
    class(INS_BoundaryFluxes_3D), intent(in) :: this

    integer :: b

    !$omp master
    if (this % part == 0) then

      write(*,*)
      write(*,'(A,ES12.5)') 'boundary fluxes at t =', this % t
      write(*,'(A)',advance='NO') '#'
      write(*,'(2X,A1,2X)',advance='NO') 'b'
      write(*,'(1X,A4,8X)',advance='NO') 'm'
      write(*,'(2X,A4,6X)',advance='NO') 'f_p1'
      write(*,'(2X,A4,6X)',advance='NO') 'f_p2'
      write(*,'(2X,A4,6X)',advance='NO') 'f_p3'
      write(*,'(2X,A4,6X)',advance='NO') 'f_d1'
      write(*,'(2X,A4,6X)',advance='NO') 'f_d2'
      write(*,'(2X,A4,6X)'             ) 'f_d3'

      do b = 1, size(this % m)
        write(*,'(I4,2X)',advance='NO') b
        write(*,'(ES10.3,2X)',advance='NO') this % m(b)
        write(*,'(ES10.3,2X)',advance='NO') this % f_p(b,1)
        write(*,'(ES10.3,2X)',advance='NO') this % f_p(b,2)
        write(*,'(ES10.3,2X)',advance='NO') this % f_p(b,3)
        write(*,'(ES10.3,2X)',advance='NO') this % f_d(b,1)
        write(*,'(ES10.3,2X)',advance='NO') this % f_d(b,2)
        write(*,'(ES10.3,2X)'             ) this % f_d(b,3)
      end do

    end if
    !$omp end master

  end subroutine PrintValues

  !=============================================================================

end module INS__Boundary_Fluxes__3D
