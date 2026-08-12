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

!> summary:  Boundary values for incompressible Navier-Stokes pressure solver
!> author:   Joerg Stiller
!> date:     2026/06/26
!===============================================================================

submodule(INS__Operator__3D) MP_GetPressureBoundaryValues
  use Mesh_Boundary__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Provision of pressure boundary values

  module subroutine GetPressureBoundaryValues(this, tau, v, bv_u, bv_p, bv_q)
    class(INS_Operator_3D), intent(in) :: this
    !< time integration method
    real(RNP), intent(in) :: tau
    !< τ, effective time step width
    real(RNP), contiguous, intent(in) :: v(:,:,:,:,:)
    !< preliminary velocity
    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values of flow variables
    class(BoundaryVariable_3D), intent(inout) :: bv_p(:)
    !< pressure boundary conditions on v-points
    class(BoundaryVariable_3D), optional, intent(inout) :: bv_q(:)
    !< pressure boundary conditions on v-points

    integer   :: b
    real(RNP) :: ct

    ct = 1 / tau

    do b = 1, this % mesh % n_bound

      select case(this % problem % bc_p(b))
      case('N')
        call BuildNeumannBV( boundary = this % mesh % boundary(b)  &
                           , n        = this % sem_u % metrics % n &
                           , ct       = ct                         &
                           , v        = v                          &
                           , vb       = bv_u(b) % val(:,:,:,1:3)   &
                           , dn_p     = bv_p(b) % val(:,:,:,1)     )
      case('D')
        call SetArray(bv_p(b) % val(:,:,:,1), bv_u(b) % val(:,:,:,4))
      end select

      if (present(bv_q)) then
        call InterpolateFaceData( A = this % iop_up % A      &
                                , p = bv_p(b) % val(:,:,:,1) &
                                , q = bv_q(b) % val(:,:,:,1) )
      end if

    end do

  end subroutine GetPressureBoundaryValues

  !-----------------------------------------------------------------------------
  !> Build pressure boundary values from normal velocity conditions

  subroutine BuildNeumannBV(boundary, ct, n, v, vb, dn_p)
    class(MeshBoundary_3D), intent(in)  :: boundary
    real(RNP),              intent(in)  :: ct
    real(RNP), contiguous,  intent(in)  :: n (0:,0:,:,:,:)
    real(RNP), contiguous,  intent(in)  :: v (0:,0:,0:,:,:)
    real(RNP), contiguous,  intent(in)  :: vb(0:,0:,:,:)
    real(RNP), contiguous,  intent(out) :: dn_p(0:,0:,:)

    integer :: e, f, i, j, k, m, po

    po = ubound(v,1)

    !$omp do
    do f = 1, boundary % n_face

      e = boundary % face(f) % element_id
      m = boundary % face(f) % element_face

      select case(m)

      case(1,2)
        i = (m - 1) * po
        do k = 0, po
        do j = 0, po
          dn_p(j,k,f) = ct * ( n(j,k,m,e,1) * (v(i,j,k,e,1) - vb(j,k,f,1)) &
                             + n(j,k,m,e,2) * (v(i,j,k,e,2) - vb(j,k,f,2)) &
                             + n(j,k,m,e,3) * (v(i,j,k,e,3) - vb(j,k,f,3)) )
        end do
        end do

      case(3,4)
        j = (m - 3) * po
        do k = 0, po
        do i = 0, po
          dn_p(i,k,f) = ct * ( n(i,k,m,e,1) * (v(i,j,k,e,1) - vb(i,k,f,1)) &
                             + n(i,k,m,e,2) * (v(i,j,k,e,2) - vb(i,k,f,2)) &
                             + n(i,k,m,e,3) * (v(i,j,k,e,3) - vb(i,k,f,3)) )
        end do
        end do

      case(5,6)
        k = (m - 5) * po
        do j = 0, po
        do i = 0, po
          dn_p(i,j,f) = ct * ( n(i,j,m,e,1) * (v(i,j,k,e,1) - vb(i,j,f,1)) &
                             + n(i,j,m,e,2) * (v(i,j,k,e,2) - vb(i,j,f,2)) &
                             + n(i,j,m,e,3) * (v(i,j,k,e,3) - vb(i,j,f,3)) )
        end do
        end do

      end select

    end do

  end subroutine BuildNeumannBV

  !-----------------------------------------------------------------------------
  !> Interpolation of face data

  subroutine InterpolateFaceData(A, p, q)
    real(RNP), intent(in)  :: A(0:,0:)   !< 1D interpolation operator
    real(RNP), intent(in)  :: p(0:,0:,:) !< variable in velocity space
    real(RNP), intent(out) :: q(0:,0:,:) !< variable in pressure space

    real(RNP), allocatable :: w(:,:)
    real(RNP) :: tmp
    integer   :: po, pq, nf
    integer   :: f, i, j, k

    po = ubound(p,1)
    pq = ubound(q,1)
    nf = ubound(q,3)

    allocate(w(0:pq,0:po))

    !$omp do
    do f = 1, nf

      ! direction 1
      do j = 0, po
      do i = 0, pq
        tmp = 0
        do k = 0, po
          tmp = tmp + A(i,k) * p(k,j,f)
        end do
        w(i,j) = tmp
      end do
      end do

      ! direction 2
      do j = 0, pq
      do i = 0, pq
        tmp = 0
        do k = 0, po
          tmp = tmp + A(j,k) * w(i,k)
        end do
        q(i,j,f) = tmp
      end do
      end do

    end do

  end subroutine InterpolateFaceData

  !=============================================================================

end submodule MP_GetPressureBoundaryValues
