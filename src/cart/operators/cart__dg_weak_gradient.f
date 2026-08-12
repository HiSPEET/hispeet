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

!> summary:  Weak gradient of a vector field
!> author:   Joerg Stiller
!> date:     2018/03/28
!===============================================================================

module CART__DG_Weak_Gradient

  use Kind_Parameters, only: RNP
  use Standard_Element_Operators__1D
  use TPO__Grad__3D_R
  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Trace_Operator

  implicit none
  private

  public :: WeakGradient

  interface WeakGradient
    module procedure WeakGradient_B
    module procedure WeakGradient_O
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Weak gradient of a scalar field on discontinuous elements - with BC

  subroutine WeakGradient_B(mesh, Ms, Ds, u, bv_u, grad_u)
    class(MeshPartition),    intent(in)  :: mesh              !< mesh partition
    real(RNP),               intent(in)  :: Ms(:)             !< std mass matrix
    real(RNP),               intent(in)  :: Ds(:,:)           !< std diff matrix
    real(RNP),               intent(in)  :: u(:,:,:,:)        !< scalar field
    class(BoundaryVariable), intent(in)  :: bv_u(:)           !< BC
    real(RNP),               intent(out) :: grad_u(:,:,:,:,:) !< weak gradient

    call WeakGradient_X(mesh, size(Ms)-1, size(u,4), Ms, Ds, u, bv_u, grad_u)

  end subroutine WeakGradient_B

  !-----------------------------------------------------------------------------
  !> Weak gradient of a scalar field on discontinuous elements - open (no BC)

  subroutine WeakGradient_O(mesh, Ms, Ds, u, grad_u)
    class(MeshPartition), intent(in)  :: mesh              !< mesh partition
    real(RNP),            intent(in)  :: Ms(:)             !< std mass matrix
    real(RNP),            intent(in)  :: Ds(:,:)           !< std diff matrix
    real(RNP),            intent(in)  :: u(:,:,:,:)        !< scalar field
    real(RNP),            intent(out) :: grad_u(:,:,:,:,:) !< weak gradient

    call WeakGradient_X(mesh, size(Ms)-1, size(u,4), Ms, Ds, u, grad_u=grad_u)

  end subroutine WeakGradient_O

  !-----------------------------------------------------------------------------
  !> Weak gradient of a scalar field on discontinuous elements - eXplicit

  subroutine WeakGradient_X(mesh, po, ne, Ms, Ds, u, bv_u, grad_u)

    ! arguments ................................................................

    !> mesh partition
    class(MeshPartition), intent(in) :: mesh
    !> polynomial order
    integer,   intent(in) :: po
    !> number of elements
    integer,   intent(in) :: ne
    !> std mass matrix
    real(RNP), intent(in) :: Ms(0:po)
    !> std diff matrix
    real(RNP), intent(in) :: Ds(0:po, 0:po)
    !> vector field
    real(RNP), intent(in) :: u(0:po, 0:po, 0:po, ne)
    !> boundary conditions for u
    class(BoundaryVariable), optional, intent(in) :: bv_u(mesh%n_boundary)
    !> weak gradient of u
    real(RNP), intent(out) :: grad_u(0:po, 0:po, 0:po, ne, 3)

    ! internal variables .......................................................

    ! shared with OpenMP
    type(TraceOperator), asynchronous, allocatable, save :: trace_op
    real(RNP), allocatable, save :: tr_u(:,:,:,:)

    real(RNP), allocatable :: Dm(:,:)
    real(RNP) :: g
    integer :: tag = 1000
    integer :: c, e, f, f1, f2, i, j, k

    ! initialization ...........................................................

    !$omp single
    allocate(trace_op)
    allocate(tr_u(0:po, 0:po, 2, mesh%nf))
    !$omp end single

    allocate(Dm(0:po,0:po))
    do j = 0, po
    do i = 0, po
      Dm(j,i) = -Ms(i) * Ds(i,j) / Ms(j)
    end do
    end do

    ! extract and start transferring traces ....................................

    if (present(bv_u)) then
      call trace_op % GetTrace_Start(mesh, u, bv_u, tr_u, tag)
    else
      call trace_op % GetTrace_Start(mesh, u, tr_u, tag)
    end if

    ! compute transposed gradient ..............................................

    call TPO_Grad_R(Ms, Dm, mesh%dx, u, v=grad_u)

    ! complete transfer ........................................................

    call trace_op % GetTrace_Finish(mesh, tr_u)

    associate(dx => mesh%dx, face => mesh%face)

      ! x1-flux contribution ...................................................

      c  = 1                      ! vector component
      g  = 1 / (dx(c) * Ms(0))    ! metric and averaging factor
      f1 = 1                      ! first face
      f2 = mesh % nf1             ! last

      !$omp do
      do f = f1, f2

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne ) then
          i = po
          do k = 0, po
          do j = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              + g * (tr_u(j,k,1,f) + tr_u(j,k,2,f))
          end do
          end do
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne ) then
          i = 0
          do k = 0, po
          do j = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              - g * (tr_u(j,k,1,f) + tr_u(j,k,2,f))
          end do
          end do
        end if

      end do

      ! x2-flux contribution ...................................................

      c  = 2                      ! vector component
      g  = 1 / (dx(c) * Ms(0))    ! metric and averaging factor
      f1 = f2 + 1                 ! first face
      f2 = f2 + mesh % nf2        ! last

      !$omp do
      do f = f1, f2

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne ) then
          j = po
          do k = 0, po
          do i = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              + g * (tr_u(i,k,1,f) + tr_u(i,k,2,f))
          end do
          end do
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne ) then
          j = 0
          do k = 0, po
          do i = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              - g * (tr_u(i,k,1,f) + tr_u(i,k,2,f))
          end do
          end do
        end if

      end do

      ! x3-flux contribution ...................................................

      c  = 3                      ! vector component
      g  = 1 / (dx(c) * Ms(0))    ! metric and averaging factor
      f1 = f2 + 1                 ! first face
      f2 = f2 + mesh % nf3        ! last

      !$omp do
      do f = f1, f2

        e = face(f) % element(1)
        if (0 < e .and. e <= mesh % ne ) then
          k = po
          do j = 0, po
          do i = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              + g * (tr_u(i,j,1,f) + tr_u(i,j,2,f))
          end do
          end do
        end if

        e = face(f) % element(2)
        if (0 < e .and. e <= mesh % ne ) then
          k = 0
          do j = 0, po
          do i = 0, po
            grad_u(i,j,k,e,c) = grad_u(i,j,k,e,c) &
                              - g * (tr_u(i,j,1,f) + tr_u(i,j,2,f))
          end do
          end do
        end if

      end do

    end associate

    ! finalization .............................................................

    !$omp single
    deallocate(tr_u)
    deallocate(trace_op)
    !$omp end single

  end subroutine WeakGradient_X

  !=============================================================================

end module CART__DG_Weak_Gradient
