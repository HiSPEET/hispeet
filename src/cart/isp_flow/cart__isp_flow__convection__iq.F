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

!> summary:  ISP flow: convective fluxes with independent Lobatto quadrature (Q≠P)
!> author:   Joerg Stiller
!> date:     2018/06/12
!===============================================================================

module CART__ISP_Flow__Convection__IQ

  use Kind_Parameters, only: RNP
  use TPO__AAA__3D
  use CART__ISP_Flow__Operators
  use CART__ISP_Flow__Convection__EQ

  implicit none
  private

  public :: WeakConvectiveFlux_IQ

contains

  !-----------------------------------------------------------------------------
  !> Weak divergence of convective fluxes for incompressible flow -- independent

  subroutine WeakConvectiveFlux_IQ(flow_op, u, div_f)

    ! arguments ................................................................

    class(FlowOperators), intent(in)  :: flow_op    !< flow operators
    real(RNP),            intent(in)  :: u          !< flow variables
    real(RNP),            intent(out) :: div_f      !< flux divergence

    dimension :: u     (:,:,:,:,:)
    dimension :: div_f (:,:,:,:,:)

    ! local variables ..........................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: uq, div_fq

    real(RNP), allocatable :: A(:,:), J(:,:)
    integer :: nc, ne, nq, pq, pu
    integer :: i, k

    ! initialization ...........................................................

    ne = size(u,4)
    nc = size(u,5)

    pu = flow_op % po_u
    pq = flow_op % po_q
    nq = pq + 1

    !$omp single
    allocate(uq(nq,nq,nq,ne,nc))
    allocate(div_fq, mold=uq)
    !$omp end single

    ! projection operator -- with diagonal (lumped) mass matrix
    allocate(A(0:pu,0:pq))
    allocate(J(0:pq,0:pu), source = flow_op % iop_uq % A)
    do k = 0, pq
    do i = 0, pu
      A(i,k) = flow_op % eop_q % w(k) * J(k,i) / flow_op % eop_u % w(i)
    end do
    end do

    ! computation ..............................................................

    ! interpolation to quadrature points
    call flow_op % iop_uq % Apply(u, uq)

    ! flux divergence at quadrature points
    call WeakConvectiveFlux_EQ( flow_op % mesh  &
                              , flow_op % eop_q &
                              , uq              &
                              , div_fq          )

    ! projection
    do i = 1, nc
      call TPO_AAA(A, div_fq(:,:,:,:,i), div_f(:,:,:,:,i))
    end do

    ! clean-up .................................................................

    !$omp barrier
    !$omp master
    deallocate(uq, div_fq)
    !$omp end master

  end subroutine WeakConvectiveFlux_IQ

  !=============================================================================

end module CART__ISP_Flow__Convection__IQ
