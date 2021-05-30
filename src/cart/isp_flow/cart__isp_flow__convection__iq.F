!> summary:  ISP flow: convective fluxes with independent GLL quadrature (Q≠P)
!> author:   Joerg Stiller
!> date:     2018/06/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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

    real(RNP), allocatable :: A(:,:)
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
    associate( wu => flow_op % eop_u  % w   &
             , wq => flow_op % eop_q  % w   &
             , J  => flow_op % iop_uq % A )
      do k = 0, pq
      do i = 0, pu
        A(i,k) = wq(k) * J(k,i) / wu(i)
      end do
      end do
    end associate

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
