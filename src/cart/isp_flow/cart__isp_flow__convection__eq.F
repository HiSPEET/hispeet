!> summary:  ISP flow: convective fluxes using embedded Lobatto quadrature (Q=P)
!> author:   Joerg Stiller
!> date:     2018/04/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Convection__EQ

  use Kind_Parameters, only: RNP
  use Constants,       only: HALF
  use Standard_Operators__1D
  use TPO__Div__3D_R
  use CART__Mesh_Partition
  use CART__Trace_Operator

  implicit none
  private

  public :: WeakConvectiveFlux_EQ

contains

  !-----------------------------------------------------------------------------
  !> Weak divergence of convective fluxes for incompressible flow -- collocation

  subroutine WeakConvectiveFlux_EQ(mesh, sop, u, div_f)

    ! arguments ................................................................

    class(MeshPartition),        intent(in)  :: mesh    !< mesh partition
    class(StandardOperators_1D),  intent(in)  :: sop     !< standard operators
    real(RNP),                   intent(in)  :: u       !< flow variables
    real(RNP),                   intent(out) :: div_f   !< flux divergence

    dimension :: u     (0:, 0:, 0:, :, :)
    dimension :: div_f (0:, 0:, 0:, :, :)

    ! internal variables .......................................................

    ! shared with OpenMP
    type(TraceOperator), asynchronous, allocatable, save :: trace_op
    real(RNP), allocatable, save :: tr_u(:,:,:,:,:)
    real(RNP), allocatable, save :: fc(:,:,:,:,:)

    integer :: tag = 1000
    integer :: po, nc
    integer :: c, e, f, f1, f2, s
    integer :: i, j, k

    real(RNP) :: g(2), vn(2)
    real(RNP), allocatable :: hn(:,:)

    ! initialization ...........................................................

    po = sop % po
    nc = size(u,5)

    !$omp single
    allocate(trace_op)
    allocate(tr_u ( 0:po, 0:po, 2   , mesh%nf, nc ))
    allocate(fc   ( 0:po, 0:po, 0:po, mesh%ne, 3  ))
    !$omp end single

    allocate(hn(0:po,0:po))

    ! extract and transfer traces ..............................................

    call trace_op % GetTrace_Start(mesh, u, tr_u, tag)

    ! elementwise flux divergence ..............................................

    do c = 1, nc

      if (c == 4) cycle  ! skip pressure

      !$omp do
      do e = 1, mesh%ne
        do k = 0, po
        do j = 0, po
        do i = 0, po
          fc(i,j,k,e,1) = u(i,j,k,e,1) * u(i,j,k,e,c)
          fc(i,j,k,e,2) = u(i,j,k,e,2) * u(i,j,k,e,c)
          fc(i,j,k,e,3) = u(i,j,k,e,3) * u(i,j,k,e,c)
        end do
        end do
        end do
      end do

      call TPO_Div(Ds=sop%D, dx=mesh%dx, u=fc, v=div_f(:,:,:,:,c))

    end do

    ! complete transfer ........................................................

    call trace_op % GetTrace_Finish(mesh, tr_u)

    ! x1-element-boundaries ....................................................

    g(1) =  2 / (mesh%dx(1) * sop%w(0))   ! =  2 / dx Ms(0)
    g(2) = -2 / (mesh%dx(1) * sop%w(po))  ! = -2 / dx Ms(P)

    f1 = 1
    f2 = mesh % nf1

    do c = 1, nc
      if (c == 4) cycle

      !$omp do private(hn,vn)
      do f = f1, f2

        do k = 0, po
        do j = 0, po
          vn = tr_u(j,k,:,f,1)
          hn(j,k) = HALF * ( vn(1) * tr_u(j,k,1,f,c)                   &
                           + vn(2) * tr_u(j,k,2,f,c)                   &
                           )                                           &
                  + max(abs(vn(1)), abs(vn(2))) * ( tr_u(j,k,1,f,c)    &
                                                  - tr_u(j,k,2,f,c)    &
                                                  )
        end do
        end do

        ! face sides
        do s = 1, 2
          e = mesh % face(f) % element(s)
          if (0 < e .and. e <= mesh % ne ) then
            i = (2 - s) * po  !  =  [0,po]  for  s = [1,2]
            do k = 0, po
            do j = 0, po
              div_f(i,j,k,e,c) = div_f(i,j,k,e,c)  &
                               + g(s) * ( hn(j,k)  &
                                        - tr_u(j,k,s,f,1) * tr_u(j,k,s,f,c) )
            end do
            end do
          end if
        end do

      end do
    end do

    ! x2-element-boundaries ....................................................

    g(1) =  2 / (mesh%dx(2) * sop%w(0))   ! =  2 / dy Ms(0)
    g(2) = -2 / (mesh%dx(2) * sop%w(po))  ! = -2 / dy Ms(P)

    f1 = f2 + 1
    f2 = f2 + mesh % nf2

    do c = 1, nc
      if (c == 4) cycle

      !$omp do private(hn,vn)
      do f = f1, f2

        do k = 0, po
        do i = 0, po
          vn = tr_u(i,k,:,f,2)
          hn(i,k) = HALF * ( vn(1) * tr_u(i,k,1,f,c)                   &
                           + vn(2) * tr_u(i,k,2,f,c)                   &
                           )                                           &
                  + max(abs(vn(1)), abs(vn(2))) * ( tr_u(i,k,1,f,c)    &
                                                  - tr_u(i,k,2,f,c)    &
                                                  )
        end do
        end do

        ! face sides
        do s = 1, 2
          e = mesh % face(f) % element(s)
          if (0 < e .and. e <= mesh % ne ) then
            j = (2 - s) * po
            do k = 0, po
            do i = 0, po
              div_f(i,j,k,e,c) = div_f(i,j,k,e,c)  &
                               + g(s) * ( hn(i,k)  &
                                        - tr_u(i,k,s,f,2) * tr_u(i,k,s,f,c) )
            end do
            end do
          end if
        end do

      end do
    end do

    ! x3-element-boundaries ....................................................

    g(1) =  2 / (mesh%dx(3) * sop%w(0))   ! =  2 / dz Ms(0)
    g(2) = -2 / (mesh%dx(3) * sop%w(po))  ! = -2 / dz Ms(P)

    f1 = f2 + 1
    f2 = f2 + mesh % nf3

    do c = 1, nc
      if (c == 4) cycle

      !$omp do private(hn,vn)
      do f = f1, f2

        do j = 0, po
        do i = 0, po
          vn = tr_u(i,j,:,f,3)
          hn(i,j) = HALF * ( vn(1) * tr_u(i,j,1,f,c)                   &
                           + vn(2) * tr_u(i,j,2,f,c)                   &
                           )                                           &
                  + max(abs(vn(1)), abs(vn(2))) * ( tr_u(i,j,1,f,c)    &
                                                  - tr_u(i,j,2,f,c)    &
                                                  )
        end do
        end do

        ! face sides
        do s = 1, 2
          e = mesh % face(f) % element(s)
          if (0 < e .and. e <= mesh % ne ) then
            k = (2 - s) * po
            do j = 0, po
            do i = 0, po
              div_f(i,j,k,e,c) = div_f(i,j,k,e,c)  &
                               + g(s) * ( hn(i,j)  &
                                        - tr_u(i,j,s,f,3) * tr_u(i,j,s,f,c) )
            end do
            end do
          end if
        end do

      end do
    end do

    ! finalization .............................................................

    !$omp barrier
    !$omp master
    deallocate(fc, tr_u)
    deallocate(trace_op)
    !$omp end master

  end subroutine WeakConvectiveFlux_EQ

  !=============================================================================

end module CART__ISP_Flow__Convection__EQ
