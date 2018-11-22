!> summary:  Application of the IP/DG elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Application of the IP/DG elliptic operator
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Apply
  use Constants, only: ZERO, ONE, HALF
  use Array_Assignments
  implicit none

  !-----------------------------------------------------------------------------
  ! private variables

  !> mapping of boundary face orientation to inner trace side
  integer, parameter :: inner_side(-3:3) = [ 2, 2, 2, 0, 1, 1, 1 ]

contains

!===============================================================================
! Common procedures

!-------------------------------------------------------------------------------
!> Elementwise computation of derivatives parallel to face normals
!>
!> Computes the normal components of grad(u) for all element boundary points.
!> Entries corresponding to interior points or tangential components are set
!> to zero.

subroutine NormalDerivatives(np, ne, Ds, dx, u, v)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RNP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RNP), intent(in)  :: dx(3)            !< element extensions
  real(RNP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RNP), intent(out) :: v(np,np,np,ne,3) !< element-wise gradient of u

  real(RNP), allocatable :: D1(:), D2(:)
  real(RNP) :: g(3), tmp1, tmp2
  integer   :: e, i, j, k, m
  integer   :: vec_len

  ! initialization .............................................................

  ! OpenACC vector length
  if (np < 8) then
    vec_len = 128
  else
    vec_len = 256
  end if

  ! transposed diff operators for first and last point
  allocate(D1, source=Ds( 1,:))
  allocate(D2, source=Ds(np,:))

  ! metric coefficients
  g = 2 / dx

  ! result
  call AssignScalar(v, ZERO, multi=.true.)

  !$acc data present(u,v) copyin(D1,D2,g) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker

  !$omp do private(e)
  do e = 1, ne

    ! v1 = du/dx1 @ face 1,2 ...................................................

    !$acc loop collapse(2) vector
    do k = 1, np
    do j = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + D1(m) * u(m,j,k,e)
        tmp2 = tmp2 + D2(m) * u(m,j,k,e)
      end do
      v( 1,j,k,e,1) = g(1) * tmp1
      v(np,j,k,e,1) = g(1) * tmp2

    end do
    end do

    ! v2 = du/dx2 @ face 3,4 ...................................................

    !$acc loop collapse(2) vector
    do k = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + D1(m) * u(i,m,k,e)
        tmp2 = tmp2 + D2(m) * u(i,m,k,e)
      end do
      v(i, 1,k,e,2) = g(2) * tmp1
      v(i,np,k,e,2) = g(2) * tmp2

    end do
    end do

    ! v3 = du/dx3 @ face 5,6 ...................................................

    !$acc loop collapse(2) vector
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + D1(m) * u(i,j,m,e)
        tmp2 = tmp2 + D2(m) * u(i,j,m,e)
      end do
      v(i,j, 1,e,3) = g(3) * tmp1
      v(i,j,np,e,3) = g(3) * tmp2

    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine NormalDerivatives

!-------------------------------------------------------------------------------
!> Modify boundary traces to yield correct contribution to the operator

subroutine ApplyBoundaryConditions(mesh, bc, tr_u, tr_dn_u)
  class(MeshPartition), intent(in)    :: mesh               !< mesh partition
  character,            intent(in)    :: bc(:)              !< boundary conds.
  real(RNP),            intent(inout) :: tr_u   (0:,0:,:,:) !< trace of u
  real(RNP),            intent(inout) :: tr_dn_u(0:,0:,:,:) !< trace of du/dn

  integer :: po
  integer :: b, f, i, k, o, p, q

  po = ubound(tr_u, 1)

  do b = 1, size(bc)
    associate(face => mesh % boundary(b) % face)

      select case(bc(b))

      case('D')
        ! Dirichlet: interior solution contributes twice to [u]
        do k = 1, size(face)
          f = face(k) % mesh_face % id          ! mesh face
          i = inner_side(face(k) % orientation) ! inner side
          o = 3 - i                             ! outer side
          do q = 0, po
          do p = 0, po
            tr_u(p,q,o,f) = -tr_u(p,q,i,f)
          end do
          end do
        end do

      case('N')
        ! Neumann: du/dn does not contribute, [u] = 0 due to extrapolation
        do k = 1, size(face)
          f = face(k) % mesh_face % id
          do q = 0, po
          do p = 0, po
            tr_dn_u(p,q,1,f) = 0
            tr_dn_u(p,q,2,f) = 0
          end do
          end do
        end do

      end select
    end associate
  end do

end subroutine ApplyBoundaryConditions

!-------------------------------------------------------------------------------
!> Compute jumps

subroutine ComputeJumps(tr_u, J_u)
  real(RNP), intent(in)  :: tr_u (0:,0:,:,:) !< trace of u
  real(RNP), intent(out) :: J_u  (0:,0:,:)   !< [u]ᵢ

  integer :: po
  integer :: f, p, q

  po = ubound(J_u, 1)

  do f = 1, size(J_u, 3)
    do q = 0, po
    do p = 0, po
      J_u(p,q,f) = tr_u(p,q,1,f) - tr_u(p,q,2,f)
    end do
    end do
  end do

end subroutine ComputeJumps

!-------------------------------------------------------------------------------
!> Compute average normal derivatives

subroutine ComputeNormalDerivatives(tr_dn_u, D_u)
  real(RNP), intent(in)  :: tr_dn_u (0:,0:,:,:) !< trace of du/dn
  real(RNP), intent(out) :: D_u     (0:,0:,:)   !< {∂ᵢu}

  integer :: po
  integer :: f, p, q

  po = ubound(D_u, 1)

  do f = 1, size(D_u, 3)
    do q = 0, po
    do p = 0, po
      D_u(p,q,f) = (tr_dn_u(p,q,1,f) - tr_dn_u(p,q,2,f)) * HALF
    end do
    end do
  end do

end subroutine ComputeNormalDerivatives

!===============================================================================

end submodule MP_Apply
