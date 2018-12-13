!> summary:  Application of the IP/DG elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Application of the IP/DG elliptic operator
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Apply
  use Execution_Control
  use CART__TPO_Elliptic_CI
  use CART__Mesh_Partition
  use CART__Trace_Operator
  use CART__Normal_Trace_Operator
  implicit none

  !> mapping of boundary face orientation to inner trace side
  integer, parameter :: inner_side(-3:3) = [ 2, 2, 2, 0, 1, 1, 1 ]

contains

!-------------------------------------------------------------------------------
!> Application of the IP/DG elliptic operator

module subroutine Apply(this, bc, u, v)
  class(EllipticOperator3D_IP), intent(in) :: this
  character, intent(in)  :: bc(:)          !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  if (allocated(this % nu_ci)) then
    call Apply_CI(this, bc, u, v)
  else
    call AssignScalar(v, ZERO)
    call Error( 'Apply', 'operator not properly initialized', &
                'MP_Apply @ CART__Elliptic_Operator_IP'       )
  end if

end subroutine Apply

!===============================================================================
! Specific procedures for constant isotropic diffusivity (CI)

!-------------------------------------------------------------------------------
!> Application of the operator with constant isotropic diffusivity

subroutine Apply_CI(this, bc, u, v)
  class(EllipticOperator3D_IP), intent(in) :: this
  character, intent(in)  :: bc(:)          !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  ! local variables ............................................................

  procedure(TPO_Elliptic_CI_Proc), pointer, save :: StiffnessOperator

  ! trace operators
  type(TraceOperator),       allocatable, save :: trace_op
  type(NormalTraceOperator), allocatable, save :: normal_trace_op

  ! trace variables
  real(RNP), allocatable, save :: tr_u(:,:,:,:)
  real(RNP), allocatable, save :: tr_dn_u(:,:,:,:)

  ! normal gradient on element faces
  real(RNP), allocatable, save :: grad_u(:,:,:,:,:)

  ! jump and average derivative in direction xᵢ // face normal
  real(RNP), allocatable, save :: J_u(:,:,:) ! [u]ᵢ
  real(RNP), allocatable, save :: D_u(:,:,:) ! {u}ᵢ

  integer :: po, ne, np = -1

  associate( mesh   => this % mesh,   eop => this % eop,  &
             lambda => this % lambda, nu  => this % nu_ci )

    ! initialization ...........................................................

    po = eop  % po
    ne = mesh % ne

    ! procedure for evaluating the element operators
    if (np /= po + 1) then
      np  = po + 1
      call TPO_Elliptic_CI_Assign(np, StiffnessOperator)
    end if

    ! workspace and operators
    !$omp single
    allocate(grad_u(0:po,0:po,0:po,ne,3))
    allocate(tr_u(0:po, 0:po, 2, mesh%nf))
    allocate(tr_dn_u, mold=tr_u)
    allocate(J_u(0:po, 0:po, mesh%nf))
    allocate(D_u, mold=J_u)
    allocate(trace_op)
    allocate(normal_trace_op)
    !$omp end single

    call AssignScalar(tr_u   , ZERO)
    call AssignScalar(tr_dn_u, ZERO)
    call AssignScalar(J_u    , ZERO)
    call AssignScalar(D_u    , ZERO)

    ! start generation of traces ...............................................

    call NormalDerivatives(np, ne, eop%D, mesh%dx, u, grad_u)

    call trace_op        % GetTrace_Start(mesh, u     , tr_u   , tag=1000)
    call normal_trace_op % GetTrace_Start(mesh, grad_u, tr_dn_u, tag=2000)

    ! apply element stiffness operator .........................................

    call StiffnessOperator(np, ne, eop%w, eop%L, lambda, nu, mesh%dx, u, v)

    ! finish generation of traces ..............................................

    call trace_op        % GetTrace_Finish(mesh, tr_u   )
    call normal_trace_op % GetTrace_Finish(mesh, tr_dn_u)

    call ApplyBoundaryConditions(mesh, bc, tr_u, tr_dn_u)

    ! jumps and average derivatives ............................................

    call ComputeJumps(tr_u, J_u)
    call ComputeNormalDerivatives(tr_dn_u, D_u)

    ! add fluxes ...............................................................

    call AddFluxes_CI(mesh, eop, nu, J_u, D_u, v)

    ! clean-up .................................................................

    !$omp single
    deallocate(grad_u, tr_u, tr_dn_u, J_u, D_u )
    deallocate(trace_op, normal_trace_op)
    !$omp end single

  end associate

end subroutine Apply_CI

!-------------------------------------------------------------------------------
!> Compute and add fluxes through element boundaries

subroutine AddFluxes_CI(mesh, eop, nu, J_u, D_u, v)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in) :: mesh !< mesh partition
  class(IP_ElementOperators1D), intent(in) :: eop  !< ID/DG element operators

  real(RNP), intent(in)    :: nu            !< diffusivity
  real(RNP), intent(in)    :: J_u(0:,0:,:)  !< [u]ᵢ
  real(RNP), intent(in)    :: D_u(0:,0:,:)  !< {∂ᵢu}
  real(RNP), intent(inout) :: v(0:,0:,0:,:) !< result

  ! local variables ............................................................

  real(RNP), dimension(:,:), allocatable :: M1, M2, M3
  real(RNP), dimension(:),   allocatable :: delta_0, delta_P
  real(RNP) :: c(3), mu(3)
  integer :: i, j, k, e, f(6)

  associate( po => eop  % po, &
             Ms => eop  % w,  &
             Ds => eop  % D,  &
             dx => mesh % dx  )

    ! auxiliaries .............................................................

    ! face mass matrices

    allocate(M1(0:po,0:po), M2(0:po,0:po), M3(0:po,0:po))

    c(1) = dx(2) * dx(3) / 4
    c(2) = dx(3) * dx(1) / 4
    c(3) = dx(1) * dx(2) / 4

    do j = 0, po
    do i = 0, po
      M1(i,j) = c(1) * Ms(i) * Ms(j)
      M2(i,j) = c(2) * Ms(i) * Ms(j)
      M3(i,j) = c(3) * Ms(i) * Ms(j)
    end do
    end do

    ! delta function
    allocate(delta_0(0:po), source = [ ONE, (ZERO, i=1,po) ])
    allocate(delta_P(0:po), source = [ (ZERO, i=1,po), ONE ])

    ! penalties
    mu(1) = eop % PenaltyFactor(dx(1))
    mu(2) = eop % PenaltyFactor(dx(2))
    mu(3) = eop % PenaltyFactor(dx(3))

    ! add fluxes ...............................................................

    c = ONE / dx

    !$omp do
    do e = 1, mesh % ne

      f = mesh % element(e) % face % id

      do k = 0, po
      do j = 0, po
      do i = 0, po


        v(i,j,k,e)                                                             &

        = v(i,j,k,e)                                                           &

        + nu *                                                                 &
          (                                                                    &
            M1(j,k) * (                                                        &
                      - (c(1)*Ds( 0,i) + mu(1)*delta_0(i)) * J_u(j,k,f(1))     &
                      +                        delta_0(i)  * D_u(j,k,f(1))     &
                      - (c(1)*Ds(po,i) - mu(1)*delta_P(i)) * J_u(j,k,f(2))     &
                      -                        delta_P(i)  * D_u(j,k,f(2))     &
                      )                                                        &

          + M2(i,k) * (                                                        &
                      - (c(2)*Ds( 0,j) + mu(2)*delta_0(j)) * J_u(i,k,f(3))     &
                      +                        delta_0(j)  * D_u(i,k,f(3))     &
                      - (c(2)*Ds(po,j) - mu(2)*delta_P(j)) * J_u(i,k,f(4))     &
                      -                        delta_P(j)  * D_u(i,k,f(4))     &
                      )                                                        &

          + M3(i,j) * (                                                        &
                      - (c(3)*Ds( 0,k) + mu(3)*delta_0(k)) * J_u(i,j,f(5))     &
                      +                        delta_0(k)  * D_u(i,j,f(5))     &
                      - (c(3)*Ds(po,k) - mu(3)*delta_P(k)) * J_u(i,j,f(6))     &
                      -                        delta_P(k)  * D_u(i,j,f(6))     &
                      )                                                        &
          )
      end do
      end do
      end do

    end do
    !$omp end do

  end associate

end subroutine AddFluxes_CI

!===============================================================================
! Specific procedures for variable isotropic diffusivity (VI)

! ... TBD ...

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
