!> summary:  Application of the IP/DG elliptic operator: const isotropic
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Application of the IP/DG elliptic operator: const isotropic diffusivity
!===============================================================================

submodule(CART__Elliptic_Operator_IP:MP_Apply) MP_Apply_CI
  use CART__TPO_Elliptic_CI
  use CART__Mesh_Partition
  use CART__Trace_Operator
  use CART__Normal_Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> Application of the operator with constant isotropic diffusivity

module subroutine Apply_CI(this, lambda, nu, bc, u, v)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)  :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)  :: nu             !< diffusivity
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

  associate(mesh => this % mesh, eop => this % eop)

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

!===============================================================================

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

    allocate(delta_0(0:po), source = ZERO)
    delta_0(0) = ONE

    allocate(delta_P(0:po), source = ZERO)
    delta_P(po) = ONE

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

end submodule MP_Apply_CI
