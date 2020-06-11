!> summary:  Application of the IP/DG elliptic operator with variable diffusivity
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Application of the IP/DG elliptic operator with constant diffusivity
!===============================================================================

submodule(CART__Elliptic_Operator_IP:MP_Apply) MP_Apply_VI

  use TPO__Elliptic_3d_RLVI
  use CART__Trace_Operator
  use CART__Normal_Trace_Operator
  implicit none

contains

!-------------------------------------------------------------------------------
!> Application of the operator with variable isotropic diffusivity

module subroutine Apply_VI(this, u, v)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  ! trace operators
  type(TraceOperator),       allocatable, save :: trace_op
  type(NormalTraceOperator), allocatable, save :: normal_trace_op

  ! trace variables
  real(RNP), allocatable, save :: tr_u (:,:,:,:)
  real(RNP), allocatable, save :: tr_qn(:,:,:,:)

  ! normal gradient on element faces
  real(RNP), allocatable, save :: q(:,:,:,:,:)

  ! jump and average derivative in direction xᵢ // face normal
  real(RNP), allocatable, save :: J_u(:,:,:) ! [u]ᵢ
  real(RNP), allocatable, save :: A_q(:,:,:) ! {q}ᵢ = {ν∇u}ᵢ

  integer :: po, ne, np 

  select type(eop => this % eop)
  class is (IP_ElementOperators1D)

    associate( mesh   => this % mesh,   &
               lambda => this % lambda, &
               nu     => this % nu_vi,  &
               nu_hat => this % nu_hat  )

      ! initialization .........................................................

      po = eop  % po
      ne = mesh % ne
      np = po + 1 

      ! workspace and operators
      !$omp single
      allocate(q(0:po,0:po,0:po,ne,3))
      allocate(tr_u(0:po, 0:po, 2, mesh%nf))
      allocate(tr_qn, mold=tr_u)
      allocate(J_u(0:po, 0:po, mesh%nf))
      allocate(A_q, mold=J_u)
      allocate(trace_op)
      allocate(normal_trace_op)
      !$omp end single

      call SetArray(tr_u , ZERO)
      call SetArray(tr_qn, ZERO)
      call SetArray(J_u  , ZERO)
      call SetArray(A_q  , ZERO)

      ! start generation of traces .............................................

      call ComputeNormalFluxes(np, ne, eop%D, mesh%dx, nu, u, q)

      call trace_op        % GetTrace_Start(mesh, u, tr_u , tag=1000)
      call normal_trace_op % GetTrace_Start(mesh, q, tr_qn, tag=2000)

      ! apply element stiffness operator .......................................

      call TPO_Elliptic_RLVI(eop%w, eop%D, lambda, nu, mesh%dx, u, v)

      ! finish generation of traces ............................................

      call trace_op        % GetTrace_Finish(mesh, tr_u   )
      call normal_trace_op % GetTrace_Finish(mesh, tr_qn)

      call ApplyBoundaryConditions(mesh, this%bc, tr_u, tr_qn)

      ! jumps and average derivatives ..........................................

      call ComputeJumps(tr_u, J_u)
      call ComputeAverages(tr_qn, A_q, normal=.true.)

      ! add fluxes .............................................................

      call AddFluxes(mesh, eop, nu, nu_hat, J_u, A_q, v)

      ! clean-up ...............................................................

      !$omp single
      deallocate(q, tr_u, tr_qn, J_u, A_q )
      deallocate(trace_op, normal_trace_op)
      !$omp end single

    end associate
  end select

end subroutine Apply_VI

!-------------------------------------------------------------------------------
!> Compute and add fluxes through element boundaries

subroutine AddFluxes(mesh, eop, nu, nu_hat, J_u, A_q, v)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in) :: mesh !< mesh partition
  class(IP_ElementOperators1D), intent(in) :: eop  !< ID/DG element operators

  real(RNP), intent(in)    :: nu(0:,0:,0:,:)  !< diffusivity
  real(RNP), intent(in)    :: nu_hat(0:,0:,:) !< diffusivity @ faces
  real(RNP), intent(in)    :: J_u(0:,0:,:)    !< [u]ᵢ
  real(RNP), intent(in)    :: A_q(0:,0:,:)    !< {q}ᵢ
  real(RNP), intent(inout) :: v(0:,0:,0:,:)   !< result

  ! local variables ............................................................

  real(RNP), dimension(:,:), allocatable :: Mf1, Mf2, Mf3
  real(RNP), dimension(:),   allocatable :: delta_0, delta_P
  real(RNP) :: g(3), mu(3)
  real(RNP) :: cd_0, cd_P, cp_0, cp_P
  integer   :: i, j, k, e, f(6), f_0, f_P

  associate( P  => eop  % po, &
             Ms => eop  % w,  &
             Ds => eop  % D,  &
             dx => mesh % dx  )

    ! auxiliaries .............................................................

    ! face mass matrices

    allocate(Mf1(0:P,0:P), Mf2(0:P,0:P), Mf3(0:P,0:P))

    g(1) = dx(2) * dx(3) / 4
    g(2) = dx(3) * dx(1) / 4
    g(3) = dx(1) * dx(2) / 4

    do j = 0, P
    do i = 0, P
      Mf1(i,j) = g(1) * Ms(i) * Ms(j)
      Mf2(i,j) = g(2) * Ms(i) * Ms(j)
      Mf3(i,j) = g(3) * Ms(i) * Ms(j)
    end do
    end do

    ! delta function
    allocate(delta_0(0:P), source = [ ONE, (ZERO, i=1,P) ])
    allocate(delta_P(0:P), source = [ (ZERO, i=1,P), ONE ])

    ! penalties
    mu(1) = eop % PenaltyFactor(dx(1))
    mu(2) = eop % PenaltyFactor(dx(2))
    mu(3) = eop % PenaltyFactor(dx(3))

    ! add fluxes ...............................................................

    g = ONE / dx

    !$omp do
    do e = 1, mesh % ne

      f = mesh % element(e) % face % id

      ! direction 1: v = v + Mf₁ (-[ϕ]₁{ν∇u}₁ - {ν∇ϕ}₁[u]₁ + μ⟨ν⟩[ϕ]₁[u]₁)
      !
      ! where, with ϕ = ℓ_ijk, at  ξ = -1
      !
      !   [ϕ]₁   = -delta_0(i)
      !   {ν∇ϕ}₁ =  nu(0,j,k,e) * g(1) *  Ds(0,i)
      !   [u]₁   =  J_u(j,k,f₁)
      !   {ν∇u}₁ =  A_q(j,k,f₁)
      !
      ! and, at  ξ = +1
      !
      !   [ϕ]₁   =  delta_P(i)
      !   {ν∇ϕ}₁ =  nu(P,j,k,e) * g(1) *  Ds(P,i)
      !   [u]₁   =  J_u(j,k,f₂)
      !   {ν∇u}₁ =  A_q(j,k,f₂)

      f_0 = f(1)
      f_P = f(2)

      do k = 0, P
      do j = 0, P

        cd_0 = -nu(0,j,k,e) * g(1)
        cd_P = -nu(P,j,k,e) * g(1)
        cp_0 = -nu_hat(j,k,f(1)) * mu(1)
        cp_P =  nu_hat(j,k,f(2)) * mu(1)

        do i = 0, P

          v(i,j,k,e) = v(i,j,k,e)                                              &
            + Mf1(j,k) * ( delta_0(i) * A_q(j,k,f_0)                           &
                         - delta_P(i) * A_q(j,k,f_P)                           &
                         + (cd_0 * Ds(0,i) + cp_0 * delta_0(i)) * J_u(j,k,f_0) &
                         + (cd_P * Ds(P,i) + cp_P * delta_P(i)) * J_u(j,k,f_P) )
        end do
      end do
      end do

      ! direction 2: v = v + Mf₂ (-[ϕ]₂{ν∇u}₂ - {ν∇ϕ}₂[u]₂ + μ⟨ν⟩[ϕ]₂[u]₂)

      f_0 = f(3)
      f_P = f(4)

      do k = 0, P
      do i = 0, P

        cd_0 = -nu(i,0,k,e) * g(2)
        cd_P = -nu(i,P,k,e) * g(2)
        cp_0 = -nu_hat(i,k,f(3)) * mu(2)
        cp_P =  nu_hat(i,k,f(4)) * mu(2)

        do j = 0, P

          v(i,j,k,e) = v(i,j,k,e)                                              &
            + Mf2(i,k) * ( delta_0(j) * A_q(i,k,f_0)                           &
                         - delta_P(j) * A_q(i,k,f_P)                           &
                         + (cd_0 * Ds(0,j) + cp_0 * delta_0(j)) * J_u(i,k,f_0) &
                         + (cd_P * Ds(P,j) + cp_P * delta_P(j)) * J_u(i,k,f_P) )
        end do
      end do
      end do

      ! direction 3: v = v + Mf₃ (-[ϕ]₃{ν∇u}₃ - {ν∇ϕ}₃[u]₃ + μ⟨ν⟩[ϕ]₃[u]₃)

      f_0 = f(5)
      f_P = f(6)

      do j = 0, P
      do i = 0, P

        cd_0 = -nu(i,j,0,e) * g(3)
        cd_P = -nu(i,j,P,e) * g(3)
        cp_0 = -nu_hat(i,j,f(5)) * mu(3)
        cp_P =  nu_hat(i,j,f(6)) * mu(3)

        do k = 0, P

          v(i,j,k,e) = v(i,j,k,e)                                              &
            + Mf3(i,j) * ( delta_0(k) * A_q(i,j,f_0)                           &
                         - delta_P(k) * A_q(i,j,f_P)                           &
                         + (cd_0 * Ds(0,k) + cp_0 * delta_0(k)) * J_u(i,j,f_0) &
                         + (cd_P * Ds(P,k) + cp_P * delta_P(k)) * J_u(i,j,f_P) )
        end do
      end do
      end do

    end do
    !$omp end do

  end associate

end subroutine AddFluxes

!-------------------------------------------------------------------------------
!> Elementwise computation of diffusive fluxes parallel to face normals
!>
!> Computes the normal components of ν grad(u) for all element boundary points.
!> Entries corresponding to interior points or tangential components are set
!> to zero.

subroutine ComputeNormalFluxes(np, ne, Ds, dx, nu, u, q)
  integer,   intent(in)  :: np               !< number of points per direction
  integer,   intent(in)  :: ne               !< number of elements
  real(RNP), intent(in)  :: Ds(np,np)        !< 1D standard diff matrix
  real(RNP), intent(in)  :: dx(3)            !< element extensions
  real(RNP), intent(in)  :: nu(np,np,np,ne)  !< diffusivity
  real(RNP), intent(in)  :: u(np,np,np,ne)   !< 3D scalar field
  real(RNP), intent(out) :: q(np,np,np,ne,3) !< element-wise gradient of u

  real(RNP), allocatable :: Ds_0(:), Ds_P(:)
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
  allocate(Ds_0, source=Ds( 1,:))
  allocate(Ds_P, source=Ds(np,:))

  ! metric coefficients
  g = 2 / dx

  ! result
  call SetArray(q, ZERO, multi=.true.)

  !$acc data present(u,q) copyin(Ds_0,Ds_P,g) async
  !$acc parallel async &
  !$acc & device_type(nvidia) num_workers(1024/vec_len) vector_length(vec_len)
  !$acc loop gang worker

  !$omp do private(e)
  do e = 1, ne

    ! q1 = ν ∂u/∂x1 @ face 1,2 .................................................

    !$acc loop collapse(2) vector
    do k = 1, np
    do j = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + Ds_0(m) * u(m,j,k,e)
        tmp2 = tmp2 + Ds_P(m) * u(m,j,k,e)
      end do
      q( 1,j,k,e,1) = g(1) * nu( 1,j,k,e) * tmp1
      q(np,j,k,e,1) = g(1) * nu(np,j,k,e) * tmp2

    end do
    end do

    ! q2 = ν ∂u/∂x2 @ face 3,4 .................................................

    !$acc loop collapse(2) vector
    do k = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + Ds_0(m) * u(i,m,k,e)
        tmp2 = tmp2 + Ds_P(m) * u(i,m,k,e)
      end do
      q(i, 1,k,e,2) = g(2) * nu(i, 1,k,e) * tmp1
      q(i,np,k,e,2) = g(2) * nu(i,np,k,e) * tmp2

    end do
    end do

    ! q3 = ν ∂u/∂x3 @ face 5,6 .................................................

    !$acc loop collapse(2) vector
    do j = 1, np
    do i = 1, np

      tmp1 = 0
      tmp2 = 0
      do m = 1, np
        tmp1 = tmp1 + Ds_0(m) * u(i,j,m,e)
        tmp2 = tmp2 + Ds_P(m) * u(i,j,m,e)
      end do
      q(i,j, 1,e,3) = g(3) * nu(i,j, 1,e) * tmp1
      q(i,j,np,e,3) = g(3) * nu(i,j,np,e) * tmp2

    end do
    end do

  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine ComputeNormalFluxes

!===============================================================================

end submodule MP_Apply_VI
