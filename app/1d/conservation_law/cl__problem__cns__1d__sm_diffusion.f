submodule(CL__Problem__CNS__1D) SM_Diffusion
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Returns the physical diffusivity matrix

  pure module function PhysicalDiffusivity(this, u) result(A_pd)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:) !< conservative variables
    real(RNP) :: A_pd(3,3)

    real(RNP) :: amv, aev, aeT, c1, cc, c3, dv1, dv2, dT1, dT2, dT3

    c1  =  1  / u(1)                 ! 1/ρ
    cc  =  c1 / this % c_v           ! 1/ρc_v

    amv =  4 * THIRD * this % eta    !  momentum diffusivity related to ∂v/∂x
    aev =  amv * c1 * u(2)           !  energy   diffusivity related to ∂v/∂x
    aeT =  this % lambda             !  energy   diffusivity related to ∂T/∂x

    dv1 = -c1 * u(2)                 ! ∂v/∂u₁
    dv2 =  c1                        ! ∂v/∂u₂

    dT1 =  cc * c1**2 * (u(2)**2 - u(1)*u(3))  ! ∂T/∂u₁
    dT2 = -cc * c1 * u(2)                      ! ∂T/∂u₂
    dT3 =  cc                                  ! ∂T/∂u₃

    A_pd(1,:) = [ ZERO                  , ZERO                  , ZERO      ]
    A_pd(2,:) = [ amv * dv1             , amv * dv2             , ZERO      ]
    A_pd(3,:) = [ aev * dv1 + aeT * dT1 , aev * dv2 + aeT * dT2 , aeT * dT3 ]

  end function PhysicalDiffusivity

  !-----------------------------------------------------------------------------
  !> Returns the streamline diffusivity matrix

  pure module function StreamlineDiffusivity(this, theta, u) result(A_d)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:)  !< conservative variables
    real(RNP), intent(in) :: theta !< streamline-diffusion time scale
    real(RNP) :: A_d(3,3)

    A_d = ConvectiveJacobian(this, u)
    A_d = (HALF * theta) * matmul(A_d, A_d)

  end function StreamlineDiffusivity

  !-----------------------------------------------------------------------------
  !> Unified CNS diffusion term combining physical and streamline contributions
  !>
  !> The composition of the diffusion term is controlled by the argument `comp`:
  !> 'P' selects physical and 'S' streamline diffusion, whereas 'T' yields the
  !> total diffusion as the sum of both.
  !>

  module subroutine GetHybridDiffusionTerm &
      (this, cl_operator, comp, theta, bv, u_0, u, r_d)

    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator !< spatial operators
    character,                intent(in)  :: comp        !< composition flag
    real(RNP),                intent(in)  :: theta       !< SD time scale
    real(RNP),                intent(in)  :: bv(:,:)     !< boundary values
    real(RNP), contiguous,    intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
    real(RNP), contiguous,    intent(in)  :: u  (0:,:,:) !< u(x,t)
    real(RNP), contiguous,    intent(out) :: r_d(0:,:,:) !< diffusion RHS

    ! local variables ..........................................................

    real(RNP), allocatable, save :: dx_u(:,:,:)
    ! element solution derivatives ∂u/∂x (0:po,1:ne,1:3)

    real(RNP), allocatable, save :: A(:,:,:,:)
    ! element diffusivity matrices, A = A_pd + A_sd (0:po,1:ne,1:3,1:3)

    real(RNP), allocatable, save :: A_hat(:,:)
    ! diagonal interface diffusivity matrices Â (0:ne,1:3)

    real(RNP), allocatable, save :: jmp_u(:,:)
    ! jump of solution across element boundaries (0:ne,1:3)

    real(RNP), allocatable, save :: avg_q(:,:)
    ! average of diffusive fluxes over element boundaries (0:ne,1:3)

    real(RNP) :: g1, mu

    logical :: has_pd ! switch for physical diffusion
    logical :: has_sd ! switch for streamline diffusion

    integer :: e, i, k

    !$omp master !!! not ready for OpenMP, enforcing sequential execution !!!

    associate( activity => cl_operator % activity &
             , eop      => cl_operator % eop      &
             , M        => cl_operator % eop % w  &
             , D        => cl_operator % eop % D  &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , dx       => cl_operator % dx       )

      ! preliminaries ..........................................................

      has_pd = scan(comp,'PT') > 0
      has_sd = scan(comp,'ST') > 0

      ! work space .............................................................

      allocate(A(0:po,ne,3,3), source = ZERO)
      allocate(A_hat(0:po,3))
      allocate(jmp_u(0:po,3))
      allocate(avg_q(0:po,3))

      ! diffusivity ............................................................

      !$omp do
      do e = 1, ne
        select case(activity(e))

        case(1) ! active element

          if (has_pd) then ! add physical diffusivity
            do i = 0, po
              A(i,e,1:3,1:3) = PhysicalDiffusivity(this, u(i,e,:))
            end do
          end if

          if (has_sd) then ! add streamline diffusivity
            do i = 0, po
              A(i,e,1:3,1:3) = A(i,e,1:3,1:3) &
                             + StreamlineDiffusivity(this, theta, u(i,e,:))
            end do
          end if

        case(0) ! frozen element, only boundary values required

          if (has_pd) then
            !BW! A( 0,e,1:3,1:3) = <A_pd>
            !BW! A(po,e,1:3,1:3) = <A_pd>
          end if

          if (has_sd) then
            !BW! ...
          end if

        end select
      end do

!### part to be adapted from DG__Elliptic_Operator__1D  % Eval_RV ##############

      ! application of interior operator .......................................

!BW! code for scalar case -- to be adapted
!BW!      !$omp do
!BW!      do e = 1, ne
!BW!        if (mask(e)) then
!BW!          dx_u(:,e) = g1 * matmul(D, u(:,e))
!BW!          r(:,e) = matmul(M * nu(:,e) * dx_u(:,e), D)
!BW!        else
!BW!          dx_u( 0,e) = g1 * dot_product(D( 0,:), u(:,e))
!BW!          dx_u(po,e) = g1 * dot_product(D(po,:), u(:,e))
!BW!          r(:,e) = ZERO
!BW!        end if
!BW!      end do

      ! interior fluxes (and Â?) ..................................

      !!! theoretical description to be derived !!!

!BW! code for scalar case -- to be adapted
!BW!      !$omp do collapse(2)
!BW!      do k = 1, 3
!BW!      do e = 1, ne-1
!BW!        jmp_u (e) = u(po,e) - u(0,e+1)
!BW!        ql = nu(po,e  ) * dx_u(po,e  )
!BW!        qr = nu( 0,e+1) * dx_u( 0,e+1)
!BW!        jmp_q(e) = (ql - qr)
!BW!        avg_q(e) = (ql + qr) * HALF
!BW!        A_hat(e,k) = max(?,?)
!BW!      end do
!BW!      end do


      ! boundary fluxes (and Â?) .............................................

      !!! theoretical description to be derived !!!

      ! left boundary
      select case(this % bc(1))
      case('P')
      case default ! nu jumps, inner fluxes
        do k = 1, 3
          jmp_u(0,k) = 0
          avg_q(0,k) = A(0,1,k,1) * dx_u(0,1,1) &
                     + A(0,1,k,2) * dx_u(0,1,2) &
                     + A(0,1,k,3) * dx_u(0,1,3)
          A_hat(0,k) = A(0,1,k,k)
        end do
      end select
!#! !!! for start consider only Neumann and periodic cases  !!!

!#! code for scalar case -- to be adapted
!#!      ! left boundary
!#!      select case(bc(1))
!#!
!#!      case('D')
!#!        nu_max(0) = nu(0,1)
!#!        if (has_bv) then
!#!          jmp_u(0) = 2 * (bv(1) - u(0,1))
!#!        else
!#!          jmp_u(0) = 2 * (      - u(0,1))
!#!        end if
!#!        avg_q(0) = nu(0,1) * dx_u(0,1)
!#!
!#!      case('N')
!#!        nu_max(0) = nu(0,1)
!#!        if (has_bv) then
!#!          avg_q(0) = nu(0,1) * bv(1)
!#!        end if
!#!
!#!      case('P')
!#!        nu_max(0) = max(nu(po,ne), nu(0,1))
!#!        jmp_u (0) = u(po,ne) - u(0,1)
!#!        ql = nu(po,ne) * dx_u(po,ne)
!#!        qr = nu( 0, 1) * dx_u( 0, 1)
!#!        jmp_q(0) = (ql - qr)
!#!        avg_q(0) = (ql + qr) * HALF
!#!      end select
!#!
!#!      ! right boundary
!#!      select case(bc(2))
!#!
!#!      case('D')
!#!        nu_max(ne) = nu(po,ne)
!#!        if (has_bv) then
!#!          jmp_u(ne) = 2 * (u(po,ne) - bv(2))
!#!        else
!#!          jmp_u(ne) = 2 * (u(po,ne)        )
!#!        end if
!#!        avg_q(ne) = nu(po,ne) * dx_u(po,ne)
!#!
!#!      case('N')
!#!        nu_max(ne) = nu(po,ne)
!#!        if (has_bv) then
!#!          avg_q(ne) = nu(po,ne) * bv(2)
!#!        end if
!#!
!#!      case('P')
!#!        nu_max(ne) = nu_max(0)
!#!        jmp_u (ne) = jmp_u (0)
!#!        jmp_q (ne) = jmp_q (0)
!#!        avg_q (ne) = avg_q (0)
!#!      end select

      ! apply ..................................................................

      !!! theoretical description to be derived !!!

!#! code for scalar case -- to be adapted
!#!      !$omp do
!#!      do e = 1, ne
!#!        if (mask(e)) then
!#!
!#!          ! - {ν∂v/∂x}[u]
!#!          cl = -g1/2 * nu( 0,e) * jmp_u(e-1)
!#!          cr = -g1/2 * nu(po,e) * jmp_u(e  )
!#!          do i = 0, po
!#!            r(i,e) = r(i,e) + cl * D(0,i) + cr * D(po,i)
!#!          end do
!#!
!#!          ! - [v]{ν∂u/∂x}
!#!          r( 0,e) = r( 0,e) + avg_q(e-1)
!#!          r(po,e) = r(po,e) - avg_q(e  )
!#!
!#!          ! + μ⟨ν⟩[v][u]
!#!          r( 0,e) = r( 0,e) - mu * nu_max(e-1) * jmp_u(e-1)
!#!          r(po,e) = r(po,e) + mu * nu_max(e  ) * jmp_u(e  )
!#!
!#!
!#!          if (has_f) then
!#!            r(:,e) = f(:,e) - r(:,e)     !??????? VORZEICHEN ???????
!#!          end if
!#!
!#!        end if
!#!      end do

!###############################################################################

      !BW! add deallocation of remaining arrays

    end associate

    !$omp end master

  end subroutine GetHybridDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Implicit CNS diffusion solver

  module subroutine DiffusionSolver &
      (this, cl_operator, dt, theta, bv, f, u_0, u, method, i_max, r_red, r_max)

    class(CL_Problem_CNS_1D), intent(in)    :: this
    class(CL_Operator_1D),    intent(in)    :: cl_operator
    real(RNP),                intent(in)    :: dt          !< ∆t = t - t₀
    real(RNP),                intent(in)    :: theta       !< SD time scale θ
    real(RNP),                intent(in)    :: bv (:,:)    !< boundary values
    real(RNP), contiguous,    intent(in)    :: f  (0:,:,:) !< sources
    real(RNP), contiguous,    intent(in)    :: u_0(0:,:,:) !< frozen solution
    real(RNP), contiguous,    intent(inout) :: u  (0:,:,:) !< approx solution
    integer,                  intent(in)    :: method      !< solution method
    integer,                  intent(in)    :: i_max       !< max num iterations
    real(RNP), optional,      intent(in)    :: r_red       !< residual reduction
    real(RNP), optional,      intent(in)    :: r_max       !< max residual

  end subroutine DiffusionSolver

  !=============================================================================

end submodule SM_Diffusion
