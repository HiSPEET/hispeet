submodule(CL__Problem__CNS__1D) SM_Diffusion
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Returns the physical diffusivity matrix

  pure module function PhysicalDiffusivity(this, u) result(A_pd)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:) !< conservative variables
    real(RNP) :: A_pd(3,3)

    real(RNP) :: amv, aev, aeT, c1, cc, dv1, dv2, dT1, dT2, dT3

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
  !> Homogeneous boundary conditions are applied in the case that `bv` omitted.

!!!!!!!!!!!!!!!!!!!!!! CHOOSE SIGN FOR TERM PLACED ON RHS !!!!!!!!!!!!!!!!!!!!!!

  module subroutine GetHybridDiffusionTerm &
      (this, cl_operator, comp, theta, bv, u_0, u, r_d)

    class(CL_Problem_CNS_1D), intent(in)  :: this
    class(CL_Operator_1D),    intent(in)  :: cl_operator !< spatial operators
    character,                intent(in)  :: comp        !< composition flag
    real(RNP),                intent(in)  :: theta       !< SD time scale
    real(RNP), optional,      intent(in)  :: bv(:,:)     !< boundary values
    real(RNP), contiguous,    intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
    real(RNP), contiguous,    intent(in)  :: u  (0:,:,:) !< u(x,t)
    real(RNP), contiguous,    intent(out) :: r_d(0:,:,:) !< diffusion RHS

    ! local variables ..........................................................

    real(RNP), allocatable, save :: A(:,:,:,:)
    ! element diffusivity matrices, A = A_pd + A_sd (0:po,1:ne,1:3,1:3)

    real(RNP), allocatable, save :: A_hat(:,:)
    ! diagonal interface diffusivity matrices Â (0:ne,1:3)

    real(RNP), allocatable, save :: jmp_u(:,:)
    ! jump of solution across element boundaries (0:ne,1:3)

    real(RNP), allocatable, save :: avg_q(:,:)
    ! average of diffusive fluxes over element boundaries (0:ne,1:3)

    real(RNP), allocatable, save :: u_l(:,:), u_r(:,:)
    ! left and right traces of solution at element interfaces (0:ne,1:3)

    real(RNP), allocatable, save :: q_l(:,:), q_r(:,:)
    ! left and right diffusive fluxes at element interfaces (0:ne,1:3)

    real(RNP), allocatable, save :: A_l(:,:), A_r(:,:)
    ! left and right traces of diagonal interface diffusivity matrices (0:ne,1:3)

    real(RNP), allocatable :: MD_t(:,:)
    ! transpose of weighted derivative matrix (MD)ᵗ

    real(RNP), allocatable :: dx_u(:,:)
    ! element solution derivatives ∂u/∂x (0:po,1:3)

    real(RNP), allocatable :: q(:,:)
    ! element diffusive flux

    logical :: has_pd ! switch for physical diffusion
    logical :: has_sd ! switch for streamline diffusion

    real(RNP) :: g, mu
    real(RNP) :: c_0, c_po

    integer :: e, i, j, k

    !$omp master !!! not ready for OpenMP, enforcing sequential execution !!!

    associate( activity => cl_operator % activity &
             , eop      => cl_operator % eop      &
             , M        => cl_operator % eop % w  &
             , D        => cl_operator % eop % D  &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , dx       => cl_operator % dx       )

      ! initialization .........................................................

      has_pd = scan(comp,'PT') > 0
      has_sd = scan(comp,'ST') > 0

      ! shared workspace
      allocate(A(0:po,ne,3,3), A_hat(0:po,3))
      allocate(A_l, A_r, u_l, u_r, q_l, q_r, jmp_u, avg_q, mold = A_hat)

      ! private workspace
      allocate(MD_t(0:po,0:po), dx_u(0:po,3), q(0:po,3))

      mu = eop % PenaltyFactor(dx)

      do i = 0, po
      do j = 0, po
        MD_t(i,j) = M(j) * D(j,i)
      end do
      end do

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
                             + StreamlineDiffusivity(this, theta, u_0(i,e,:))
            end do
          end if

          do k = 1, 3
            A_r(e-1,k) = A( 0,e,k,k)
            A_l(e  ,k) = A(po,e,k,k)
          end do

        case(0) ! frozen element, only boundary values required

          if (has_pd) then
            A( 0,e,1:3,1:3) = PhysicalDiffusivity(this, u(0 ,e,:))
            A(po,e,1:3,1:3) = PhysicalDiffusivity(this, u(po,e,:))
          end if

          if (has_sd) then
            A( 0,e,1:3,1:3) = A( 0,e,1:3,1:3) &
                            + StreamlineDiffusivity(this, theta, u_0(0 ,e,:))
            A(po,e,1:3,1:3) = A(po,e,1:3,1:3) &
                            + StreamlineDiffusivity(this, theta, u_0(po,e,:))
          end if

          do k = 1, 3
            A_r(e-1,k) = A( 0,e,k,k)
            A_l(e  ,k) = A(po,e,k,k)
          end do

        case default

          A( 0,e,1:3,1:3) = 0
          A(po,e,1:3,1:3) = 0

          A_r(e-1,1:3) = 0
          A_l(e  ,1:3) = 0

        end select
      end do

      ! application of interior operator .......................................

      g = 2 / dx

      !$omp do
      do e = 1, ne
        select case(activity(e))

        case(1) ! active element

          dx_u = g * matmul(D, u(0:po,e,1:3))

          do k = 1, 3
            q(:,k) = A(:,e,k,1) * dx_u(:,1) &
                   + A(:,e,k,2) * dx_u(:,2) &
                   + A(:,e,k,3) * dx_u(:,3)
          end do

          r_d(0:po,e,1:3) = -matmul(MD_t, q)

          ! save traces for computation of fluxes
          u_r(e-1,1:3) = u( 0, e, 1:3)
          u_l(e  ,1:3) = u(po, e, 1:3)
          q_r(e-1,1:3) = q( 0,    1:3)
          q_l(e  ,1:3) = q(po,    1:3)

        case(0) ! frozen element, only traces required

          r_d(0:po,e,1:3) = 0

          u_r(e-1,1:3) = u( 0, e, 1:3)
          u_l(e  ,1:3) = u(po, e, 1:3)

          dx_u( 0,1:3) = g * matmul(D( 0,:), u(0:po,e,1:3))
          dx_u(po,1:3) = g * matmul(D(po,:), u(0:po,e,1:3))

          q_r(e-1,1:3) = A( 0,e,1:3,1) * dx_u( 0,1) &
                       + A( 0,e,1:3,2) * dx_u( 0,2) &
                       + A( 0,e,1:3,3) * dx_u( 0,3)

          q_l(e  ,1:3) = A(po,e,1:3,1) * dx_u(po,1) &
                       + A(po,e,1:3,2) * dx_u(po,2) &
                       + A(po,e,1:3,3) * dx_u(po,3)

        case default

          r_d(0:po,e,1:3) = 0

          u_r(e-1,1:3) = 0
          u_l(e  ,1:3) = 0
          q_r(e-1,1:3) = 0
          q_l(e  ,1:3) = 0

        end select
      end do

      ! application of boundary conditions .....................................

      ! left boundary
      select case(this % bc(1))

      case('P') ! periodic
        u_l(0,1:3) = u_l(ne,1:3)
        q_l(0,1:3) = q_l(ne,1:3)
        A_l(0,1:3) = A_l(ne,1:3)

      case default ! extrapolation
        u_l(0,1:3) = u_r(0,1:3)
        q_l(0,1:3) = q_r(0,1:3)
        A_l(0,1:3) = A_r(0,1:3)

      end select

      ! right boundary
      select case(this % bc(2))

      case('P') ! periodic
        u_r(ne,1:3) = u_r(0,1:3)
        q_r(ne,1:3) = q_r(0,1:3)
        A_r(ne,1:3) = A_r(0,1:3)

      case default ! extrapolation
        u_r(ne,1:3) = u_l(ne,1:3)
        q_r(ne,1:3) = q_l(ne,1:3)
        A_r(ne,1:3) = A_l(ne,1:3)

      end select

      ! element interface values ...............................................

      !$omp do collapse(2)
      do k = 1, 3
      do e = 0, ne
        jmp_u(e,k) =  u_l(e,k) - u_r(e,k)
        avg_q(e,k) = (q_l(e,k) + q_r(e,k)) * HALF
        A_hat(e,k) =  max(A_l(e,k), A_r(e,k))
      end do
      end do

      ! apply fluxes ...........................................................

      g = 1 / dx

      !$omp do
      do e = 1, ne
        if (activity(e) < 1) cycle

        do k = 1, 3

          ! r_d += {v'A}[u]  . . . . . . . . . . . . . . . . . . . . . . . . . .

          c_0  = g * ( A( 0,e,k,1) * jmp_u(e-1,1) &
                     + A( 0,e,k,2) * jmp_u(e-1,2) &
                     + A( 0,e,k,3) * jmp_u(e-1,3) )

          c_po = g * ( A(po,e,k,1) * jmp_u(e  ,1) &
                     + A(po,e,k,2) * jmp_u(e  ,2) &
                     + A(po,e,k,3) * jmp_u(e  ,3) )

          r_d(:,e,k) = r_d(:,e,k) + D(0,i) * c_0 + D(po,:) * c_po

          ! r_d += [v]{Au'}  . . . . . . . . . . . . . . . . . . . . . . . . . .

          r_d( 0,e,k) = r_d( 0,e,k) - avg_q(e-1,k)
          r_d(po,e,k) = r_d(po,e,k) + avg_q(  e,k)

          ! r_d -= μ[v]Â[u]  . . . . . . . . . . . . . . . . . . . . . . . . . .

          r_d( 0,e,k) = r_d( 0,e,k) + mu * A_hat(e-1,k) * jmp_u(e-1,k)
          r_d(po,e,k) = r_d(po,e,k) - mu * A_hat(  e,k) * jmp_u(e  ,k)

        end do
      end do

      ! finalization ...........................................................

      ! free shared workspace
      deallocate(A, A_hat, A_l, A_r, u_l, u_r, q_l, q_r, jmp_u, avg_q)

    end associate

    !$omp end master

  end subroutine GetHybridDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Implicit CNS diffusion solver
  !>
  !> FGMRES: Van der Vorst, Fig. 6.4

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

    ! internal variables .......................................................

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3

    real(RNP), allocatable, save :: h(:,:)     ! Hessenberg matrix
    real(RNP), allocatable, save :: r(:,:,:)   ! residual
    real(RNP), allocatable, save :: v(:,:,:,:) !
    real(RNP), allocatable, save :: w(:,:,:)   ! ...
    real(RNP), allocatable, save :: y(:)       ! ...
    real(RNP), allocatable, save :: z(:,:,:,:) ! ...
    real(RNP) :: beta

    ! auxiliary variables for solving the least-squares problem
    real(RNP), allocatable, save :: g(:,:) ! r       in VDV03
    real(RNP), allocatable, save :: b(:)   ! \hat{b} in VDV03
    real(RNP), allocatable, save :: c(:)   ! c       in VDV03
    real(RNP), allocatable, save :: s(:)   ! s       in VDV03
    real(RNP) :: delta, gamma, rho

    real(RNP), save :: r_term
    logical  , save :: converged

    character :: comp
    logical   :: check_convergence
    integer   :: e, i, j, k

    !$omp master !!! so far

    associate( nc => this % nc              &
             , nk => this % n_krylov        &
             , ne => cl_operator % ne       &
             , po => cl_operator % eop % po &
             , Me => cl_operator % Me       )

      ! initialization .........................................................

      if (theta > 0) then
        ! physical and streamline (total) diffusion
        comp = 'T'
      else
        ! physical diffusion only
        comp = 'P'
      end if

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      !$omp master
      allocate(h(nk+1,nk), source = ZERO)
      allocate(g(nk  ,nk), source = ZERO)
      allocate(b(nk+1)   , source = ZERO)
      allocate(c(nk), s(nk), y(nk))
      allocate(r, w, mold = u)
      allocate(v(0:po,ne,nc,nk))
      allocate(z, mold = v)
      !$omp end master
      !$omp barrier

      OUTER_ITERATION: do k = 1, i_max

        ! initial residual, r = f - Au .........................................

        call GetHybridDiffusionTerm( this, cl_operator, comp, theta &
                                   , bv, u_0, u, r                  )

        !$omp do collapse(2)
        do i = 1, nc
        do e = 1, ne
          r(:,e,i) = Me * (f(:,e,i) - 1/dt * u(:,e,i)) + r(:,e,i)
        end do
        end do

        beta = sqrt(ScalarProduct(r, r))
        b(1) = beta

        ! convergence check ....................................................

        if (check_convergence) then

          !$omp single
          if (k == 1) then
            ! set terminal condition
            r_term = huge(r_term)
            if (present(r_max)) then
              if (r_max > 0) r_term = r_max
            end if
            if (present(r_red)) then
              if (r_red > 0) r_term = min(r_term, beta * r_red)
            end if
          end if
          converged = beta <= r_term
          !$omp end single

          if (converged) exit OUTER_ITERATION

        end if

        ! first Krylov vector ..................................................

        call MergeArrays(ZERO, v(:,:,:,1), 1/max(beta,eps), r, multi=.true.)

        INNER_ITERATION: do j = 1, nk
          associate(zj => z(:,:,:,j))

            ! preconditioning, z(j) = K⁻¹v(j) ..................................

            ! using identity, so far
            zj = v(:,:,:,j)

            ! application of homogeneous operator, w = A z(j) ..................

            call GetHybridDiffusionTerm( this, cl_operator, comp, theta &
                                       , u_0 = u_0, u = zj, r_d = w     )

            call MergeArrays(-ONE, w, 1/dt, zj, multi=.true.)

            !$omp do collapse(2)
            do i = 1, nc
            do e = 1, ne
              w(:,e,i) = 1/dt * Me * u(:,e,i) - w(:,e,i)
            end do
            end do

            ! computation of new Krylov vector .................................

            ! orthogonalization against old Krylov vectors
            do i = 1, j
              h(i,j) = ScalarProduct(w, v(:,:,:,i))
              call MergeArrays(ONE, w, h(i,j), v(:,:,:,i), multi=.true.)
            end do

            ! normalization, v(j+1) = w / ‖w‖
            h(j+1,j) = sqrt(ScalarProduct(w, w))
            call MergeArrays(ZERO, v(:,:,:,j+1), 1/h(j+1,j), w, multi=.true.)

            ! Givens rotation transforming h to upper triagonal matrix g .......

            g(1,j) = h(1,j)

            do i = 2, j
              gamma    =  c(i-1) * g(i-1,j) + s(i-1) * h(i,j)
              g(i,j)   = -s(i-1) * g(i-1,j) + c(i-1) * h(i,j)
              g(i-1,j) = gamma
            end do

            delta  =  max(sqrt(g(j,j)**2 + h(j+1,j)**2), eps)
            c(j)   =  g(j  ,j) / delta
            s(j)   =  h(j+1,j) / delta
            g(j,j) =  c(j) * g(j,j) + s(j) * h(j+1,j)
            b(j+1) = -s(j) * b(j)
            b(j)   =  c(j) * b(j)

            ! convergence test .................................................

            ! residual norm if inner iterations were exited now
            rho = abs(b(j+1))

            if (check_convergence) then
              converged = rho <= r_term
            end if

            if (converged) exit INNER_ITERATION

          end associate
        end do INNER_ITERATION

        ! intermediate solution ................................................

        j = min(j, nk)

        ! solve least-squares problem for y using backward substitution
        y(j) = b(j) / g(j,j)
        do i = j-1, 1, -1
          y(i) = (b(i) - dot_product(g(i,i+1:j), y(i+1:j))) / g(i,i)
        end do

        ! improved approximate solution
        do i = 1, j
          call MergeArrays(ONE, u, y(i), z(:,:,:,i), multi=.true.)
        end do

      end do OUTER_ITERATION

      ! finalization ...........................................................

      !$omp master
      deallocate(b, c, g, h, r, s, v, w, y, z)
      !$omp end master

    end associate

    !$omp end master

  end subroutine DiffusionSolver

  !=============================================================================

end submodule SM_Diffusion
