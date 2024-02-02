submodule(CL__Problem__CNS__1D) SM_Diffusion
  use Array_Assignments
  use Array_Reductions
! use, intrinsic :: IEEE_Arithmetic
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Returns the physical diffusivity matrix
  !>
  !> If present, `eta` overrides `this%eta`. Specification of `prandtl` implies
  !> the recomputation of `lamdda`as  well.

  pure module function PhysicalDiffusivity(this, u, eta, prandtl) result(A_d)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP),           intent(in) :: u(:)    !< conservative variables
    real(RNP), optional, intent(in) :: eta     !< dynamic viscosity
    real(RNP), optional, intent(in) :: prandtl !< prandtl number
    real(RNP) :: A_d(3,3)

    real(RNP) :: eta_, lambda_
    real(RNP) :: amv, aev, aeT, c1, cc, dv1, dv2, dT1, dT2, dT3

    if (present(eta)) then
      eta_ = eta
    else
      eta_ = this % eta
    end if

    if (present(prandtl)) then
      lambda_ = eta_ * this%c_p / prandtl
    else
      lambda_ = this % lambda
    end if

    c1  =  1  / u(1)          ! 1/ρ
    cc  =  c1 / this % c_v    ! 1/ρc_v

    amv =  4 * THIRD * eta_   !  momentum diffusivity related to ∂v/∂x
    aev =  amv * c1 * u(2)    !  energy   diffusivity related to ∂v/∂x
    aeT =  lambda_            !  energy   diffusivity related to ∂T/∂x

    dv1 = -c1 * u(2)          ! ∂v/∂u₁
    dv2 =  c1                 ! ∂v/∂u₂

    dT1 =  cc * c1**2 * (u(2)**2 - u(1)*u(3))  ! ∂T/∂u₁
    dT2 = -cc * c1 * u(2)                      ! ∂T/∂u₂
    dT3 =  cc                                  ! ∂T/∂u₃

    A_d(1,:) = [ ZERO                  , ZERO                  , ZERO      ]
    A_d(2,:) = [ amv * dv1             , amv * dv2             , ZERO      ]
    A_d(3,:) = [ aev * dv1 + aeT * dT1 , aev * dv2 + aeT * dT2 , aeT * dT3 ]

  end function PhysicalDiffusivity

  !-----------------------------------------------------------------------------
  !> Unified CNS diffusion term combining physical and streamline contributions
  !>
  !> The composition of the diffusion term is controlled by the flags contained
  !> in the argument `comp`:
  !>
  !>   - `P`  physical
  !>   - `D`  discontinuity capturing
  !>   - `S`  streamline
  !>   - `T`  total
  !>
  !> Several flags can be specified, e.g. 'PS' combines physical and streamline
  !> diffusion. Flag `T` yields the total diffusion resulting from adding all
  !> contributions.
  !>
  !> Homogeneous boundary conditions are applied in the case that `bv` omitted.

  module subroutine GetHybridDiffusionTerm &
      (this, cl_operator, comp, theta, bv, u_0, u, r_d)

    class(CL_Problem_CNS_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator !< spatial operators
    character(len=*),      intent(in)  :: comp        !< composition flag
    real(RNP),             intent(in)  :: theta       !< SD time scale
    real(RNP), optional,   intent(in)  :: bv(:,:)     !< boundary values
    real(RNP), contiguous, intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
    real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
    real(RNP), contiguous, intent(out) :: r_d(0:,:,:) !< diffusion RHS

    ! local variables ..........................................................

    real(RNP), allocatable, save :: A(:,:,:,:)
    ! element diffusivity matrices, A(0:po,1:ne,1:3,1:3)

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
    logical :: has_dc ! switch for physical diffusion
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

      has_pd = scan(comp,'PT') > 0 .and. 0 < this % eta
      has_dc = scan(comp,'DT') > 0 .and. 0 < this % dc_diffusivity
      has_sd = scan(comp,'ST') > 0 .and. 0 < theta

      ! shared workspace
      allocate(A(0:po,ne,3,3), source = ZERO)
      allocate(A_hat(0:ne,3),  source = ZERO)
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

      call GetHybridDiffusity(this, cl_operator, comp, theta, u_0, A)

      !$omp do
      do e = 1, ne

        select case(activity(e))

        case(0:)
          do k = 1, 3
            A_r(e-1,k) = A( 0,e,k,k)
            A_l(e  ,k) = A(po,e,k,k)
          end do

        case default
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

      case('D','S')
        u_l(0,1:3) = -u_r(0,1:3)
        q_l(0,1:3) =  q_r(0,1:3)
        A_l(0,1:3) =  A_r(0,1:3)
        if (present(bv)) then
          u_l(0,1:3) = u_l(0,1:3) + 2*bv(:,1)
        end if

      case default ! extrapolation
        u_l(0,1:3) =  u_r(0,1:3)
        q_l(0,1:3) = -q_r(0,1:3)
        A_l(0,1:3) =  A_r(0,1:3)

      end select

      ! right boundary
      select case(this % bc(2))

      case('P') ! periodic
        u_r(ne,1:3) = u_r(0,1:3)
        q_r(ne,1:3) = q_r(0,1:3)
        A_r(ne,1:3) = A_r(0,1:3)

      case('D','S')
        u_r(ne,1:3) = -u_l(ne,1:3)
        q_r(ne,1:3) =  q_l(ne,1:3)
        A_r(ne,1:3) =  A_l(ne,1:3)
        if (present(bv)) then
          u_r(ne,1:3) = u_r(ne,1:3) + 2*bv(:,2)
        end if

      case default ! extrapolation
        u_r(ne,1:3) =  u_l(ne,1:3)
        q_r(ne,1:3) = -q_l(ne,1:3)
        A_r(ne,1:3) =  A_l(ne,1:3)

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

          r_d(:,e,k) = r_d(:,e,k) + D(0,:) * c_0 + D(po,:) * c_po

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
  !> Composition of hybrid diffusivity

  subroutine GetHybridDiffusity(this, cl_operator, comp, theta, u, A_d)
    class(CL_Problem_CNS_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator   !< spatial operators
    character(len=*),      intent(in)  :: comp          !< composition flags
    real(RNP),             intent(in)  :: theta         !< SD time scale
    real(RNP), contiguous, intent(in)  :: u(0:,:,:)     !< u(x,t)
    real(RNP), contiguous, intent(out) :: A_d(0:,:,:,:) !< diffusivity tensor

    logical :: has_pd ! switch for physical diffusion
    logical :: has_dc ! switch for discontinuity capturing diffusion
    logical :: has_sd ! switch for streamline diffusion
    logical :: has_ad ! switch for artificial diffusion

    ! discontinuity capturing
    real(RNP), allocatable :: V_inv_dc(:,:) ! inverse Vandermonde matrix
    real(RNP) :: c_dc, d_dc, s_dc           ! parameters

    ! artificial diffusivity filtering
    real(RNP), allocatable :: Q_ad(:,:)     ! filtering operator

    ! element variables
    real(RNP), allocatable :: A_de(:,:,:), u_e(:,:), v_e(:), eta_e(:)

    real(RNP) :: A_c(3,3)
    integer   :: po_cut
    integer   :: e, i, j, k

    associate( dx => cl_operator % dx       &
             , po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      ! initialization .........................................................

      has_pd = scan(comp,'PT') > 0 .and. 0 < this % eta
      has_dc = scan(comp,'DT') > 0 .and. 0 < this % dc_diffusivity
      has_sd = scan(comp,'ST') > 0 .and. 0 < theta
      has_ad = has_dc .or. has_sd

      ! discontinuity capturing
      if (has_dc) then

        c_dc = this % dc_scaling_coeff * dx / po
        d_dc = this % dc_sensor_delta
        s_dc = this % dc_sensor_coeff * log10(real(po,RNP)) &
             + this % dc_sensor_const

        select case(this % dc_diffusivity)
        case(1)
          allocate(V_inv_dc(0:po,0:po))
          call cl_operator % eop % Get_Inverse_Legendre_VDM(V_inv_dc)
        end select

      end if

      ! artificial diffusivity filtering
      if (has_ad) then

        select case(this % ad_filter_method)
        case(1)
          po_cut = this % ad_filter_degree
        case(2)
          po_cut = po / 2
        case default
          po_cut = -1
        end select

        if (po_cut >= 0)  then
          allocate(Q_ad(0:po,0:po), source = ZERO)
          select case(this % ad_filter_basis)
          case(1)
            call cl_operator % eop % Get_Legendre_CutoffFilter(po_cut, Q_ad)
          case(2)
            call cl_operator % eop % Get_Bubble_CutoffFilter(po_cut, Q_ad)
          end select
        end if

      end if

      ! element variables
      allocate(A_de(0:po,3,3), u_e(0:po,3), v_e(0:po), eta_e(0:po))

      ! composition of hybrid diffusivity ......................................

      do e = 1, ne

        A_de = 0

        if (cl_operator % activity(e) >= 0) then

          u_e = u(0:po,e,1:3)

          ! streamline diffusivity
          if (has_sd) then

            ! streamline diffusivity tensor
            do i = 0, po
              A_c = ConvectiveJacobian(this, u(i,e,1:3))
              A_de(i,1:3,1:3) = theta/2 * matmul(A_c, A_c)
            end do

            ! filtering
            select case(this % ad_filter_method)
            case(1:2)
              do k = 1, 3
                A_de(:,:,k) = matmul(Q_ad, A_de(:,:,k))
              end do
            case(3)
              do k = 1, 3
              do j = 1, 3
                if (j == k) then
                  A_de(:,j,k) = maxval(A_de(:,j,k))
                else
                  A_de(:,j,k) = 0
                end if
              end do
              end do
            end select

          end if

          ! discontinuity capturing diffusivity
          if (has_dc) then

            ! sensor variable
            select case(this % dc_sensor_var)
            case(1)
              ! density
              v_e = u_e(:,1)
            case(2)
              ! entropy
              do i = 0, po
                call ConservativeToPrimitive(this, u_e(i,:), s = v_e(i))
              end do
            end select

            ! dynamic viscosity
            select case(this % dc_diffusivity)
            case(1)
              eta_e = PerssonViscosity(this, c_dc, d_dc, s_dc, v_e, u_e, V_inv_dc)
            case(2)
              ! ...
              ! filtering
              select case(this % ad_filter_method)
              case(1:2)
                eta_e = matmul(Q_ad, eta_e)
              case(3)
                eta_e = maxval(eta_e)
              end select
            end select

            ! diffusivity tensor
            do i = 0, po
              A_de(i,1:3,1:3) = A_de(i,1:3,1:3) &
                              + PhysicalDiffusivity(this, u_e(i,1:3), eta_e(i))
            end do

          end if

          ! physical diffusivity
          if (has_pd) then
            do i = 0, po
              A_de(i,1:3,1:3) = A_de(i,1:3,1:3) &
                              + PhysicalDiffusivity(this, u_e(i,1:3))
            end do
          end if

        end if

        A_d(0:po,e,1:3,1:3) = A_de

      end do

    end associate

  end subroutine GetHybridDiffusity

  !-----------------------------------------------------------------------------
  !> Discontinuity capturing viscosity in spirit of Persson & Preraire (2006)


  pure function PerssonViscosity(this, c_max, delta_s, s_ref, v, u, V_inv) &
      result(eta)

    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: c_max        !< max diffusivity scaling factor
    real(RNP), intent(in) :: delta_s      !< sensor half width
    real(RNP), intent(in) :: s_ref        !< sensor reference value
    real(RNP), intent(in) :: v(0:)        !< sensor variable
    real(RNP), intent(in) :: u(0:,:)      !< conservation variables
    real(RNP), intent(in) :: V_inv(0:,0:) !< inverse Vandermonde matrix
    real(RNP) :: eta

    real(RNP), parameter :: eps = epsilon(ONE)
    real(RNP) :: v_hat(0:ubound(u,1))
    real(RNP) :: lambda(3)
    real(RNP) :: norm_v, eta_max, s
    integer   :: i, po

    po = ubound(u,1)

    ! sensoring
    v_hat  = matmul(V_inv, v)
    norm_v = 0
    do i = 0, po
      norm_v = norm_v + v_hat(i)**2 / (2*i + 1)
    end do
    s = (v_hat(po)**2 / (2*po + 1)) / max(norm_v, eps)
    s = log10(max(s, eps))

    ! evaluation
    if (s <= s_ref - delta_s) then

      eta = 0

    else

      ! maximum artificial viscosity: η_max = c_max max(ρ|Λ(u)|)
      eta_max = 0
      do i = 0, po
        call ConvectiveEigensystem(this, u(i,:), lambda)
        eta_max = max(eta_max, u(i,1) * maxval(abs(lambda)))
      end do
      eta_max = c_max * eta_max

      if (s > s_ref + delta_s) then
        eta = eta_max
      else
        eta = eta_max/2 * (1 + sin(PI * (s - s_ref) / (2*delta_s)))
      end if

    end if

  end function PerssonViscosity

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
    real(RNP), allocatable, save :: v(:,:,:,:)
    real(RNP), allocatable, save :: w(:,:,:)
    real(RNP), allocatable, save :: y(:)
    real(RNP), allocatable, save :: z(:,:,:,:)
    real(RNP) :: beta

    ! auxiliary variables for solving the least-squares problem
    real(RNP), allocatable, save :: r(:,:) ! r       in VDV03
    real(RNP), allocatable, save :: b(:)   ! \hat{b} in VDV03
    real(RNP), allocatable, save :: c(:)   ! c       in VDV03
    real(RNP), allocatable, save :: s(:)   ! s       in VDV03
    real(RNP) :: delta, gamma, rho

    real(RNP), save :: r_term
    logical  , save :: converged

    logical   :: check_convergence
    integer   :: e, i, j, k
    integer   :: ic

    !$omp master !!! so far

    associate( nc => this % nc              &
             , nk => this % n_krylov        &
             , ne => cl_operator % ne       &
             , po => cl_operator % eop % po &
             , Me => cl_operator % Me       )

      ! initialization .........................................................

      if (theta > 0) then
        ! full diffusivity tensor, including first component
        ic   =  1
      else
        ! no difffusion in first component (density)
        ic   =  2
      end if

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      !$omp master
      allocate(v(0:po,ne,nc,nk+1)           , source = ZERO)
      allocate(w(0:po,ne,nc)                , source = ZERO)
      allocate(z(0:po,ne,nc,nk)             , source = ZERO)
      allocate(b(nk+1), c(nk), s(nk), y(nk) , source = ZERO)
      allocate(r(nk,nk), h(nk+1,nk)         , source = ZERO)
      !$omp end master
      !$omp barrier

      OUTER_ITERATION: do k = 1, i_max

        ! initial residual, v₁ = f - Au ........................................

        associate(v1 => v(:,:,:,1))

          call GetHybridDiffusionTerm( this, cl_operator, 'T', theta &
                                     , bv, u_0, u, r_d = v1          )

          !$omp do collapse(2)
          do i = ic, nc
          do e = 1, ne
            v1(:,e,i) = 1/dt * Me * (f(:,e,i) - u(:,e,i)) + v1(:,e,i)
          end do
          end do

          beta = sqrt(ScalarProduct(v1, v1))
          b(1) = beta

          if (log_level_outer_iteration > 0) then
            write(*,'(2X,A,I4,A,ES12.5)') &
              'CNS DiffusionSolver, outer iteration k =',k,': beta =', beta
          end if

          ! convergence check ..................................................

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

          ! first Krylov vector ................................................

          call ScaleArray(v1, 1/max(beta,eps), multi=.true.)

        end associate

        INNER_ITERATION: do j = 1, nk
          associate(zj => z(:,:,:,j))

            ! preconditioning, z(j) = K⁻¹v(j) ..................................

            ! using identity, so far
            zj = v(:,:,:,j)

            ! application of homogeneous operator, w = A z(j) ..................

            call GetHybridDiffusionTerm( this, cl_operator, 'T', theta &
                                       , u_0 = u_0, u = zj, r_d = w    )

            !$omp do collapse(2)
            do i = ic, nc
            do e = 1, ne
              w(:,e,i) = 1/dt * Me * zj(:,e,i) - w(:,e,i)
            end do
            end do

            ! computation of new Krylov vector .................................

            ! orthogonalization against old Krylov vectors
            do i = 1, j
              h(i,j) = ScalarProduct(v(:,:,:,i), w)
              call MergeArrays(ONE, w, -h(i,j), v(:,:,:,i), multi=.true.)
            end do

            ! normalization, v(j+1) = w / ‖w‖
            h(j+1,j) = sqrt(ScalarProduct(w, w))
            call MergeArrays(ZERO, v(:,:,:,j+1), 1/h(j+1,j), w, multi=.true.)

            ! Givens rotation transforming h to upper triagonal matrix r .......

            r(1,j) = h(1,j)

            do i = 2, j
              gamma    =  c(i-1) * r(i-1,j) + s(i-1) * h(i,j)
              r(i,j)   = -s(i-1) * r(i-1,j) + c(i-1) * h(i,j)
              r(i-1,j) = gamma
            end do

            delta  =  max(sqrt(r(j,j)**2 + h(j+1,j)**2), eps)
            c(j)   =  r(j  ,j) / delta
            s(j)   =  h(j+1,j) / delta
            r(j,j) =  c(j) * r(j,j) + s(j) * h(j+1,j)
            b(j+1) = -s(j) * b(j)
            b(j)   =  c(j) * b(j)

            ! convergence test .................................................

            ! residual norm if inner iterations were exited now
            rho = abs(b(j+1))

            if (log_level_inner_iteration > 0) then
              write(*,'(2X,A,I4,A,ES12.5)') &
                'CNS DiffusionSolver, inner iteration j =',j,': rho  =', rho
            end if

            if (check_convergence) then
              converged = rho <= r_term
            end if

            if (converged) exit INNER_ITERATION

          end associate
        end do INNER_ITERATION

        ! intermediate solution ................................................

        j = min(j, nk)

        ! solve least-squares problem for y using backward substitution
        y(j) = b(j) / r(j,j)
        do i = j-1, 1, -1
          y(i) = (b(i) - dot_product(r(i,i+1:j), y(i+1:j))) / r(i,i)
        end do

        ! improved approximate solution
        do i = 1, j
          call MergeArrays(ONE, u, y(i), z(:,:,:,i), multi=.true.)
        end do

      end do OUTER_ITERATION

      ! finalization ...........................................................

      !$omp master
      deallocate(b, c, h, r, s, v, w, y, z)
      !$omp end master

    end associate

    !$omp end master

  end subroutine DiffusionSolver

  !=============================================================================

end submodule SM_Diffusion
