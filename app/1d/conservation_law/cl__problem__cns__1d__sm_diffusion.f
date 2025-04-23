submodule(CL__Problem__CNS__1D) SM_Diffusion
  use Array_Assignments
  use Array_Reductions
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Returns the physical diffusivity matrix
  !>
  !> If present, `eta` overrides `this%eta`. Specification of `prandtl` implies
  !> the recomputation of `lamdda`as  well.

  pure function PhysicalDiffusivity(this, u, eta, prandtl) result(A_d)
    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP),           intent(in) :: u(:)    !< conservative variables
    real(RNP), optional, intent(in) :: eta     !< dynamic viscosity
    real(RNP), optional, intent(in) :: prandtl !< prandtl number
    real(RNP) :: A_d(3,3)

    real(RNP) :: eta_, lambda_
    real(RNP) :: amv, aev, aeT, c1, cc, dv1, dv2, dT1, dT2, dT3, v

    if (present(eta)) then
      eta_ = eta
    else
      eta_ = this % eta
    end if

    if (present(prandtl)) then
      lambda_ = eta_ * this%c_p / prandtl
    else
      lambda_ = eta_ * this%c_p / this % prandtl
    end if

    c1  =  1  / u(1)              ! 1/ρ
    cc  =  c1 / this % c_v        ! 1/ρc_v
    v   =  c1 * u(2)              ! v

    amv =  4 * THIRD * eta_       ! momentum diffusivity related to ∂v/∂x
    aev =  amv * v                ! energy   diffusivity related to ∂v/∂x
    aeT =  lambda_                ! energy   diffusivity related to ∂T/∂x

    dv1 = -c1 * v                 ! ∂v/∂u₁
    dv2 =  c1                     ! ∂v/∂u₂

    dT1 =  cc * (v**2 - c1*u(3))  ! ∂T/∂u₁
    dT2 = -cc *  v                ! ∂T/∂u₂
    dT3 =  cc                     ! ∂T/∂u₃

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

    real(RNP), allocatable, save :: A_hat(:,:,:)
    ! interface diffusivity matrices Â (0:ne,1:3.1:3)

    real(RNP), allocatable, save :: jmp_u(:,:)
    ! jump of solution across element boundaries (0:ne,1:3)

    real(RNP), allocatable, save :: avg_q(:,:)
    ! average of diffusive fluxes over element boundaries (0:ne,1:3)

    real(RNP), allocatable, save :: u_l(:,:), u_r(:,:)
    ! left and right traces of solution at element interfaces (0:ne,1:3)

    real(RNP), allocatable, save :: q_l(:,:), q_r(:,:)
    ! left and right diffusive fluxes at element interfaces (0:ne,1:3)

    real(RNP), allocatable, save :: A_l(:,:,:), A_r(:,:,:)
    ! left and right traces of diffusivity matrices (0:ne,1:3,1:3)

    real(RNP), allocatable :: MD_t(:,:)
    ! transpose of weighted derivative matrix (MD)ᵗ

    real(RNP), allocatable :: dx_u(:,:)
    ! element solution derivatives ∂u/∂x (0:po,1:3)

    real(RNP), allocatable :: q(:,:)
    ! element diffusive flux

    real(RNP) :: g, mu
    real(RNP) :: c_0, c_po

    integer :: e, i, j, k

    !!$omp master !!! not ready for OpenMP, enforcing sequential execution !!!

    associate( activity => cl_operator % activity &
             , eop      => cl_operator % eop      &
             , M        => cl_operator % eop % w  &
             , D        => cl_operator % eop % D  &
             , po       => cl_operator % eop % po &
             , ne       => cl_operator % ne       &
             , dx       => cl_operator % dx       )

      ! initialization .........................................................

      ! shared workspace
      allocate(A(0:po,ne,3,3),  source = ZERO)
      allocate(A_hat(0:ne,3,3), source = ZERO)
      allocate(A_l, A_r, mold = A_hat)

      allocate(avg_q(0:ne,3), source = ZERO)
      allocate(u_l, u_r, q_l, q_r, jmp_u, mold = avg_q)

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

      !!$omp do
      do e = 1, ne

        select case(activity(e))

        case(0:)
          A_r(e-1,1:3,1:3) = A( 0,e,:,:)
          A_l(e  ,1:3,1:3) = A(po,e,:,:)

        case default
          A_r(e-1,1:3,1:3) = 0
          A_l(e  ,1:3,1:3) = 0

        end select
      end do

      ! application of interior operator .......................................

      g = 2 / dx

      !!$omp do
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
        u_l(0,1:3)     = u_l(ne,1:3)
        q_l(0,1:3)     = q_l(ne,1:3)
        A_l(0,1:3,1:3) = A_l(ne,1:3,1:3)

      case('D','S')
        u_l(0,1:3)     = -u_r(0,1:3)
        q_l(0,1:3)     =  q_r(0,1:3)
        A_l(0,1:3,1:3) =  A_r(0,1:3,1:3)
        if (present(bv)) then
          u_l(0,1:3) = u_l(0,1:3) + 2*bv(:,1)
        end if

      case default ! extrapolation
        u_l(0,1:3)     =  u_r(0,1:3)
        q_l(0,1:3)     = -q_r(0,1:3)
        A_l(0,1:3,1:3) =  A_r(0,1:3,1:3)

      end select

      ! right boundary
      select case(this % bc(2))

      case('P') ! periodic
        u_r(ne,1:3)     = u_r(0,1:3)
        q_r(ne,1:3)     = q_r(0,1:3)
        A_r(ne,1:3,1:3) = A_r(0,1:3,1:3)

      case('D','S')
        u_r(ne,1:3)     = -u_l(ne,1:3)
        q_r(ne,1:3)     =  q_l(ne,1:3)
        A_r(ne,1:3,1:3) =  A_l(ne,1:3,1:3)
        if (present(bv)) then
          u_r(ne,1:3) = u_r(ne,1:3) + 2*bv(:,2)
        end if

      case default ! extrapolation
        u_r(ne,1:3)     =  u_l(ne,1:3)
        q_r(ne,1:3)     = -q_l(ne,1:3)
        A_r(ne,1:3,1:3) =  A_l(ne,1:3,1:3)

      end select

      ! element interface values ...............................................

      !!$omp do collapse(2)
      do e = 0, ne
        do k = 1, 3
          jmp_u(e,k) =  u_l(e,k) - u_r(e,k)
          avg_q(e,k) = (q_l(e,k) + q_r(e,k)) * HALF
          do i = 1, 3
            if (i == k) then
              A_hat(e,i,k) =  max(A_l(e,i,k), A_r(e,i,k))
            else
              A_hat(e,i,k) =  (A_l(e,i,k) + A_r(e,i,k)) * HALF
            end if
          end do
        end do
      end do

      ! apply fluxes ...........................................................

      g = 1 / dx

      !!$omp do
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

          r_d( 0,e,k) = r_d( 0,e,k)                          &
                      + mu * ( A_hat(e-1,k,1) * jmp_u(e-1,1) &
                             + A_hat(e-1,k,2) * jmp_u(e-1,2) &
                             + A_hat(e-1,k,3) * jmp_u(e-1,3) )

          r_d(po,e,k) = r_d(po,e,k)                      &
                      - mu * ( A_hat(e,k,1) * jmp_u(e,1) &
                             + A_hat(e,k,2) * jmp_u(e,2) &
                             + A_hat(e,k,3) * jmp_u(e,3) )

        end do
      end do

      ! finalization ...........................................................

      ! free shared workspace
      deallocate(A, A_hat, A_l, A_r, u_l, u_r, q_l, q_r, jmp_u, avg_q)

    end associate

    !!$omp end master

  end subroutine GetHybridDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Composition of hybrid diffusivity

  subroutine GetHybridDiffusity &
      (this, cl_operator, comp, theta, u, A_d, lmb_d, vr_d, vl_d)

    ! arguments ................................................................

    class(CL_Problem_CNS_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator   !< spatial operators
    character(len=*),      intent(in)  :: comp          !< composition flags
    real(RNP),             intent(in)  :: theta         !< SD time scale
    real(RNP), contiguous, intent(in)  :: u(0:,:,:)     !< u(x,t)
    real(RNP), contiguous, intent(out) :: A_d(0:,:,:,:) !< diffusivity tensor
    real(RNP), contiguous, intent(out) :: lmb_d(:,:)    !< eigenvalues
    real(RNP), contiguous, intent(out) :: vr_d(:,:,:)   !< right eigenvectors
    real(RNP), contiguous, intent(out) :: vl_d(:,:,:)   !< left eigenvectors

    optional :: A_d, lmb_d, vr_d, vl_d

    ! internal variables .......................................................

    logical :: has_pd ! switch for physical diffusion
    logical :: has_dc ! switch for discontinuity capturing diffusion
    logical :: has_sd ! switch for streamline diffusion
    logical :: has_ad ! switch for artificial diffusion

    integer :: precon

    real(RNP), allocatable :: VL_inv(:,:) ! inverse Legendre Vandermonde matrix
    real(RNP), allocatable :: Q_ad(:,:)   ! artificial diffusivity filter

    ! discontinuity capturing parameters
    real(RNP) :: c_dc, d_dc, s_dc

    ! element variables
    real(RNP), allocatable :: A_de(:,:,:), u_e(:,:), q_e(:)

    real(RNP) :: A_c(3,3)
    real(RNP) :: u_0(3), lmb_0(3), vl_0(3,3), vr_0(3,3)
    real(RNP) :: nu_0, kappa_0, T_0, v_0, vv_0
    real(RNP) :: a31_0, r31_0
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

      ! identify preconditioner
      if (present(lmb_d) .and. present(vr_d) .and. present(vl_d)) then
        precon = this % precon
      else
        precon = 0
      end if

      ! Vandermonde matrix
      allocate(VL_inv(0:po,0:po))
      call cl_operator % eop % Get_Inverse_Legendre_VDM(VL_inv)

      ! discontinuity capturing
      if (has_dc) then
        c_dc = this % dc_scaling_coeff * dx / po
        d_dc = this % dc_sensor_delta
        s_dc = this % dc_sensor_coeff * log10(real(po,RNP)) &
             + this % dc_sensor_const
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
      allocate(A_de(0:po,3,3), u_e(0:po,3), q_e(0:po))

      ! composition of hybrid diffusivity ......................................

      do e = 1, ne

        A_de = 0

        if (precon > 0) then
          lmb_0 = 0
          vr_0 = reshape([ 1,0,0, 0,1,0, 0,0,1 ], [3,3])
          vl_0 = vr_0
        end if

        if (cl_operator % activity(e) >= 0) then

          u_e = u(0:po,e,1:3)             ! element variables
          u_0 = matmul(VL_inv(0,:), u_e)  ! mean solution values

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

            ! contribution to preconditioning
            select case(precon)
            case(1)
              ! mean diagonal values
              do k = 1, 3
                lmb_0(k) = dot_product(VL_inv(0,:), A_de(:,k,k))
              end do
            case(2)
              ! physical diffusion eigensystem ← average of mean diagonal values
              lmb_0 = THIRD * dot_product(VL_inv(0,:), A_de(:,1,1) &
                                                     + A_de(:,2,2) &
                                                     + A_de(:,3,3) )
            case(3)
              ! streamline diffusion eigensystem
              call ConvectiveEigensystem(this, u_0, lmb_0, vr_0, vl_0)
              lmb_0 = theta/2 * lmb_0**2
            end select

          end if

          ! discontinuity capturing diffusivity
          if (has_dc) then

            ! sensor variable
            select case(this % dc_sensor_var)
            case(1)
              ! density
              q_e = u_e(:,1)
            case(2)
              ! entropy
              do i = 0, po
                call ConservativeToPrimitive(this, u_e(i,:), s = q_e(i))
              end do
            end select

            ! diffusivity tensor
            select case(this % dc_diffusivity)
            case(1)
              nu_0 = PerssonDiffusivity(this, c_dc, d_dc, s_dc, q_e, u_e, VL_inv)
              do i = 0, po
                A_de(i,1,1) = A_de(i,1,1) + nu_0
                A_de(i,2,2) = A_de(i,2,2) + nu_0
                A_de(i,3,3) = A_de(i,3,3) + nu_0
              end do
            end select

            ! contribution to preconditioning
            if (precon > 0) then
              lmb_0 = lmb_0 + nu_0
            end if

          end if

          ! physical diffusivity
          if (has_pd) then

            do i = 0, po
              A_de(i,1:3,1:3) = A_de(i,1:3,1:3) &
                              + PhysicalDiffusivity(this, u_e(i,1:3))
            end do

            ! contribution to preconditioning
            nu_0    = this%eta * 4 / (3 * u_0(1))
            kappa_0 = this%eta / (u_0(1) * this%prandtl)
            select case(precon)
            case(1)
              ! using mean diagonal values
              lmb_0 = lmb_0 + [ ZERO, nu_0, kappa_0 ]
            case(2)
              ! using physical diffusion eigensystem
              call ConservativeToPrimitive(this, u_0, v=v_0, T=T_0)
              ! auxiliary
              vv_0  = v_0 ** 2
              a31_0 = -nu_0 * vv_0 + kappa_0 * (vv_0/2 - this%c_v * T_0)
              r31_0 = -(a31_0 + (nu_0 - kappa_0) * vv_0) / kappa_0
              ! eigenvalues
              lmb_0 = lmb_0 + [ ZERO, nu_0, kappa_0 ]
              ! right eigenvectors
              vr_0(1,:) = [   ONE, ZERO, ZERO ]
              vr_0(2,:) = [   v_0,  ONE, ZERO ]
              vr_0(3,:) = [ r31_0,  v_0,  ONE ]
              ! left eigenvectors
              vl_0(1,:) = [          ONE, ZERO, ZERO ]
              vl_0(2,:) = [         -v_0,  ONE, ZERO ]
              vl_0(3,:) = [ vv_0 - r31_0, -v_0,  ONE ]
            case(3)
              ! using streamline diffusion eigensystem
              lmb_0 = lmb_0 + nu_0
            end select

          end if

        end if

        if (present(A_d)) then
          A_d(0:po,e,1:3,1:3) = A_de
        end if

        if (precon > 0) then
          lmb_d(e,1:3) = lmb_0
          vr_d(e,1:3,1:3) = vr_0
          vl_d(e,1:3,1:3) = vl_0
        end if

      end do

    end associate

  end subroutine GetHybridDiffusity

  !-----------------------------------------------------------------------------
  !> Discontinuity capturing viscosity in spirit of Persson & Preraire (2006)

  pure function PerssonDiffusivity(this, c_max, delta_s, s_ref, q, u, V_inv) &
      result(nu)

    class(CL_Problem_CNS_1D), intent(in) :: this
    real(RNP), intent(in) :: c_max        !< max diffusivity scaling factor
    real(RNP), intent(in) :: delta_s      !< sensor half width
    real(RNP), intent(in) :: s_ref        !< sensor reference value
    real(RNP), intent(in) :: q(0:)        !< sensor variable
    real(RNP), intent(in) :: u(0:,:)      !< conservation variables
    real(RNP), intent(in) :: V_inv(0:,0:) !< inverse Vandermonde matrix
    real(RNP) :: nu

    real(RNP), parameter :: eps = epsilon(ONE)
    real(RNP) :: q_hat(0:ubound(q,1))
    real(RNP) :: lambda(3)
    real(RNP) :: norm_q, nu_max, s
    integer   :: i, po

    po = ubound(u,1)

    ! sensoring
    q_hat  = matmul(V_inv, q)
    norm_q = 0
    do i = 0, po
      norm_q = norm_q + q_hat(i)**2 / (2*i + 1)
    end do
    s = (q_hat(po)**2 / (2*po + 1)) / max(norm_q, eps)
    s = log10(max(s, eps))

    ! evaluation
    if (s <= s_ref - delta_s) then

      nu = 0

    else

      ! maximum artificial viscosity: η_max = c_max max(ρ|Λ(u)|)
      nu_max = 0
      do i = 0, po
        call ConvectiveEigensystem(this, u(i,:), lambda)
        nu_max = max(nu_max, maxval(abs(lambda)))
      end do
      nu_max = c_max * nu_max

      if (s > s_ref + delta_s) then
        nu = nu_max
      else
        nu = nu_max/2 * (1 + sin(PI * (s - s_ref) / (2*delta_s)))
      end if

    end if

  end function PerssonDiffusivity

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

    ! auxiliary variables for preconditioning
    logical,   allocatable, save :: mask(:)     ! element mask
    integer,   allocatable, save :: cfg(:)      ! Schwarz domain configuration
    real(RNP), allocatable, save :: lmb_d(:,:)  ! diffusivity tensor eigenvalues
    real(RNP), allocatable, save :: vr_d(:,:,:) ! right eigenvectors
    real(RNP), allocatable, save :: vl_d(:,:,:) ! left eigenvectors

    real(RNP), save :: r_term
    logical  , save :: converged

    character :: bc(2)
    logical   :: check_convergence
    integer   :: e, i, j, k, ni, no
    integer   :: ic = 1

    !!$omp master !!! so far

    associate( nc => this % nc              &
             , ne => cl_operator % ne       &
             , po => cl_operator % eop % po &
             , Me => cl_operator % Me       )

      ! initialization .........................................................

      ! number of inner and outer iterations
      if (i_max > 0) then
        ni = min(this%n_krylov, i_max)
        no = i_max/ni + min(mod(i_max,ni), 1)
      else
        return
      end if

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      !!$omp master
      allocate(v(0:po,ne,nc,ni+1)           , source = ZERO)
      allocate(w(0:po,ne,nc)                , source = ZERO)
      allocate(z(0:po,ne,nc,ni)             , source = ZERO)
      allocate(b(ni+1), c(ni), s(ni), y(ni) , source = ZERO)
      allocate(r(ni,ni), h(ni+1,ni)         , source = ZERO)
      if (this % precon > 0) then
        allocate(mask(ne), source = cl_operator % activity > 0)
        allocate(cfg(ne), lmb_d(ne,nc), vr_d(ne,nc,nc), vl_d(ne,nc,nc))
      end if
      !!$omp end master
      !!$omp barrier

      ! setup of Schwarz preconditioner
      if (this % precon > 0) then

        ! boundary conditions
        do i = 1, 2
          select case(this % bc(i))
          case('P')
            bc(i) = 'P'
          case('D','S')
            bc(i) = 'D'
          case default
            bc(i) = 'N'
          end select
        end do

        ! domain configurations
        call cl_operator%elliptic_op%schwarz%ConfigureSubdomains(bc, cfg, mask)

      end if

      OUTER_ITERATION: do k = 1, no

        ! initial residual, v₁ = f - Au ........................................

        associate(v1 => v(:,:,:,1))

          call GetHybridDiffusionTerm( this, cl_operator, 'T', theta &
                                     , bv, u_0, u, r_d = v1          )

          if (this % precon > 0) then
            call GetHybridDiffusity( this, cl_operator, 'T', theta, u_0      &
                                   , lmb_d = lmb_d, vr_d = vr_d, vl_d = vl_d )
          end if

          !!$omp do collapse(2)
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

            !!$omp single
            if (k == 1) then
              ! set terminal condition
              r_term = huge(r_term)
              if (present(r_max)) then
                if (r_max > 0) r_term = r_max
              end if
              if (present(r_red)) then
                if (r_red > 0) r_term = max(r_term, beta * r_red)
              end if
            end if

            if (log_level > 0) then
              if (present(r_red)) print '(99(G0,X))', 'r_red  = ', r_red
              if (present(r_max)) print '(99(G0,X))', 'r_max  = ', r_max
              print '(99(G0,X))', 'r_term = ', r_term
            end if

            converged = beta <= r_term
            !!$omp end single

            if (converged) exit OUTER_ITERATION

          end if

          ! first Krylov vector ................................................

          call ScaleArray(v1, 1/max(beta,eps), multi=.true.)

        end associate

        INNER_ITERATION: do j = 1, ni
          associate(vj => v(:,:,:,j), zj => z(:,:,:,j))

            ! preconditioning, z(j) = K⁻¹v(j) ..................................

            select case(this % precon)
            case(1)
              call DiagonalPreconditioner( cl_operator, bc, cfg, dt  &
                                         , lmb_d, vj, zj             )
            case(2:3)
              call SpectralPreconditioner( cl_operator, bc, cfg, dt  &
                                         , lmb_d, vr_d, vl_d, vj, zj )
            case default
              call SetArray(zj, vj, multi = .true.)
            end select

            ! application of homogeneous operator, w = A z(j) ..................

            call GetHybridDiffusionTerm( this, cl_operator, 'T', theta &
                                       , u_0 = u_0, u = zj, r_d = w    )

            !!$omp do collapse(2)
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

        j = min(j, ni)

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

      !!$omp master
      deallocate(b, c, h, r, s, v, w, y, z)
      if (allocated(mask )) deallocate(mask)
      if (allocated(cfg  )) deallocate(cfg)
      if (allocated(lmb_d)) deallocate(lmb_d)
      if (allocated(vr_d )) deallocate(vr_d)
      if (allocated(vl_d )) deallocate(vl_d)
      !!$omp end master

    end associate

    !!$omp end master

  end subroutine DiffusionSolver

  !-----------------------------------------------------------------------------
  !>

  subroutine DiagonalPreconditioner(cl_operator, bc, cfg, dt, nu, r, z)
    class(CL_Operator_1D), intent(in)  :: cl_operator
    character,             intent(in)  :: bc(2)
    integer,               intent(in)  :: cfg(:)
    real(RNP),             intent(in)  :: dt
    real(RNP), contiguous, intent(in)  :: nu(:,:)
    real(RNP), contiguous, intent(in)  :: r(:,:,:)
    real(RNP), contiguous, intent(out) :: z(:,:,:)

    real(RNP), dimension(:,:), allocatable, save :: rs, zs

    real(RNP) :: lambda
    logical   :: periodic
    integer   :: nc, ne, np, ns
    integer   :: c

    associate( schwarz => cl_operator % elliptic_op % schwarz &
             , dx      => cl_operator % dx                    )

      ne = size(nu,1)
      nc = size(nu,2)
      np = schwarz % po + 1
      ns = schwarz % no * 2 + np

      allocate(rs(ns,ne))
      allocate(zs(ns,ne))

      lambda   = 1 / dt
      periodic = all(bc == 'P')

      do c = 1, nc
        call SetArray(z(:,:,c), ZERO)
        call schwarz % RestrictResidual(periodic, r(:,:,c), rs)
        call schwarz % Apply(cfg, dx, lambda, nu(:,c), rs, zs)
        call schwarz % MergeCorrections(periodic, zs, z(:,:,c))
      end do

      deallocate(rs, zs)

    end associate

  end subroutine DiagonalPreconditioner

  !-----------------------------------------------------------------------------
  !>

  subroutine SpectralPreconditioner(cl_operator, bc, cfg, dt, nu, vr, vl, r, z)
    class(CL_Operator_1D), intent(in)  :: cl_operator
    character,             intent(in)  :: bc(2)
    integer,               intent(in)  :: cfg(:)
    real(RNP),             intent(in)  :: dt
    real(RNP), contiguous, intent(in)  :: nu(:,:)
    real(RNP), contiguous, intent(in)  :: vr(:,:,:)
    real(RNP), contiguous, intent(in)  :: vl(:,:,:)
    real(RNP), contiguous, intent(in)  :: r(:,:,:)
    real(RNP), contiguous, intent(out) :: z(:,:,:)

    real(RNP), dimension(:,:,:), allocatable, save :: rs, zs

    real(RNP) :: lambda
    logical   :: periodic
    integer   :: nc, ne, np, ns
    integer   :: c, e

    associate( schwarz => cl_operator % elliptic_op % schwarz &
             , dx      => cl_operator % dx                    )

      ne = size(nu,1)
      nc = size(nu,2)
      np = schwarz % po + 1
      ns = schwarz % no * 2 + np

      allocate(rs(ns,ne,nc))
      allocate(zs(ns,ne,nc))

      lambda   = 1 / dt
      periodic = all(bc == 'P')

      ! restrict residual to subdomains
      do c = 1, nc
        call schwarz % RestrictResidual(periodic, r(:,:,c), rs(:,:,c))
      end do

      ! transform residual to eigenspace
      do e = 1, ne
        if (cl_operator % activity(e) > 0) then
          rs(:,e,:) = matmul(rs(:,e,:), transpose(vl(e,:,:)))
        end if
      end do

      ! apply Schwarz operator
      do c = 1, nc
        call schwarz % Apply(cfg, dx, lambda, nu(:,c), rs(:,:,c), zs(:,:,c))
      end do

      ! transform subdomain results back from eigenspace
      do e = 1, ne
        if (cl_operator % activity(e) > 0) then
          zs(:,e,:) = matmul(zs(:,e,:), transpose(vr(e,:,:)))
        end if
      end do

      ! merge results
      do c = 1, nc
        call SetArray(z(:,:,c), ZERO)
        call schwarz % MergeCorrections(periodic, zs(:,:,c), z(:,:,c))
      end do

      deallocate(rs, zs)

    end associate

  end subroutine SpectralPreconditioner

  !=============================================================================

end submodule SM_Diffusion
