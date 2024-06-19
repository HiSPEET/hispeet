!> summary:  Base class for Burgers problems
!> author:   Joerg Stiller
!> date:     2023/05/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__Burgers__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE, PI
  use Execution_Control
  use Array_Assignments

  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_Burgers_1D
  public :: CL_Problem_Burgers_Options_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D Burgers problems
  !>
  !>       ∂u/∂t + ∂(u²/2)/∂x = ν∂²u/∂x² + f_s(x,t)
  !>
  !> with constant diffusivity ν. Available boundary conditions are
  !>
  !>   - Dirichlet with `bc = 'D'` and `bv = u`
  !>   - Neumann   with `bc = 'N'` and `bv = ∂u/∂x`
  !>   - Periodic  with `bc = 'P'`
  !>   - Static    with `bc = 'S'`

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_Burgers_1D

    real(RNP) :: nu  !< viscosity

  contains

    procedure :: HasDiffusion
    procedure :: ConvectiveFlux
    procedure :: ConvectiveJacobian
    procedure :: ConvectiveEigensystem
    procedure :: GetConvectionTerm
    procedure :: GetHybridDiffusionTerm
    procedure :: GetDiffusionTerm
    procedure :: GetSDTerm
    procedure :: DiffusionSolver
    procedure :: GetMaxVelocity
    procedure :: GetMaxDiffusivity

    procedure :: Init_CL_Problem_Burgers_1D

  end type CL_Problem_Burgers_1D

  !-----------------------------------------------------------------------------
  !> 1D Burgers options

  type, extends(CL_Problem_Options_1D) :: CL_Problem_Burgers_Options_1D
    real(RNP) :: nu = 0.1 !< viscosity
  end type CL_Problem_Burgers_Options_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the base type

  subroutine Init_CL_Problem_Burgers_1D(this, opt)
    class(CL_Problem_Burgers_1D),         intent(inout) :: this
    class(CL_Problem_Burgers_Options_1D), intent(in)    :: opt

    call this % Init_CL_Problem_1D(opt)

    this % nc = 1
    this % nu = opt % nu

  end subroutine Init_CL_Problem_Burgers_1D

  !-----------------------------------------------------------------------------
  !> Query whether problem has non vanishing physical diffusion

  logical function HasDiffusion(this)
    class(CL_Problem_Burgers_1D), intent(in) :: this

    HasDiffusion = this % nu > 0 .or. this % dc_diffusivity > 0

  end function HasDiffusion

  !-----------------------------------------------------------------------------
  !> Calculates the convective flux with given conservative variables

  pure function ConvectiveFlux(this, u) result(f_c)
    class(CL_Problem_Burgers_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:)
    real(RNP) :: f_c(this%nc)

    f_c(1) = HALF * u(1)**2

  end function ConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Calculates the Jacobian ot the Convective Term

  pure function ConvectiveJacobian(this, u) result(A_c)
    class(CL_Problem_Burgers_1D), intent(in) :: this
    real(RNP), intent(in) :: u(:)
    real(RNP) :: A_c(this%nc,this%nc)

    A_c(1,1) = u(1)

  end function ConvectiveJacobian

  !-----------------------------------------------------------------------------
  !> Calculates the Eigendecomposition of the Convective Jacobian

  pure subroutine ConvectiveEigensystem(this, u, lambda, R, L)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    real(RNP),                    intent(in)  :: u(:)
    real(RNP), optional,          intent(out) :: lambda(:)
    real(RNP), optional,          intent(out) :: R(:,:)
    real(RNP), optional,          intent(out) :: L(:,:)

    if(present(lambda)) then
      lambda(1) = u(1)
    end if

    if(present(R)) then
      R(1,1) = ONE
    end if

    if(present(L)) then
      L(1,1) = ONE
    end if

  end subroutine ConvectiveEigensystem

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

  subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_c(0:,:,:)

    real(RNP), allocatable :: MD_t(:,:), f_q(:), u_q(:), h_c(:)
    logical :: use_interpolation
    integer :: e, i, j, j0

    associate( eop      => cl_operator % eop      &
             , qop      => cl_operator % qop      &
             , iop      => cl_operator % iop_uq   &
             , po       => cl_operator % eop % po &
             , qo       => cl_operator % qop % po &
             , ne       => cl_operator % ne       &
             , activity => cl_operator % activity )

      ! initialization .........................................................

      !$omp master
      allocate(h_c(0:ne))

      allocate(MD_t(0:po, 0:qo), f_q(0:qo), u_q(0:qo))

      use_interpolation = qo /= po

      if (use_interpolation) then
        j0 = lbound(iop%A,1)
        do i = 0, po
        do j = 0, qo
          MD_t(i,j) = qop%w(j) * dot_product(iop%A(j0+j,:), eop%D(:,i))
        end do
        end do
      else
        do i = 0, po
        do j = 0, po
          MD_t(i,j) = eop%w(j) * eop%D(j,i)
        end do
        end do
      end if

      ! element integrals ......................................................

      !$omp do
      do e = 1, ne
        if (activity(e) > 0) then

          ! interpolate solution to quadrature points
          if (use_interpolation) then
            u_q = matmul(iop%A, u(:,e,1))
          else
            u_q = u(:,e,1)
          end if

          ! compute convective fluxes
          f_q = HALF * u_q**2

          ! apply mass-weighted transposed diff-matrix
          r_c(:,e,1) = matmul(MD_t, f_q)

        else
          r_c(:,e,1) = 0
        end if
      end do

      ! element boundary fluxes ................................................

      call GetNumericalConvectiveFlux(this, bv, u(:,:,1), h_c)

      do e = 1, ne
        if (activity(e) > 0) then
          r_c( 0,e,1) = r_c( 0,e,1) + h_c(e-1)
          r_c(po,e,1) = r_c(po,e,1) - h_c(e)
        end if
      end do

      ! finalization ...........................................................

      deallocate(h_c)
      !$omp end master

    end associate

  end subroutine GetConvectionTerm

  !-----------------------------------------------------------------------------
  !> Get numerical convective flux

  subroutine GetNumericalConvectiveFlux(this, bv, u, h_c)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    real(RNP),                    intent(in)  :: bv(:,:)
    real(RNP), contiguous,        intent(in)  :: u(0:,:)
    real(RNP), contiguous,        intent(out) :: h_c(0:)

    real(RNP), allocatable :: u_l(:), u_r(:)

    integer :: po, ne

    po = ubound(u,1)
    ne = ubound(u,2)

    allocate(u_l(0:ne), u_r(0:ne))

    ! inner traces
    u_l(1:ne)   = u(po, 1:ne)
    u_r(0:ne-1) = u( 0, 1:ne)

    ! left boundary
    select case(this % bc(1))
    case('D','S')
      u_l(0) = 2 * bv(1,1) - u(0,1)
    case('P')
      u_l(0) = u(po,ne)
    case default
      u_l(0) = u(0,1)
    end select

    ! right boundary
    select case(this % bc(2))
    case('D','S')
      u_r(ne) = 2 * bv(1,2) - u(po,ne)
    case('P')
      u_r(ne) = u(0,1)
    case default
      u_r(ne) = u(po,ne)
    end select

    h_c = RiemannFlux(u_l, u_r)
!   h_c = LocalLaxFriedrichsFlux(u_l, u_r)

  end subroutine GetNumericalConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Numerical convective flux based on Riemann solver

  elemental function RiemannFlux(u_l, u_r) result(h_c)
    real(RNP), intent(in) :: u_l
    real(RNP), intent(in) :: u_r
    real(RNP) :: h_c

    if (u_l < ZERO .and. u_r > ZERO) then
      h_c = ZERO
    else if (u_l + u_r > ZERO) then
      h_c = HALF * u_l ** 2
    else
      h_c = HALF * u_r ** 2
    end if

  end function RiemannFlux

  !-----------------------------------------------------------------------------
  !> Numerical convective flux based on local Lax-Friedrichs method

! elemental function LocalLaxFriedrichsFlux(u_l, u_r) result(h_c)
!   real(RNP), intent(in) :: u_l
!   real(RNP), intent(in) :: u_r
!   real(RNP) :: h_c
!
!   h_c = HALF * ( (u_l**2 + u_r**2) * HALF             &
!                + max(abs(u_l),abs(u_r)) * (u_l - u_r) )
!
! end function LocalLaxFriedrichsFlux

  !-----------------------------------------------------------------------------
  !> Unified diffusion term combining physical and streamline contributions
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

  subroutine GetHybridDiffusionTerm &
      (this, cl_operator, comp, theta, bv, u_0, u, r_d)

    class(CL_Problem_Burgers_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator !< spatial operators
    character(len=*),      intent(in)  :: comp        !< composition flags
    real(RNP),             intent(in)  :: theta       !< SD time scale
    real(RNP), optional,   intent(in)  :: bv(:,:)     !< boundary values
    real(RNP), contiguous, intent(in)  :: u_0(0:,:,:) !< u₀(x,t)
    real(RNP), contiguous, intent(in)  :: u  (0:,:,:) !< u(x,t)
    real(RNP), contiguous, intent(out) :: r_d(0:,:,:) !< diffusion RHS

    logical,   allocatable, save :: mask(:) ! active element mask
    real(RNP), allocatable, save :: nu(:,:) ! total diffusivity
    real(RNP), allocatable, save :: f(:,:)  ! source distribution (= 0)

    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)

    logical :: has_pd ! switch for physical diffusion
    logical :: has_dc ! switch for physical diffusion
    logical :: has_sd ! switch for streamline diffusion

    integer :: i

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      ! initialization .........................................................

      has_pd = scan(comp,'PT') > 0 .and. 0 < this % nu
      has_dc = scan(comp,'DT') > 0 .and. 0 < this % dc_diffusivity
      has_sd = scan(comp,'ST') > 0 .and. 0 < theta

      ! boundary conditions
      do i = 1, 2
        select case(this % bc(i))
        case('D','S')
          elliptic_bc(i) = 'D'
        case('N')
          elliptic_bc(i) = 'N'
        case('P')
          elliptic_bc(i) = 'P'
        end select
      end do

      ! boundary values
      if (present(bv)) then
        elliptic_bv = bv(1,1:2)
      else
        elliptic_bv = 0
      end if

      ! shared workspace
      allocate(mask(ne),    source = cl_operator % activity > 0)
      allocate(nu(0:po,ne), source = ZERO)
      allocate(f(0:po,ne),  source = ZERO)

      if (has_dc .or. has_sd) then

        ! variable total diffusivity ...........................................

        call GetHybridDiffusity(this, cl_operator, comp, theta, u_0(:,:,1), nu)

        call cl_operator % elliptic_op % Residual( elliptic_bc               &
                                                 , elliptic_bv               &
                                                 , dx     = cl_operator % dx &
                                                 , lambda = ZERO             &
                                                 , nu     = nu               &
                                                 , f      = f                &
                                                 , u      = u  (:,:,1)       &
                                                 , r      = r_d(:,:,1)       &
                                                 , mask   = mask             )

      else if (has_pd .and. this % nu > 0) then

        ! constant physical diffusivity ........................................

        call cl_operator % elliptic_op % Residual( elliptic_bc               &
                                                 , elliptic_bv               &
                                                 , dx     = cl_operator % dx &
                                                 , lambda = ZERO             &
                                                 , nu     = this % nu        &
                                                 , f      = f                &
                                                 , u      = u  (:,:,1)       &
                                                 , r      = r_d(:,:,1)       &
                                                 , mask   = mask             )

      else

        ! no diffusivity .......................................................

        call SetArray(r_d(:,:,1), ZERO)

      end if

      ! finalization ...........................................................

      ! free shared workspace
      deallocate(mask, nu, f)

    end associate

  end subroutine GetHybridDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Composition of hybrid diffusivity

  subroutine GetHybridDiffusity(this, cl_operator, comp, theta, u, nu)
    class(CL_Problem_Burgers_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator !< spatial operators
    character(len=*),      intent(in)  :: comp        !< composition flags
    real(RNP),             intent(in)  :: theta       !< SD time scale
    real(RNP), contiguous, intent(in)  :: u (0:,:)    !< u(x,t)
    real(RNP), contiguous, intent(out) :: nu(0:,:)    !< diffusivity

    logical :: has_pd ! switch for physical diffusion
    logical :: has_dc ! switch for discontinuity capturing diffusion
    logical :: has_sd ! switch for streamline diffusion
    logical :: has_ad ! switch for artificial diffusion

    ! discontinuity capturing
    real(RNP), allocatable :: V_inv_dc(:,:) ! inverse Vandermonde matrix
    real(RNP) :: c_dc, d_dc, s_dc           ! parameters

    ! artificial diffusivity filtering
    real(RNP), allocatable :: Q_ad(:,:)     ! filtering operator

    integer :: e, po_cut

    associate( dx => cl_operator % dx       &
             , po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      ! initialization .........................................................

      has_pd = scan(comp,'PT') > 0 .and. 0 < this % nu
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

      ! composition of hybrid diffusivity ......................................

      do e = 1, ne

        nu(:,e) = 0

        if (cl_operator % activity(e) < 0) cycle

        ! artificial diffusivity
        if (has_ad) then

          ! discontinuity capturing diffusivity
          if (has_dc) then
            select case(this % dc_diffusivity)
            case(1)
              nu(:,e) = PerssonDiffusivity(c_dc, d_dc, s_dc, u(:,e), V_inv_dc)
            end select
          end if

          ! streamline diffusivity
          if (has_sd) then
            nu(:,e) = nu(:,e) + theta/2 * u(:,e)**2
          end if

          ! filtering
          select case(this % ad_filter_method)
          case(1:2)
            nu(:,e) = matmul(Q_ad, nu(:,e))
          case(3)
            nu(:,e) = maxval(nu(:,e))
          end select

        end if

        ! physical diffusivity
        if (has_pd) then
          nu(:,e) = nu(:,e) + this % nu
        end if

      end do

    end associate

  end subroutine GetHybridDiffusity

  !-----------------------------------------------------------------------------
  !> Discontinuity capturing diffusivity of Persson & Preraire (2006)

  pure function PerssonDiffusivity(c_max, delta_s, s_ref, u, V_inv) result(nu)
    real(RNP), intent(in) :: c_max
    real(RNP), intent(in) :: delta_s
    real(RNP), intent(in) :: s_ref
    real(RNP), intent(in) :: u(0:)
    real(RNP), intent(in) :: V_inv(0:,0:)
    real(RNP) :: nu

    real(RNP), parameter :: eps = epsilon(ONE)
    real(RNP) :: u_hat(0:ubound(u,1))
    real(RNP) :: norm_u, nu_max, s
    integer   :: i, po

    po = ubound(u,1)

    u_hat  = matmul(V_inv, u)
    nu_max = c_max * maxval(abs(u))
    norm_u = 0
    do i = 0, po
      norm_u = norm_u + u_hat(i)**2 / (2*i + 1)
    end do

    s = (u_hat(po)**2 / (2*po + 1)) / max(norm_u, eps)
    s = log10(max(s, eps))

    if (s > s_ref + delta_s) then
      nu = nu_max
    else if (s > s_ref - delta_s) then
      nu = c_max/2 * maxval(abs(u)) * (1 + sin(PI * (s - s_ref) / (2*delta_s)))
    else
      nu = 0
    end if

  end function PerssonDiffusivity

  !-----------------------------------------------------------------------------
  !> Diffusive contribution to RHS of DG-SEM formulation

  subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_d(0:,:,:)

    call GetHybridDiffusionTerm(this, cl_operator, 'PD', ZERO, bv, u, u, r_d)

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Streamline-diffusion contribution to RHS of DG-SEM formulation

  subroutine GetSDTerm(this, cl_operator, theta, bv, u_0, u, r_sd)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: theta
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u_0 (0:,:,:)
    real(RNP), contiguous,        intent(in)  :: u   (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_sd(0:,:,:)

    call GetHybridDiffusionTerm(this, cl_operator, 'S', theta, bv, u_0, u, r_sd)

  end subroutine GetSDTerm

  !---------------------------------------------------------------------------
  !> Implicit diffusion solver for Burgers problems
  !>
  !> Implicit method for solving or relaxing the diffusion subproblem
  !>
  !>       u = f + ∆t [r_d(bv,u) + r_ds(bv,θ,u₀,u)]
  !>
  !> The streamline-diffusion term `r_ds` is evaluated with `u₀` and included
  !> only if θ > 0.
  !> At present, the following solution methods are available:
  !>
  !> 1. Direct hybrid solver
  !>      - requires `cl_operator%eop%hybrid = .true.`
  !>      - not suited for cases with frozen elements
  !>      - not suited for streamline-diffusion unless f_c'(u) is constant
  !>
  !> 2. Conjugate gradient method
  !>      - basic solver
  !>      - bad smoother
  !>
  !> 3. Schwarz method
  !>      - requires proper initialization of `cl_operator%eop%schwarz`
  !>      - bad solver
  !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`
  !>
  !> 4. Inexact preconditioned conjugate gradient method (IPCG)
  !>      - requires proper initialization of `cl_operator%eop%schwarz`
  !>      - best iterative solver when used with no overlap, i.e.,
  !>        `schwarz % delta = -1`, `schwarz % no_min = -1`
  !>      - good smoother when used with overlap `schwarz % delta ≈ 0.25`

  subroutine DiffusionSolver( this, cl_operator, dt, theta, bv, f, u_0, u &
                            , method, i_max, r_red, r_max                 )
    class(CL_Problem_Burgers_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP),             intent(in)    :: dt          !< ∆t = t - t₀
    real(RNP),             intent(in)    :: theta       !< SD time scale θ
    real(RNP),             intent(in)    :: bv (:,:)    !< boundary values
    real(RNP), contiguous, intent(in)    :: f  (0:,:,:) !< sources
    real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< frozen solution
    real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< approx solution
    integer,               intent(in)    :: method      !< solution method
    integer,               intent(in)    :: i_max       !< max num iterations
    real(RNP), optional,   intent(in)    :: r_red       !< residual reduction
    real(RNP), optional,   intent(in)    :: r_max       !< max residual

    real(RNP), allocatable, save :: g(:,:), nu(:,:)
    logical,   allocatable, save :: mask(:)

    logical   :: has_variable_nu
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)
    real(RNP) :: lambda
    integer   :: e, i

    associate( po          => cl_operator % eop % po    &
             , ne          => cl_operator % ne          &
             , dx          => cl_operator % dx          &
             , Me          => cl_operator % Me          &
             , activity    => cl_operator % activity    &
             , elliptic_op => cl_operator % elliptic_op )

      ! check for variable diffusivity
      has_variable_nu = 0 < this%dc_diffusivity .or. 0 < theta

      ! elliptic boundary conditions
      do i = 1, 2
        select case(this % bc(i))
        case('P')
          elliptic_bc(i) = 'P'
        case('D','S')
          elliptic_bc(i) = 'D'
        case default
          elliptic_bc(i) = 'N'
        end select
      end do

      ! boundary values
      elliptic_bv = bv(1,:)

      ! Helmholtz parameter
      lambda = 1 / dt

      ! activity mask
      allocate(mask(ne), source = activity > 0)

      ! sources
      allocate(g(0:po,ne))
      do e = 1, ne
        if (mask(e)) then
          g(:,e) = lambda * Me * f(:,e,1)
        else
          g(:,e) = 0
        end if
      end do

      if (has_variable_nu) then

        allocate(nu(0:po,ne))
        call GetHybridDiffusity(this, cl_operator, 'T', theta, u_0(:,:,1), nu)

        select case(method)
        case(1)
          call Error( 'DiffusionSolver'             &
                    , 'method 1 not suited for ISD' &
                    , 'CL__Problem__Burgers__1D'    )
        case(2)
          call elliptic_op % CG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(3)
          call elliptic_op % Schwarz_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(4)
          call elliptic_op % SchwarzPCG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        end select

      else

        select case(method)
        case(1)
          call elliptic_op % HybridSolver &
                  ( elliptic_bc, elliptic_bv, dx, lambda, this%nu, g, u(:,:,1) )
        case(2)
          call elliptic_op % CG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, this%nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(3)
          call elliptic_op % Schwarz_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, this%nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(4)
          call elliptic_op % SchwarzPCG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, this%nu, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        end select

      end if

      deallocate(g, mask)
      if (allocated(nu)) deallocate(nu)

    end associate

  end subroutine DiffusionSolver

  !-----------------------------------------------------------------------------
  !> Provides the maximum velocity based on eigenvalues of advective Jacobian

  subroutine GetMaxVelocity(this, cl_operator, u, v_max)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
    real(RNP),                    intent(out) :: v_max     !< maximum velocity

    integer :: e

    v_max = 0

    do e = 1, cl_operator % ne
      if (cl_operator % activity(e) > 0) then
        v_max = max(v_max, maxval(abs(u(:,e,1))))
      end if
    end do

    ! avoid compiler warning
    if (this % nc > 0) return

  end subroutine GetMaxVelocity

  !-----------------------------------------------------------------------------
  !> Provides the maximum diffusivity

  subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
    real(RNP),                    intent(out) :: nu_max     !< maximum velocity

    nu_max = this % nu

    ! avoid compiler warnings
    if (cl_operator % ne > 1 .or. size(u) > 0) return

  end subroutine GetMaxDiffusivity

  !=============================================================================

end module CL__Problem__Burgers__1D
