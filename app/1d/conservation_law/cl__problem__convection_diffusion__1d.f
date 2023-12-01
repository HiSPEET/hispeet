!> summary:  Base class for convection-diffusion problems
!> author:   Robin Fraenzel, Joerg Stiller
!> date:     2023/09/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__Convection_Diffusion__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF
  use Execution_Control
  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_ConvectionDiffusion_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D convection-diffusion problems
  !>
  !>     ∂u/∂t + v ∂u/∂x = nu ∂²u/∂x² + f_s(x,t)
  !>
  !> with constant velocity v and diffusivity ν. Available boundary conditions
  !> are
  !>
  !>   - Dirichlet with `bc = 'D'` and `bv = u`
  !>   - Neumann   with `bc = 'N'` and `bv = ∂u/∂x`
  !>   - Periodic  with `bc = 'P'`

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_ConvectionDiffusion_1D

    real(RNP) :: v  = 1  !< velocity
    real(RNP) :: nu = 0  !< viscosity

  contains

    procedure :: GetConvectionTerm
    procedure :: GetDiffusionTerm
    procedure :: GetSDTerm
    procedure :: DiffusionSolver
    procedure :: GetMaxVelocity
    procedure :: GetMaxDiffusivity

  end type CL_Problem_ConvectionDiffusion_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Convective contribution to RHS of DG-SEM formulation

  subroutine GetConvectionTerm(this, cl_operator, bv, u, r_c)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
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

          ! compute fluxes in quadrature points
          if (use_interpolation) then
            u_q = matmul(iop%A, u(:,e,1))
          else
            u_q = u(:,e,1)
          end if
          f_q = ConvectiveFlux(this%v, u_q)

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
  !> Convective flux fc(v,u)

  elemental function ConvectiveFlux(v, u) result(f_c)
    real(RNP), intent(in) :: v
    real(RNP), intent(in) :: u
    real(RNP) :: f_c

    f_c = v * u

  end function ConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Get numerical convective flux

  subroutine GetNumericalConvectiveFlux(this, bv, u, h_c)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
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
    case('D')
      u_l(0) = 2 * bv(1,1) - u(0,1)
    case('P')
      u_l(0) = u(po,ne)
    case default
      u_l(0) = u(0,1)
    end select

    ! right boundary
    select case(this % bc(2))
    case('D')
      u_r(ne) = 2 * bv(1,2) - u(po,ne)
    case('P')
      u_r(ne) = u(0,1)
    case default
      u_r(ne) = u(po,ne)
    end select

    h_c = RiemannFlux(this % v, u_l, u_r)

  end subroutine GetNumericalConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Numerical convective flux based on Riemann solver

  elemental function RiemannFlux(v, u_l, u_r) result(h_c)
    real(RNP), intent(in) :: v
    real(RNP), intent(in) :: u_l
    real(RNP), intent(in) :: u_r
    real(RNP) :: h_c

    if (v >= ZERO) then
      h_c = v * u_l
    else
      h_c = v * u_r
    end if

  end function RiemannFlux

  !-----------------------------------------------------------------------------
  !> Diffusive contribution to RHS of DG-SEM formulation

  subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_d(0:,:,:)

    real(RNP), allocatable, save :: f(:,:)
    logical,   allocatable, save :: mask(:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      !$omp master
      allocate(mask(ne), source = cl_operator % activity > 0)
      allocate(f(0:po,ne), source = ZERO)

      elliptic_bc = this % bc(:)(1:1)
      elliptic_bv = bv(1,:)

      call cl_operator % elliptic_op % Residual( elliptic_bc               &
                                               , elliptic_bv               &
                                               , dx     = cl_operator % dx &
                                               , lambda = ZERO             &
                                               , nu     = this % nu        &
                                               , f      = f                &
                                               , u      = u  (:,:,1)       &
                                               , r      = r_d(:,:,1)       &
                                               , mask   = mask             )

      deallocate(f, mask)
      !$omp end master

    end associate

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Streamline-diffusion contribution to RHS of DG-SEM formulation

  subroutine GetSDTerm(this, cl_operator, theta, bv, u_0, u, r_sd)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: theta
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u_0 (0:,:,:)
    real(RNP), contiguous,        intent(in)  :: u   (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_sd(0:,:,:)

    real(RNP), allocatable, save :: f(:,:), nu_sd(:,:)
    logical,   allocatable, save :: mask(:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      !$omp master
      allocate(mask(ne), source = cl_operator % activity > 0)
      allocate(f(0:po,ne), source = ZERO)
      allocate(nu_sd(0:po,ne), source = theta/2 * this%v**2)

      elliptic_bc = this % bc(:)(1:1)
      elliptic_bv = bv(1,:)

      call cl_operator % elliptic_op % Residual( elliptic_bc               &
                                               , elliptic_bv               &
                                               , dx     = cl_operator % dx &
                                               , lambda = ZERO             &
                                               , nu     = nu_sd            &
                                               , f      = f                &
                                               , u      = u   (:,:,1)      &
                                               , r      = r_sd(:,:,1)      &
                                               , mask   = mask             )

      deallocate(f, mask, nu_sd)
      !$omp end master

    end associate

    ! avoid compiler warning
    if (size(u_0) > 0) return

  end subroutine GetSDTerm

  !---------------------------------------------------------------------------
  !> Implicit diffusion solver for convection-diffusion problems
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
                            , method, i_max, r_red, r_max               )
    class(CL_Problem_ConvectionDiffusion_1D), intent(in) :: this
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

    real(RNP), allocatable, save :: g(:,:), nu_tot(:,:)
    logical,   allocatable, save :: mask(:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)
    real(RNP) :: lambda
    integer   :: e

    associate( nu          => this % nu                 &
             , po          => cl_operator % eop % po    &
             , ne          => cl_operator % ne          &
             , dx          => cl_operator % dx          &
             , Me          => cl_operator % Me          &
             , activity    => cl_operator % activity    &
             , elliptic_op => cl_operator % elliptic_op )

      !$omp master

      ! elliptic boundary conditions and boundary values
      elliptic_bc = this % bc(:)(1:1)
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

      if (theta > ZERO) then

        allocate(nu_tot(0:po,ne), source = this%nu + theta/2 * this%v**2)

        select case(method)
        case(1)
          call Error( 'DiffusionSolver'             &
                    , 'method 1 not suited for ISD' &
                    , 'CL__Problem__Convection_Diffusion__1D'    )
        case(2)
          call elliptic_op % CG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(3)
          call elliptic_op % Schwarz_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(4)
          call elliptic_op % SchwarzPCG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, g, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        end select

      else

        select case(method)
        case(1)
          call elliptic_op % HybridSolver &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, g, u(:,:,1) )
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

      end if

      deallocate(g, mask)
      if (allocated(nu_tot)) deallocate(nu_tot)
      !$omp end master

    end associate

    ! avoid compiler warning
    if (size(u_0) > 0) return

  end subroutine DiffusionSolver

  !-----------------------------------------------------------------------------
  !> Provides the maximum velocity based on eigenvalues of advective Jacobian

  subroutine GetMaxVelocity(this, cl_operator, u, v_max)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
    real(RNP),                    intent(out) :: v_max     !< maximum velocity

    v_max = abs(this % v)

    ! avoid compiler warning
    if (this % nc > 0 .or. cl_operator % ne > 0 .or. size(u) > 0) return

  end subroutine GetMaxVelocity

  !-----------------------------------------------------------------------------
  !> Provides the maximum diffusivity

  subroutine GetMaxDiffusivity(this, cl_operator, u, nu_max)
    class(CL_Problem_ConvectionDiffusion_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP), contiguous,        intent(in)  :: u(0:,:,:) !< solution variable
    real(RNP),                    intent(out) :: nu_max     !< maximum velocity

    nu_max = this % nu

    ! avoid compiler warnings
    if (this % nc > 0 .or. cl_operator % ne > 1 .or. size(u) > 0) return

  end subroutine GetMaxDiffusivity

  !=============================================================================

end module CL__Problem__Convection_Diffusion__1D
