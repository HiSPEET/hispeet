!> summary:  Base class for Burgers problems
!> author:   Joerg Stiller
!> date:     2023/05/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   Consistent approach to filtering of streamline-diffusivity
!===============================================================================

module CL__Problem__Burgers__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF
  use Execution_Control
  use CL__Problem__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_Burgers_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D Burgers problems
  !>
  !>       ∂u/∂t + ∂(u²/2)/∂x = ν∂²u/∂x² + f_s(x,t)
  !>
  !> with constant diffusivity ν. Available boundary conditions are
  !>
  !>   - Dirichlet with `bc = 'D'` and `bv = u`
  !>   - Neumann   with `bc = 'N'` and `bv = q = ν∂u/∂x`
  !>   - Periodic  with `bc = 'P'`

  type, abstract, extends(CL_Problem_1D) :: CL_Problem_Burgers_1D

    real(RNP) :: nu = 0  !< viscosity

  contains

    procedure :: GetConvectionTerm
    procedure :: GetDiffusionTerm
    procedure :: GetSDTerm
    procedure :: DiffusionSolver

  end type CL_Problem_Burgers_1D

  !=============================================================================

contains

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

    associate( eop  => cl_operator % eop      &
             , qop  => cl_operator % qop      &
             , iop  => cl_operator % iop_uq   &
             , mask => cl_operator % mask     &
             , po   => cl_operator % eop % po &
             , qo   => cl_operator % qop % po &
             , ne   => cl_operator % ne       )

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
        if (mask(e)) then

          ! compute fluxes in quadrature points
          if (use_interpolation) then
            u_q = matmul(iop%A, u(:,e,1))
          else
            u_q = u(:,e,1)
          end if
          f_q = ConvectiveFlux(u_q)

          ! apply mass-weighted transposed diff-matrix
          r_c(:,e,1) = matmul(MD_t, f_q)

        else
          r_c(:,e,1) = 0
        end if
      end do

      ! element boundary fluxes ................................................

      call GetNumericalConvectiveFlux(this, bv, u(:,:,1), h_c)

      do e = 1, ne
        if (mask(e)) then
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
  !> Convective flux fc(u)

  elemental function ConvectiveFlux(u) result(f_c)
    real(RNP), intent(in) :: u
    real(RNP) :: f_c

    f_c = HALF * u ** 2

  end function ConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Get numerical convective flux

  subroutine GetNumericalConvectiveFlux(this, bv, u, h_c)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    real(RNP),                    intent(in)  :: bv(:,:)
    real(RNP), contiguous,        intent(in)  :: u(0:,:)
    real(RNP), contiguous,        intent(out) :: h_c(0:)

    real(RNP), allocatable :: ul(:), ur(:)

    integer :: po, ne

    po = ubound(u,1)
    ne = ubound(u,2)

    allocate(ul(0:ne), ur(0:ne))

    ! inner traces
    ul(1:ne)   = u(po, 1:ne)
    ur(0:ne-1) = u( 0, 1:ne)

    ! left boundary
    select case(this % bc(1))
    case('D')
      ul(0) = 2 * bv(1,1) - u(0,1)
    case('P')
      ul(0) = u(po,ne)
    case default
      ul(0) = u(0,1)
    end select

    ! right boundary
    select case(this % bc(2))
    case('D')
      ur(ne) = 2 * bv(1,2)- u(po,ne)
    case('P')
      ur(ne) = u(0,1)
    case default
      ur(ne) = u(po,ne)
    end select

    h_c = RiemannFlux(ul, ur)
!   h_c = LLF_Flux(ul, ur)

  end subroutine GetNumericalConvectiveFlux

  !-----------------------------------------------------------------------------
  !> Numerical convective flux hc(ul,ur) based on Riemann solver

  elemental function RiemannFlux(ul, ur) result(h_c)
    real(RNP), intent(in) :: ul
    real(RNP), intent(in) :: ur
    real(RNP) :: h_c

    if (ul < ZERO .and. ur > ZERO) then
      h_c = ZERO
    else if (ul + ur > ZERO) then
      h_c = HALF * ul ** 2
    else
      h_c = HALF * ur ** 2
    end if

  end function RiemannFlux

  !-----------------------------------------------------------------------------
  !> Numerical convective flux hc(ul,ur) based on Riemann solver

  elemental function LLF_Flux(ul, ur) result(h_c)
    real(RNP), intent(in) :: ul
    real(RNP), intent(in) :: ur
    real(RNP) :: h_c

    h_c = (ul**2 + ur**2)/4 + max(abs(ul),abs(ur)) * (ul - ur)

  end function LLF_Flux

  !-----------------------------------------------------------------------------
  !> Diffusive contribution to RHS of DG-SEM formulation

  subroutine GetDiffusionTerm(this, cl_operator, bv, u, r_d)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u  (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_d(0:,:,:)

    real(RNP), allocatable, save :: f(:,:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      !$omp master
      allocate(f(0:po,ne), source = ZERO)

      elliptic_bc = this % bc(:)(1:1)
      elliptic_bv = bv(1,:)

      call cl_operator % elliptic_op % Residual( elliptic_bc                 &
                                               , elliptic_bv                 &
                                               , dx     = cl_operator % dx   &
                                               , lambda = ZERO               &
                                               , nu     = this % nu          &
                                               , f      = f                  &
                                               , u      = u  (:,:,1)         &
                                               , r      = r_d(:,:,1)         &
                                               , mask   = cl_operator % mask )

      deallocate(f)
      !$omp end master

    end associate

  end subroutine GetDiffusionTerm

  !-----------------------------------------------------------------------------
  !> Streamline-diffusion contribution to RHS of DG-SEM formulation

  subroutine GetSDTerm(this, cl_operator, tau, bv, u, r_sd)
    class(CL_Problem_Burgers_1D), intent(in)  :: this
    class(CL_Operator_1D),        intent(in)  :: cl_operator
    real(RNP),                    intent(in)  :: tau
    real(RNP),                    intent(in)  :: bv (:,:)
    real(RNP), contiguous,        intent(in)  :: u   (0:,:,:)
    real(RNP), contiguous,        intent(out) :: r_sd(0:,:,:)

    real(RNP), allocatable, save :: f(:,:), nu_sd(:,:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)

    associate( po => cl_operator % eop % po &
             , ne => cl_operator % ne       )

      !$omp master
      allocate(f    (0:po,ne), source = ZERO)
      allocate(nu_sd(0:po,ne))
      call GetStreamlineDiffusivity(this, tau, u(:,:,1), nu_sd)

      elliptic_bc = this % bc(:)(1:1)
      where(elliptic_bc == 'D')
        elliptic_bv = bv(1,:)
      else where
        elliptic_bv = 0
      end where

      call cl_operator % elliptic_op % Residual( elliptic_bc                 &
                                               , elliptic_bv                 &
                                               , dx     = cl_operator % dx   &
                                               , lambda = ZERO               &
                                               , nu     = nu_sd              &
                                               , f      = f                  &
                                               , u      = u   (:,:,1)        &
                                               , r      = r_sd(:,:,1)        &
                                               , mask   = cl_operator % mask )

      deallocate(f, nu_sd)
      !$omp end master

    end associate

  end subroutine GetSDTerm

  !---------------------------------------------------------------------------
  !> Implicit diffusion solver for Burgers problems
  !>
  !> Implicit method for solving or relaxing the diffusion subproblem
  !>
  !>       u = u₀ + ∆t [r_d(bv,u) + r_ds(bv,τ,u)]
  !>
  !> The streamline-diffusion term `f_ds` is included only if τ > 0.
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
  !>      - good smoother when used with overlap `schwarz%no ≈ po/4`
  !>
  !> 4. Inexact preconditioned conjugate gradient method (IPCG)
  !>      - requires proper initialization of `cl_operator%eop%schwarz`
  !>      - best iterative solver when used with no overlap, `schwarz%no = 0`
  !>      - good smoother when used with overlap `schwarz%no ≈ po/4`

  subroutine DiffusionSolver( this, cl_operator, dt, tau, bv, u_0, u &
                            , method, i_max, r_red, r_max            )
    class(CL_Problem_Burgers_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)    :: cl_operator
    real(RNP),             intent(in)    :: dt          !< ∆t = t - t₀
    real(RNP),             intent(in)    :: tau         !< SD time scale τ
    real(RNP),             intent(in)    :: bv (:,:)    !< boundary values
    real(RNP), contiguous, intent(in)    :: u_0(0:,:,:) !< initial value u₀
    real(RNP), contiguous, intent(inout) :: u  (0:,:,:) !< approx solution u
    integer,               intent(in)    :: method      !< solution method
    integer,               intent(in)    :: i_max       !< max num iterations
    real(RNP), optional,   intent(in)    :: r_red       !< residual reduction
    real(RNP), optional,   intent(in)    :: r_max       !< max residual

    real(RNP), allocatable, save :: f(:,:), nu_tot(:,:)
    character :: elliptic_bc(2)
    real(RNP) :: elliptic_bv(2)
    real(RNP) :: lambda
    integer   :: e

    associate( nu          => this % nu                 &
             , po          => cl_operator % eop % po    &
             , ne          => cl_operator % ne          &
             , dx          => cl_operator % dx          &
             , mask        => cl_operator % mask        &
             , Me          => cl_operator % Me          &
             , elliptic_op => cl_operator % elliptic_op )

      !$omp master

      ! elliptic boundary conditions and boundary values
      elliptic_bc = this % bc(:)(1:1)
      elliptic_bv = bv(1,:)

      ! Helmholtz parameter
      lambda = 1 / dt

      ! sources
      allocate(f(0:po,ne))
      do e = 1, ne
        if (mask(e)) then
          f(:,e) = lambda * Me * u_0(:,e,1)
        else
          f(:,e) = 0
        end if
      end do

      if (tau > ZERO) then

        allocate(nu_tot(0:po,ne))
        call GetStreamlineDiffusivity(this, tau, u(:,:,1), nu_tot)
        nu_tot = nu_tot + this % nu

        select case(method)
        case(1)
          call Error( 'DiffusionSolver'             &
                    , 'method 1 not suited for ISD' &
                    , 'CL__Problem__Burgers__1D'    )
        case(2)
          call elliptic_op % CG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(3)
          call elliptic_op % Schwarz_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(4)
          call elliptic_op % SchwarzPCG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu_tot, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        end select

      else

        select case(method)
        case(1)
          call elliptic_op % HybridSolver &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, f, u(:,:,1) )
        case(2)
          call elliptic_op % CG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(3)
          call elliptic_op % Schwarz_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        case(4)
          call elliptic_op % SchwarzPCG_Method &
                  ( elliptic_bc, elliptic_bv, dx, lambda, nu, f, u(:,:,1) &
                  , i_max, r_red, r_max, mask )
        end select

      end if

      deallocate(f)
      if (allocated(nu_tot)) deallocate(nu_tot)
      !$omp end master

    end associate

  end subroutine DiffusionSolver

  !---------------------------------------------------------------------------
  !> Computation of streamline diffusivity

  subroutine GetStreamlineDiffusivity(cl_problem, tau, u, nu_sd)
    class(CL_Problem_Burgers_1D), intent(in) :: cl_problem
    real(RNP),             intent(in)  :: tau         !< SD time scale τ
    real(RNP), contiguous, intent(in)  :: u    (0:,:) !< approx solution u
    real(RNP), contiguous, intent(out) :: nu_sd(0:,:) !< approx solution u

    integer :: e

    nu_sd = tau/2 * u**2

    select case(cl_problem % nu_sd_filter)
    case(1)
      ! element maximum
      do e = 1, size(nu_sd,2)
        nu_sd(:,e) = maxval(nu_sd(:,e))
      end do
    end select

  end subroutine GetStreamlineDiffusivity

  !=============================================================================

end module CL__Problem__Burgers__1D

