!> summary:  VARK-SDC method with CG-SEM for 1D convection-diffusion
!> author:   Joerg Stiller
!> date:     2019/03/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Variable additive RK SDC method with CG-SEM for 1D convection-diffusion
!>
!> **TBD:** Rework description, extend to SDC
!>
!> This module provides the type `ConvDiff_VARK` which extends the variable
!> additive Runge-Kutta methods defined in `VARK_Method` for advancing the
!> solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> in time. Spatial discretization is based on the continuous spectral-element
!> method using nodal base functions along with GLL quadrature.
!>
!===============================================================================

module CG_ConvDiff_1D__VARK_SDC
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use CG_ConvDiff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Solver_1D
  use VA_Runge_Kutta_Method
  use Harmonic_Wave_Package
  use Spectral_Deferred_Correction
  implicit none
  private

  public :: VARK_SDC_Method1D
  public :: VARK_SDC_Options1D

  !-----------------------------------------------------------------------------
  !> Propagator for VARK-based SDC for 1D convection-diffusion with SEM

  type, extends(VARK_Method) :: VARK_Propagator1D

    real(RNP), allocatable :: u      (:,:,:) !< solution
    real(RNP), allocatable :: F_im   (:,:,:) !< new implicit RHS
    real(RNP), allocatable :: F_im_o (:,:,:) !< old implicit RHS
    real(RNP), allocatable :: F_ex   (:,:,:) !< new explicit RHS
    real(RNP), allocatable :: F_ex_o (:,:,:) !< old explicit RHS

    real(RNP), allocatable :: w      (:,:)   !< stage quadrature weights
    real(RNP), allocatable :: S      (:,:,:) !< stage integrals ∫Fdt

  contains

    procedure :: Init_VARK_Propagator1D

  end type VARK_Propagator1D

  ! constructor
  interface VARK_Propagator1D
    module procedure New_Propagator1D
  end interface

  !-----------------------------------------------------------------------------
  !> Implementation of the VARK-SDC method for 1D convection-diffusion
  !>
  !> One VARK per time step, for start

  type, extends(SDC_Method) :: VARK_SDC_Method1D

    integer :: po  = -1 !< polynomial order
    integer :: ne  = -1 !< number of elements

    type(VARK_Propagator1D) :: vark !< VARK propagator

  contains

    procedure :: Init_VARK_SDC_Method1D => Init_VARK_SDC
    procedure :: TimeStep

  end type VARK_SDC_Method1D

  ! constructors
  interface VARK_SDC_Method1D
    module procedure New_VARK_SDC
  end interface

  !-----------------------------------------------------------------------------
  !> Type bundling VARK-SDC options

  type, extends(SDC_Options) :: VARK_SDC_Options1D
    integer   :: so_p =  1 !< order of principal stages
    real(RNP) :: r_ex = -1 !< last coefficient of R_ex(z)
  end type VARK_SDC_Options1D

contains

!===============================================================================
! VARK_Propagator1D

!-------------------------------------------------------------------------------
!> New 1D VARK propagator

function New_Propagator1D(po, ne, t, set, range, so_p, r_ex) result(this)
  integer,           intent(in) :: po       !< polynomial order
  integer,           intent(in) :: ne       !< number of elements
  real(RNP),         intent(in) :: t(0:)    !< SDC points
  integer,           intent(in) :: set      !< point set: equidistant (1) GLL (2)
  integer, optional, intent(in) :: range(2) !< point range, default: all
  integer,           intent(in) :: so_p     !< VARK stage order at SDC  points
  real(RNP),         intent(in) :: r_ex     !< last coefficient of R_ex(z)
  type(VARK_Propagator1D)       :: this

  call Init_VARK_Propagator1D(this, po, ne, t, set, range, so_p, r_ex)

end function New_Propagator1D

!-------------------------------------------------------------------------------
!> Initialize 1D VARK propagator

subroutine Init_VARK_Propagator1D(this, po, ne, t, set, range, so_p, r_ex)
  class(VARK_Propagator1D), intent(inout) :: this
  integer,           intent(in) :: po       !< polynomial order
  integer,           intent(in) :: ne       !< number of elements
  real(RNP),         intent(in) :: t(0:)    !< SDC points
  integer,           intent(in) :: set      !< point set: equistant (1) GLL (2)
  integer, optional, intent(in) :: range(2) !< point range, default: all
  integer,           intent(in) :: so_p     !< VARK stage order at SDC  points
  real(RNP),         intent(in) :: r_ex     !< last coefficient of R_ex(z)

  real(RNP), allocatable :: x(:), w(:)
  real(RNP) :: tk, yk, w_ji
  integer   :: i, j, k, n_sub

  ! safeguard ..................................................................

  if (allocated(this % u     )) deallocate(this % u     )
  if (allocated(this % F_im  )) deallocate(this % F_im  )
  if (allocated(this % F_im_o)) deallocate(this % F_im_o)
  if (allocated(this % F_ex  )) deallocate(this % F_ex  )
  if (allocated(this % F_ex_o)) deallocate(this % F_ex_o)
  if (allocated(this % w     )) deallocate(this % w     )
  if (allocated(this % S     )) deallocate(this % S     )

  ! preliminaries ..............................................................

  ! number of SDC subintervals
  n_sub = ubound(t,1)

  ! selected point range
  if (present(range)) then
    i = range(1)
    j = range(2)
  else
    i = 0
    j = n_sub
  end if

  ! VARK method ................................................................

  call this % Init_VARK_Method(t(i:j), so_p, r_ex)
  call this % Write()

  ! workspace ..................................................................

  allocate(this % u(0:po, ne, this%n_stage))

  allocate(this % F_im,   mold = this % u)
  allocate(this % F_im_o, mold = this % u)
  allocate(this % F_ex  , mold = this % u)
  allocate(this % F_ex_o, mold = this % u)

  allocate(this % w(0:n_sub,  2:this%n_stage))
  allocate(this % S(0:po, ne, 2:this%n_stage))

  ! quadrature weights .........................................................

  ! GLL points and weights in [-1,1]
  allocate(x(0:n_sub), source = GLL_Points(n_sub))
  allocate(w(0:n_sub), source = GLL_Weights(x))

  do i = 2, this % n_stage
  do j = 0, n_sub
    w_ji = 0
    do k = 0, n_sub
      ! tk = k-th GLL point in [0,cᵢ] mapped to [0,1]
      tk  = this % c(i) * HALF * (x(k) + 1)
      ! yk = value of j-th Lagrange polynomial at tk
      select case(set)
      case(1) ! use equidistant Lagrange polynomial in [0,1]
        yk = LagrangePolynomial(j, t, tk)
      case default ! use GLL Lagrange polynomial in [-1,1]
        yk = GLL_Polynomial(j, x, 2*tk-1)
      end select
      w_ji = w_ji + w(k) * yk
    end do
    ! weight of j-th SDC point in stage i scaled from [-1,1] to [0,cᵢ]
    this % w(j,i) = this % c(i) * HALF * w_ji
  end do
  end do

end subroutine Init_VARK_Propagator1D

!-------------------------------------------------------------------------------
!> VARK propagator

subroutine Propagator(vark, eop, dx, dt, M, wave, v, nu, bc, x, t0)
  class(VARK_Propagator1D),     intent(inout) :: vark
  class(CG_ElementOperators1D), intent(in)    :: eop       !< element operators
  real(RNP),                    intent(in)    :: dx        !< element length
  real(RNP),                    intent(in)    :: dt        !< time step size
  real(RNP),                    intent(in)    :: M(0:,:)   !< global mass matrix
  class(HarmonicWavePackage),   intent(in)    :: wave      !< exact wave solution
  real(RNP),                    intent(in)    :: v         !< convection velicity
  real(RNP),                    intent(in)    :: nu        !< diffusivity
  character,                    intent(in)    :: bc(:)     !< boundary conditions
  real(RNP),                    intent(in)    :: x(0:,:)   !< mesh points
  real(RNP),                    intent(in)    :: t0        !< time t₀

  real(RNP), allocatable :: f(:,:)
  real(RNP) :: c, t
  integer   :: i, j

  ! initialization .............................................................

  allocate(f, mold = M)

  associate( ns     => vark % n_stage , b_ex   => vark % b_ex    &
           , a_im   => vark % a_im    , a_ex   => vark % a_ex    &
           , F_im   => vark % F_im    , F_ex   => vark % F_ex    &
           , F_im_o => vark % F_im_o  , F_ex_o => vark % F_ex_o  &
           , u      => vark % u       , S      => vark % S       )

    ! stage 1 ..................................................................

    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t0, u(:,:,1), F_im(:,:,1))
    call GetLinearConvectionTerm(eop, v, bc, u(:,:,1), F_ex(:,:,1))

    ! stages 2 to s-1 ..........................................................

    do i = 2, ns-1

      t = t0 + vark%c(i) * dt

      f = M * u(:,:,1)
      do j = 1, i-1
        f = f + dt * a_im(i,j) * (F_im(:,:,j) - F_im_o(:,:,j))  &
              + dt * a_ex(i,j) * (F_ex(:,:,j) - F_ex_o(:,:,j))
      end do
      f = f - dt * a_im(i,i) * F_im_o(:,:,i) + S(:,:,i)

      if (nu > 0 .and. a_im(i,i) /= 0) then
        c = 1 / (dt * a_im(i,i))
        f = c * f
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u(:,:,i), f)
        call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u(:,:,i))
      else
        u(:,:,i) = f / M
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u(:,:,i))
      end if

      call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t, u(:,:,i), F_im(:,:,i))
      call GetLinearConvectionTerm(eop, v, bc, u(:,:,i), F_ex(:,:,i))

    end do

    ! last stage ...............................................................

    t = t0 + dt

    f = M * u(:,:,1)
    do j = 1, ns-1
      f = f + dt * a_im(ns,j) * (F_im(:,:,j) - F_im_o(:,:,j))  &
            + dt * b_ex(   j) * (F_ex(:,:,j) - F_ex_o(:,:,j))
    end do
    f = f - dt * a_im(ns,ns) * F_im_o(:,:,ns) + S(:,:,ns)

    ! implicit part
    if (nu > 0) then
      c = 1 / (dt * a_im(ns,ns))
      f = c * f
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u(:,:,ns), f)
      call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u(:,:,ns))
    else
      u(:,:,ns) = f / M
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u(:,:,ns))
    end if

    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t, u(:,:,ns), F_im(:,:,ns))
    call GetLinearConvectionTerm(eop, v, bc, u(:,:,ns), F_ex(:,:,ns))

  end associate

end subroutine Propagator

!===============================================================================
! VARK_SDC_Method1D

!-------------------------------------------------------------------------------
!> New 1D VARK-SDC method

type(VARK_SDC_Method1D) function New_VARK_SDC(opt, po, ne) result(this)
  class(VARK_SDC_Options1D), intent(in) :: opt !< VARK-SDC options
  integer,                   intent(in) :: po  !< polynomial order
  integer,                   intent(in) :: ne  !< number of elements

  call Init_VARK_SDC(this, opt, po, ne)

end function New_VARK_SDC

!-------------------------------------------------------------------------------
!> Initialize 1D VARK-SDC method

subroutine Init_VARK_SDC(this, opt, po, ne)
  class(VARK_SDC_Method1D),  intent(inout) :: this
  class(VARK_SDC_Options1D), intent(in)    :: opt !< VARK-SDC options
  integer,                   intent(in)    :: po  !< polynomial order
  integer,                   intent(in)    :: ne  !< number of elements

  call this % Init_SDC_Method(opt)

  this % vark = VARK_Propagator1D( po, ne, this%t, this%set &
                                 , so_p = opt % so_p        &
                                 , r_ex = opt % r_ex        )

end subroutine Init_VARK_SDC

!-------------------------------------------------------------------------------
!> IMEX-Euler SDC method with CG-SEM for 1D convection-diffusion equation

subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
  class(VARK_SDC_Method1D),     intent(inout) :: this     !< VARK SDC method
  class(CG_ElementOperators1D), intent(in)    :: eop      !< element operators
  real(RNP),                    intent(in)    :: dx       !< element length
  real(RNP),                    intent(in)    :: dt       !< time step size ∆t
  real(RNP),                    intent(in)    :: M(0:,:)  !< global mass matrix
  class(HarmonicWavePackage),   intent(in)    :: wave     !< exact wave solution
  real(RNP),                    intent(in)    :: v        !< convection velicity
  real(RNP),                    intent(in)    :: nu       !< diffusivity
  character,                    intent(in)    :: bc(:)    !< boundary conditions
  real(RNP),                    intent(in)    :: x(0:,:)  !< mesh points
  real(RNP),                    intent(in)    :: t0       !< time t₀ = t - ∆t
  real(RNP),                    intent(in)    :: u0(0:,:) !< solution u(t₀)
  real(RNP),                    intent(out)   :: u (0:,:) !< solution u(t)

  ! local variables  ...........................................................

  real(RNP), dimension(:,:,:), allocatable :: F

  integer :: po, ne, n_sub
  integer :: k, s1, s2

  ! initialization ...........................................................

  po = ubound(u, 1)
  ne = ubound(u, 2)
  n_sub = this % n_sub

  allocate(F(0:po, 1:ne, 0:n_sub), source = ZERO)

  associate(vark => this % vark)

    s1 = 1                    ! 1-st principal stage
    s2 = 2 + vark % helpers   ! 2-nd principal stage (skipping helpers)

    ! predictor ................................................................

    vark % u(:,:,1) = u0
    vark % F_im_o   = 0
    vark % F_ex_o   = 0
    vark % S        = 0

    call Propagator(vark, eop, dx, dt, M, wave, v, nu, bc, x, t0)

    ! corrector sweeps .........................................................

    do k = 1, this % n_sweep

      vark % F_im_o = vark % F_im
      vark % F_ex_o = vark % F_ex

      ! time derivative at SDC points
      F(:,:,0 ) = vark % F_im(:,:,s1 ) + vark % F_ex(:,:,s1 )
      F(:,:,1:) = vark % F_im(:,:,s2:) + vark % F_ex(:,:,s2:)

      ! integrals of time derivative over subintervals [0,tᵢ]
      call GetSubintegrals(dt, vark%w, F, vark%S)

      call Propagator(vark, eop, dx, dt, M, wave, v, nu, bc, x, t0)

    end do

    ! finalization .............................................................

    u = vark % u(:,:,vark % n_stage)

  end associate

end subroutine TimeStep

!-------------------------------------------------------------------------------
!> Computes the integrals of the time derivative F over subintervals [0,tᵢ]

subroutine GetSubintegrals(dt, w, F, S)
  real(RNP), intent(in)  :: dt        !< time step width
  real(RNP), intent(in)  :: w(:,:)    !< weights
  real(RNP), intent(in)  :: F(:,:,:)  !< time derivatives
  real(RNP), intent(out) :: S(:,:,:)  !< subinterval integrals

  integer :: nm, ni

  nm = size(F,1) * size(F,2)
  ni = size(F,3) ! = size(w,1)

  S = reshape(matmul(reshape(F,[nm,ni]), dt*w), shape(S))

end subroutine GetSubintegrals

!===============================================================================

end module CG_ConvDiff_1D__VARK_SDC
