!> summary:  Variable additive RK method with CG-SEM for 1D convection-diffusion
!> author:   Joerg Stiller
!> date:     2019/03/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Variable additive RK method with CG-SEM for 1D convection-diffusion
!>
!> This module provides the type `VARK_Method1D` which extends the variable
!> additive Runge-Kutta methods defined in `VARK_Method` for advancing the
!> solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> in time. Spatial discretization is based on the continuous spectral-element
!> method using nodal base functions along with GLL quadrature.
!>
!> Typical usage:
!>
!>     type(VARK_Method1D) :: vark
!>
!>     vark = VARK_Method1D(po, ne, t, so_p [,r_ex])
!>     ! po   :  polynomial order and
!>     ! ne   :  number of elements
!>     ! t    :  equidistant (1) or GLL (2) points
!>     ! so_p :  principal stage order
!>     ! r_ex :  last coefficient of the stability function R_ex(z) of the
!>     !         explicit part
!>
!>     call vark % TimeStep(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
!>     ! for description of arguments see below
!>
!===============================================================================

module CG_ConvDiff_1D__VARK
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, HALF
  use Gauss_Jacobi
  use CG_ConvDiff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Solver_1D
  use VA_Runge_Kutta_Method
  use Harmonic_Wave_Package
  implicit none
  private

  public :: VARK_Method1D

  !-----------------------------------------------------------------------------
  !> Implementation of the VARK method for 1D convection-diffusion

  type, extends(VARK_Method) :: VARK_Method1D
    real(RNP), allocatable :: F_ex(:,:,:) !< explicit RHS per stage
    real(RNP), allocatable :: F_im(:,:,:) !< implicit RHS per stage
  contains
    generic :: Init_ConvDiff_VARK => Init_VARK_t, Init_VARK_s
    procedure, private :: Init_VARK_t
    procedure, private :: Init_VARK_s
    procedure :: TimeStep
  end type VARK_Method1D

  ! constructors
  interface VARK_Method1D
    module procedure New_VARK_t
    module procedure New_VARK_s
  end interface

contains

!-------------------------------------------------------------------------------
!> New ConvDiff VARK method from given points

type(VARK_Method1D) function New_VARK_t(po, ne, t, o_ps, r_ex) result(this)
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  real(RNP), intent(in) :: t(:)  !< principal nodes
  integer,   intent(in) :: o_ps  !< stage order at principal nodes
  real(RNP), intent(in) :: r_ex  !< last coefficient of R_ex(z)
  optional :: r_ex

  call Init_VARK_t(this, po, ne, t, o_ps, r_ex)

end function New_VARK_t

!-------------------------------------------------------------------------------
!> New ConvDiff VARK method from specified point set

type(VARK_Method1D) function New_VARK_s(po, ne, set, np, o_ps, r_ex) &
    result(this)
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  integer,   intent(in) :: set   !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: np    !< number of principal nodes
  integer,   intent(in) :: o_ps  !< stage order at principal nodes
  real(RNP), intent(in) :: r_ex  !< last coefficient of R_ex(z)
  optional :: r_ex

  call Init_VARK_s(this, po, ne, set, np, o_ps, r_ex)

end function New_VARK_s

!-------------------------------------------------------------------------------
!> Initialize ConvDiff VARK method for given points

subroutine Init_VARK_t(this, po, ne, t, o_ps, r_ex)
  class(VARK_Method1D), intent(inout) :: this
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  real(RNP), intent(in) :: t(:)  !< principal nodes
  integer,   intent(in) :: o_ps  !< stage order at principal nodes
  real(RNP), intent(in) :: r_ex  !< last coefficient of R_ex(z)
  optional :: r_ex

  ! initialize VARK
  call this % Init_VARK_Method(t, o_ps, r_ex)

  ! workspace
  allocate(this % F_ex(0:po, ne, this%n_stage-1))
  allocate(this % F_im(0:po, ne, this%n_stage  ))

end subroutine Init_VARK_t

!-------------------------------------------------------------------------------
!> Initialize ConvDiff VARK method for specified point set

subroutine Init_VARK_s(this, po, ne, set, n_ps, o_ps, r_ex)
  class(VARK_Method1D), intent(inout) :: this
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  integer,   intent(in) :: set   !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: n_ps  !< number of principal nodes
  integer,   intent(in) :: o_ps  !< stage order at principal nodes
  real(RNP), intent(in) :: r_ex  !< last coefficient of R_ex(z)
  optional :: r_ex

  real(RNP) :: t(n_ps)
  integer   :: i, ni

  ! number of intervals
  ni = n_ps - 1

  select case(set)
  case(1) ! equidistant
    t = [ ((i-1)*ONE/ni, i = 1,n_ps) ]
  case default ! GLL
    t = HALF * (GLL_Points(ni) + ONE)
  end select

  call Init_VARK_t(this, po, ne, t, o_ps, r_ex)

end subroutine Init_VARK_s

!-------------------------------------------------------------------------------
!> Performs a single VARK time step

subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
  class(VARK_Method1D),         intent(inout) :: this
  class(CG_ElementOperators1D), intent(in)    :: eop      !< element operators
  real(RNP),                    intent(in)    :: dx       !< element length
  real(RNP),                    intent(in)    :: dt       !< time step size
  real(RNP),                    intent(in)    :: M(0:,:)  !< global mass matrix
  class(HarmonicWavePackage),   intent(in)    :: wave     !< exact wave solution
  real(RNP),                    intent(in)    :: v        !< convection velicity
  real(RNP),                    intent(in)    :: nu       !< diffusivity
  character,                    intent(in)    :: bc(:)    !< boundary conditions
  real(RNP),                    intent(in)    :: x(0:,:)  !< mesh points
  real(RNP),                    intent(in)    :: t0       !< time t₀
  real(RNP),                    intent(in)    :: u0(0:,:) !< solution u(t₀)
  real(RNP),                    intent(out)   :: u (0:,:) !< solution u(t₀+∆t)

  real(RNP), allocatable :: f(:,:)
  real(RNP) :: c, t
  integer   :: i, j

  ! initialization .............................................................

  allocate(f, mold = u)

  associate( ns    => this % n_stage, b_ex => this % b_ex  &
           , a_im  => this % a_im   , a_ex => this % a_ex  &
           , F_im  => this % F_im   , F_ex => this % F_ex  )

    ! stage 1 ..................................................................

    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t0, u0, F_im(:,:,1))
    call GetLinearConvectionTerm(eop, v, bc, u0, F_ex(:,:,1))

    ! stages 2 to s-1 ..........................................................

    do i = 2, ns-1

      t = t0 + this%c(i) * dt

      f = M * u0
      do j = 1, i-1
        f = f + dt * a_im(i,j) * F_im(:,:,j) + dt * a_ex(i,j) * F_ex(:,:,j)
      end do

      if (nu > 0 .and. a_im(i,i) /= 0) then
        c = 1 / (dt * a_im(i,i))
        f = c * f
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
        call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
      else
        u = f / M
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
      end if

      call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t, u, F_im(:,:,i))
      call GetLinearConvectionTerm(eop, v, bc, u, F_ex(:,:,i))

    end do

    ! last stage ...............................................................

    t = t0 + dt

    f = M * u0
    do i = 1, ns-1
      f = f + dt * a_im(ns,i) * F_im(:,:,i) + dt * b_ex(i) * F_ex(:,:,i)
    end do

    ! implicit part
    if (nu > 0) then
      c = 1 / (dt * a_im(ns,ns))
      f = c * f
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
      call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
    else
      u = f / M
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
    end if

  end associate

end subroutine TimeStep

!===============================================================================

end module CG_ConvDiff_1D__VARK
