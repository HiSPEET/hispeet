!> summary:  Variable additive RK method with CG-SEM for 1D convection-diffusion
!> author:   Joerg Stiller
!> date:     2019/03/24
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Variable additive RK method with CG-SEM for 1D convection-diffusion
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
!> Typical usage:
!>
!>     type(ConvDiff_VARK) :: vark
!>
!>     vark = ConvDiff_VARK(po, ne, t, sh [,r])
!>     ! po :  polynomial order and
!>     ! ne :  number of elements
!>     ! t  :  equidistant (1) or GLL (2) points
!>     ! sh :  number of helper stages
!>     ! r  :  last coefficient of the stability function R_ex(z) of the
!>     !       explicit part
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

  public :: ConvDiff_VARK

  !-----------------------------------------------------------------------------
  !> Implementation of the VARK method for 1D convection-diffusion

  type, extends(VARK_Method) :: ConvDiff_VARK
    real(RNP), allocatable :: f_ex(:,:,:) !< explicit RHS per stage
    real(RNP), allocatable :: f_im(:,:,:) !< implicit RHS per stage
  contains
    generic :: Init_ConvDiff_VARK => Init_VARK_t, Init_VARK_s
    procedure, private :: Init_VARK_t
    procedure, private :: Init_VARK_s
    procedure :: TimeStep
  end type ConvDiff_VARK

  ! constructors
  interface ConvDiff_VARK
    module procedure New_VARK_t
    module procedure New_VARK_s
  end interface

contains

!-------------------------------------------------------------------------------
!> New ConvDiff VARK method from given points

type(ConvDiff_VARK) function New_VARK_t(po, ne, t, sh, r) result(this)
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  real(RNP), intent(in) :: t(0:) !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: sh    !< number of helper stages
  real(RNP), intent(in) :: r     !< last coefficient of R_ex(z)
  optional :: r

  call Init_VARK_t(this, po, ne, t, sh, r)

end function New_VARK_t

!-------------------------------------------------------------------------------
!> New ConvDiff VARK method from specified point set

type(ConvDiff_VARK) function New_VARK_s(po, ne, set, sp, sh, r) result(this)
  integer,   intent(in) :: po  !< polynomial order
  integer,   intent(in) :: ne  !< number of elements
  integer,   intent(in) :: set !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: sp  !< number of principal stages
  integer,   intent(in) :: sh  !< number of helper stages
  real(RNP), intent(in) :: r   !< last coefficient of R_ex(z)
  optional :: r

  call Init_VARK_s(this, po, ne, set, sp, sh, r)

end function New_VARK_s

!-------------------------------------------------------------------------------
!> Initialize ConvDiff VARK method for given points

subroutine Init_VARK_t(this, po, ne, t, sh, r)
  class(ConvDiff_VARK), intent(inout) :: this
  integer,   intent(in) :: po    !< polynomial order
  integer,   intent(in) :: ne    !< number of elements
  real(RNP), intent(in) :: t(0:) !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: sh    !< number of helper stages
  real(RNP), intent(in) :: r     !< last coefficient of R_ex(z)
  optional :: r

  ! initialize VARK
  call this % Init_VARK_Method(t, sh, r)

  ! workspace
  allocate(this % f_ex(0:po, ne, this%s-1))
  allocate(this % f_im(0:po, ne, this%s  ))

end subroutine Init_VARK_t

!-------------------------------------------------------------------------------
!> Initialize ConvDiff VARK method for specified point set

subroutine Init_VARK_s(this, po, ne, set, sp, sh, r)
  class(ConvDiff_VARK), intent(inout) :: this
  integer,   intent(in) :: po  !< polynomial order
  integer,   intent(in) :: ne  !< number of elements
  integer,   intent(in) :: set !< equidistant (1) or GLL (2) points
  integer,   intent(in) :: sp  !< number of principal stages
  integer,   intent(in) :: sh  !< number of helper stages
  real(RNP), intent(in) :: r   !< last coefficient of R_ex(z)
  optional :: r

  real(RNP) :: t(0:sp-1)
  integer   :: i, n

  n = ubound(t,1)

  select case(set)
  case(1) ! equidistant
    t = [ (i*ONE/n, i = 0,n) ]
  case default ! GLL
    t = HALF * (GLL_Points(n) + ONE)
  end select

  call Init_VARK_t(this, po, ne, t, sh, r)

end subroutine Init_VARK_s

!-------------------------------------------------------------------------------
!> Performs a single VARRK time step

subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
  class(ConvDiff_VARK),         intent(inout) :: this
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

  associate( s    => this % s    , b_ex => this % b_ex  &
           , a_im => this % a_im , a_ex => this % a_ex  &
           , f_im => this % f_im , f_ex => this % f_ex  )

    ! stage 1 ..................................................................

    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t0, u0, f_im(:,:,1))
    call GetLinearConvectionTerm(eop, v, bc, u0, f_ex(:,:,1))

    ! stages 2 to s-1 ..........................................................

    do i = 2, s-1

      t = t0 + this%c(i) * dt

      f = M * u0
      do j = 1, i-1
        f = f + dt * a_im(i,j) * f_im(:,:,j) + dt * a_ex(i,j) * f_ex(:,:,j)
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

      call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t, u, f_im(:,:,i))
      call GetLinearConvectionTerm(eop, v, bc, u, f_ex(:,:,i))

    end do

    ! last stage ...............................................................

    t = t0 + dt

    f = M * u0
    do i = 1, s-1
      f = f + dt * a_im(s,i) * f_im(:,:,i) + dt * b_ex(i) * f_ex(:,:,i)
    end do

    ! implicit part
    if (nu > 0) then
      c = 1 / (dt * a_im(s,s))
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
