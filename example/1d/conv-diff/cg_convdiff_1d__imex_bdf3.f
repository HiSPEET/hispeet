!> summary:  IMEX BDF3 with CG-SEM for 1D convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/10/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX BDF3 with CG-SEM for 1D convection-diffusion equation
!>
!> Advances the solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> according to the implicit-explicit second-order backward-difference method
!>
!>     u = 18/11 u₀ - 9/11 u₁ + 2/11 u₂
!>       + 18/11 ∆t C(u₀) - 18/11 ∆t C(u₁) + 6/11∆t C(u₂)
!>       +  6/11 ∆t D(u)
!>
!> where
!>
!>     u₀ = u(t₀)
!>     u₁ = u(t₀ - ∆t)
!>     u₂ = u(t₀ - 2 ∆t)
!>     u  = u(t₀ + ∆t)
!>
!> C` and `D` are discretized using continuous spectral elements.
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_BDF3
  use Kind_Parameters, only: RNP
  use Constants, only: ONE
  use CG_ConvDiff_1D__Utils
  use CG_Elliptic_Operator_1D
  use CG_Condensed_Solver_1D
  use Harmonic_Wave_Package
  implicit none
  private

  public :: IMEX_BDF3

contains

!-------------------------------------------------------------------------------
!> IMEX BDF3 method with CG-SEM for 1D convection-diffusion equation

subroutine IMEX_BDF3(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u1, u2, u)
  class(CG_EllipticOperator1D), intent(in)  :: eop      !< element operators
  real(RNP),                    intent(in)  :: dx       !< element length
  real(RNP),                    intent(in)  :: dt       !< time step size
  real(RNP),                    intent(in)  :: M(0:,:)  !< global mass matrix
  class(HarmonicWavePackage),   intent(in)  :: wave     !< exact wave solution
  real(RNP),                    intent(in)  :: v        !< convection velicity
  real(RNP),                    intent(in)  :: nu       !< diffusivity
  character,                    intent(in)  :: bc(:)    !< boundary conditions
  real(RNP),                    intent(in)  :: x(0:,:)  !< mesh points
  real(RNP),                    intent(in)  :: t0       !< time t₀
  real(RNP),                    intent(in)  :: u0(0:,:) !< solution u(t₀)
  real(RNP),                    intent(in)  :: u1(0:,:) !< solution u(t₀-∆t)
  real(RNP),                    intent(in)  :: u2(0:,:) !< solution u(t₀-2∆t)
  real(RNP),                    intent(out) :: u (0:,:) !< solution u(t₀+∆t)

  real(RNP), allocatable :: ux(:,:), f(:,:)
  real(RNP) :: t, a, b0, b1, b2, c0, c1, c2, c

  t  =  t0 + dt

  a  =  ONE / 11

  b0 =  18 * a
  b1 =  -9 * a
  b2 =   2 * a

  c0 =  18 * a
  c1 = -18 * a
  c2 =   6 * a

  allocate(f , mold = u)
  allocate(ux, mold = u)

  if (abs(v) > 0) then
    ux = c0*u0 + c1*u1 + c2*u2
    call GetLinearConvectionTerm(eop, v, bc, ux, fc=f)
  else
    f = 0
  end if

  if (nu > 0) then
    c = ONE / (6 * a * dt)
    f = c * (M * (b0*u0 + b1*u1 + b2*u2) + dt * f)
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
    call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
  else
    u = b0*u0 + b1*u1 + b2*u2 + dt * f / M
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
  end if

end subroutine IMEX_BDF3

!===============================================================================

end module CG_ConvDiff_1D__IMEX_BDF3
