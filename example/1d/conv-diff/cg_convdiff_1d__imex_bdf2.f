!> summary:  IMEX BDF2 with CG-SEM for 1D convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/09/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX BDF2 with CG-SEM for 1D convection-diffusion equation
!>
!> Advances the solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> according to the implicit-explicit second-order backward-difference method
!>
!>     u = 4/3 u₀ - 1/3 u₁ + 4/3 ∆t C(u₀) - 2/3 ∆t C(u₁) + 2/3 ∆t D(u)
!>
!> where
!>
!>     u₀ = u(t₀)
!>     u₁ = u(t₀ - ∆t)
!>     u  = u(t₀ + ∆t)
!>
!> C` and `D` are discretized using continuous spectral elements.
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_BDF2
  use Kind_Parameters, only: RNP
  use Constants, only: ONE
  use CG_ConvDiff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Solver_1D
  use Harmonic_Wave_Package
  implicit none
  private

  public :: IMEX_BDF2

contains

!-------------------------------------------------------------------------------
!> IMEX BDF2 method with CG-SEM for 1D convection-diffusion equation

subroutine IMEX_BDF2(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u1, u)
  class(CG_ElementOperators1D), intent(in)  :: eop      !< element operators
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
  real(RNP),                    intent(out) :: u (0:,:) !< solution u(t₀+∆t)

  real(RNP), allocatable :: ux(:,:), f(:,:)
  real(RNP) :: t, a, b0, b1, c0, c1, c

  t  =  t0 + dt

  a  =  ONE / 3
  b0 =  4 * a
  b1 = -1 * a
  c0 =  4 * a
  c1 = -2 * a

  allocate(f , mold = u)
  allocate(ux, mold = u)

  if (abs(v) > 0) then
    ux = c0*u0 + c1*u1
    call GetLinearConvectionTerm(eop, v, bc, ux, fc=f)
  else
    f = 0
  end if

  if (nu > 0) then
    c = ONE / (2 * a * dt)
    f = c * (M * (b0*u0 + b1*u1) + dt * f)
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
    call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
  else
    u = b0*u0 + b1*u1 + dt * f / M
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
  end if

end subroutine IMEX_BDF2

!===============================================================================

end module CG_ConvDiff_1D__IMEX_BDF2
