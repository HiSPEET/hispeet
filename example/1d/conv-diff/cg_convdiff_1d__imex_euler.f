!> summary:  IMEX Euler with CG-SEM for 1D convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX Euler with CG-SEM for 1D convection-diffusion equation
!>
!> Advances the solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> according to the implicit-explicit Euler method
!>
!>     u = u₀ + ∆t C(u₀) - ∆t D(u)
!>
!> where
!>
!>     u₀ = u(t₀)
!>     u  = u(t₀ + ∆t)
!>
!> C` and `D` are discretized using continuous spectral elements.
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_Euler
  use Kind_Parameters, only: RNP
  use CG_ConvDiff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Solver_1D
  use Harmonic_Wave_Package
  implicit none
  private

  public :: IMEX_Euler

contains

!-------------------------------------------------------------------------------
!> IMEX Euler method with CG-SEM for 1D convection-diffusion equation

subroutine IMEX_Euler(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
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
  real(RNP),                    intent(out) :: u (0:,:) !< solution u(t₀+∆t)

  real(RNP), allocatable :: f(:,:)
  real(RNP) :: c, t

  t = t0 + dt

  allocate(f, mold = u)
  if (abs(v) > 0) then
    call GetLinearConvectionTerm(eop, v, bc, u0, f)
  else
    f = 0
  end if

  if (nu > 0) then
    c = 1 / dt
    f = f + c * M * u0
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
    call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u, standby = .true.)
  else
    u = u0 + dt * f / M
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
  end if

end subroutine IMEX_Euler

!===============================================================================

end module CG_ConvDiff_1D__IMEX_Euler
