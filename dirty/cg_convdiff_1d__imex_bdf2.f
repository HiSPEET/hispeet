module CG_ConvDiff_1D__IMEX_BDF2
  use Kind_Parameters, only: RNP
  use Constants, only: HALF
  use CG_Conv_Diff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Elliptic_Solver_1D
  use Harmonic_Wave_Package
  implicit none
  private

  public :: IMEX_BDF2

contains

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
  real(RNP) :: b0, b1, c, t

  t  =  t0 + dt

  b0 =  4 * HALF / dt
  b1 = -1 * HALF / dt
  c  =  3 * HALF / dt

  allocate(ux, mold = u)
  allocate(f , mold = u)
  ux = 2*u0 - u1
  call GetLinearConvectionTerm(eop, v, bc, ux, f)

  if (nu > 0) then
    f = f + M * (b0*u0 + b1 *u1)
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
    call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u)
  else
    u = 1/c * (b0*u0 + b1 *u1 + f / M)
  end if

end subroutine IMEX_BDF2


end module CG_ConvDiff_1D__IMEX_BDF2
