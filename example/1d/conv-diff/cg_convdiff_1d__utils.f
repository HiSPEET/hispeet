!> summary:  Linear convection and diffusion terms for 1D CG-SEM
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Linear convection and diffusion terms for 1D CG-SEM
!>
!> The modules provides procedures for evaluating the convection and diffusion
!> terms for the 1D convection-diffusion equation
!>
!>     ∂u/∂t + v ∂u/∂x = nu  ∂²u/∂x²
!>
!> with constant velocity `v` and diffusivity `nu`. Discretization is based on
!> the continuous Galerkin spectral-element method.
!>
!> @note
!> Implementation restricted to GLL Langrange basis
!> @endnote
!===============================================================================

module CG_ConvDiff_1D__Utils
  use Kind_Parameters, only: RNP
  use Standard_Operators_1D
  use Harmonic_Wave_Package
  use CG_Utilities_1D, only: Assembly
  implicit none
  private

  public :: GetLinearConvectionTerm
  public :: GetDiffusionTerm
  public :: ApplyBoundaryConditions

contains

!-------------------------------------------------------------------------------
!> Linear convection term for constant velocity
!>
!> Computes for all global basis functions φᵢ
!>
!>     fcᵢ = -∫ φᵢ v ∂u/∂x dx
!>
!> Periodic boundary conditions are assumed if `all(bc ='P')` and taken into
!> account when assembling the element integrals, correspondingly.
!> The result is returned as an element variable.

subroutine GetLinearConvectionTerm(sop, v, bc, u, fc)
  class(StandardOperators1D), intent(in)  :: sop      !< standard operators
  real(RNP),                  intent(in)  :: v        !< velocity
  character,                  intent(in)  :: bc(:)    !< left/right BC
  real(RNP),                  intent(in)  :: u (0:,:) !< solution
  real(RNP),                  intent(out) :: fc(0:,:) !< convection term

  real(RNP), allocatable :: Ce(:,:)
  integer :: j

  associate(po => sop%po, w => sop%w, Ds => sop%D)

    allocate(Ce(0:po,0:po))
    do j = 0, po
      Ce(:,j) = -v * w * Ds(:,j)
    end do

    fc = matmul(Ce, u)

    call Assembly(bc, fc)

  end associate

end subroutine GetLinearConvectionTerm

!-------------------------------------------------------------------------------
!> Diffusion term
!>
!> Computes for all global basis functions φᵢ
!>
!>     fdᵢ = -∫ ∂φᵢ/∂x nu ∂u/∂x dx

subroutine GetDiffusionTerm(sop, dx, wave, v, nu, bc, x, t, u, fd)
  class(StandardOperators1D), intent(in)  :: sop      !< standard operators
  real(RNP),                  intent(in)  :: dx       !< element length
  class(HarmonicWavePackage), intent(in)  :: wave     !< exact wave solution
  real(RNP),                  intent(in)  :: v        !< velocity
  real(RNP),                  intent(in)  :: nu       !< diffusivity
  character,                  intent(in)  :: bc(:)    !< left/right BC
  real(RNP),                  intent(in)  :: x(0:,:)  !< mesh points
  real(RNP),                  intent(in)  :: t        !< time
  real(RNP),                  intent(in)  :: u (0:,:) !< solution
  real(RNP),                  intent(out) :: fd(0:,:) !< diffusion term

  if (nu > 0) then

    fd = -nu * 2/dx * matmul(sop%L, u)
    call Assembly(bc, fd)
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, f=fd)

  else

    fd = 0

  end if

end subroutine GetDiffusionTerm

!-------------------------------------------------------------------------------
!> Apply boundary conditions to solution `u` and RHS `f`
!>
!> Injects Dirichlet conditions into the solution vector `u` and adds Neumann
!> conditions `nu ∂u/∂x` to the right hand side `f`.

subroutine ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
  class(HarmonicWavePackage), intent(in)    :: wave     !< exact wave solution
  real(RNP),                  intent(in)    :: v        !< velocity
  real(RNP),                  intent(in)    :: nu       !< diffusivity
  character,                  intent(in)    :: bc(:)    !< left/right BC
  real(RNP),                  intent(in)    :: x(0:,:)  !< mesh points
  real(RNP),                  intent(in)    :: t        !< time
  real(RNP),        optional, intent(inout) :: u(0:,:)  !< solution
  real(RNP),        optional, intent(inout) :: f(0:,:)  !< RHS

  real(RNP) :: dx_u
  integer   :: po, ne

  po = ubound(x,1)
  ne = ubound(x,2)

  if (present(u)) then
    if (bc(1) == 'D') call wave % GetAmplitude(v, nu, x(0,1), t, u(0,1))
    if (bc(2) == 'D') call wave % GetAmplitude(v, nu, x(po,ne), t, u(po,ne))
  end if

  if (present(f)) then
    if (bc(1) == 'N') then
      call wave % GetDerivative(v, nu, x(0,1), t, dx_u)
      f(0,1) = f(0,1) - nu * dx_u
    end if
    if (bc(2) == 'N') then
      call wave % GetDerivative(v, nu, x(po,ne), t, dx_u)
      f(po,ne) = f(po,ne) + nu * dx_u
    end if
  end if

end subroutine ApplyBoundaryConditions

!===============================================================================

end module CG_ConvDiff_1D__Utils
