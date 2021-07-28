!> summary:  IMEX-Euler SDC with CG-SEM for 1D convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX-Euler SDC with CG-SEM for 1D convection-diffusion equation
!>
!> **TBD:** Rework description, extend to SDC
!>
!> Advances the solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ -C(u) + D(u)
!>
!> according to the implicit-explicit Euler method
!>
!>     u = u₀ - ∆t C(u₀) + ∆t D(u)
!>
!> where
!>
!>     u₀ = u(t₀)
!>     u  = u(t₀ + ∆t)
!>
!> `C` and `D` are discretized using continuous spectral elements.
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_Euler_SDC
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use CG__Element_Operators__1D
  use CG_Condensed_Solver_1D
  use CG_ConvDiff_1D__Utils
  use Harmonic_Wave_Package
  use Spectral_Deferred_Correction
  implicit none
  private

  public :: IMEX_Euler_SDC

contains

!-------------------------------------------------------------------------------
!> IMEX-Euler SDC method with CG-SEM for 1D convection-diffusion equation

subroutine IMEX_Euler_SDC(sdc, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
  class(SDC_Method),            intent(in)    :: sdc      !< SDC parameters
  class(CG_ElementOperators_1D), intent(in)    :: eop      !< element operators
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

  real(RNP), dimension(:),     allocatable :: ts, dts
  real(RNP), dimension(:,:,:), allocatable :: us, F, H, S

  integer :: po, ne, ns
  integer :: i, k

  ! initialization ...........................................................

  po = ubound(u, 1)
  ne = ubound(u, 2)
  ns = sdc % n_sub

  allocate(  ts (0:ns)) ! {tᵢ}
  allocate( dts (1:ns)) ! {∆tᵢ}
  ts (0:ns) = sdc % IntermediateTimes(t0, dt)
  dts(1:ns) = ts(1:) - ts(0:ns-1)

  allocate(  us (0:po, 1:ne, 0:ns), source = ZERO) ! {uᵢ}
  allocate(  F  (0:po, 1:ne, 0:ns), source = ZERO) ! {Fᵢ}
  allocate(  H  (0:po, 1:ne, 1:ns), source = ZERO) ! {Hᵢ}
  allocate(  S  (0:po, 1:ne, 1:ns), source = ZERO) ! {Sᵢ}

  ! start values and time derivatives at t₀
  us(:,:,0) = u0
  call GetTimeDerivative(eop, dx, M, wave, v, nu, bc, x, t0, u0, F(:,:,0))

  ! predictor ..................................................................

  do i = 1, ns

    call Propagator( eop, dx, M, wave, v, nu, bc, x &
                   , t0 = ts(i-1)                   &
                   , dt = dts(i)                    &
                   , u0 = us(:,:,i-1)               &
                   , u  = us(:,:,i)                 &
                   , H  = H(:,:,i)                  )
  end do

  ! corrector sweeps ...........................................................

  do k = 1, sdc % n_sweep

    ! time derivatives: Fᵢ = Fᵏ⁻¹(tᵢ) = F(tᵢ, uᵏ⁻¹(tᵢ))
    do i = 1, ns
      call GetTimeDerivative( eop, dx, M, wave, v, nu, bc, x, &
                              ts(i), us(:,:,i), F(:,:,i)      )
    end do

    ! subinterval integrals: Sᵏ⁻¹(tᵢ) = F · wᵢ
    call GetSubintegrals(dt, sdc%w_sub, F, S)

    ! correction
    do i = 1, ns
      call Propagator( eop, dx, M, wave, v, nu, bc, x &
                     , t0 = ts(i-1)                   &
                     , dt = dts(i)                    &
                     , u0 = us(:,:,i-1)               &
                     , u  = us(:,:,i)                 &
                     , H  = H(:,:,i)                  &
                     , S  = S(:,:,i)                  )
    end do

  end do

  ! finalization ...............................................................

  u = us(:,:,ns)

end subroutine IMEX_Euler_SDC

!-------------------------------------------------------------------------------
!> Computes the integrals of the time derivative F over subintervals [tᵢ₋₁,tᵢ]

subroutine GetSubintegrals(dt, w_sub, F, S)
  real(RNP), intent(in)  :: dt         !< time step width
  real(RNP), intent(in)  :: w_sub(:,:) !< weights
  real(RNP), intent(in)  :: F(:,:,:)   !< time derivatives
  real(RNP), intent(out) :: S(:,:,:)   !< subinterval integrals

  integer :: nm, ni

  nm = size(F,1) * size(F,2)
  ni = size(F,3) ! = size(w,1)

  S = reshape(matmul(reshape(F,[nm,ni]), dt*w_sub), shape(S))

end subroutine GetSubintegrals


!-------------------------------------------------------------------------------
!> IMEX-Euler SDC propagator with CG-SEM for 1D convection-diffusion equation
!>
!> Note that with SDC
!>
!>     t0  =  t₀  = tᵢ₋₁
!>     dt  =  ∆tᵢ
!>     t   =  tᵢ
!>     u0  =  uᵏ(tᵢ₋₁)
!>     u   =  uᵏ(tᵢ)

subroutine Propagator(eop, dx, M, wave, v, nu, bc, x, t0, dt, u0, u, H, S)
  class(CG_ElementOperators_1D), intent(in)    :: eop      !< element operators
  real(RNP),                    intent(in)    :: dx       !< element length
  real(RNP),                    intent(in)    :: M(0:,:)  !< global mass matrix
  class(HarmonicWavePackage),   intent(in)    :: wave     !< exact wave solution
  real(RNP),                    intent(in)    :: v        !< convection velicity
  real(RNP),                    intent(in)    :: nu       !< diffusivity
  character,                    intent(in)    :: bc(:)    !< boundary conditions
  real(RNP),                    intent(in)    :: x(0:,:)  !< mesh points
  real(RNP),                    intent(in)    :: t0       !< time t₀ = tᵢ₋₁
  real(RNP),                    intent(in)    :: dt       !< time step size ∆tᵢ
  real(RNP),                    intent(in)    :: u0(0:,:) !< solution uᵏ(tᵢ₋₁)
  real(RNP),                    intent(out)   :: u (0:,:) !< solution uᵏ(tᵢ)
  real(RNP),                    intent(inout) :: H (0:,:) !< Hᵏ⁻¹(tᵢ) → Hᵏ(tᵢ)
  real(RNP),          optional, intent(in)    :: S (0:,:) !< Sᵏ⁻¹(tᵢ)

  real(RNP), allocatable :: w(:,:)
  real(RNP) :: c, t

  allocate(w, mold = u)
  t = t0 + dt
  c = 1 / dt

  ! w = M ũ/∆tᵢ = M (uᵏ(tᵢ₋₁) + Hᵏ⁻¹ + Sᵏ⁻¹)/∆tᵢ - C uᵏ(tᵢ₋₁)  .................

  ! w = -C uᵏ(tᵢ₋₁)
  if (abs(v) > 0) then
    call GetLinearConvectionTerm(eop, v, bc, u0, w)
  else
    w = 0
  end if

  ! w += M (uᵏ(tᵢ₋₁) + Hᵏ⁻¹ + Sᵏ⁻¹)/∆tᵢ
  if (present(S)) then ! corrector
    w = w + c * M * (u0 - H + S)
  else ! predictor
    w = w + c * M * u0
  end if

  ! uᵏ(tᵢ) ← Solve(∆tᵢ⁻¹ Mu + ν Lu = w)  .......................................

  if (nu > 0) then
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, w)
    call CondensedEllipticSolver(eop, dx, c, nu, bc, w, u)
  else
    u = dt * w / M
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
  end if

  ! Hᵏ(tᵢ) = uᵏ(tᵢ) - uᵏ(tᵢ₋₁) + Hᵏ⁻¹(tᵢ) - Sᵏ⁻¹(tᵢ)  ..........................

  if (present(S)) then
    H = u - u0 + H - S
  else
    H = u - u0
  end if

end subroutine Propagator

!===============================================================================

end module CG_ConvDiff_1D__IMEX_Euler_SDC
