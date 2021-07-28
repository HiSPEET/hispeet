!> summary:  IMEX-Euler Trapezoidal Rule for 1D convection-diffusion equation
!> author:
!> date:
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX-Trapezoidal Rule for 1D convection-diffusion equation, in two Variants
!>
!> **TBD:** Rework description, extend to SDC
!>
!> Advances the solution of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ -C(u) + D(u)
!>
!>
!>
!> where
!>
!>     u₀ = u(t₀)
!>     u  = u(t₀ + ∆t)
!>
!> `C` and `D` are discretized using continuous spectral elements.
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_TR_SDC
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use CG__Element_Operators__1D
  use CG_Condensed_Solver_1D
  use CG_ConvDiff_1D__Utils
  use Harmonic_Wave_Package
  use Spectral_Deferred_Correction
  implicit none
  private

  public :: IMEX_TR_SDC

contains

  !-------------------------------------------------------------------------------
  !> IMEX-TR SDC method with CG-SEM for 1D convection-diffusion equation

  subroutine IMEX_TR_SDC(sdc, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
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

    real(RNP), dimension(:),     allocatable :: ts       !< time sdc-subintervals
    real(RNP), dimension(:),     allocatable :: dts      !< time-step widths
    real(RNP), dimension(:,:,:), allocatable :: us       !< solution u(ts)
    real(RNP), dimension(:,:,:), allocatable :: S        !< subinterval integrals
    real(RNP), dimension(:,:,:), allocatable :: Fex_old  !< old time derivatives at t_{i-1}, explicit
    real(RNP), dimension(:,:,:), allocatable :: Fim_old  !< old time derivatives at t_{i-1}, implicit
    real(RNP), dimension(:,:,:), allocatable :: Fex      !< before: F^ex(t_{i}  , u_{i}^{k-1}), after: F^ex(t_{i}  , u_{i}^{k})
    real(RNP), dimension(:,:,:), allocatable :: Fim      !< before: F^im(t_{i}  , u_{i}^{k-1}), after: F^im(t_{i}  , u_{i}^{k})
    real(RNP), dimension(:,:,:), allocatable :: Fex_0_old    !< F^ex(t_{i-i}  , u_{i-1}^{k-1})
    real(RNP), dimension(:,:,:), allocatable :: Fim_0_old    !< F^im(t_{i-i}  , u_{i-1}^{k-1})

    integer :: po, ne, ns
    integer :: i, k

    ! initialization ...........................................................

    po = ubound(u, 1)
    ne = ubound(u, 2)
    ns = sdc % n_sub

    allocate(  ts (0:ns)) ! {tᵢ}
    allocate( dts (1:ns)) ! {∆tᵢ}
    ts  = sdc % IntermediateTimes(t0, dt)
    dts = ts(1:) - ts(0:ns-1)

    allocate(  us        (0:po, 1:ne, 0:ns), source = ZERO) ! {uᵢ}
    allocate(  S         (0:po, 1:ne, 1:ns), source = ZERO) ! {Sᵢ}
    allocate(  Fex_0_old (0:po, 1:ne, 0:ns), source = ZERO)
    allocate(  Fim_0_old (0:po, 1:ne, 0:ns), source = ZERO)
    allocate(  Fex      (0:po, 1:ne, 0:ns), source = ZERO)
    allocate(  Fim      (0:po, 1:ne, 0:ns), source = ZERO)

    ! predictor ..................................................................

    us(:,:,0) = u0

    ! F_ex(t0,u0) i.e. at i = 0
    !call ApplyBoundaryConditions(wave, v, nu, bc, x, t0, us(:,:,0)) ! us instead of u0
    call GetLinearConvectionTerm(eop, v, bc, u0, Fex(:,:,0))
    Fex(:,:,0) = Fex(:,:,0) / M

    ! F_im(t0,u0) i.e. at i = 0
    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t0, u0, Fim(:,:,0))
    Fim(:,:,0) = Fim(:,:,0) / M

    do i = 1, ns

     call Propagator( eop, dx, M, wave, v, nu, bc, x   &
                    , t0         =  ts(i-1)            &
                    , dt         =  dts(i)             &
                    , u0         =  us(:,:,i-1)        &
                    , u          =  us(:,:,i)          &
                    , Fex_0      =  Fex(:,:,i-1)       &
                    , Fim_0      =  Fim(:,:,i-1)       &
                    , Fex_0_old  =  Fex_old(:,:,i-1)   & ! ZERO
                    , Fim_0_old  =  Fim_old(:,:,i-1)   & ! ZERO
                    , Fex        =  Fex(:,:,i)         &
                    , Fim        =  Fim(:,:,i)         )

    end do

    ! correction sweeps ........................................................

    do k = 1, sdc % n_sweep

      ! subinterval integrals: Sᵏ⁻¹(tᵢ) = F · wᵢ
      call GetSubintegrals(dt, sdc%w_sub, Fex + Fim, S)

      ! save old time derivatives
      Fex_old = Fex
      Fim_old = Fim

      ! correction
      do i = 1, ns

       call Propagator( eop, dx, M, wave, v, nu, bc, x   &
                      , t0         =  ts(i-1)            &
                      , dt         =  dts(i)             &
                      , u0         =  us(:,:,i-1)        &
                      , u          =  us(:,:,i)          &
                      , Fex_0      =  Fex(:,:,i-1)       &
                      , Fim_0      =  Fim(:,:,i-1)       &
                      , Fex_0_old  =  Fex_old(:,:,i-1)   &
                      , Fim_0_old  =  Fim_old(:,:,i-1)   &
                      , Fex        =  Fex(:,:,i)         &
                      , Fim        =  Fim(:,:,i)         &
                      , S          =  S(:,:,i)           )

      end do

    end do

    ! finalization ...............................................................

    u = us(:,:,ns)

  end subroutine IMEX_TR_SDC


  !-----------------------------------------------------------------------------
  !> ...
  !>
  !> <documentation TBD>
  !>
  !> Acts as the
  !>   * SDC predictor, if S is absent
  !>   * SDC corrector, if S is present
  !>

  subroutine Propagator( eop, dx, M, wave, v, nu, bc, x, t0, dt, u0, u,  &
                         Fex_0, Fim_0, Fex_0_old, Fim_0_old, Fex, Fim, S )


  !!! IMPORTANT NOTE
  !!!   - correct is:  u          = u_{i}^{k-1}                   on input
  !!!                  u          = u_{i}^{k}                     on output
  !!!                  Fex_0      = F^ex(t_{i-1}, u_{i-1}^{k})    on input
  !!!                  Fex_0_old  = F^ex(t_{i-1}, u_{i-1}^{k-1})  on input
  !!!
  !!! etc.
  !!!
  !!! FIX:             k    -->   k-1
  !!!                  k+1  -->   k

    class(CG_ElementOperators_1D), intent(in)    :: eop          !< element operators
    class(HarmonicWavePackage),   intent(in)    :: wave         !< exact wave solution

    real(RNP), intent(in)    :: dx              !< element length
    real(RNP), intent(in)    :: M(0:,:)         !< global mass matrix
    real(RNP), intent(in)    :: v               !< convection velocity
    real(RNP), intent(in)    :: nu              !< diffusivity ν
    character, intent(in)    :: bc(:)           !< boundary conditions
    real(RNP), intent(in)    :: x(0:,:)         !< mesh points
    real(RNP), intent(in)    :: t0              !< time t₀ = t - ∆t = t_{i-1}
    real(RNP), intent(in)    :: dt              !< time step size ∆t
    real(RNP), intent(in)    :: u0(0:,:)        !< solution u(t₀) = u_{i-1}^{k}
    real(RNP), intent(out)   :: u (0:,:)        !< u_{i}^{k}
    real(RNP), intent(in)    :: Fex_0(0:,:)     !< F^ex(t_{i-1}, u_{i-1}^{k })
    real(RNP), intent(in)    :: Fim_0(0:,:)     !< F^im(t_{i-1}, u_{i-1}^{k })
    real(RNP), intent(in)    :: Fex_0_old(0:,:) !< F^ex(t_{i-1}, u_{i-1}^{k-1})
    real(RNP), intent(in)    :: Fim_0_old(0:,:) !< F^im(t_{i-1}, u_{i-1}^{k-1})
    real(RNP), intent(inout) :: Fex(0:,:)       !< before: F^ex(t_{i}  , u_{i}^{k-1}), after: F^ex(t_{i}  , u_{i}^{k})
    real(RNP), intent(inout) :: Fim(0:,:)       !< before: F^im(t_{i}  , u_{i}^{k-1}), after: F^im(t_{i}  , u_{i}^{k})
    real(RNP), optional, intent(in) :: S(0:,:)  !< subinterval integrals


    ! local variables  ...........................................................

    logical   :: corrector
    real(RNP) :: c, t2, t3
    real(RNP), dimension(:,:), allocatable :: u1, u2, Fex_2, w

    integer :: po, ne


    ! initialization ...........................................................

    corrector = present(S)

    po = ubound(u, 1)
    ne = ubound(u, 2)

    !t1 = t0    ! t1 not needed
    t2 = t0 + dt
    t3 = t2

    c = 2 / dt

    allocate(  u1    (0:po, 1:ne), source = ZERO)
    allocate(  u2    (0:po, 1:ne), source = ZERO)
    allocate(  Fex_2 (0:po, 1:ne), source = ZERO)
    allocate(  w     (0:po, 1:ne), source = ZERO)

    ! stage 1 ..................................................................

    u1 = u0

    ! stage 2 ..................................................................
    ! (variant 1)

    u2 = u0 + dt * (Fex_0 + Fim_0)

    if (corrector) then
      u2 = u2 - dt * (Fex_0_old + Fim_0_old) + S
    end if

    call ApplyBoundaryConditions(wave, v, nu, bc, x, t2, u2)

    ! F^ex(t2,u2)
    call GetLinearConvectionTerm(eop, v, bc, u2, Fex_2)
    Fex_2 = Fex_2 / M

    ! stage 3 ....................................................................
    ! u = u3

    ! w = sum of known contributions
    w = u0 + dt/2 * (Fex_0 + Fim_0 + Fex_2)
    if (corrector) then
      w  =  w  -  dt/2 * (Fex_0_old + Fim_0_old + Fex + Fim) + S
    end if

    ! solve elliptic problem (requires special treatment if nu = 0)
    if (nu > 0) then ! solve diffusive part

      ! scale w with M/dt to act as RHS for elliptic solver
      w = c * M * w

      ! add Neumann BC to RHS w, and apply Dirichlet BC in u
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t3, u, w)

      ! solve  (2/∆t)Mu + νLu = w  using static condensation -> u = u3
      call CondensedEllipticSolver(eop, dx, c, nu, bc, w, u)

    else ! no diffusion

      u = w
      ! enforce Dirichlet BC in u
      call ApplyBoundaryConditions(wave, v, nu, bc, x, t3, u)

    end if

    ! Fex = F^ex(t3,u3)
    call GetLinearConvectionTerm(eop, v, bc, u, Fex)
    Fex = Fex / M

    ! Fim = F^im(t3,u3)
    call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t3, u, Fim)
    Fim = Fim / M


  end subroutine Propagator

!-------------------------------------------------------------------------------
!> Computes the integrals of the time derivative F over subintervals [tᵢ₋₁,tᵢ]

subroutine GetSubintegrals(dt, w_sub, F, S)
  real(RNP), intent(in)  :: dt         !< time step width
  real(RNP), intent(in)  :: w_sub(:,:) !< weights
  real(RNP), intent(in)  :: F(:,:,:)   !< time derivatives
  real(RNP), intent(out) :: S(:,:,:)   !< subinterval integrals

  integer :: nm, ni

  nm = size(F,1) * size(F,2)
  ni = size(F,3)

  S = reshape(matmul(reshape(F,[nm,ni]), dt*w_sub), shape(S))

end subroutine GetSubintegrals

!===============================================================================

end module CG_ConvDiff_1D__IMEX_TR_SDC
