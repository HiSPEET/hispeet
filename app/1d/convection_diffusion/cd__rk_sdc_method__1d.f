!> summary:  IMEX Runge-Kutta method with CG-SEM for 1D convection-diffusion
!> author:   Joerg Stiller
!> date:     2018/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX Runge-Kutta method with CG-SEM for 1D convection-diffusion
!>
!> This module provides the type `CD_IMEX_RK_Method_1D` which extends the IMEX
!> Runge-Kutta methods defined in `CD_IMEX_RK_Method_1D` for advancing the solution
!> of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> in time. Spatial discretization is based on the continuous spectral-element
!> method using nodal base functions along with Lobatto quadrature.
!>
!> Typical usage:
!>
!>     type(CD_IMEX_RK_Method_1D) :: imex_rk
!>
!>     imex_rk = CD_IMEX_RK_Method_1D(ns, method, po, ne)
!>     ! ns     :  number of stages
!>     ! method :  method, if several with ns stages exist (optional)
!>     ! po     :  polynomial order and
!>     ! ne     :  number of elements
!>
!>     call imex_rk % TimeStep(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
!>     ! for description of arguments see below
!>
!===============================================================================

module CD__RK_SDC_Method__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use CD__Utils__1D
  use CG__Element_Operators__1D
  use CG__Condensed_Solver__1D
  use IMEX_Runge_Kutta_Method
  use Harmonic_Wave_Package
  use Spectral_Deferred_Correction
  implicit none
  private

  public :: CD_RK_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Implementation of the IMEX-RK-SDC method for 1D convection-diffusion
  !>
  !>### Dimensions
  !>
  !>#### Solution and RHS
  !>
  !>      u(0:po, 1:ne, 1:n_stage, 1:n_sub)
  !>
  !> where
  !>   - `po`       is the polynomial order,
  !>   - `ne`       is the number of elements
  !>   - `n_stage`  is the number of RK stages
  !>   - `n_sub`    is the number of SDC subintervals
  !>
  !> All the `F_*` are dimensioned equivalently.
  !>
  !>####  SDC subinterval integrals
  !>
  !> The integrals over full SDC subintervals are defined as follows
  !>
  !>   \[
  !>      \int_{t_{m-1}}^{t_m} (F^{im} + F^{ex}) dt
  !>      \approx
  !>      \Delta t \sum_{j=0}^{M} (F^{im}(t_j) + F^{ex}(t_j)) w^s_{j,m}
  !>   \]
  !>
  !> where
  !>   - \( \Delta t = t^{n+1} - t^{n} \)   is the full time step
  !>   - \( \{ t_j       \}            \)   are the SDC points
  !>   - \( \{ w^s_{j,m} \}            \)   are the weights for SDC point \{t_j\}
  !>                                        for subinterval \( (t_{m-1}, t_m) \)
  !>
  !> Given an object of the type `SDC_Method`, the weights correspond to its
  !>  `w_sub` component.
  !>
  !> The result is stored in `S(0:po, 1:ne, 1:n_sub)`.
  !>
  !>#### RK subinterval integrals
  !>
  !> For RK-SDC one needs additionally the integrals
  !>
  !>   \[
  !>      \int_{t_{m-1}}^{\tilde{t}_i} (F^{im} + F^{ex}) dt
  !>      \approx
  !>      \Delta t \sum_{j=0}^{M} (F^{im}(t_j) + F^{ex}(t_j)) w^{RK}_{j,i,m}
  !>   \]
  !>
  !> where
  !>   - \( \Delta t = t^{n+1} - t^{n} \)   is the full time step
  !>   - \( \{ t_{m-1}        \}       \)   is the starting point of the SDC interval
  !>   - \( \{ \tilde{t}_i    \}       \)   are the RK points within SDC interval `m`
  !>   - \( \{ w^{RK}_{j,i,m} \}       \)   are the weights for SDC point \{t_j\} for
  !>                                        subinterval \( (t_{m-1}, \tilde{t}_i) \)
  !>
  !> The quadrature weights are stored in
  !>
  !>      w_rk(0:n_sub, 1:n_stage, 1:n_sub)
  !>

  type, extends(SDC_Method) :: CD_RK_SDC_Method_1D

    integer :: po = -1 !< polynomial order
    integer :: ne = -1 !< number of elements

    type(IMEX_RK_Method)   :: imex_rk     !< underlying IMEX Runge-Kutta method
    real(RNP), allocatable :: w_rk(:,:,:) !< quadrature weights at RK points
    real(RNP), allocatable :: l_rk(:,:,:) !< interpolation weights at RK points

  contains

    procedure :: Init_CD_RK_SDC_Method_1D
    procedure :: TimeStep

  end type CD_RK_SDC_Method_1D

  ! constructor
  interface CD_RK_SDC_Method_1D
    module procedure New_CD_RK_SDC_Method_1D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor

  function New_CD_RK_SDC_Method_1D(po, ne, sdc_opt, n_stage, method) result(this)
    integer,           intent(in) :: po      !< polynomial order
    integer,           intent(in) :: ne      !< number of elements
    type(SDC_Options), intent(in) :: sdc_opt !< SDC options
    integer,           intent(in) :: n_stage !< number of RK stages
    integer, optional, intent(in) :: method  !< RK scheme [1]

    type(CD_RK_SDC_Method_1D) :: this

    call Init_CD_RK_SDC_Method_1D(this, po, ne, sdc_opt, n_stage, method)

  end function New_CD_RK_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Init IMEX RK SDC method for 1D CG-SE convection diffusion solver

  subroutine Init_CD_RK_SDC_Method_1D(this, po, ne, sdc_opt, n_stage, method)
    class(CD_RK_SDC_Method_1D), intent(inout) :: this
    integer,              intent(in)  :: po      !< polynomial order
    integer,              intent(in)  :: ne      !< number of elements
    type(SDC_Options),    intent(in)  :: sdc_opt !< SDC options
    integer,              intent(in)  :: n_stage !< number of RK stages
    integer,    optional, intent(in)  :: method  !< RK method [1]

    ! local variables ..........................................................

    real(RNP) :: t0, t1, ti
    integer   :: i, j, m

    ! mesh dimensions ..........................................................

    this % po = po
    this % ne = ne

    ! initialize SDC method, i.e. the parent type ..............................

    call this % Init_SDC_Method(sdc_opt)

    ! initialize the IMEX RK method ............................................

    call this % imex_rk % Init_IMEX_RK_Method(n_stage, method)

    associate(n_sub => this % n_sub)

      ! RK components ..........................................................

      if (allocated(this % w_rk)) deallocate(this % w_rk)
      if (allocated(this % l_rk)) deallocate(this % l_rk)

      allocate(this % w_rk(0:n_sub, 1:n_stage, 1:n_sub), source = ZERO)
      allocate(this % l_rk(0:n_sub, 1:n_stage, 1:n_sub), source = ZERO)

      ! weights for RK subinterval integrals ...................................

      associate(imex_rk => this % imex_rk)

        subintervals: do m = 1, n_sub

          ! starting and end point of the SDC interval
          t0 = this % t(m-1)
          t1 = this % t(m)

          rk_points: do i = 1, n_stage

            ! position the current RK point
            ti = t0 + (t1 - t0) * imex_rk % c(i)

            ! quadrature weights
            !!! These are the weights for numerical integration over interval
            !!! (t0,ti) where, t0 is the start and ti is the time of RK stage i
            !!! both in subinterval m
            this % w_rk(:,i,m) = this % SubintervalWeights(t0, ti)

            ! interpolation weights
            !!! 'l_rk(j,i,m)' is the Lagrange polynomial to SDC point 'j'
            !!! evaluated at RK point 'i' in subinterval 'm'
            do j = 0, n_sub
              this % l_rk(j,i,m) = LagrangePolynomial(j, this % t, ti)
            end do

          end do rk_points
        end do subintervals

      end associate

    end associate

  end subroutine Init_CD_RK_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Performs a single RK-SDC time step

  subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
    class(CD_RK_SDC_Method_1D),    intent(inout) :: this
    class(CG_ElementOperators_1D), intent(in)    :: eop      !< element operators
    real(RNP),                     intent(in)    :: dx       !< element length
    real(RNP),                     intent(in)    :: dt       !< time step size
    real(RNP),                     intent(in)    :: M(0:,:)  !< global mass matrix
    class(HarmonicWavePackage),    intent(in)    :: wave     !< exact wave solution
    real(RNP),                     intent(in)    :: v        !< convection velicity
    real(RNP),                     intent(in)    :: nu       !< diffusivity
    character,                     intent(in)    :: bc(:)    !< boundary conditions
    real(RNP),                     intent(in)    :: x(0:,:)  !< mesh points
    real(RNP),                     intent(in)    :: t0       !< time t₀
    real(RNP),                     intent(in)    :: u0(0:,:) !< solution u(t₀)
    real(RNP),                     intent(out)   :: u (0:,:) !< solution u(t₀+∆t)

    ! local variables  .........................................................

    ! variables related to SDC
    real(RNP), allocatable :: t_sdc(:)         ! SDC points {tᵢ} in [t0, t0+dt]
    real(RNP), allocatable :: dt_sdc(:)        ! SDC subinterval lengths {∆tᵢ}
    real(RNP), allocatable :: u_sdc(:,:,:)     ! solution at SDC points
    real(RNP), allocatable :: F_im_sdc(:,:,:)  ! implicit RHS at SDC points
    real(RNP), allocatable :: F_ex_sdc(:,:,:)  ! explicit RHS at SDC points
    real(RNP), allocatable :: S_sdc(:,:,:)     ! SDC subinterval integrals

    ! variables related to the Runge-Kutta propagator
    ! all these quantities refer to RK stages for each SDC subinterval
    real(RNP), allocatable :: F_im_rk(:,:,:,:) ! implicit RHS at RK points
    real(RNP), allocatable :: F_ex_rk(:,:,:,:) ! explicit RHS at RK points
    real(RNP), allocatable :: S_rk(:,:,:,:)    ! Runge-Kutta stage integrals

    integer :: po, ne, n_sub, n_stage
    integer :: i, k

    ! initialization ...........................................................

    po      = ubound(u, 1)             ! polynomial order
    ne      = ubound(u, 2)             ! number of elements
    n_sub   = this % n_sub             ! number of SDC subintervals
    n_stage = this % imex_rk % n_stage ! number of RK stages

    allocate(  t_sdc (0:n_sub) ) ! {tᵢ}
    allocate( dt_sdc (1:n_sub) ) ! {∆tᵢ}

    allocate(  u_sdc    (0:po, 1:ne, 0:n_sub), source = ZERO)
    allocate(  F_im_sdc (0:po, 1:ne, 0:n_sub), source = ZERO)
    allocate(  F_ex_sdc (0:po, 1:ne, 0:n_sub), source = ZERO)
    allocate(  S_sdc    (0:po, 1:ne, 1:n_sub), source = ZERO)

    allocate(  F_im_rk  (0:po, 1:ne, 1:n_stage, 1:n_sub), source = ZERO)
    allocate(  F_ex_rk  (0:po, 1:ne, 1:n_stage, 1:n_sub), source = ZERO)
    allocate(  S_rk     (0:po, 1:ne, 1:n_stage, 1:n_sub), source = ZERO)

    t_sdc = this % IntermediateTimes(t0, dt)
    dt_sdc = t_sdc(1:) - t_sdc(0:n_sub-1)

    ! first SDC point coincides with starting time t0
    u_sdc(:,:,0) = u0

    ! predictor ................................................................

    do i = 1, n_sub

      call Propagator( this % imex_rk                  &
                     , wave, eop, dx, M, v, nu, bc, x  &
                     , t0    = t_sdc(i-1)              &
                     , dt    = dt_sdc(i)               &
                     , u0    = u_sdc(:,:,i-1)          &
                     , u     = u_sdc(:,:,i)            &
                     , F_im  = F_im_rk(:,:,:,i)        &
                     , F_ex  = F_ex_rk(:,:,:,i)        )

    end do

    ! corrector sweeps .........................................................

     do k = 1, this % n_sweep

      ! compute "old" RHS at SDC points
      !
      ! remarks
      ! 1) recomputes first point though unchanged
      ! 2) can be skipped, when reusing F_im_rk, F_ex_rk from previous sweep
      call GetRHS(wave, eop, dx, v, nu, bc, x, t_sdc, u_sdc, F_im_sdc, F_ex_sdc)

      ! compute SDC subinterval integrals
      call Get_SDC_Subintegrals(dt, this%w_sub, F_im_sdc, F_ex_sdc, S_sdc)

      ! for all subintervals: compute RK stage integrals
      call Get_RK_Subintegrals(dt, this%w_rk, F_im_sdc, F_ex_sdc, S_rk)

      ! for all subintervals: provide RHS at RK points
      ! options
      !   1. interpolate F_im, F_ex
      !   2. interpolate u and then evaluate F_im, F_ex (later?)
      call Interpolate_RHS(this%l_rk, F_im_sdc, F_ex_sdc, F_im_rk, F_ex_rk)

      ! correction
      do i = 1, n_sub

        call Propagator( this % imex_rk                  &
                       , wave, eop, dx, M, v, nu, bc, x  &
                       , t0    = t_sdc(i-1)              &
                       , dt    = dt_sdc(i)               &
                       , u0    = u_sdc(:,:,i-1)          &
                       , u     = u_sdc(:,:,i)            &
                       , F_im  = F_im_rk(:,:,:,i)        &
                       , F_ex  = F_ex_rk(:,:,:,i)        &
                       , S     = S_rk(:,:,:,i)           &
                       , S_sdc = S_sdc(:,:,i)            )

       end do

     end do

    ! finalization .............................................................

    u = u_sdc(:,:,n_sub)

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> Compose RHS contributions

  subroutine GetRHS(wave, eop, dx, v, nu, bc, x, t, u, F_im, F_ex)

    class(HarmonicWavePackage),    intent(in) :: wave !< exact wave solution
    class(CG_ElementOperators_1D), intent(in) :: eop  !< element operators

    real(RNP), intent(in)  :: dx           !< element length
    real(RNP), intent(in)  :: v            !< convection velocity
    real(RNP), intent(in)  :: nu           !< diffusivity
    character, intent(in)  :: bc(:)        !< boundary conditions
    real(RNP), intent(in)  :: x(:,:)       !< mesh points
    real(RNP), intent(in)  :: t(:)         !< time points
    real(RNP), intent(in)  :: u(:,:,:)     !< solution at different times t(:)
    real(RNP), intent(out) :: F_im(:,:,:)  !< implicit RHS at times t(:)
    real(RNP), intent(out) :: F_ex(:,:,:)  !< explicit RHS at times t(:)

    integer :: i

    do i = 1, size(t)
      call CD_GetDiffusionTerm_1D(eop, dx, wave, v, nu, bc, x,  &
                            t(i),u(:,:,i), F_im(:,:,i)    )
      call CD_GetLinearConvectionTerm_1D(eop, v, bc, u(:,:,i), F_ex(:,:,i))
    end do

  end subroutine GetRHS

  !-----------------------------------------------------------------------------
  !> Computes `S(:,:,i) = dt * sum_j [(F_im(:,:,j) + F_ex(:,:,j)) * w_sub(j,i)]`

  subroutine Get_SDC_Subintegrals(dt, w_sub, F_im, F_ex, S)
    real(RNP), intent(in)  :: dt           !< time step width
    real(RNP), intent(in)  :: w_sub(0:,:)  !< weights
    real(RNP), intent(in)  :: F_im(:,:,0:) !< implicit RHS at SDC points
    real(RNP), intent(in)  :: F_ex(:,:,0:) !< explicit RHS at SDC points
    real(RNP), intent(out) :: S(:,:,:)     !< subinterval integrals

    integer :: n_sub
    integer :: j, m

    n_sub = size(w_sub, 2)

    SDC_Subintervals: do m = 1, n_sub

      S(:,:,m) = 0
      do j = 0, n_sub
        S(:,:,m) = S(:,:,m) + dt * (F_im(:,:,j) + F_ex(:,:,j)) * w_sub(j,m)
      end do

    end do SDC_Subintervals

  end subroutine Get_SDC_Subintegrals

  !-----------------------------------------------------------------------------
  !> Compute RK stage integrals for each SDC subinterval

  subroutine Get_RK_Subintegrals(dt, w_rk, F_im, F_ex, S_rk)
    real(RNP), intent(in)  :: dt            !< time step width
    real(RNP), intent(in)  :: w_rk(0:,:,:)  !< quadrature weights
    real(RNP), intent(in)  :: F_im(:,:,0:)  !< implicit RHS at SDC points
    real(RNP), intent(in)  :: F_ex(:,:,0:)  !< explicit RHS at SDC points
    real(RNP), intent(out) :: S_rk(:,:,:,:) !< RK stage integrals in SDC intervals

    integer :: n_stage, n_sub
    integer :: i, j, m

    n_stage = size(w_rk, 2)
    n_sub   = size(w_rk, 3)

    SDC_Subintervals: do m = 1, n_sub

      RK_Stage_Integrals: do i = 1, n_stage

        S_rk(:,:,i,m) = 0
        do j = 0, n_sub
          S_rk(:,:,i,m) = S_rk(:,:,i,m) &
                        + dt * (F_im(:,:,j) + F_ex(:,:,j)) * w_rk(j,i,m)
        end do

      end do RK_Stage_Integrals

    end do SDC_Subintervals

  end subroutine Get_RK_Subintegrals

  !-----------------------------------------------------------------------------
  !> Interpolates F_im/ex given at SDC points to RK points in all subintervals
  !> (instead of first interpolating u from SDC to RK points and then eval F(u))

  subroutine Interpolate_RHS(l_rk, F_im_sdc, F_ex_sdc, F_im_rk, F_ex_rk)
    real(RNP), intent(in)  :: l_rk(0:,:,:)     !< interpolation weights
    real(RNP), intent(in)  :: F_im_sdc(:,:,0:) !< implicit RHS at SDC points
    real(RNP), intent(in)  :: F_ex_sdc(:,:,0:) !< explicit RHS at SDC points
    real(RNP), intent(out) :: F_im_rk(:,:,:,:) !< implicit RHS at RK points
    real(RNP), intent(out) :: F_ex_rk(:,:,:,:) !< explicit RHS at RK points

    integer :: n_stage, n_sub
    integer :: i, j, m

    n_stage = size(l_rk, 2)
    n_sub   = size(l_rk, 3)

    SDC_Subintervals: do m = 1, n_sub

      RK_Stages: do i = 1, n_stage

        F_im_rk(:,:,i,m) = 0
        F_ex_rk(:,:,i,m) = 0
        do j = 0, n_sub
          !!! note that 'l_rk(j,i,m)' is the Lagrange polynomial to SDC point 'j'
          !!! evaluated at RK point 'i' in subinterval 'm'
          !!! i.e. just normal Lagrange interpolation
          F_im_rk(:,:,i,m) = F_im_rk(:,:,i,m) + F_im_sdc(:,:,j) * l_rk(j,i,m)
          F_ex_rk(:,:,i,m) = F_ex_rk(:,:,i,m) + F_ex_sdc(:,:,j) * l_rk(j,i,m)
        end do

      end do RK_Stages

    end do SDC_Subintervals

  end subroutine Interpolate_RHS

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX RK time step

  subroutine Propagator( imex_rk, wave, eop, dx, M, v, nu, bc, x, t0, dt, &
                         u0, u, F_im, F_ex, S, S_sdc                      )

    ! arguments ................................................................

    class(IMEX_RK_Method),         intent(inout) :: imex_rk !< IMEX RK method
    class(HarmonicWavePackage),    intent(in)    :: wave    !< exact wave solution
    class(CG_ElementOperators_1D), intent(in)    :: eop     !< element operators

    real(RNP), intent(in)  :: dx      !< element length
    real(RNP), intent(in)  :: M(0:,:) !< global mass matrix
    real(RNP), intent(in)  :: v       !< convection velocity
    real(RNP), intent(in)  :: nu      !< diffusivity
    character, intent(in)  :: bc(:)   !< boundary conditions
    real(RNP), intent(in)  :: x(0:,:) !< mesh points
    real(RNP), intent(in)  :: t0      !< start  of current SDC subinterval
    real(RNP), intent(in)  :: dt      !< length of current SDC subinterval

    real(RNP), intent(in)  :: u0(0:,:) !< solution u(t₀)
    real(RNP), intent(out) :: u (0:,:) !< solution u(t₀+∆t)

    !!! These are the implicit/explicit F from previous SDC iteration
    !!! interpolated to RK stage times, e.g., F_im(:,:,2) = F_im(t2, Iu^k(t2))
    real(RNP), intent(inout) :: F_im(0:,:,:) !< F_im(0:po, 1:ne, 1:n_stage)
    real(RNP), intent(inout) :: F_ex(0:,:,:) !< F_ex(0:po, 1:ne, 1:n_stage)

    ! subinterval integrals, needed only when operating as a corrector
    real(RNP), optional, intent(in) :: S(0:,:,:)   !< S for RK stages
    real(RNP), optional, intent(in) :: S_sdc(0:,:) !< S for SDC subinterval

    ! local variables ..........................................................

    real(RNP), allocatable :: G_im(:,:,:)  !
    real(RNP), allocatable :: G_ex(:,:,:)  ! ...
    real(RNP), allocatable :: w(:,:)       ! workspace

    real(RNP) :: c, t
    integer   :: i, j

    ! initialization ...........................................................

    allocate( w   , mold = u    )
    allocate( G_im, mold = F_im )
    allocate( G_ex, mold = F_ex )

    G_im = -F_im   !!! here, e.g: G_im(:,:,2) = -F_im(t2, Iu^k(t2))
    G_ex = -F_ex

    associate( ns   => imex_rk % n_stage                            &
             , b_im => imex_rk % b_im    , b_ex => imex_rk % b_ex   &
             , a_im => imex_rk % a_im    , a_ex => imex_rk % a_ex   )


      ! stage 1 ................................................................

      call CD_GetDiffusionTerm_1D(eop, dx, wave, v, nu, bc, x, t0, u0, F_im(:,:,1))
      call CD_GetLinearConvectionTerm_1D(eop, v, bc, u0, F_ex(:,:,1))

      !!! now, e.g: G_im(:,:,1) = F_im(t1,u^{k+1}(t1)) - F_im(t1, Iu^k(t1))
      G_im(:,:,1) = G_im(:,:,1) + F_im(:,:,1)
      G_ex(:,:,1) = G_ex(:,:,1) + F_ex(:,:,1)

      ! stages 2 to ns .........................................................

      do i = 2, ns

        t = t0 + imex_rk%c(i) * dt
        c = 1 / (dt * a_im(i,i))

        w = M * u0
        do j = 1, i-1
          w = w + dt * (a_im(i,j) * G_im(:,:,j) + a_ex(i,j) * G_ex(:,:,j))
        end do
        w = w + dt * a_im(i,i) * G_im(:,:,i)

        if (present(S)) then
          w = w + S(:,:,i)
        end if

        if (nu > 0) then
          w = c * w
          call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t, u, w)
          call CG_CondensedEllipticSolver_1D(eop, dx, c, nu, bc, w, u)
        else
          u = w / M
          call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t, u)
        end if

        call CD_GetDiffusionTerm_1D(eop, dx, wave, v, nu, bc, x, t, u, F_im(:,:,i))
        call CD_GetLinearConvectionTerm_1D(eop, v, bc, u, F_ex(:,:,i))

        !!! now, e.g: G_im(:,:,i) = F_im(ti,u^{k+1}(ti)) - F_im(ti, Iu^k(ti))
        G_im(:,:,i) = G_im(:,:,i) + F_im(:,:,i)
        G_ex(:,:,i) = G_ex(:,:,i) + F_ex(:,:,i)

      end do

      ! result .................................................................

      w = 0
      do i = 1, ns
        w = w + dt * (b_im(i) * G_im(:,:,i) + b_ex(i) * G_ex(:,:,i))
      end do

      if (present(S_sdc)) then
        !!! now add the integral S
        w = w + S_sdc
      end if

      u = u0 + w / M
      call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t0 + dt, u)

    end associate

  end subroutine Propagator

  !=============================================================================

end module CD__RK_SDC_Method__1D
