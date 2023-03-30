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

module CD__IMEX_RK_Method__1D
  use Kind_Parameters, only: RNP
  use CD__Utils__1D
  use CG__Element_Operators__1D
  use CG__Condensed_Solver__1D
  use IMEX_Runge_Kutta_Method
  use Harmonic_Wave_Package
  implicit none
  private

  public :: CD_IMEX_RK_Method_1D

  !-----------------------------------------------------------------------------
  !> Implementation of the IMEX RK method for 1D convection-diffusion

  type, extends(IMEX_RK_Method) :: CD_IMEX_RK_Method_1D
    real(RNP), allocatable :: F_ex(:,:,:) !< explicit RHS per stage
    real(RNP), allocatable :: F_im(:,:,:) !< implicit RHS per stage
  contains
    procedure :: Init_CD_IMEX_RK_Method_1D
    procedure :: TimeStep
  end type CD_IMEX_RK_Method_1D

  ! constructor
  interface CD_IMEX_RK_Method_1D
    module procedure New_CD_IMEX_RK_Method_1D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor

  type(CD_IMEX_RK_Method_1D) function New_CD_IMEX_RK_Method_1D &
                                               (po, ne, ns, method) result(this)

    integer,           intent(in) :: po     !< polynomial order
    integer,           intent(in) :: ne     !< number of elements
    integer,           intent(in) :: ns     !< number of stages
    integer, optional, intent(in) :: method !< RK scheme [1]

    call Init_CD_IMEX_RK_Method_1D(this, po, ne, ns, method)

  end function New_CD_IMEX_RK_Method_1D

  !-----------------------------------------------------------------------------
  !> Init IMEX RK method for 1D CG-SE convection diffusion solver

  subroutine Init_CD_IMEX_RK_Method_1D(this, po, ne, ns, method)
    class(CD_IMEX_RK_Method_1D), intent(inout) :: this
    integer,                 intent(in)    :: po     !< polynomial order
    integer,                 intent(in)    :: ne     !< number of elements
    integer,                 intent(in)    :: ns     !< number of stages
    integer,       optional, intent(in)    :: method !< RK method [1]

    call this % Init_IMEX_RK_Method(ns, method)

    if (allocated(this % F_im)) deallocate(this % F_im)
    if (allocated(this % F_ex)) deallocate(this % F_ex)

    allocate(this % F_ex(0:po, ne, ns))
    allocate(this % F_im(0:po, ne, ns))

  end subroutine Init_CD_IMEX_RK_Method_1D

  !-----------------------------------------------------------------------------
  !> Performs a single IMEX RK time step

  subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
    class(CD_IMEX_RK_Method_1D),   intent(inout) :: this
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

    real(RNP), allocatable :: f(:,:)
    real(RNP) :: c, t
    logical   :: first = .true.
    integer   :: i, j

    ! initialization ...........................................................

    allocate(f, mold = u)

    associate( ns   => this % n_stage                         &
             , b_im => this % b_im    , b_ex => this % b_ex   &
             , a_im => this % a_im    , a_ex => this % a_ex   &
             , F_im => this % F_im    , F_ex => this % F_ex  )

      ! stage 1 ................................................................

      if (first .or. this%c(ns) /= 1) then
        call CD_GetDiffusionTerm_1D(eop, dx, wave, v, nu, bc, x, t0, u0, F_im(:,:,1))
        call CD_GetLinearConvectionTerm_1D(eop, v, bc, u0, F_ex(:,:,1))
      else
        ! assume FSAL scheme
        F_im(:,:,1) = F_im(:,:,ns)
        F_ex(:,:,1) = F_ex(:,:,ns)
      end if

      ! stages 2 to ns .........................................................

      do i = 2, ns

        t = t0 + this%c(i) * dt
        c = 1 / (dt * a_im(i,i))

        f = c * M * u0
        do j = 1, i-1
          f = f + a_im(i,j) / a_im(i,i) * F_im(:,:,j)  &
                + a_ex(i,j) / a_im(i,i) * F_ex(:,:,j)
        end do

        if (nu > 0) then
          call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t, u, f)
          call CG_CondensedEllipticSolver_1D(eop, dx, c, nu, bc, f, u)
        else
          u = f / (c * M)
          call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t, u)
        end if

        call CD_GetDiffusionTerm_1D(eop, dx, wave, v, nu, bc, x, t, u, F_im(:,:,i))
        call CD_GetLinearConvectionTerm_1D(eop, v, bc, u, F_ex(:,:,i))
        ! NOTE that the diffusion term could be obtained cheaper from
        ! F_im(:,:,i) = c * M * u - f
        ! HOWEVER need to check, if BC are treated correctly that way

      end do

      ! result .................................................................

      f = 0
      do i = 1, ns
        f = f + b_im(i) * F_im(:,:,i) + b_ex(i) * F_ex(:,:,i)
      end do

      u = u0 + dt * f / M
      call CD_ApplyBoundaryConditions_1D(wave, v, nu, bc, x, t0 + dt, u)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CD__IMEX_RK_Method__1D
