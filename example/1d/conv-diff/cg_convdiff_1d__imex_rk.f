!> summary:  IMEX Runge-Kutta method with CG-SEM for 1D convection-diffusion
!> author:   Joerg Stiller
!> date:     2018/09/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### IMEX Runge-Kutta method with CG-SEM for 1D convection-diffusion
!>
!> This module provides the type `ConvDiff_IMEX_RK` which extends the IMEX
!> Runge-Kutta methods defined in `ConvDiff_IMEX_RK` for advancing the solution
!> of the semi-discrete 1D convection-diffusion equation
!>
!>     ∂u/∂t = -v ∂u/∂v + nu ∂²u/∂u² ≡ C(u) + D(u)
!>
!> in time. Spatial discretization is based on the continuous spectral-element
!> method using nodal base functions along with GLL quadrature.
!>
!> Typical usage:
!>
!>     type(ConvDiff_IMEX_RK) :: imex_rk
!>
!>     call imex_rk % New(s, m, po, ne)
!>     ! s  :  number of stages
!>     ! m  :  method (optional)
!>     ! po :  polynomial order and
!>     ! ne :  number of elements
!>
!>     call imex_rk % TimeStep(eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
!>     ! for description of arguments see below
!>
!===============================================================================

module CG_ConvDiff_1D__IMEX_RK
  use Kind_Parameters, only: RNP
  use CG_ConvDiff_1D__Utils
  use CG_Element_Operators_1D
  use CG_Condensed_Solver_1D
  use IMEX_Runge_Kutta_Method
  use Harmonic_Wave_Package
  implicit none
  private

  public :: ConvDiff_IMEX_RK

  !-----------------------------------------------------------------------------
  !> Implementation of the IMEX RK method for 1D convection-diffusion

  type, extends(IMEX_RK_Method) :: ConvDiff_IMEX_RK
    real(RNP), allocatable :: f_ex(:,:,:) !< explicit RHS per stage
    real(RNP), allocatable :: f_im(:,:,:) !< implicit RHS per stage
  contains
    generic,   public  :: New => New_ConvDiff_IMEX_RK
    procedure, private :: New_ConvDiff_IMEX_RK
    procedure, public  :: TimeStep
    final              :: Delete_ConvDiff_IMEX_RK
  end type ConvDiff_IMEX_RK

contains

!-------------------------------------------------------------------------------
!> New IMEX RK method for 1D CG-SE convection diffusion solver

subroutine New_ConvDiff_IMEX_RK(this, s, m, po, ne)
  class(ConvDiff_IMEX_RK), intent(inout) :: this
  integer,                 intent(in)    :: s   !< number of stages
  integer,       optional, intent(in)    :: m   !< RK scheme [1]
  integer,                 intent(in)    :: po  !< polynomial order
  integer,                 intent(in)    :: ne  !< number of elements

  ! initialize IMEX_RK_Method
  call this % New(s, m)

  allocate(this % f_ex(0:po, ne, s))
  allocate(this % f_im, mold = this%f_ex)

end subroutine New_ConvDiff_IMEX_RK

!-------------------------------------------------------------------------------
!> Performs a single IMEX RK time step

subroutine TimeStep(this, eop, dx, dt, M, wave, v, nu, bc, x, t0, u0, u)
  class(ConvDiff_IMEX_RK),      intent(inout) :: this
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
  logical   :: first = .true.
  integer   :: i, j

  ! initialization .............................................................

  allocate(f, mold = u)

  associate( s    => this % s    , b    => this % b     &
           , a_im => this % a_im , a_ex => this % a_ex  &
           , f_im => this % f_im , f_ex => this % f_ex  )

    ! stage 1 ..................................................................

    if (first) then
      call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t0, u0, f_im(:,:,1))
      call GetLinearConvectionTerm(eop, v, bc, u0, f_ex(:,:,1))
    else
      ! assume FSAL scheme
      f_im(:,:,1) = f_im(:,:,s)
      f_ex(:,:,1) = f_ex(:,:,s)
    end if

    ! stages 2 to s ............................................................

    do i = 2, s

      t = t0 + this%c(i) * dt
      c = 1 / (dt * a_im(i,i))

      f = c * M * u0
      do j = 1, i-1
        f = f + a_im(i,j) / a_im(i,i) * f_im(:,:,j)  &
              + a_ex(i,j) / a_im(i,i) * f_ex(:,:,j)
      end do

      if (nu > 0) then
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u, f)
        call CondensedEllipticSolver(eop, dx, c, nu, bc, f, u, standby = .true.)
      else
        u = f / (c * M)
        call ApplyBoundaryConditions(wave, v, nu, bc, x, t, u)
      end if

      call GetDiffusionTerm(eop, dx, wave, v, nu, bc, x, t, u, f_im(:,:,i))
      call GetLinearConvectionTerm(eop, v, bc, u, f_ex(:,:,i))
      ! NOTE that the diffusion term could be obtained cheaper from
      ! f_im(:,:,i) = c * M * u - f
      ! HOWEVER need to check, if BC are treated correctly that way

    end do

    ! result ...................................................................

    f = 0
    do i = 1, s
      f = f + b(i) * (f_im(:,:,i) + f_ex(:,:,i))
    end do

    u = u0 + dt * f / M
    call ApplyBoundaryConditions(wave, v, nu, bc, x, t0 + dt, u)

  end associate

end subroutine TimeStep

!-------------------------------------------------------------------------------
!> Deletes the given `ConvDiff_IMEX_RK` object.

subroutine Delete_ConvDiff_IMEX_RK(this)
  type(ConvDiff_IMEX_RK), intent(inout) :: this

  if (allocated(this % f_im)) deallocate(this % f_im)
  if (allocated(this % f_ex)) deallocate(this % f_ex)

end subroutine Delete_ConvDiff_IMEX_RK

!===============================================================================

end module CG_ConvDiff_1D__IMEX_RK
