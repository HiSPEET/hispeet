!> summary:  Spectral deferred correction base type
!> author:   Joerg Stiller
!> date:     2019/03/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>###   Spectral deferred correction base type
!===============================================================================

module Spectral_Deferred_Correction

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use XMPI

  implicit none
  private

  public :: SDC_Method
  public :: SDC_Options

  !-----------------------------------------------------------------------------
  !> Spectral deferred correction parameters and procedures

  type SDC_Method

    integer :: n_sub     = -1  !< number of subintervals (M)
    integer :: n_sweep   = -1  !< max num correction sweeps (K)
    integer :: point_set = -1  !< equidistant (1) or GLL (2) points

    real(RNP), allocatable :: t(:)    !< nodes tᵢ in [0,1]
    real(RNP), allocatable :: w(:)    !< quadrature weights for [0, 1]
    real(RNP), allocatable :: wt(:,:) !< quadrature weights for [0, tᵢ]
    real(RNP), allocatable :: ws(:,:) !< quadrature weights for [tᵢ₋₁,tᵢ]

  contains

    procedure :: Init_SpectralDeferredCorrection  =>  Init_SDC
    procedure :: HasEquidistantPoints
    procedure :: HasLobattoPoints
    procedure :: IntermediateTimes

  end type SDC_Method

  ! constructor
  interface SDC_Method
    module procedure New_SDC
  end interface

  !-----------------------------------------------------------------------------
  !> Type bundling spectral deferred correction options

  type SDC_Options
    integer :: n_sub     = 1  !< number of subintervals
    integer :: n_sweep   = 0  !< max number of correction sweeps
    integer :: point_set = 2  !< equidistant (1) or GLL (2) points
  end type SDC_Options

contains

!===============================================================================
! SDC_Method: type-bound procedures

!-------------------------------------------------------------------------------
!> Constructor

type(SDC_Method) function New_SDC(opt) result(this)
  type(SDC_Options), intent(in) :: opt  !< SDC options

  call Init_SDC(this, opt)

end function New_SDC

!-------------------------------------------------------------------------------
!> Initialization of spectral deferred correction

subroutine Init_SDC(this, opt)
  class(SDC_Method),  intent(inout) :: this
  class(SDC_Options), intent(in)    :: opt  !< SDC options

  ! local variables ............................................................

  real(RNP), allocatable :: x(:), w(:)
  real(RNP) :: ts, ys, tmp
  integer   :: i, j, k, ns

  ! initialization .............................................................

  ns = opt % n_sub

  if (this % n_sub > 0 .and. this % n_sub /= ns) then
    deallocate(this % t )
    deallocate(this % w )
    deallocate(this % wt)
    deallocate(this % ws)
  end if
  if (.not. allocated(this % t )) allocate(this % t  (0:ns)     )
  if (.not. allocated(this % w )) allocate(this % w  (0:ns)     )
  if (.not. allocated(this % wt)) allocate(this % wt (0:ns, ns) )
  if (.not. allocated(this % ws)) allocate(this % ws (0:ns, ns) )

  this % n_sub     = ns
  this % n_sweep   = max(0, opt % n_sweep)
  this % point_set = max(1, min(2, opt % point_set))

  ! GLL points and weights in [-1,1]
  allocate(x(0:ns), source = GLL_Points(ns))
  allocate(w(0:ns), source = GLL_Weights(x))

  ! points and quadrature weights in [0,1] .....................................

  select case(this % point_set)
  case(1) ! equidistant
    this % t(0:ns) = [ ZERO, (i*ONE/ns, i = 1,ns-1), ONE ]
    this % w(0:ns) = GaussLagrangeWeights(this% t )
  case default ! GLL
    this % t(0:ns) = HALF * (x + ONE)
    this % w(0:ns) = HALF * w
  end select

  associate(t => this % t )

    ! quadrature weights in [0, tᵢ] ............................................

    do i = 1, ns
    do j = 0, ns
      tmp = 0
      do k = 0, ns
        ! ts = k-th GLL point in [0,tᵢ] mapped to [0,1]
        ts  = t(i) * HALF * (x(k) + 1)
        ! ys = value of j-th Lagrange polynomial at ts
        select case(this % point_set)
        case(1) ! use equidistant Lagrange polynomial in [0,1]
          ys = LagrangePolynomial(j, t, ts)
        case default ! use GLL Lagrange polynomial in [-1,1]
          ys = GLL_Polynomial(j, x, 2*ts-1)
        end select
        tmp = tmp + w(k) * ys
      end do
      this % wt(j,i) = HALF * tmp
    end do
    end do

    ! quadrature weights in [tᵢ₋₁,tᵢ] ..........................................

    do i = 1, ns
    do j = 0, ns
      tmp = 0
      do k = 0, ns
        ! ts = k-th GLL point in [tᵢ₋₁,tᵢ] mapped to [0,1]
        ts  = t(i-1) + (t(i) - t(i-1)) * HALF * (x(k) + 1)
        ! ys = value of j-th Lagrange polynomial at ts
        select case(this % point_set)
        case(1) ! use equidistant Lagrange polynomial in [0,1]
          ys = LagrangePolynomial(j, t, ts)
        case default! use GLL Lagrange polynomial in [-1,1]
          ys = GLL_Polynomial(j, x, 2*ts-1)
        end select
        tmp = tmp + w(k) * ys
      end do
      this % ws(j,i) = HALF * tmp
    end do
    end do

  end associate

end subroutine Init_SDC

!-------------------------------------------------------------------------------
!> Query if point set is equidistant

pure logical function HasEquidistantPoints(this)
  class(SDC_Method), intent(in) :: this

  HasEquidistantPoints = this % point_set == 1

end function HasEquidistantPoints

!-------------------------------------------------------------------------------
!> Query if point set is based on Lobatto (GLL) points

pure logical function HasLobattoPoints(this)
  class(SDC_Method), intent(in) :: this

  HasLobattoPoints = this % point_set == 2

end function HasLobattoPoints

!-------------------------------------------------------------------------------
!> Returns the intermediate times within a given time interval

pure function IntermediateTimes(this, t0, dt) result(t)
  class(SDC_Method), intent(in) :: this
  real(RNP), intent(in)  :: t0              !< start of the time interval
  real(RNP), intent(in)  :: dt              !< length of the time interval
  real(RNP)              :: t(0:this%n_sub) !< intermediate times

  t = t0 + dt * this % t

end function IntermediateTimes

!===============================================================================

end module Spectral_Deferred_Correction
