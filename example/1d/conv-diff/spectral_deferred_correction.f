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

    integer :: n_sub   = -1  !< number of subintervals (M)
    integer :: n_sweep = -1  !< max num correction sweeps (K)
    integer :: set     = -1  !< equidistant (1) or GLL (2) points

    real(RNP), allocatable :: t(:)    !< nodes tᵢ in [0,1]
    real(RNP), allocatable :: w(:)    !< quadrature weights for [0, 1]
    real(RNP), allocatable :: ws(:,:) !< quadrature weights for [tᵢ₋₁,tᵢ]

  contains

    procedure :: Init_SDC_Method  =>  Init_SDC
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
    integer :: n_sub   = 1  !< number of subintervals
    integer :: n_sweep = 0  !< max number of correction sweeps
    integer :: set     = 2  !< equidistant (1) or GLL (2) points
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
  real(RNP) :: delta, tk, yk, ws_ji
  integer   :: i, j, k, n_sub

  ! initialization .............................................................

  n_sub = opt % n_sub

  if (this % n_sub > 0 .and. this % n_sub /= n_sub) then
    deallocate(this % t )
    deallocate(this % w )
    deallocate(this % ws)
  end if
  if (.not. allocated(this % t )) allocate(this % t  (0:n_sub)        )
  if (.not. allocated(this % w )) allocate(this % w  (0:n_sub)        )
  if (.not. allocated(this % ws)) allocate(this % ws (0:n_sub, n_sub) )

  this % n_sub   = n_sub
  this % n_sweep = max(0, opt % n_sweep)
  this % set     = max(1, min(2, opt % set))

  ! GLL points and weights in [-1,1]
  allocate(x(0:n_sub), source = GLL_Points(n_sub))
  allocate(w(0:n_sub), source = GLL_Weights(x))

  ! points and quadrature weights in [0,1] .....................................

  select case(this % set)
  case(1) ! equidistant
    this % t(0:n_sub) = [ ZERO, (i*ONE/n_sub, i = 1,n_sub-1), ONE ]
    this % w(0:n_sub) = GaussLagrangeWeights(this% t)
  case default ! GLL
    this % t(0:n_sub) = HALF * (x + ONE)
    this % w(0:n_sub) = HALF * w
  end select

  ! quadrature weights in [tᵢ₋₁,tᵢ] ............................................

  ! The weights are obtained as follows
  ! 1) map the GLL points x(0:n_sub) to [tᵢ₋₁,tᵢ]
  ! 2) evaluate the Lagrange polynomials defined in [0,1] in these points
  ! 3) obtain the weight of the j-th Langrange polynomial by
  !    3.1) computing the sum of GLL weights scaled with the values of
  !         the polynomial determined in step 2
  !    3.2) scale the result to match the interval [tᵢ₋₁,tᵢ]

  associate(t => this % t )

    do i = 1, n_sub
      delta = (t(i) - t(i-1)) * HALF
      do j = 0, n_sub
        ws_ji = 0
        do k = 0, n_sub
          ! tk = k-th GLL point in [tᵢ₋₁,tᵢ]
          tk  = t(i-1) + delta * (x(k) + 1)
          ! yk = value of j-th Lagrange polynomial at tk
          select case(this % set)
          case(1) ! use equidistant Lagrange polynomial in [0,1]
            yk = LagrangePolynomial(j, t, tk)
          case default ! use GLL Lagrange polynomial in [-1,1]
            yk = GLL_Polynomial(j, x, 2*tk-1)
          end select
          ! add contribution of k-th GLL point
          ws_ji = ws_ji + w(k) * yk
        end do
        ! scale the weight to match the length of interval [tᵢ₋₁,tᵢ]
        this % ws(j,i) = delta * ws_ji
      end do
    end do

  end associate

end subroutine Init_SDC

!-------------------------------------------------------------------------------
!> Query if point set is equidistant

pure logical function HasEquidistantPoints(this)
  class(SDC_Method), intent(in) :: this

  HasEquidistantPoints = this % set == 1

end function HasEquidistantPoints

!-------------------------------------------------------------------------------
!> Query if point set is based on Lobatto (GLL) points

pure logical function HasLobattoPoints(this)
  class(SDC_Method), intent(in) :: this

  HasLobattoPoints = this % set == 2

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
