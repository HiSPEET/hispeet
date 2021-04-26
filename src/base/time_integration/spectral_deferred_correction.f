!> summary:  Spectral deferred correction base type
!> author:   Joerg Stiller, Martina Grotteschi
!> date:     2019/03/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Deferred_Correction

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE, HALF
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use XMPI

  implicit none
  private

  public :: SDC_Options
  public :: SDC_Method

  !-----------------------------------------------------------------------------
  !> Type bundling spectral deferred correction options

  type SDC_Options
    integer :: n_sub     = 1  !< number of subintervals
    integer :: n_sweep   = 0  !< max number of correction sweeps
    integer :: point_set = 2  !< equidistant (1) or GLL (2) points
  contains
    procedure :: Bcast => Bcast_SDC_Options
  end type SDC_Options

  !-----------------------------------------------------------------------------
  !> Spectral deferred correction parameters and procedures
  !>
  !> This type defines a subdivision of the reference interval `[0,1]` into M
  !> subintervals `[τᵢ₋₁,τᵢ]`. Two choices exist for the  point set `{τᵢ}`:
  !>
  !>   1. the equidistant partition of `[0,1]`, or
  !>   2. the Gauss-Legendre-Lobatto (GLL) points mapped to `[0,1]`.
  !>
  !> Within the type, `t(i) = τᵢ` represents the i-th point, `w(i) = wᵢ` the
  !> corresponding quadrature weight and `n_sub = M` the number of subintervals.
  !> The integral of a function f over the reference interval is approximated by
  !>
  !>   \[ \int_{0}^{1} f d\tau \approx \sum_{i=0}^{M} w_i\, f(\tau_i) \]
  !>
  !> The quadrature will be exact for polynomials of degree `M` with equidistant
  !> points and degree `2M-1` with GLL points.
  !>
  !> Similarly, integrals over the subintervals `[τᵢ₋₁,τᵢ]` can be evaluated by
  !>
  !>   \[
  !>      \int_{\tau_{i-1}}^{\tau_i} f d\tau
  !>      \approx
  !>      \sum_{j=0}^{M} w^s_{j,i}\, f(\tau_i)
  !>   \]
  !>
  !> where the weights \(w^s_{j,i}\), denoted `w_sub(j,i)` in Fortran, are
  !> obtained by application of the GLL quadrature with `M+1` points to the
  !> Lagrange interpolant constructed from `f(τᵢ)`.

  type SDC_Method

    integer :: n_sub     = -1  !< number of subintervals (M)
    integer :: n_sweep   = -1  !< max num correction sweeps (K)
    integer :: point_set = -1  !< equidistant (1) or GLL (2) points

    real(RNP), allocatable :: t(:)       !< nodes τᵢ in [0,1]
    real(RNP), allocatable :: w(:)       !< quadrature weights for [0, 1]
    real(RNP), allocatable :: w_sub(:,:) !< quadrature weights for [τᵢ₋₁,τᵢ]

    real(RNP), allocatable, private :: x_gll(:) !< GLL nodes in [-1,1]
    real(RNP), allocatable, private :: w_gll(:) !< GLL weights to x

  contains

    procedure :: Init_SDC_Method  =>  Init_SDC
    procedure :: Show             =>  Show_SDC_Method
    procedure :: HasEquidistantPoints
    procedure :: HasLobattoPoints
    procedure :: IntermediateTimes
    procedure :: SubintervalWeights

  end type SDC_Method

  ! constructor
  interface SDC_Method
    module procedure New_SDC
  end interface

contains

  !=============================================================================
  ! SDC_Options: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Broadcasting the SDC options

  subroutine Bcast_SDC_Options(this, root, comm)
    class(SDC_Options), intent(inout) :: this
    integer,            intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),     intent(in)    :: comm !< MPI communicator

    type(MPI_Request)  :: request(3)
    integer :: n

    n = 1
    call XMPI_Ibcast( this % n_sub    , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % n_sweep  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % point_set, root, comm, request(n) )

    call MPI_Waitall(n, request, MPI_STATUSES_IGNORE)

  end subroutine Bcast_SDC_Options

  !=============================================================================
  ! SDC_Method: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor

  type(SDC_Method) function New_SDC(opt) result(this)
    type(SDC_Options), intent(in) :: opt  !< SDC options

    call Init_SDC(this, opt)

  end function New_SDC

  !-----------------------------------------------------------------------------
  !> Initialization of spectral deferred correction

  subroutine Init_SDC(this, opt)
    class(SDC_Method),  intent(inout) :: this
    class(SDC_Options), intent(in)    :: opt  !< SDC options

    ! local variables ..........................................................

    integer :: i, n_sub

    ! initialization ...........................................................

    n_sub = opt % n_sub

    if (this % n_sub > 0 .and. this % n_sub /= n_sub) then
      deallocate(this % t    )
      deallocate(this % w    )
      deallocate(this % w_sub)
    end if
    if (.not. allocated(this % t ))    allocate(this % t     (0:n_sub)      )
    if (.not. allocated(this % w ))    allocate(this % w     (0:n_sub)      )
    if (.not. allocated(this % w_sub)) allocate(this % w_sub (0:n_sub, n_sub) )

    this % n_sub     = n_sub
    this % n_sweep   = max(0, opt % n_sweep)
    this % point_set = max(1, min(2, opt % point_set))

    ! GLL points and weights in [-1,1]
    allocate(this % x_gll(0:n_sub), source = GLL_Points(n_sub))
    allocate(this % w_gll(0:n_sub), source = GLL_Weights(this % x_gll))

    ! points and quadrature weights in [0,1] ...................................

    select case(this % point_set)
    case(1) ! equidistant
      this % t(0:n_sub) = [ ZERO, (i*ONE/n_sub, i = 1,n_sub-1), ONE ]
      this % w(0:n_sub) = GaussLagrangeWeights(this% t)
    case default ! GLL
      this % t(0:n_sub) = HALF * (this % x_gll + ONE)
      this % w(0:n_sub) = HALF *  this % w_gll
    end select

    ! quadrature weights in [τᵢ₋₁,τᵢ] ..........................................

    do i = 1, n_sub
      this % w_sub(:,i) = this % SubintervalWeights(this%t(i-1), this%t(i))
    end do

  end subroutine Init_SDC

  !-----------------------------------------------------------------------------
  !> Output of SDC settings

  subroutine Show_SDC_Method(this, unit)
    class(SDC_Method), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'SDC_Method settings'
    write(io,'(A,/)') repeat('≡',80)
    write(io,'(2X,A,T15,I0)') 'n_sub'      , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweeps:'  , this % n_sweep
    write(io,'(2X,A,T15,I0)') 'point_set:' , this % point_set

  end subroutine Show_SDC_Method

  !-----------------------------------------------------------------------------
  !> Query if point set is equidistant

  pure logical function HasEquidistantPoints(this)
    class(SDC_Method), intent(in) :: this

    HasEquidistantPoints = this % point_set == 1

  end function HasEquidistantPoints

  !-----------------------------------------------------------------------------
  !> Query if point set is based on Lobatto (GLL) points

  pure logical function HasLobattoPoints(this)
    class(SDC_Method), intent(in) :: this

    HasLobattoPoints = this % point_set == 2

  end function HasLobattoPoints

  !-----------------------------------------------------------------------------
  !> Returns the intermediate times within a given time interval

  pure function IntermediateTimes(this, t0, dt) result(t)
    class(SDC_Method), intent(in) :: this
    real(RNP), intent(in)  :: t0              !< start of the time interval
    real(RNP), intent(in)  :: dt              !< length of the time interval
    real(RNP)              :: t(0:this%n_sub) !< intermediate times

    t = t0 + dt * this % t

  end function IntermediateTimes

  !-----------------------------------------------------------------------------
  !> Returns the quadrature weigths for an arbitrary subinterval of [0,1]
  !>
  !>
  !> Given the intervall $$[\tau_a, \tau_b]$$ the weights $$w^s$$ are computed
  !> such that
  !>
  !>   \[
  !>      \int_{\tau_a}^{\tau_b} f d\tau
  !>      \approx
  !>      \sum_{j=0}^{M} w^s_{j}\, f(\tau_j)
  !>   \]
  !>
  !> where $$\tau_j$$ are the SDC points.

  function SubintervalWeights(this, ta, tb) result(ws)
    class(SDC_Method), intent(in) :: this
    real(RNP), intent(in) :: ta               !< start of the subinterval
    real(RNP), intent(in) :: tb               !< end of the subinterval
    real(RNP)             :: ws(0:this%n_sub) !< weights

    real(RNP) :: delta, tk, yk
    integer   :: j, k

    integer :: n_quad

    associate(n_sub => this%n_sub, t => this%t, x => this%x_gll, w => this%w_gll)

      n_quad = ubound(x,1)
      delta  = (tb - ta) * HALF
      do j = 0, n_sub
        ws(j) = 0
        do k = 0, n_quad
          ! tk = τ(x(k)) = k-th quadrature point mapped to [τa,τb]
          tk  = ta + delta * (x(k) + 1)
          ! yk = value of j-th Lagrange polynomial to SDC points t(:) at tk
          yk = LagrangePolynomial(j, t, tk)
          ! add contribution of k-th GLL point
          ws(j) = ws(j) + w(k) * yk
        end do
        ! scale the weight to match the length of interval [τa,τb]
        ws(j) = delta * ws(j)
      end do

    end associate

  end function SubintervalWeights

  !=============================================================================

end module Spectral_Deferred_Correction
