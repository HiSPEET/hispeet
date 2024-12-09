!> summary:  Spectral deferred correction base type
!> author:   Joerg Stiller, Martina Grotteschi
!> date:     2019/03/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Spectral_Deferred_Correction

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF
  use Execution_Control, only: Error
  use Gauss_Jacobi
  use Lagrange_Interpolation
  use XMPI

  implicit none
  private

  public :: SDC_Method
  public :: SDC_Options

  !-----------------------------------------------------------------------------
  !> Spectral deferred correction parameters and procedures
  !>
  !> This type defines a subdivision of the reference interval `[0,1]` into M
  !> subintervals `[τᵢ₋₁,τᵢ]`. The following choices exist for the point set
  !> `{τᵢ}`:
  !>
  !>   - `'E'`   equidistant partition of `[0,1]`,
  !>   - `'L'`   Lobatto nodes mapped to `[0,1]`,
  !>   - `'RR'`  right-sided Radau nodes mapped to `[0,1]`.
  !>
  !> Within the type, `t(i) = τᵢ` represents the i-th point, `w(i) = wᵢ` the
  !> corresponding quadrature weight and `n_sub = M` the number of subintervals.
  !> For a unified approach, the leftmost point is always set to `τ₀ = 0`.
  !>
  !> The integral of a function f over the reference interval is approximated by
  !>
  !>   \[ \int_{0}^{1} f d\tau \approx \sum_{i=0}^{M} w_i\, f(\tau_i) \]
  !>
  !> The quadrature will be exact for polynomials of degree `M` with equidistant
  !> points, `2M-1` with Lobatto nodes and `2M-2` with right-sided Radau nodes.
  !> Note that the latter are less accurate because they do not include `τ₀`.
  !>
  !> Similarly, integrals over the subintervals `[τᵢ₋₁,τᵢ]` can be evaluated by
  !>
  !>   \[
  !>      \int_{\tau_{i-1}}^{\tau_i} f d\tau
  !>      \approx
  !>      \sum_{j=0}^{M} w^s_{j,i}\, f(\tau_i)
  !>   \]
  !>
  !> where the weights \(w^s_{j,i}\), denoted `w_nn(j,i)` in Fortran, are
  !> obtained by application of the Lobatto quadrature with `M+1` points to
  !> the Lagrange interpolant constructed from `f(τᵢ)`.
  !>
  !> Additionally, `w_0n(:,i)` provides the quadrature weights for the interval
  !> [0,τᵢ]. The transpose of `w_0n` represents the coefficients of the related
  !> collocation method.

  type SDC_Method

    character(2) :: nodes   !< type of collocation points (nodes) {'E','L','RR'}
    integer      :: p_col   !< polynomial degree of collocation points
    integer      :: n_col   !< number of collocation points
    integer      :: n_sub   !< number of subintervals (M)
    integer      :: n_sweep !< max num correction sweeps (K)

    real(RNP), allocatable :: t(:)      !< subinterval points τᵢ in [0,1]
    real(RNP), allocatable :: w(:)      !< quadrature weights for [0, 1]
    real(RNP), allocatable :: w_nn(:,:) !< quadrature weights for [τᵢ₋₁,τᵢ]
    real(RNP), allocatable :: w_0n(:,:) !< quadrature weights for [0,τᵢ]

    real(RNP), allocatable, private :: x_quad(:) !< quadrature nodes   in [-1,1]
    real(RNP), allocatable, private :: w_quad(:) !< quadrature weights in [-1,1]

  contains

    procedure :: Init_SDC_Method  =>  Init_SDC
    procedure :: Show             =>  Show_SDC_Method
    procedure :: CollocationPoints
    procedure :: CollocationWeights
    procedure :: SubintervalPoints
    procedure :: SubintervalWeights

  end type SDC_Method

  ! constructor
  interface SDC_Method
    module procedure New_SDC
  end interface

  !-----------------------------------------------------------------------------
  !> Type bundling spectral deferred correction options

  type SDC_Options
    character(2) :: nodes   = 'RR' !< type of collocation points {'E','L','RR'}
    integer      :: n_col   =  1   !< number of collocation points
    integer      :: n_sweep =  0   !< max number of correction sweeps
  contains
    procedure :: Bcast => Bcast_SDC_Options
  end type SDC_Options

contains

  !=============================================================================
  ! SDC_Options: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Broadcasting the SDC options

  subroutine Bcast_SDC_Options(this, root, comm)
    class(SDC_Options), intent(inout) :: this
    integer,            intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),     intent(in)    :: comm !< MPI communicator

    type(MPI_Request) :: request(3)
    integer :: n

    n = 1
    call XMPI_Ibcast( this % nodes  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % n_col  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % n_sweep, root, comm, request(n) )

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

    integer :: i

    ! initialization ...........................................................

    if (allocated(this % t     ))  deallocate(this % t     )
    if (allocated(this % w     ))  deallocate(this % w     )
    if (allocated(this % w_nn  ))  deallocate(this % w_nn  )
    if (allocated(this % w_0n  ))  deallocate(this % w_0n  )
    if (allocated(this % w_quad))  deallocate(this % w_quad)
    if (allocated(this % w_quad))  deallocate(this % w_quad)

    this % nodes   = opt % nodes
    this % n_col   = opt % n_col
    this % n_sweep = max(0, opt % n_sweep)

    select case(this % nodes)
    case('E','L')
      this % n_sub = this % n_col - 1
    case('RR')
      this % n_sub = this % n_col
    case default
      call Error( 'Init_SDC' &
                , 'node type "' // trim(this%nodes) // '" not supported' &
                , 'Spectral_Deferred_Correction')
    end select

    ! polynomial degree of the interpolation polynomial
    this % p_col = this % n_col - 1

    associate(n_sub => this % n_sub)

      ! auxiliary quadrature points and weights in [-1,1] ......................

      allocate(this % x_quad(0:n_sub), source = GaussPoints(n_sub))
      allocate(this % w_quad(0:n_sub), source = GaussWeights(this % x_quad))

      ! points and quadrature weights in [0,1] .................................

      allocate(this % t(0:n_sub), source = ZERO)
      allocate(this % w(0:n_sub), source = ZERO)

      select case(this % nodes)
      case('E')
        ! equidistant
        this % t(0:n_sub) = [ ZERO, (i*ONE/n_sub, i = 1,this%p_col-1), ONE ]
        this % w(0:n_sub) = GaussLagrangeWeights(this % t)
      case('L')
        ! Lobatto points and weights in [-1,1]
        this % t(0:n_sub) = LobattoPoints(this % p_col)
        this % w(0:n_sub) = LobattoWeights(this % t)
        ! transform to [0,1]
        this % t(0:n_sub) = HALF * (this % t + ONE)
        this % w(0:n_sub) = HALF * this % w
      case('RR') ! Radau right
        ! right-sided Radau points and weights in [-1,1]
        this % t(1:n_sub) = RadauPoints(this % p_col, right = .true.)
        this % w(1:n_sub) = RadauWeights(this % t(1:n_sub))
        ! transform to [0,1]
        this % t(1:n_sub) = HALF * (this % t(1:n_sub) + ONE)
        this % w(1:n_sub) = HALF *  this % w(1:n_sub)
        ! set leftmost entries
        this % t(0) = 0
        this % w(0) = 0
      end select

      ! quadrature weights in [τᵢ₋₁,τᵢ] ........................................

      allocate(this % w_nn (n_sub, 0:n_sub) )

      do i = 1, n_sub
        this % w_nn(i,:) = this % SubintervalWeights(this%t(i-1), this%t(i))
      end do

      ! quadrature weights in [0,τᵢ] ...........................................

      allocate(this % w_0n (n_sub, 0:n_sub) )

      this % w_0n(1,:) = this % w_nn(1,:)
      do i = 2, n_sub
        this % w_0n(i,:) = this % w_0n(i-1,:) + this % w_nn(i,:)
      end do

    end associate

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
    write(io,'(2X,A,T15,A )') 'nodes:'   , this % nodes
    write(io,'(2X,A,T15,I0)') 'n_sub:'   , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweeps:', this % n_sweep

  end subroutine Show_SDC_Method

  !-----------------------------------------------------------------------------
  !> Returns the collocation points transformed to [t0, t0+dt]

  pure function CollocationPoints(this, t0, dt) result(t)
    class(SDC_Method), intent(in) :: this
    real(RNP), intent(in)  :: t0              !< start of the time interval
    real(RNP), intent(in)  :: dt              !< length of the time interval
    real(RNP)              :: t(0:this%p_col) !< intermediate times

    select case(this % nodes)
    case('RR')
      t = t0 + dt * this % t(1:this%n_sub)
    case default
      t = t0 + dt * this % t
    end select

  end function CollocationPoints

  !-----------------------------------------------------------------------------
  !> Returns the quadrature weigths for the collocation points1

  pure function CollocationWeights(this) result(w)
    class(SDC_Method), intent(in) :: this
    real(RNP) :: w(0:this%p_col) !< collocation weigths

    select case(this % nodes)
    case('RR')
      w = this % w(1:this%n_sub)
    case default
      w = this % w
    end select

  end function CollocationWeights

  !-----------------------------------------------------------------------------
  !> Returns the subinterval points transformed to [t0, t0+dt]

  pure function SubintervalPoints(this, t0, dt) result(t)
    class(SDC_Method), intent(in) :: this
    real(RNP), intent(in)  :: t0              !< start of the time interval
    real(RNP), intent(in)  :: dt              !< length of the time interval
    real(RNP)              :: t(0:this%n_sub) !< intermediate times

    t = t0 + dt * this % t

  end function SubintervalPoints

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
    integer   :: j, k, n0, nq

    associate(t => this % t, xq => this % x_quad, wq => this % w_quad )

      nq = ubound(xq,1)

      ! identify first SDC point
      select case(this % nodes)
      case('RR')
        n0 = 1
      case default
        n0 = 0
      end select

      ! metric factor
      delta = HALF * (tb - ta)

      do j = 0, this % n_sub
        ws(j) = 0
        if (j >= n0) then
          do k = 0, nq
            ! tk = τ(xq(k)) = k-th quadrature point mapped to [τa,τb]
            tk  = ta + delta * (xq(k) + 1)
            ! yk = value of j-th Lagrange polynomial to SDC points t(:) at tk
            yk = LagrangePolynomial(j-n0, t(n0:), tk)
            ! add contribution of k-th collocation point
            ws(j) = ws(j) + wq(k) * yk
          end do
          ! scale the weight to match the length of interval [τa,τb]
          ws(j) = delta * ws(j)
        end if
      end do

    end associate

  end function SubintervalWeights

  !=============================================================================

end module Spectral_Deferred_Correction
