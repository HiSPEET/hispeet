module DQ__SDC__Method__Diag

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters,  only: RNP, RHP
  use Constants,        only: ONE, ZERO
  use Execution_Control
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_Diag
  public :: DQ_SDC_Options_Diag

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_Diag
    integer :: variant !< 1: MIN-SR-NS, 2: MIN-SR-FLEX
  contains
    procedure :: Init_DQ_SDC_Method_Diag
    procedure :: Show => Show_DQ_SDC_Method_Diag
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_Diag

  ! overloading the constructor
  interface DQ_SDC_Method_Diag
    module procedure New_DQ_SDC_Method_Diag
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Euler options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_Diag
    integer :: variant = 2 !< 1: MIN-SR-NS, 2: MIN-SR-FLEX
  end type DQ_SDC_Options_Diag

contains

  !=============================================================================
  ! SDC_Corrector_Diag: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Diag

  function New_DQ_SDC_Method_Diag(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Diag),       intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_Diag) :: this

    call Init_DQ_SDC_Method_Diag(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_Diag

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Diag object

  subroutine Init_DQ_SDC_Method_Diag(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_Diag),        intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Diag),       intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % variant = sdc_opt % variant

    select case(this%variant)
    case(1)
      this % corrector_name = 'Diag MIN-SR-NS'
    case(2)
      this % corrector_name = 'Diag MIN-SR-FLEX'
    end select

  end subroutine Init_DQ_SDC_Method_Diag

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Diag settings

  subroutine Show_DQ_SDC_Method_Diag(this, unit)
    class(DQ_SDC_Method_Diag), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_DQ_SDC_Method_Diag

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_Diag), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    select case (this % impl)

    case(0) ! explicit
      F_im = 0
      F_ex = lambda * u

    case(1) ! implicit
      F_im = lambda * u
      F_ex = 0

    case(2) ! IMEX
      F_im =     lambda % re * u
      F_ex = i * lambda % im * u

    case(3) ! IMEXa la Lunet
      F_im = lambda % re * u
      F_ex = 0

    end select

    if (dt > 0) return  ! just to avoid compiler warnings !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, k, t, u , F   &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_Diag), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda       !< λ
    integer     , intent(in)    :: m            !< current SDC interval index
    integer     , intent(in)    :: k            !< current corrector sweep
    real   (RNP), intent(in)    :: t(0:)        !< SDC time nodes
    complex(RNP), intent(inout) :: u(0:)        !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
    complex(RNP), intent(in)    :: F(0:)        !< Fᵏ
    complex(RNP), intent(in)    :: F_ex(0:)     !< F_exᵏ
    complex(RNP), intent(in)    :: F_im(0:)     !< F_imᵏ
    complex(RNP), intent(inout) :: F_ex_new(0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    complex(RNP), intent(inout) :: F_im_new(0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

    ! auxiliary variables .....................................................

    complex(RNP), parameter :: i = (ZERO, ONE)

    complex(RNP) :: u1, u2, S
    real(RNP)    :: dt_sub, dt_step, q_del
    integer      :: j

    associate( n_sub => this % n_sub &
             , n_col => this % n_col &
             , w_col => this % w_col &
             , tau   => this % t     )

      ! initialization ........................................................

      dt_sub  = t(m) - t(m-1)
      dt_step = t(n_sub) - t(0)

      select case(this % variant)
      case(1)
        q_del = tau(m) / n_col
      case(2)
        q_del = tau(m) / min(k, n_col)
      end select

      ! SDC quadrature ........................................................

      S = 0
      do j = 0, n_sub
        S = S + dt_step * F(j) * w_col(j,m)
      end do

      u1 = u(0) + S

     ! correction ............................................................

      select case(this % impl)

      case(0) ! explicit
        u2 = u1 + dt_step * q_del * (F_ex_new(m-1) - F_ex(m-1))

      case(1) ! implicit
        u2 = (u1 - dt_step * q_del * F_im(m)) &
           / (ONE - dt_step * q_del * lambda)

      case(2) ! IMEX
        u2 = (u1 + dt_step * q_del * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))) &
           / (ONE - dt_step * q_del * lambda%re)

      case(3) ! IMEX a la Lunet
        u2 = (u1 - dt_step * q_del * F_im(m)) &
           / (ONE - dt_step * q_del * lambda%re)

      end select

      u(m) = u2

      ! update RHS .............................................................

      call this % CorrectorRHS(lambda, dt_sub, u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__Diag
