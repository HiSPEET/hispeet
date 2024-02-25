module DQ__SDC__Method__ISD

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_ISD
  public :: DQ_SDC_Options_ISD

  !-----------------------------------------------------------------------------
  !> IMEX ISD SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_ISD
    integer :: method !< method selector
  contains
    procedure :: Init_DQ_SDC_Method_ISD
    procedure :: Show => Show_DQ_SDC_Method_ISD
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_ISD

  ! overloading the constructor
  interface DQ_SDC_Method_ISD
    module procedure New_DQ_SDC_Method_ISD
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-ISD options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_ISD
    integer :: method = 0 !< method, 1: one-stage, default: two-stage
  end type DQ_SDC_Options_ISD

contains

  !=============================================================================
  ! SDC_Corrector_ISD: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_ISD

  function New_DQ_SDC_Method_ISD(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),       intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_ISD) :: this

    call Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_ISD object

  subroutine Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_ISD),        intent(inout) :: this
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),       intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % method = sdc_opt % method

    select case(this % method)
    case(2)
      this % corrector_name = 'ISD method of order 1 with two stages'
    case default
      this % corrector_name = 'ISD method of order 1 with one stage'
    end select

  end subroutine Init_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_DQ_SDC_Method_ISD(this, unit)
    class(DQ_SDC_Method_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_ISD), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    F_im = (lambda % re - HALF * dt * lambda%im ** 2) * u
    F_ex = i * lambda % im * u

    if (this % impl == 0) return ! just to avoid compiler warning !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, t, u , F      &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_ISD), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda       !< λ
    integer     , intent(in)    :: m            !< current SDC interval index
    real   (RNP), intent(in)    :: t(0:)        !< SDC time nodes
    complex(RNP), intent(inout) :: u(0:)        !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
    complex(RNP), intent(in)    :: F(0:)        !< Fᵏ
    complex(RNP), intent(in)    :: F_ex(0:)     !< F_exᵏ
    complex(RNP), intent(in)    :: F_im(0:)     !< F_imᵏ
    complex(RNP), intent(inout) :: F_ex_new(0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    complex(RNP), intent(inout) :: F_im_new(0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

    ! auxiliary variables .....................................................

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP) :: ui, u1, u2, S
    real(RNP)    :: t0, t1, dt
    real(RNP)    :: a_inv, delta
    integer      :: j

    ! initialization ..........................................................

    t0 = t(m-1)
    t1 = t(m)
    dt = t1 - t0    ! subinterval

    ! SDC quadrature ..........................................................

    associate(n_sub => this % n_sub, w_sub => this % w_sub)

      delta = t(n_sub) - t(0)
      S = 0
      do j = 0, n_sub
        S = S + delta * F(j) * w_sub(j,m)
      end do

      ! u' = u₀ + Sᵏ
      ui = u(m-1) + S

    end associate

    ! correction ..............................................................

    a_inv = ONE / (ONE - dt * (lambda%re - dt/2 * lambda%im**2))

    select case(this % method)
    case(2)
      u1 = ui + dt * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
      u1 = u1 * a_inv
      u2 = ui + dt * (i * lambda%im * u1 - F_ex(m) - F_im(m))
      u2 = u2 * a_inv
      u(m) = u2
    case default
      u1 = ui + dt * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
      u1 = u1 * a_inv
      u(m) = u1
    end select

    ! update RHS
    call this % CorrectorRHS(lambda, dt, u(m), F_ex_new(m), F_im_new(m))

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__ISD
