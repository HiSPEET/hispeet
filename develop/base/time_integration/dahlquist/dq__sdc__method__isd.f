module DQ__SDC__Method__ISD
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_ISD
  public :: DQ_SDC_Options_ISD

  !-----------------------------------------------------------------------------
  !> IMEX ISD SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_ISD
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

    this % corrector_name = 'ISD IMEX method of order 1'

  end subroutine Init_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_DQ_SDC_Method_ISD(this, unit)
    class(DQ_SDC_Method_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    call this % Show_DQ_SDC_Method(unit)

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

    F_im = (lambda % re - dt * lambda % im ** 2) * u
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

    complex(RNP) :: S
    real(RNP)    :: t0, t1, dt
    real(RNP)    :: delta
    integer      :: i

    ! initialization ..........................................................

    t0 = t(m-1)
    t1 = t(m)
    dt = t1 - t0    ! subinterval

    ! SDC quadrature ..........................................................

    associate(n_sub => this % n_sub, w_sub => this % w_sub)

      delta = t(n_sub) - t(0)
      S = 0
      do i = 0, n_sub
        S = S + delta * F(i) * w_sub(i,m)
      end do

      ! u' = u₀ + Sᵏ
      u(m) = u(m-1) + S

    end associate

    ! correction ..............................................................

    u(m) = u(m) + dt * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
    u(m) = u(m) / (ONE - dt * lambda%re + (dt * lambda%im)**2)

    ! update RHS
    call this % CorrectorRHS(lambda, dt, u(m), F_ex_new(m), F_im_new(m))

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__ISD
