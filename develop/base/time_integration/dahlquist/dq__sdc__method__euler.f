module DQ__SDC__Method__Euler
  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_Euler
  public :: DQ_SDC_Options_Euler

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_Euler
  contains
    procedure :: Init_DQ_SDC_Method_Euler
    procedure :: Show => Show_DQ_SDC_Method_Euler
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_Euler

  ! overloading the constructor
  interface DQ_SDC_Method_Euler
    module procedure New_DQ_SDC_Method_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Euler options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_Euler
  end type DQ_SDC_Options_Euler

contains

  !=============================================================================
  ! SDC_Corrector_Euler: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Euler

  function New_DQ_SDC_Method_Euler(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Euler),     intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_Euler) :: this

    call Init_DQ_SDC_Method_Euler(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Euler object

  subroutine Init_DQ_SDC_Method_Euler(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_Euler),      intent(inout) :: this
    class(DQ_TimeIntegratorOptions), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Euler),     intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % corrector_name = 'IMEX Euler method'

  end subroutine Init_DQ_SDC_Method_Euler

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Euler settings

  subroutine Show_DQ_SDC_Method_Euler(this, unit)
    class(DQ_SDC_Method_Euler), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    call this % Show_DQ_SDC_Method(unit)

  end subroutine Show_DQ_SDC_Method_Euler

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals

  elemental subroutine CorrectorRHS(this, lambda, z, F_ex, F_im)
    class(DQ_SDC_Method_Euler), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    complex(RNP), intent(in)  :: z      !< z
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    select case (this % impl)

    case(0) ! explicit
      F_im = 0
      F_ex = lambda * z

    case(2) ! implicit
      F_im = lambda * z
      F_ex = 0

    case default ! IMEX
      F_im =     lambda % re * z
      F_ex = i * lambda % im * z

    end select

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, t, z &
                          , F_ex_old, F_ex        &
                          , F_im_old, F_im        )

   class(DQ_SDC_Method_Euler), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda        !< λ
    integer     , intent(in)    :: m             !< current SDC interval index
    real   (RNP), intent(in)    :: t(0:)         !< SDC time nodes
    complex(RNP), intent(inout) :: z(0:)         !< z(0:m-1) = zᵏ⁺¹, z(m)=zᵏ → zᵏ⁺¹ , z(m+1:) = zᵏ
    complex(RNP), intent(in)    :: F_ex_old(0:)  !< F^exᵏ
    complex(RNP), intent(inout) :: F_ex(0:)      !< F^exᵏ⁺¹, inout: F_ex(0:m-1), out: F_ex(m)
    complex(RNP), intent(in)    :: F_im_old(0:)  !< F^imᵏ
    complex(RNP), intent(inout) :: F_im(0:)      !< F^imᵏ⁺¹, like F^exᵏ⁺¹

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
        S = S + delta * ( F_im_old(i) + F_ex_old(i) ) * w_sub(i,m)
      end do

      ! z' = z₀ + Sᵏ
      z(m) = z(m-1) + S

    end associate

    ! correction ..............................................................

    select case(this % impl)

    case(0) ! explicit
      z(m) = z(m) + dt * (F_ex(m-1) - F_ex_old(m-1))
    case(2) ! implicit
      z(m) = (z(m) - dt * F_im_old(m)) / (ONE - dt * lambda)
    case default ! IMEX
      z(m) = z(m) + dt * (F_ex(m-1) - F_ex_old(m-1) - F_im_old(m))
      z(m) = z(m) / (ONE - dt * lambda%re)
    end select

    ! update RHS
    call this % CorrectorRHS( lambda           &
                            , z    = z   (m:m) &
                            , F_ex = F_ex(m:m) &
                            , F_im = F_im(m:m) )

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__Euler
