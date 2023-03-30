module CL__SDC__Method__ISD__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO

  use CL__Problem__Scalar__1D
  use CL__Time_Integrator__1D
  use CL__SDC__Method__1D

  implicit none
  private

  public :: CL_SDC_Method_ISD_1D
  public :: CL_SDC_Options_ISD_1D

  !-----------------------------------------------------------------------------
  !> IMEX ISD SDC ...

  type, extends(CL_SDC_Method_1D) :: CL_SDC_Method_ISD_1D
  contains
    procedure :: Init_CL_SDC_Method_ISD_1D
    procedure :: Show => Show_CL_SDC_Method_ISD_1D
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type CL_SDC_Method_ISD_1D

  ! overloading the constructor
  interface CL_SDC_Method_ISD_1D
    module procedure New_CL_SDC_Method_ISD_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-ISD options

  type, extends(CL_SDC_Options_1D) :: CL_SDC_Options_ISD_1D
  end type CL_SDC_Options_ISD_1D

contains

  !=============================================================================
  ! SDC_Corrector_ISD: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_ISD

  function New_CL_SDC_Method_ISD_1D(pre_opt, sdc_opt) result(this)
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor options
    class(CL_SDC_Options_ISD_1D),       intent(in) :: sdc_opt !< SDC options
    type(CL_SDC_Method_ISD_1D) :: this

    call Init_CL_SDC_Method_ISD_1D(this, pre_opt, sdc_opt)

  end function New_CL_SDC_Method_ISD_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_ISD object

  subroutine Init_CL_SDC_Method_ISD_1D(this, pre_opt, sdc_opt)
    class(CL_SDC_Method_ISD_1D),        intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor options
    class(CL_SDC_Options_ISD_1D),       intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_C_SDC_Method(pre_opt, sdc_opt)

    this % corrector_name = 'ISD IMEX method of order 1'

  end subroutine Init_CL_SDC_Method_ISD_1D

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_CL_SDC_Method_ISD_1D(this, unit)
    class(CL_SDC_Method_ISD_1D), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_C_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_CL_SDC_Method_ISD_1D

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(CL_SDC_Method_ISD_1D), intent(in) :: this
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

    class(CL_SDC_Method_ISD_1D), intent(inout) :: this
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

    complex(RNP) :: S, u1, u2
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
      u1 = u(m-1) + S

    end associate

    ! correction ..............................................................

    u1 = u1 + dt * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
    u2 = u1 / (ONE - dt * lambda%re + (dt * lambda%im)**2)

    u(m) = u2

    ! update RHS
    call this % CorrectorRHS(lambda, dt, u(m), F_ex_new(m), F_im_new(m))

  end subroutine CorrectorStep

  !=============================================================================

end module CL__SDC__Method__ISD__1D
