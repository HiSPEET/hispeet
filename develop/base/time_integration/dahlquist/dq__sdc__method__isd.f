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
    integer   :: n_stage !< number of stages
    real(RNP) :: c_isd   !< artificial diffusivity factor
    integer   :: scaling !< artificial diffusivity scaling: 0/1/2 none/k/double
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
    integer   :: n_stage = 2  !< number of stages
    real(RNP) :: c_isd   = 1  !< AD amplitude factor
    integer   :: scaling = 0  !< AD scaling: 0/1/2 none/k/double
  end type DQ_SDC_Options_ISD

contains

  !=============================================================================
  ! SDC_Corrector_ISD: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_ISD

  function New_DQ_SDC_Method_ISD(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),        intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_ISD) :: this

    call Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_ISD object

  subroutine Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_ISD),         intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),        intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % n_stage = sdc_opt % n_stage
    this % c_isd   = sdc_opt % c_isd
    this % scaling = sdc_opt % scaling

    write(this%corrector_name,'(A,G0,A)') &
        'ISD method of order 1 with ', this%n_stage, ' stage(s)'

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

    write(io,'(2X,A,T15,G0)') 'name:',    this % corrector_name
    write(io,'(2X,A,T15,G0)') 'c_isd:',   this % c_isd
    write(io,'(2X,A,T15,G0)') 'scaling:', this % scaling

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

    F_im = (lambda % re - this%c_isd * HALF * dt * lambda%im ** 2) * u
    F_ex = i * lambda % im * u

    if (this % impl == 0) return ! just to avoid compiler warning !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, k, t, u , F   &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_ISD), intent(inout) :: this
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
    complex(RNP) :: ui, uj, S
    real(RNP)    :: a_inv,  c_im, dt_isd, dt_step, dt_sub
    integer      :: j

    associate( n_sub => this % n_sub &
             , w_sub => this % w_sub &
             , c_isd => this % c_isd )

      ! initialization .........................................................

      dt_step = t(n_sub) - t(0  )
      dt_sub  = t(m    ) - t(m-1)
      dt_isd  = c_isd * HALF * dt_sub

      c_im = 1
      if (this % scaling > 0) then
        c_im = lambda % re - dt_isd * lambda%im**2
        if (abs(c_im) > 1000 * tiny(ONE)) then
          select case(this % scaling)
          case(1)
            c_im = (lambda % re - k * dt_isd * lambda%im**2) / c_im
          case(2)
            c_im = (lambda % re - 2**(k-1) * dt_isd * lambda%im**2) / c_im
          end select
        end if
      end if

      a_inv = ONE / (ONE - dt_sub * c_im * (lambda%re - dt_isd * lambda%im**2))

      ! SDC quadrature .........................................................

      S = 0
      do j = 0, n_sub
        S = S + dt_step * F(j) * w_sub(j,m)
      end do

      ! u' = u₀ + Sᵏ
      ui = u(m-1) + S

      ! correction .............................................................

      uj = (ui + dt_sub * (F_ex_new(m-1) - F_ex(m-1) - c_im*F_im(m))) * a_inv
      do j = 2, this%n_stage
        uj = (ui + dt_sub * (i * lambda%im * uj - F_ex(m) - c_im * F_im(m))) &
           * a_inv
      end do
      u(m) = uj

      ! update RHS .............................................................

      call this % CorrectorRHS(lambda, dt_sub, u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__ISD
