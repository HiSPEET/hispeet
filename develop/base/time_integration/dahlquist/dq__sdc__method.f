!> summary:  SDC method
!> author:   Joerg Stiller
!> date:     2020/05/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module DQ__SDC__Method

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Spectral_Deferred_Correction

  use Gauss_Jacobi
  use Lagrange_Interpolation

  use DQ__Time_Integrator
  use DQ__Time_Integrator__Euler
  use DQ__Time_Integrator__TR
  use DQ__Time_Integrator__RK

  implicit none
  private

  public :: DQ_SDC_Options
  public :: DQ_SDC_Method

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options

  type, extends(SDC_Options) :: DQ_SDC_Options
    integer :: impl = 2  !< 0: explicit, 1: implicit, 2/default: IMEX
  end type DQ_SDC_Options

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type, abstract, extends(SDC_Method) :: DQ_SDC_Method

    class(DQ_TimeIntegrator), allocatable :: predictor !< predictor method

    character(len=80) :: corrector_name = ''
    integer :: impl !< switch to explicit/implicit/IMEX corrector (0/1/2)

  contains

    procedure, non_overridable :: Init_DQ_SDC_Method
    procedure, non_overridable :: Show_DQ_SDC_Method
    procedure :: TimeStep
    procedure(CorrectorRHS ), deferred :: CorrectorRHS
    procedure(CorrectorStep), deferred :: CorrectorStep

  end type DQ_SDC_Method

  abstract interface

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector and F for subintegrals

  elemental subroutine CorrectorRHS(this, lambda, tau, u, F_ex, F_im)
    import

    class(DQ_SDC_Method), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: tau    !< step width
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

   end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step
  subroutine CorrectorStep( this, lambda, m, t, u  &
                          , F_ex_old, F_ex         &
                          , F_im_old, F_im         )
  import

    class(DQ_SDC_Method), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda        !< λ
    integer     , intent(in)    :: m             !< current SDC interval index
    real   (RNP), intent(in)    :: t(0:)         !< SDC time nodes
    complex(RNP), intent(inout) :: u(0:)         !< u(0:m-1) = uᵏ⁺¹, u(m)=uᵏ → uᵏ⁺¹ , u(m+1:) = uᵏ
    complex(RNP), intent(in)    :: F_ex_old(0:)  !< F^exᵏ
    complex(RNP), intent(in)    :: F_im_old(0:)  !< F^imᵏ
    complex(RNP), intent(inout) :: F_ex(0:)      !< F^exᵏ⁺¹, inout: F_ex(0:m-1), out: F_ex(m)
    complex(RNP), intent(inout) :: F_im(0:)      !< F^imᵏ⁺¹, like F^exᵏ⁺¹

   end subroutine CorrectorStep

  end interface

contains


  !=============================================================================
  ! DQ_SDC_Method: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of DQ_SDC_Method object

  subroutine Init_DQ_SDC_Method(this, pre_opt, sdc_opt)

    ! arguments ................................................................

    class(DQ_SDC_Method), intent(inout) :: this

    class(DQ_TimeIntegratorOptions),  intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options),  intent(in) :: sdc_opt !< SDC options

    ! parent type initialization ...............................................

    call this % SDC_Method % Init_SDC_Method(sdc_opt)

    ! predictor ................................................................

    select type(pre_opt)
    class is (DQ_TimeIntegrator_Options_Euler)
      this % predictor = DQ_TimeIntegrator_Euler(pre_opt)
    class is (DQ_TimeIntegrator_Options_TR)
      this % predictor = DQ_TimeIntegrator_TR(pre_opt)
    class is (DQ_TimeIntegrator_Options_RK)
      this % predictor = DQ_TimeIntegrator_RK(pre_opt)
    end select

    ! SDC ......................................................................

    this % impl = sdc_opt % impl

  end subroutine Init_DQ_SDC_Method

  !-----------------------------------------------------------------------------
  !> Output of DQ_SDC_Method settings

  subroutine Show_DQ_SDC_Method(this, unit)
    class(DQ_SDC_Method), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'DQ_SDC_Method settings'
    write(io,'(A,/)') repeat('≡',80)
    write(io,'(2X,A,T15,I0)') 'n_sub:'     , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweeps:'  , this % n_sweep
    write(io,'(2X,A,T15,I0)') 'point_set:' , this % point_set
    write(io,'(2X,A,T15,I0)') 'impl:'      , this % impl

    call this % predictor % Show(unit)

    write(*,'(/,A,T15,A)')  'Corrector:' , trim(this % corrector_name)

  end subroutine Show_DQ_SDC_Method

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, lambda, dt, u)

    ! arguments ................................................................

    class(DQ_SDC_Method), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    ! local variables  .........................................................

    real(RNP)   , dimension(:), allocatable :: t_     ! [tᵢ]
    real(RNP)   , dimension(:), allocatable :: dt_    ! [∆tᵢ]
    complex(RNP), dimension(:), allocatable :: u_     ! [uᵢ]
    complex(RNP), dimension(:), allocatable :: F_ex_  ! [F_exᵢ]
    complex(RNP), dimension(:), allocatable :: F_im_  ! [F_imᵢ]
    complex(RNP), dimension(:), allocatable :: F_ex_old
    complex(RNP), dimension(:), allocatable :: F_im_old

    integer :: n_sub, n_sweep
    integer :: i, n

    !---------------------------------------------------------------------------
    ! initialization

    n_sub   = this % n_sub
    n_sweep = this % n_sweep

    allocate(t_    (0:n_sub))
    allocate(dt_   (0:n_sub))
    allocate(u_    (0:n_sub))
    allocate(F_ex_ (0:n_sub))
    allocate(F_im_ (0:n_sub))

    allocate(F_ex_old (0:n_sub))
    allocate(F_im_old (0:n_sub))

    t_ = this % IntermediateTimes(ZERO, dt)
    dt_(1:n_sub) = t_(1:n_sub) - t_(0:n_sub-1)
    dt_(0) = dt_(1) ! required only for CorrectorRHS

    !---------------------------------------------------------------------------
    ! predictor

    ! u⁰(t_0) = u(t_0)
    u_(0) = u

    do i = 1, n_sub
      u_(i) = u_(i-1)
      call this % predictor % TimeStep(lambda, dt_(i), u_(i))
    end do

    !---------------------------------------------------------------------------
    ! corrector

    Corrector: if (n_sweep > 0) then

      ! prerequisites ..........................................................

      ! RHS for corrector and subintegrals
      call this % CorrectorRHS( lambda             &
                              , tau     = dt_      &
                              , u       = u_       &
                              , F_ex    = F_ex_old &
                              , F_im    = F_im_old )

      F_ex_(0) = F_ex_old(0)
      F_im_(0) = F_im_old(0)

      ! correction sweeps ......................................................

      Sweeps: do n = 1, n_sweep

        do i = 1, n_sub

          call this % CorrectorStep( lambda                  &
                                   , m          = i          &
                                   , t          = t_         &
                                   , u          = u_         &
                                   , F_ex_old   = F_ex_old   &
                                   , F_ex       = F_ex_      &
                                   , F_im_old   = F_im_old   &
                                   , F_im       = F_im_      )
        end do

        if (n == n_sweep) exit

        F_ex_old(1:n_sub) = F_ex_(1:n_sub)
        F_im_old(1:n_sub) = F_im_(1:n_sub)

      end do Sweeps

    end if Corrector

    ! result ...................................................................

    u = u_(n_sub)

  end subroutine TimeStep

  !=============================================================================

end module DQ__SDC__Method
