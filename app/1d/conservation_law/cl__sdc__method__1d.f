!> summary:  SDC method
!> author:   Joerg Stiller
!> date:     2020/05/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__SDC__Method__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Spectral_Deferred_Correction

  use Gauss_Jacobi
  use Lagrange_Interpolation

  use CL__Problem__Scalar__1D           ! For datatype CL_Problem_Scalar_1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
! use CL__Time_Integrator__ISD__1D
  use CL__Time_Integrator__TR__1D
  use CL__Time_Integrator__RK__1D
  use CL__Time_Integrator__TVDRK__1D

  implicit none
  private

  public :: CL_SDC_Options_1D
  public :: CL_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options

  type, extends(SDC_Options) :: CL_SDC_Options_1D
    integer :: impl = 0  !< 0: explicit, 1: implicit, 2/default: IMEX
  end type CL_SDC_Options_1D

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type, abstract, extends(SDC_Method) :: CL_SDC_Method_1D

    class(CL_TimeIntegrator_1D), allocatable :: predictor !< predictor method

    character(len=80) :: corrector_name = ''
    integer :: impl !< switch to explicit/implicit/IMEX corrector (0/1/2)

  contains

    procedure, non_overridable :: Init_C_SDC_Method
    procedure, non_overridable :: Show_C_SDC_Method
    procedure :: TimeStep
    procedure(CorrectorRHS ), deferred :: CorrectorRHS
    procedure(CorrectorRHSSweep ), deferred :: CorrectorRHSSweep !used as long as I do not understand function overloading
    procedure(CorrectorStep), deferred :: CorrectorStep

  end type CL_SDC_Method_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Computes F_ex and F_im as defined in the corrector

    subroutine CorrectorRHS(this, problem, t, dt, u, F_ex, F_im)
      import

      class(CL_SDC_Method_1D), intent(in)    :: this
      class(CL_Problem_Scalar_1D), intent(in) :: problem
      real(RNP), intent(in)    :: t
      real(RNP), intent(in)    :: dt            !< step size ∆t
      real(RNP), intent(inout) :: u(:,:,:)      !< u(t) → u(t+ ∆t)
      real(RNP), intent(out)   :: F_ex(:,:,:)   !< explicit RHS for corrector
      real(RNP), intent(out)   :: F_im(:,:,:)   !< implicit RHS for corrector

     end subroutine CorrectorRHS

    !---------------------------------------------------------------------------
    !> Overload of CorrectorRHS(...)

    subroutine CorrectorRHSSweep(this, problem, t, dt, u, F_ex, F_im)
      import

      class(CL_SDC_Method_1D), intent(in)    :: this
      class(CL_Problem_Scalar_1D), intent(in) :: problem
      real(RNP), intent(in)    :: t(:)
      real(RNP), intent(in)    :: dt(:)           !< step size ∆t
      real(RNP), intent(inout) :: u(:,:,:,:)      !< u(t) → u(t+ ∆t)
      real(RNP), intent(out)   :: F_ex(:,:,:,:)   !< explicit RHS for corrector
      real(RNP), intent(out)   :: F_im(:,:,:,:)   !< implicit RHS for corrector

     end subroutine CorrectorRHSSweep

    !---------------------------------------------------------------------------
    !> Execution of a single correction step

    subroutine CorrectorStep( this, problem, m, t, u, M_inv, F      &
                            , F_ex, F_im, F_ex_new, F_im_new )
      import

      class(CL_SDC_Method_1D), intent(inout) :: this
      class(CL_Problem_Scalar_1D), intent(in) :: problem
      real(RNP), intent(in)    :: t(0:)              !< SDC time nodes
      integer  , intent(in)    :: m                  !< current SDC interval index
      real(RNP), intent(inout) :: u(0:,:,:,0:)           !< uᵏ⁺¹(:m-1),uᵏ→uᵏ⁺¹(m),uᵏ(m+1:)
      real(RNP), intent(in)    :: M_inv(:,:,:)
      real(RNP), intent(in)    :: F(:,:,:,0:)        !< Fᵏ
      real(RNP), intent(in)    :: F_ex(:,:,:,0:)     !< F_exᵏ
      real(RNP), intent(in)    :: F_im(:,:,:,0:)     !< F_imᵏ
      real(RNP), intent(inout) :: F_ex_new(:,:,:,0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
      real(RNP), intent(inout) :: F_im_new(:,:,:,0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)

     end subroutine CorrectorStep

  end interface

contains


  !=============================================================================
  ! CL_SDC_Method_1D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of CL_SDC_Method_1D object

  subroutine Init_C_SDC_Method(this, pre_opt, sdc_opt)

    ! arguments ................................................................

    class(CL_SDC_Method_1D), intent(inout) :: this

    class(CL_TimeIntegrator_Options_1D),  intent(in) :: pre_opt !< predictor options
    class(CL_SDC_Options_1D),  intent(in)           :: sdc_opt !< SDC options

    ! parent type initialization ...............................................

    call this % SDC_Method % Init_SDC_Method(sdc_opt)

    ! predictor ................................................................

    select type(pre_opt)
    class is (CL_TimeIntegrator_Options_Euler_1D)
      this % predictor = CL_TimeIntegrator_Euler_1D(pre_opt)
    !class is (CL_TimeIntegrator_Options_ISD_1D)
      !this % predictor = CL_TimeIntegrator_ISD_1D(pre_opt)
    class is (CL_TimeIntegrator_Options_RK_1D)
      this % predictor = CL_TimeIntegrator_RK_1D(pre_opt)
    class is (CL_TimeIntegrator_Options_TVDRK_1D)
      this % predictor = CL_TimeIntegrator_TVDRK_1D(pre_opt)
    end select

    ! SDC ......................................................................

    this % impl = sdc_opt % impl

  end subroutine Init_C_SDC_Method

  !-----------------------------------------------------------------------------
  !> Output of CL_SDC_Method_1D settings

  subroutine Show_C_SDC_Method(this, unit)
    class(CL_SDC_Method_1D), intent(in) :: this
    integer, optional, intent(in)   :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    !write(io,'(/,A)') 'CL_SDC_Method_1D settings'
    !write(io,'(A,/)') repeat('≡',80)
    write(io,'(2X,A,T15,I0)') 'n_sub:'     , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweeps:'  , this % n_sweep
    write(io,'(2X,A,T15,I0)') 'point_set:' , this % point_set
    write(io,'(2X,A,T15,I0)') 'impl:'      , this % impl

    call this % predictor % Show(unit)

    write(io,'(/,A)') 'Corrector settings'
    write(io,'(A,/)') repeat('=',80)

  end subroutine Show_C_SDC_Method

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, problem, t, dt, u, M_inv)
    ! arguments ................................................................

    class(CL_SDC_Method_1D), intent(inout) :: this
    class(CL_Problem_Scalar_1D), intent(in) :: problem
    real(RNP),    intent(inout) :: t
    real(RNP),    intent(in)    :: dt              !< step size ∆t
    real(RNP),    intent(inout) :: u(0:,:,:)       !< u(t) → u(t+ ∆t)
    real(RNP),    intent(in)    :: M_inv(0:,:,:)    !< transformation matrix

    ! local variables  .........................................................

    real(RNP), allocatable :: t_(:)              ! [tᵢ]  between t and t+dt
    real(RNP), allocatable :: dt_(:)             ! [∆tᵢ] ! dimension(:): seems to be needed because of overgiving it to functions
    real(RNP), allocatable :: u_(:,:,:,:)        ! [uᵢ]
    real(RNP), allocatable :: F_(:,:,:,:)        ! [Fᵢ]ᵏ
    real(RNP), allocatable :: F_ex_(:,:,:,:)     ! [F_exᵢ]ᵏ
    real(RNP), allocatable :: F_im_(:,:,:,:)     ! [F_imᵢ]ᵏ
    real(RNP), allocatable :: F_ex_new(:,:,:,:)  ! [F_exᵢ]ᵏ⁺¹
    real(RNP), allocatable :: F_im_new(:,:,:,:)  ! [F_imᵢ]ᵏ⁺¹

    integer :: n_sub, n_sweep
    integer :: i, j, n

    !---------------------------------------------------------------------------
    ! initialization

    associate( n_sub   => this % n_sub       &
             , n_sweep => this % n_sweep     &
    !bounds ....................................................................
             , po      => problem % eop % po &
             , ne      => problem % ne       &
             , nc      => problem % nc       )

      allocate(t_       (0:n_sub))
      allocate(dt_      (0:n_sub))
      allocate(u_       (0:po,ne,nc,0:n_sub))
      allocate(F_       (0:po,ne,nc,0:n_sub))
      allocate(F_ex_    (0:po,ne,nc,0:n_sub))
      allocate(F_im_    (0:po,ne,nc,0:n_sub))
      allocate(F_ex_new (0:po,ne,nc,0:n_sub))
      allocate(F_im_new (0:po,ne,nc,0:n_sub))


      t_  = this % IntermediateTimes(t, dt)
      dt_(0)       = dt ! never used
      dt_(1:n_sub) = t_(1:n_sub) - t_(0:n_sub-1)

      !---------------------------------------------------------------------------
      ! predictor

      ! u⁰(t_0) = u(t_0)
      u_(:,:,:,0) = u
      do i = 1, n_sub
        u_(:,:,:,i) = u_(:,:,:,i-1)
        call this % predictor % TimeStep(problem, t_(i-1), dt_(i), u_(:,:,:,i), M_inv) !Timestep from TimeIntegrator needs t!

      end do

      !---------------------------------------------------------------------------
      ! corrector

      Corrector: if (n_sweep > 0) then

        ! prerequisites ..........................................................

        ! RHS for high-order quadrature
        do j = 0, n_sub
          F_(:,:,:,j) = problem % RHS_Convection(t_(j),u_(:,:,:,j)) &
                      + problem % RHS_Diffusion (t_(j),u_(:,:,:,j))
        end do

        ! RHS for corrector
        call this % CorrectorRHSSweep(problem, t_, dt_, u_, F_ex_, F_im_) ! F_ex_ = f( u^0_m )
        F_ex_new(:,:,:,0) = F_ex_(:,:,:,0) ! RHS on first sub-timestep
        F_im_new(:,:,:,0) = F_im_(:,:,:,0)

        ! correction sweeps ......................................................

        Sweeps: do n = 1, n_sweep ! in christlieb n->k

          do i = 1, n_sub ! Loop for de-jure-computation of δ^k_m and alternating η^k+1_m over subinterval de-facto one uses η^k_m and η^k+1_m because δ^k+1_m + η^k_m = η^k+1_m
            call this % CorrectorStep( problem             & ! performs a corrector timestep on sub-timesteps in Christlieb they are indicated with m
                                     , m        = i        &
                                     , t        = t_       & ! could cause problems with inout <- actually does not but could
                                     , u        = u_       &
                                     , M_inv    = M_inv    &
                                     , F        = F_       &
                                     , F_ex     = F_ex_    & ! -> F_ex_    = f( t_(i-1), u^(k-1)_i )  i in (1, ..., n_sub)
                                     , F_im     = F_im_    &
                                     , F_ex_new = F_ex_new & ! -> F_ex_new = f( t_(m-1), u^k_m )
                                     , F_im_new = F_im_new )
          end do          ! k <- k+1

          if (n == n_sweep) exit

          do j = n, n_sub
            F_(:,:,:,j) = problem % RHS_Convection(t_(j),u_(:,:,:,j)) &
                        + problem % RHS_Diffusion (t_(j),u_(:,:,:,j)) ! same stuff as in correctorRHS is done
          end do

          F_ex_(:,:,:,1:n_sub) = F_ex_new(:,:,:,1:n_sub)     ! F_ex_ = f( t_(i-1), u^k_i )  i in (1, ..., n_sub)
          F_im_(:,:,:,1:n_sub) = F_im_new(:,:,:,1:n_sub)

        end do Sweeps

      end if Corrector

      ! result ...................................................................

      u = u_(:,:,:,n_sub)

    end associate

  end subroutine TimeStep

  !=============================================================================

end module CL__SDC__Method__1D
