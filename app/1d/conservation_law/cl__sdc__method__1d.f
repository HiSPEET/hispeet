!> summary:  SDC method for 1D conservation laws
!> author:   Robin Fraenzel, Joerg Stiller
!> date:     2023/06/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__SDC__Method__1D

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE
  use Array_Assignments
  use Spectral_Deferred_Correction

  use Gauss_Jacobi
  use Lagrange_Interpolation

  use CL__Problem__1D
  use CL__Operator__1D

  use CL__Time_Integrator__1D
  use CL__Time_Integrator__Euler__1D
  use CL__Time_Integrator__ISD1__1D
  use CL__Time_Integrator__ISD2__1D
  use CL__Time_Integrator__RK__1D

  implicit none
  private

  public :: CL_SDC_Options_1D
  public :: CL_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Type for providing SDC options

  type, extends(SDC_Options) :: CL_SDC_Options_1D
    integer   :: diffusion_method = 1       !< implicit diffusion method
    integer   :: diffusion_i_max  = 10      !< max num iterations
    real(RNP) :: diffusion_r_red  = 1e-10   !< residual reduction
    real(RNP) :: diffusion_r_max  = 1e-12   !< max residual
    logical   :: final_assembly   = .false. !< assembly using collocation method
  end type CL_SDC_Options_1D

  !-----------------------------------------------------------------------------
  !> Modular SDC method with flexible choice of predictor and corrector

  type, abstract, extends(SDC_Method) :: CL_SDC_Method_1D

    class(CL_TimeIntegrator_1D), allocatable :: predictor !< predictor method

    character(len=80) :: corrector_name = ''

    integer   :: diffusion_method  !< implicit diffusion method
    integer   :: diffusion_i_max   !< max num iterations
    real(RNP) :: diffusion_r_red   !< residual reduction
    real(RNP) :: diffusion_r_max   !< max residual
    logical   :: final_assembly    !< assembly using collocation method

  contains

    procedure :: Init_CL_SDC_Method_1D
    procedure :: Show_CL_SDC_Method_1D
    procedure :: TimeStep
    procedure, nopass :: GetHighOrderRHS

    procedure(GetCorrectorRHS), deferred :: GetCorrectorRHS
    procedure(CorrectorStep),   deferred :: CorrectorStep

  end type CL_SDC_Method_1D

  abstract interface

    !---------------------------------------------------------------------------
    !> Computes F_ex and F_im as defined in the corrector
    !>
    !> The initial values `u₀` are used for evaluating nonlinear coefficients,
    !> such as streamline diffusivity.

    subroutine GetCorrectorRHS( this, cl_problem, cl_operator &
                              , t, dt, u_0, u, F_ex, F_im     )
      import
      class(CL_SDC_Method_1D), intent(in)  :: this
      class(CL_Problem_1D),    intent(in)  :: cl_problem
      class(CL_Operator_1D),   intent(in)  :: cl_operator
      real(RNP),               intent(in)  :: t           !< time
      real(RNP),               intent(in)  :: dt          !< time step size
      real(RNP), contiguous,   intent(in)  :: u_0 (:,:,:) !< u₀(x,t)
      real(RNP), contiguous,   intent(in)  :: u   (:,:,:) !< u(x,t)
      real(RNP), contiguous,   intent(out) :: F_ex(:,:,:) !< explicit RHS part
      real(RNP), contiguous,   intent(out) :: F_im(:,:,:) !< implicit RHS part

    end subroutine GetCorrectorRHS

    !---------------------------------------------------------------------------
    !> Execution of a single correction step

    subroutine CorrectorStep( this, cl_problem, cl_operator, m, t, u &
                            , F, F_ex, F_im, F_ex_new, F_im_new, G   )
      import

      class(CL_SDC_Method_1D), intent(in) :: this
      class(CL_Problem_1D),    intent(in) :: cl_problem
      class(CL_Operator_1D),   intent(in) :: cl_operator

      real(RNP), intent(in) :: t(0:)
        !< SDC time nodes
      integer, intent(in) :: m
        !< current SDC interval
      real(RNP), intent(inout) :: u(0:,:,:,0:)
        !< solution, in: [uᵏ⁺¹(:m-1),uᵏ(m:)], out: [uᵏ⁺¹(:m),uᵏ(m+1:)]
      real(RNP), contiguous, intent(in) :: F(0:,:,:,0:)
        !< RHS for high-order quadrature, Fᵏ
      real(RNP), contiguous, intent(in) :: F_ex(0:,:,:,0:)
        !< explicit part of low-order RHS, F_exᵏ
      real(RNP), contiguous, intent(in) :: F_im(0:,:,:,0:)
        !< implicit part of low-order RHS, F_imᵏ
      real(RNP), contiguous, intent(inout) :: F_ex_new(0:,:,:,0:)
        !< explicit of new low-order RHS, in: F_exᵏ⁺¹(0:m-1), out: F_exᵏ⁺¹(0:m)
      real(RNP), contiguous, intent(inout) :: F_im_new(0:,:,:,0:)
        !< implicit of new low-order RHS, in: F_imᵏ⁺¹(0:m-1), out: F_imᵏ⁺¹(0:m)
      real(RNP), contiguous, optional, intent(in) :: G(0:,:,:,0:)
        !< FAS defect correction term

    end subroutine CorrectorStep

  end interface

contains

  !=============================================================================
  ! CL_SDC_Method_1D: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Initialization of CL_SDC_Method_1D object

  subroutine Init_CL_SDC_Method_1D(this, pre_opt, sdc_opt)

    class(CL_SDC_Method_1D),             intent(inout) :: this
    class(CL_TimeIntegrator_Options_1D), intent(in) :: pre_opt !< predictor opts
    class(CL_SDC_Options_1D),            intent(in) :: sdc_opt !< SDC options

    ! parent type initialization ...............................................

    call this % SDC_Method % Init_SDC_Method(sdc_opt)

    ! predictor ................................................................

    select type(pre_opt)
    class is (CL_TimeIntegrator_Options_Euler_1D)
      this % predictor = CL_TimeIntegrator_Euler_1D(pre_opt)
    class is (CL_TimeIntegrator_Options_ISD1_1D)
      this % predictor = CL_TimeIntegrator_ISD1_1D(pre_opt)
    class is (CL_TimeIntegrator_Options_ISD2_1D)
      this % predictor = CL_TimeIntegrator_ISD2_1D(pre_opt)
    class is (CL_TimeIntegrator_Options_RK_1D)
      this % predictor = CL_TimeIntegrator_RK_1D(pre_opt)
    end select

    ! SDC ......................................................................

    this % diffusion_method = sdc_opt % diffusion_method
    this % diffusion_i_max  = sdc_opt % diffusion_i_max
    this % diffusion_r_red  = sdc_opt % diffusion_r_red
    this % diffusion_r_max  = sdc_opt % diffusion_r_max
    this % final_assembly   = sdc_opt % final_assembly

  end subroutine Init_CL_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> Output of CL_SDC_Method_1D settings

  subroutine Show_CL_SDC_Method_1D(this, unit)
    class(CL_SDC_Method_1D), intent(in) :: this
    integer, optional, intent(in)   :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    write(io,'(/,A)') 'CL_SDC_Method_1D settings'
    write(io,'(A,/)') repeat('≡',80)
    write(io,'(2X,A,T15,A )') 'point_set:' , this % point_set
    write(io,'(2X,A,T15,I0)') 'n_col:'     , this % n_col
    write(io,'(2X,A,T15,I0)') 'n_sub:'     , this % n_sub
    write(io,'(2X,A,T15,I0)') 'n_sweep:'   , this % n_sweep

    call this % predictor % Show(unit)

    write(io,'(/,A)') 'Corrector settings'
    write(io,'(A,/)') repeat('=',80)

    write(io,'(2X,2A)') 'name: ', this % corrector_name

    write(io,'(/,A)') 'Solver'
    write(io,'(A,/)') repeat('-',80)
    write(io,'(2X,A,T22,I0)')     'diffusion_method:', this % diffusion_method
    write(io,'(2X,A,T22,I0)')     'diffusion_i_max:' , this % diffusion_i_max
    write(io,'(2X,A,T21,ES12.5)') 'diffusion_r_red:' , this % diffusion_r_red
    write(io,'(2X,A,T21,ES12.5)') 'diffusion_r_max:' , this % diffusion_r_max
    write(io,'(2X,A,T22,L0)')     'final_assembly:'  , this % final_assembly

  end subroutine Show_CL_SDC_Method_1D

  !-----------------------------------------------------------------------------
  !> SDC time step

  subroutine TimeStep(this, cl_problem, cl_operator, dt, t_0, u_0, u)
    class(CL_SDC_Method_1D), intent(in)    :: this
    class(CL_Problem_1D),    intent(in)    :: cl_problem
    class(CL_Operator_1D),   intent(in)    :: cl_operator
    real(RNP),               intent(in)    :: dt          !< step size ∆t
    real(RNP),               intent(in)    :: t_0         !< initial time
    real(RNP), contiguous,   intent(in)    :: u_0(0:,:,:) !< u(t₀)
    real(RNP), contiguous,   intent(inout) :: u  (0:,:,:) !< u(t₀+∆t)

    ! local variables  .........................................................

    real(RNP), allocatable, save :: t_(:)              ! [tᵢ]  SDC nodes
    real(RNP), allocatable, save :: dt_(:)             ! [∆tᵢ]
    real(RNP), allocatable, save :: u_(:,:,:,:)        ! [uᵢ]
    real(RNP), allocatable, save :: F_(:,:,:,:)        ! [Fᵢ]ᵏ
    real(RNP), allocatable, save :: F_ex_(:,:,:,:)     ! [F_exᵢ]ᵏ
    real(RNP), allocatable, save :: F_im_(:,:,:,:)     ! [F_imᵢ]ᵏ
    real(RNP), allocatable, save :: F_ex_new(:,:,:,:)  ! [F_exᵢ]ᵏ⁺¹
    real(RNP), allocatable, save :: F_im_new(:,:,:,:)  ! [F_imᵢ]ᵏ⁺¹

    integer :: i, k

    associate( n_sub   => this % n_sub           &
             , n_sweep => this % n_sweep         &
             , eop     => cl_operator % eop      &
             , po      => cl_operator % eop % po &
             , ne      => cl_operator % ne       &
             , nc      => cl_problem  % nc       )

      ! initialization .........................................................

      allocate(t_       (0:n_sub))
      allocate(dt_      (1:n_sub))
      allocate(u_       (0:po,ne,nc,0:n_sub))
      allocate(F_       (0:po,ne,nc,0:n_sub))
      allocate(F_ex_    (0:po,ne,nc,0:n_sub))
      allocate(F_im_    (0:po,ne,nc,0:n_sub))
      allocate(F_ex_new (0:po,ne,nc,0:n_sub))
      allocate(F_im_new (0:po,ne,nc,0:n_sub))

      t_ (0:) = this % IntermediateTimes(t_0, dt)
      dt_(1:) = t_(1:n_sub) - t_(0:n_sub-1)

      call SetArray(u_(:,:,:,0), u_0, multi = .true.)

      ! predictor ..............................................................

      do i = 1, n_sub
        call this % predictor % TimeStep( cl_problem          &
                                        , cl_operator         &
                                        , dt  = dt_(i)        &
                                        , t_0 = t_(i-1)       &
                                        , u_0 = u_(:,:,:,i-1) &
                                        , u   = u_(:,:,:,i)   )
      end do

      ! corrector ..............................................................

      Corrector: if (n_sweep > 0) then

        do i = 0, n_sub

          call this % GetHighOrderRHS( cl_problem       &
                                     , cl_operator      &
                                     , t  = t_(i)       &
                                     , u  = u_(:,:,:,i) &
                                     , F  = F_(:,:,:,i) )

          call this % GetCorrectorRHS( cl_problem                      &
                                     , cl_operator                     &
                                     , t    = t_    (i)                &
                                     , dt   = dt_   (max(i,1))         &
                                     , u_0  = u_    (:,:,:,max(i-1,0)) &
                                     , u    = u_    (:,:,:,i)          &
                                     , F_ex = F_ex_ (:,:,:,i)          &
                                     , F_im = F_im_ (:,:,:,i)          )

        end do

        ! initialization of new corrector RHS
        call SetArray(F_ex_new(:,:,:,0), F_ex_(:,:,:,0), multi = .true.)
        call SetArray(F_im_new(:,:,:,0), F_im_(:,:,:,0), multi = .true.)

        Sweeps: do k = 1, n_sweep

          do i = 1, n_sub
            call this % CorrectorStep( cl_problem          &
                                     , cl_operator         &
                                     , m        = i        &
                                     , t        = t_       &
                                     , u        = u_       &
                                     , F        = F_       &
                                     , F_ex     = F_ex_    &
                                     , F_im     = F_im_    &
                                     , F_ex_new = F_ex_new &
                                     , F_im_new = F_im_new )
          end do

          if (k == n_sweep .and. .not. this % final_assembly) exit

          do i = 1, n_sub

            call this % GetHighOrderRHS( cl_problem       &
                                       , cl_operator      &
                                       , t  = t_(i)       &
                                       , u  = u_(:,:,:,i) &
                                       , F  = F_(:,:,:,i) )

            call SetArray(F_ex_(:,:,:,i), F_ex_new(:,:,:,i), multi = .true.)
            call SetArray(F_im_(:,:,:,i), F_im_new(:,:,:,i), multi = .true.)

          end do

        end do Sweeps

      end if Corrector

      ! result .................................................................

      if (this % final_assembly) then
        call SetArray(u, u_0, multi = .true.)
        do i = 0, n_sub
          associate(w_i => this%w_col(i,n_sub), u_i => u_(:,:,:,i))
            if (w_i /= 0) then
              call MergeArrays(ONE, u, dt*w_i, u_i, multi=.true.)
            end if
          end associate
        end do
      else
        call SetArray(u, u_(:,:,:,n_sub), multi = .true.)
      end if

      ! clean-up ...............................................................

      deallocate(t_, dt_, u_, F_, F_ex_, F_im_, F_ex_new, F_im_new)

    end associate

  end subroutine TimeStep

  !-----------------------------------------------------------------------------
  !> RHS for high-order quadrature

  subroutine GetHighOrderRHS(cl_problem, cl_operator, t, u, F)
    class(CL_Problem_1D),    intent(in)  :: cl_problem
    class(CL_Operator_1D),   intent(in)  :: cl_operator
    real(RNP),               intent(in)  :: t        !< time
    real(RNP), contiguous,   intent(in)  :: u(:,:,:) !< u(x,t)
    real(RNP), contiguous,   intent(out) :: F(:,:,:) !< RHS

    real(RNP), allocatable, save :: r_c(:,:,:)
    real(RNP), allocatable, save :: r_d(:,:,:)
    real(RNP), allocatable, save :: f_s(:,:,:)
    real(RNP), allocatable, save :: bv(:,:)

    real(RNP), allocatable :: Me_inv(:)

    integer :: e, j

    associate( eop => cl_operator % eop      &
             , dx  => cl_operator % dx       &
             , po  => cl_operator % eop % po &
             , ne  => cl_operator % ne       &
             , nc  => cl_problem  % nc       )

      allocate(r_c, mold = u)
      allocate(r_d, mold = u)
      allocate(f_s, mold = u)
      allocate(bv(nc,2))
      allocate(Me_inv(0:po), source = ONE/(dx/2 * eop%w))

      call cl_problem % GetBoundaryValues (t, bv)
      call cl_problem % GetConvectionTerm (cl_operator, bv, u, r_c)
      call cl_problem % GetDiffusionTerm  (cl_operator, bv, u, r_d)
      call cl_problem % GetSources        (cl_operator, t , u, f_s)

      do j = 1, nc
      do e = 1, ne
        F(:,e,j) = Me_inv * (r_c(:,e,j) + r_d(:,e,j)) + f_s(:,e,j)
      end do
      end do

      deallocate(r_c, r_d, f_s, bv)

    end associate

  end subroutine GetHighOrderRHS

  !=============================================================================

end module CL__SDC__Method__1D
