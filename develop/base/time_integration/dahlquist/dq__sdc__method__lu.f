!> summary:  SDC method for consercation laws based on LU decomposition
!> author:   Joerg Stiller
!> date:     2024/07/15
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This module implements the implicit SDC method of
!> [Weiser 2015](https://link.springer.com/article/10.1007/s10543-014-0540-y)
!> and provides corresponding explicit and IMEX versions.
!===============================================================================

module DQ__SDC__Method__LU

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters,  only: RNP, RHP
  use Constants,        only: ONE, ZERO
  use Linear_Equations, only: LU_Decomposition
  use Execution_Control
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_LU
  public :: DQ_SDC_Options_LU

  !-----------------------------------------------------------------------------
  !> SDC corrector using an LU decomposition of the quadrature matrix

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_LU
    real(RNP), allocatable :: q_del(:,:) !< LU-based iteration matrix
  contains
    procedure :: Init_DQ_SDC_Method_LU
    procedure :: Show => Show_DQ_SDC_Method_LU
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_LU

  ! overloading the constructor
  interface DQ_SDC_Method_LU
    module procedure New_DQ_SDC_Method_LU
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Euler options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_LU
  end type DQ_SDC_Options_LU

contains

  !=============================================================================
  ! SDC_Corrector_LU: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_LU

  function New_DQ_SDC_Method_LU(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_LU),         intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_LU) :: this

    call Init_DQ_SDC_Method_LU(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_LU

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_LU object

  subroutine Init_DQ_SDC_Method_LU(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_LU),          intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_LU),         intent(in) :: sdc_opt !< SDC options

    real(RNP), allocatable :: St(:,:)
    real(RHP), allocatable :: LU(:,:)
    integer :: i, j

    ! sanity check .............................................................

    if (sdc_opt % nodes /= 'RR') then
      call Error('Init_DQ_SDC_Method_LU','LU corrector requires RR points')
    end if

    ! intialize parent type ....................................................

    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    ! initialize specific components ...........................................

    this % corrector_name = 'LU'

    allocate(this % q_del(this%n_sub,this%n_sub), source = ZERO)

    associate( tau   => this % t     &
             , n_sub => this % n_sub &
             , w_nn  => this % w_nn  &
             , q_del => this % q_del )

      allocate(St(n_sub,n_sub))
      allocate(LU(n_sub,n_sub))

      do j = 1, n_sub
        St(1:n_sub,j) = w_nn(j,1:n_sub) / (tau(j) - tau(j-1))
      end do

      call LU_Decomposition(St, LU)

      do i = 1, n_sub
      do j = 1, i
        q_del(i,j) = real(LU(j,i), RNP)
      end do
      end do

    end associate

  end subroutine Init_DQ_SDC_Method_LU

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_LU settings

  subroutine Show_DQ_SDC_Method_LU(this, unit)
    class(DQ_SDC_Method_LU), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_DQ_SDC_Method_LU

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_LU), intent(in) :: this
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

    end select

    if (dt > 0) return  ! just to avoid compiler warnings !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, k, t, u , F   &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_LU), intent(inout) :: this
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
    real(RNP)    :: dt_sub, dt_step
    integer      :: j

    associate( n_sub => this % n_sub &
             , w_nn  => this % w_nn  &
             , q_del => this % q_del )

      ! initialization ........................................................

      dt_sub  = t(m) - t(m-1)
      dt_step = t(n_sub) - t(0)

      ! SDC quadrature ........................................................

      S = 0
      do j = 0, n_sub
        S = S + dt_step * F(j) * w_nn(m,j)
      end do

      ! u' = u₀ + Sᵏ
      u1 = u(m-1) + S

     ! correction ............................................................

      select case(this % impl)

      case(0) ! explicit
        u2 = u1
        do j = 2, m
          u2 = u2 + dt_sub * q_del(m,j) * (F_ex_new(j-1) - F_ex(j-1))
        end do

      case(1) ! implicit
        u2 = u1
        do j = 1, m-1
          u2 = u2 + dt_sub * q_del(m,j) * (F_im_new(j) - F_im(j))
        end do
        u2 = (u2 - dt_sub * q_del(m,m) * F_im(m)) &
           / (ONE - dt_sub * q_del(m,m) * lambda)

      case(2)  ! IMEX
        u2 = u1
        do j = 1, m-1
          u2 = u2 + dt_sub * q_del(m,j) * ( F_ex_new(j-1) - F_ex(j-1) &
                                          + F_im_new(j  ) - F_im(j  ) )
        end do
        u2 = u2 + dt_sub * q_del(m,m) * ( F_ex_new(m-1) - F_ex(m-1) &
                                                        - F_im(m  ) )
        u2 = u2 / (ONE - dt_sub * q_del(m,m) * lambda%re)

     !case(3) ! IMEX partitioned -- TBD

      end select

      u(m) = u2

      ! update RHS .............................................................

      call this % CorrectorRHS(lambda, dt_sub, u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__LU
