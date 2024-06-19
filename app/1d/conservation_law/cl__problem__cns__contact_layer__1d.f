!> summary:  Shock Tube Problem for the compressible Navier-Stokes equations
!> author:   Benedikt Wex, Joerg Stiller
!> date:     2023/12/14
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__CNS__Contact_Layer__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: PI, ZERO
  use Execution_Control
  use CL__Problem__CNS__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_ContactLayer_1D
  public :: CL_Problem_CNS_ContactLayer_Options_1D

  !-----------------------------------------------------------------------------
  !> Type defining a harmonic acoustic wave

  type, extends(CL_Problem_CNS_1D) :: CL_Problem_CNS_ContactLayer_1D

    real(RNP) :: p_0    !< pressure
    real(RNP) :: T_0    !< mean temperature
    real(RNP) :: dT     !< temperature difference
    real(RNP) :: delta  !< thickness of the layer

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues
    procedure :: GetSources
    procedure :: GetExactSolution

  end type CL_Problem_CNS_ContactLayer_1D

  !-----------------------------------------------------------------------------
  !> Acoustic wave options

  type, extends(CL_Problem_CNS_Options_1D) :: &
      CL_Problem_CNS_ContactLayer_Options_1D

    real(RNP) :: p_0   = 1000        !< pressure
    real(RNP) :: T_0   =   15        !< mean temperature
    real(RNP) :: dT    =    5        !< temperature difference
    real(RNP) :: delta =    1E-1_RNP !< thickness of the layer

  end type CL_Problem_CNS_ContactLayer_Options_1D

contains

  !---------------------------------------------------------------------------
  !> Initialization of the flow problem

  subroutine SetProblem(this, file)
    class(CL_Problem_CNS_ContactLayer_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    type(CL_Problem_CNS_ContactLayer_Options_1D) :: opt
    namelist/cns_contact_layer_prm/ opt

    logical :: exists, opened
    integer :: prm

    ! preset options ...........................................................

    opt % xb1 = -1
    opt % xb2 =  1
    opt % bc  = 'S'

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Error( 'SetProblem'                                  &
                  , 'Input file "'//trim(file)//'.prm" not found' &
                  , 'CL__Problem__CNS__Contact_Layer__1D'         )
      end if

      ! read parameters ........................................................

      read(prm, nml=cns_contact_layer_prm)
      if (.not. opened) close(prm)

    end if

    ! set parameters ...........................................................

    ! initialize base type
    call this % Init_CL_Problem_CNS_1D(opt)

    ! specific parameters
    this % p_0   = opt % p_0
    this % T_0   = opt % T_0
    this % dT    = opt % dT
    this % delta = opt % delta

    this % has_exact_solution = .true.

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_CNS_ContactLayer_1D), intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:ne,1:nc)

    real(RNP) :: T
    integer   :: e, i

    associate( x     => cl_operator % x        &
             , ne    => cl_operator % ne       &
             , po    => cl_operator % eop % po &
             , p_0   => this % p_0             &
             , T_0   => this % T_0             &
             , dT    => this % dT              &
             , delta => this % delta           )

      do e = 1, ne
      do i = 0, po
        T = T_0 + dT * tanh(x(i,e) / delta)
        call this % PrimitiveToConservative(ZERO, T, p_0, u(i,e,:))
      end do
      end do

    end associate

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides boundary values as Values from inside of the domain

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_CNS_ContactLayer_1D), intent(in)  :: this
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: bv(:,:)   !< boundary values

    real(RNP) :: T_b

    if (any(this % bc /= 'S')) then
      call Error( 'GetBoundaryValues'                   &
                , 'This problem requires constant BC'   &
                , 'CL__Problem__CNS__ContactLayer__1D' )
    else

      associate( xb1   => this % xb1   &
               , xb2   => this % xb2   &
               , p_0   => this % p_0   &
               , T_0   => this % T_0   &
               , dT    => this % dT    &
               , delta => this % delta )

         ! left
         T_b = T_0 + dT * tanh(xb1 / delta)
         call this % PrimitiveToConservative(ZERO, T_b, p_0, bv(:,1))

         ! right
         T_b = T_0 + dT * tanh(xb2 / delta)
         call this % PrimitiveToConservative(ZERO, T_b, p_0, bv(:,2))

      end associate
    end if

    ! avoid compiler warning
    if (t > 0) return

  end subroutine GetBoundaryValues

  !-----------------------------------------------------------------------------
  !> Provides the sources

  subroutine GetSources(this, cl_operator, t, u, f_s)
    class(CL_Problem_CNS_ContactLayer_1D),  intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator !< provides x
    real(RNP),             intent(in)  :: t           !< time
    real(RNP), contiguous, intent(in)  :: u(0:,:,:)   !< u(x,t)
    real(RNP), contiguous, intent(out) :: f_s(0:,:,:) !< sources

    integer :: e, i

    associate( x      => cl_operator % x        &
             , ne     => cl_operator % ne       &
             , po     => cl_operator % eop % po &
             , lambda => this % lambda          &
             , dT     => this % dT              &
             , delta  => this % delta           )

      do e = 1, ne
      do i = 0, po
        f_s(i,e,1) = 0
        f_s(i,e,2) = 0
        f_s(i,e,3) = lambda * dT / delta**2 &
                   * (2 * tanh((x(i,e))/delta) / cosh((x(i,e))/delta)**2)
      end do
      end do

    end associate

    ! avoid compiler warning
    if (t > 0 .or. size(u) > 0) return

  end subroutine GetSources

  !-----------------------------------------------------------------------------
  !> Provides the exact solution u(x,t)

  subroutine GetExactSolution(this, cl_operator, t, u)
    class(CL_Problem_CNS_ContactLayer_1D), intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP),             intent(in)  :: t
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u(x,t) at mesh points

    call GetInitialValues(this, cl_operator, u)

    ! just to avoid compiler warnings ;)
    if ( t > 0) return

  end subroutine GetExactSolution

  !=============================================================================

end module CL__Problem__CNS__Contact_Layer__1D
