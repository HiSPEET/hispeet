!> summary:  Acoustic wave problem for the compressible Navier-Stokes equations
!> author:   Joerg Stiller
!> date:     2023/12/08
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__Problem__CNS__Acoustic_Wave__1D
  use Kind_Parameters, only: RNP
  use Constants,       only: PI
  use Execution_Control
  use CL__Problem__CNS__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_CNS_AcousticWave_1D
  public :: CL_Problem_CNS_AcousticWave_Options_1D

  !-----------------------------------------------------------------------------
  !> Type defining a harmonic acoustic wave

  type, extends(CL_Problem_CNS_1D) :: CL_Problem_CNS_AcousticWave_1D

    real(RNP) :: mach !< mean flow Mach number
    real(RNP) :: p_0  !< mean pressure
    real(RNP) :: T_0  !< mean temperature
    real(RNP) :: c_w  !< wave amplitude (v'/a)
    integer   :: k_w  !< wave number (wave length / domain length)

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetBoundaryValues

  end type CL_Problem_CNS_AcousticWave_1D

  !-----------------------------------------------------------------------------
  !> Acoustic wave options

  type, extends(CL_Problem_CNS_Options_1D) :: &
      CL_Problem_CNS_AcousticWave_Options_1D

    real(RNP) :: xb1  = 0    !< left boundary position
    real(RNP) :: xb2  = 1    !< right boundary position
    real(RNP) :: mach = 0.1  !< mean flow Mach number
    real(RNP) :: p_0  = 1E3  !< mean pressure
    real(RNP) :: T_0  = 3E2  !< mean temperature
    real(RNP) :: c_w  = 1E-3 !< wave amplitude
    integer   :: k_w  = 1    !< wave number

  end type CL_Problem_CNS_AcousticWave_Options_1D

contains

  !---------------------------------------------------------------------------
  !> Initialization of the flow problem

  subroutine SetProblem(this, file)
    class(CL_Problem_CNS_AcousticWave_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    type(CL_Problem_CNS_AcousticWave_Options_1D) :: acoustic_wave_opt

    namelist/cns_acoustic_wave_prm/ acoustic_wave_opt

    logical :: exists, opened
    integer :: prm

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
                  , 'CL__Problem__CNS__Acoustic_Wave__1D'         )
      end if

      ! read parameters ........................................................

      read(prm, nml=cns_acoustic_wave_prm)
      if (.not. opened) close(prm)

    end if

    ! set parameters ...........................................................

    ! initialize base type
    call this % Init_CL_Problem_CNS_1D(acoustic_wave_opt)

    ! enforce periodicity
    this % bc = 'P'

    ! specific parameters
    this % xb1  = acoustic_wave_opt % xb1
    this % xb2  = acoustic_wave_opt % xb2
    this % mach = acoustic_wave_opt % mach
    this % p_0  = acoustic_wave_opt % p_0
    this % T_0  = acoustic_wave_opt % T_0
    this % c_w  = acoustic_wave_opt % c_w
    this % k_w  = acoustic_wave_opt % k_w

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_CNS_AcousticWave_1D), intent(in)  :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:ne,1:nc)

    real(RNP) :: a, c, r, s, v, w, z(3)
    real(RNP) :: kappa
    integer   :: e, i

    associate( r_gas => this % r_gas           &
             , gamma => this % gamma           &
             , p_0   => this % p_0             &
             , T_0   => this % T_0             &
             , c_p   => this % c_p             &
             , po    => cl_operator % eop % po &
             , ne    => cl_operator % ne       &
             , x     => cl_operator % x        )

      kappa = 2 * PI * this%k_w / (this%xb2 - this%xb1)

      a = sqrt(gamma * r_gas * T_0)
      r = a * 2 / (gamma - 1)
      c = a * this % c_w / 2
      v = a * this % mach
      s = c_p * log(T_0) - r_gas * log(p_0)

      do e = 1, ne
      do i = 0, po
        w    = c * (1 - cos(kappa * x(i,e)))
        z(1) = r - v
        z(2) = s
        z(3) = r + v + w
        call this % CharacteristicToConservative(z, u(i,e,:))
      end do
      end do

    end associate

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Provides NO boundary values, as the problem is always periodic

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_CNS_AcousticWave_1D), intent(in)  :: this
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: bv(:,:) !< boundary values

    if (any(this % bc /= 'P')) then
      call Error( 'GetBoundaryValues'                   &
                , 'This problem requires periodic BC'   &
                , 'CL__Problem__CNS__Acoustic_Wave__1D' )
    else
      bv = 0
    end if

  end subroutine GetBoundaryValues

  !=============================================================================

end module CL__Problem__CNS__Acoustic_Wave__1D
