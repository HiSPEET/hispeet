!> summary:  Wave package  problem for the convection-diffusion equation
!> author:   Joerg Stiller
!> date:     2023/09/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> This problem is based on the wave package
!>
!>     u(x,t) = ∑ᵢ aᵢ sin[κᵢ(x - sᵢ - vt)] exp[-(κᵢ)² nu t]
!>
!> with wave numbers
!>
!>     κᵢ = 2πkᵢ/lw
!>
!> which constitutes an exact solution of the convection diffusion equation.
!> The problem is stated in the domain [0,1] with periodic or any combination
!> of Dirichlet and Neumann conditions.
!>
!> Initialization requires the namelists
!>
!>     - `convection_diffusion_wave_package_prm`
!>     - `wave_package_dim`
!>     - `wave_package_prm`
!>
!> must be specified in a common input file.
!>
!> For more details on the wave package it is referred to the description in the
!> module `Harmonic_Wave_Package`.
!>
!===============================================================================

module CL__Problem__Convection_Diffusion__Wave_Package__1D

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, HALF, ZERO
  use Execution_Control
  use Harmonic_Wave_Package
  use CL__Problem__Convection_Diffusion__1D
  use CL__Operator__1D

  implicit none
  private

  public :: CL_Problem_ConvectionDiffusion_WavePackage_1D
  public :: CL_Problem_ConvectionDiffusion_WavePackage_Options_1D

  !-----------------------------------------------------------------------------
  !> Type for defining and handling 1D convection-diffusion problem

  type, extends(CL_Problem_ConvectionDiffusion_1D) :: &
      CL_Problem_ConvectionDiffusion_WavePackage_1D

    type(HarmonicWavePackage) :: wave !< exact wave package solution

  contains

    procedure :: SetProblem
    procedure :: GetInitialValues
    procedure :: GetExactSolution
    procedure :: GetBoundaryValues

  end type CL_Problem_ConvectionDiffusion_WavePackage_1D

  !-----------------------------------------------------------------------------
  !> 1D convection-diffusion wave package options

  type, extends(CL_Problem_ConvectionDiffusion_Options_1D) :: &
      CL_Problem_ConvectionDiffusion_WavePackage_Options_1D
  end type CL_Problem_ConvectionDiffusion_WavePackage_Options_1D

  !=============================================================================

contains

  !-----------------------------------------------------------------------------
  !> Initialization of the convection-diffusion wave package problem

  subroutine SetProblem(this, file)
    class(CL_Problem_ConvectionDiffusion_WavePackage_1D), intent(inout) :: this
    character(len=*), optional, intent(in) :: file !< (*.prm)

    type(CL_Problem_ConvectionDiffusion_WavePackage_Options_1D) :: opt
    namelist /convection_diffusion_wave_package_prm/ opt

    logical :: exists, opened
    integer :: prm

    ! preset options ...........................................................

    opt % xb1 =  0
    opt % xb2 =  1
    opt % v   =  1
    opt % nu  =  0.1
    opt % bc  = ['D','N']

    ! check for input file .....................................................

    if (present(file)) then

      inquire(file=trim(file)//'.prm', exist=exists, opened=opened, number=prm)

      if (exists .and. .not. opened) then
        open(newunit=prm, file=trim(file)//'.prm', action='READ')
      else if (opened) then
        rewind(prm)
      else
        call Error( 'SetProblem'                                           &
                  , 'Input file "'//trim(file)//'.prm" not found'          &
                  ,  'CL__Problem__Convection_Diffusion__Wave_Package__1D' )
      end if

    else

      call Error( 'SetProblem'                                           &
                , 'Wave package problem requires input file'             &
                ,  'CL__Problem__Convection_Diffusion__Wave_Package__1D' )

    end if

    ! read and set parameters  .................................................

    read(prm, nml=convection_diffusion_wave_package_prm)

    call this % Init_CL_Problem_ConvectionDiffusion_1D(opt)

    this % wave = HarmonicWavePackage(trim(file)//'.prm')
    this % has_exact_solution = .true.

    if (.not. opened) close(prm)

  end subroutine SetProblem

  !-----------------------------------------------------------------------------
  !> Provides the initial values u(x,0)

  subroutine GetInitialValues(this, cl_operator, u)
    class(CL_Problem_ConvectionDiffusion_WavePackage_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u⁰(0:po,1:n1,1:nc)

    associate(x => cl_operator % x)
      call this%wave % Get_Amplitude(this%v, this%nu, x, ZERO, u(:,:,1))
    end associate

  end subroutine GetInitialValues

  !-----------------------------------------------------------------------------
  !> Exact solution u(x,t)

  subroutine GetExactSolution(this, cl_operator, t, u)
    class(CL_Problem_ConvectionDiffusion_WavePackage_1D), intent(in) :: this
    class(CL_Operator_1D), intent(in)  :: cl_operator
    real(RNP),             intent(in)  :: t
    real(RNP), contiguous, intent(out) :: u(0:,:,:) !< u(x,t) at mesh points

    associate(x => cl_operator % x)
      call this%wave % Get_Amplitude(this%v, this%nu, x, t, u(:,:,1))
    end associate

  end subroutine GetExactSolution

  !---------------------------------------------------------------------------
  !> Provides the left and right boundary values for time t

  subroutine GetBoundaryValues(this, t, bv)
    class(CL_Problem_ConvectionDiffusion_WavePackage_1D), intent(in) :: this
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: bv(:,:) !< boundary values

    associate( wave => this % wave &
             , v    => this % v    &
             , nu   => this % nu   )

      ! left boundary
      select case(this % bc(1))
      case('D')
        call wave % Get_Amplitude(v, nu, this%xb1, t, bv(1,1))
      case('N')
        call wave % Get_SpatialDerivative(v, nu, this%xb1, t, dx_u = bv(1,1))
      case default
        bv(1,1) = 0
      end select

      ! right boundary
      select case(this % bc(2))
      case('D')
        call wave % Get_Amplitude(v, nu, this%xb2, t, bv(1,2))
      case('N')
        call wave % Get_SpatialDerivative(v, nu, this%xb2, t, dx_u = bv(1,2))
      case default
        bv(1,2) = 0
      end select

    end associate

  end subroutine GetBoundaryValues

  !=============================================================================

end module CL__Problem__Convection_Diffusion__Wave_Package__1D

