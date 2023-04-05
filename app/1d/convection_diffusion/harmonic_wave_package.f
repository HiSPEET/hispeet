!> summary:  1D harmonic wave packages
!> author:   Joerg Stiller
!> date:     2018/09/26
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### 1D harmonic wave packages
!>
!> This module provides an exact wave solution for the convection-diffusion
!> equation with constant velocity `v` and diffusivity `nu`. The wave is defined
!> as a package of superimposed sine waves with different wave numbers and phase
!> shifts. For details see the description of the type `HarmonicWavePackage`.
!>
!> @todo
!>   *  Extend Get_FirstDerivative to second and time derivative using optional
!>      arguments as needed
!> @endtodo
!===============================================================================

module Harmonic_Wave_Package
  use Kind_Parameters,   only: RNP
  use Constants,         only: ONE, ZERO, PI
  use Execution_Control, only: Error
  implicit none
  private

  public :: HarmonicWavePackage

  !-----------------------------------------------------------------------------
  !> Provides a harmonic wave package
  !>
  !> The package is defined as a superposition of sine waves in the form
  !>
  !>     u(x,t) = ∑ᵢ aᵢ sin[2πkᵢ/lw (x - sᵢ - vt)] exp[-(2πkᵢ)² nu t]
  !>
  !> The velocity `v` and diffusivity `nu` are considered as external
  !> parameters and, hence, not included in the definition of the wave
  !> package.

  type HarmonicWavePackage

    private
    integer   :: nw = 0            !< number of waves
    real(RNP) :: lw = 1            !< length of waves with k=1
    integer,   allocatable :: k(:) !< integer wave numbers kᵢ
    real(RNP), allocatable :: a(:) !< maximum amplitudes aᵢ
    real(RNP), allocatable :: s(:) !< start positions sᵢ

  contains

    generic,   public  :: Get_Amplitude =>   &
                            Get_Amplitude_0, &
                            Get_Amplitude_1, &
                            Get_Amplitude_2
    procedure, private :: Get_Amplitude_0
    procedure, private :: Get_Amplitude_1
    procedure, private :: Get_Amplitude_2

    generic,   public  :: Get_FirstDerivative =>   &
                            Get_FirstDerivative_0, &
                            Get_FirstDerivative_1, &
                            Get_FirstDerivative_2
    procedure, private :: Get_FirstDerivative_0
    procedure, private :: Get_FirstDerivative_1
    procedure, private :: Get_FirstDerivative_2

  end type HarmonicWavePackage

  ! overload constructors
  interface HarmonicWavePackage
    module procedure :: New_from_Args
    module procedure :: New_from_File
  end interface

contains

  !=============================================================================
  ! Constructors

  !-----------------------------------------------------------------------------
  !> Generate new wave package from arguments

  type(HarmonicWavePackage) function New_from_Args(nw, lw, k, a, s) result(this)
    integer,   intent(in) :: nw    !< number of waves
    real(RNP), intent(in) :: lw    !< length of waves with k=1
    integer,   intent(in) :: k(nw) !< integer wave numbers kᵢ
    real(RNP), intent(in) :: a(nw) !< start amplitudes aᵢ
    real(RNP), intent(in) :: s(nw) !< start positions sᵢ

    this % nw = nw
    this % lw = lw

    allocate(this % k, source = k)
    allocate(this % a, source = a)
    allocate(this % s, source = s)

  end function New_from_Args

  !-----------------------------------------------------------------------------
  !> Generate new wave package from namelist input specified in a file
  !>
  !> The `wave_file` must contain sections for the following two namelists
  !>
  !>     namelist /wave_package_dim/ nw, lw
  !>     namelist /wave_package_prm/ k, a, s
  !>
  !> e.g.,
  !>
  !>     &wave_package_dim
  !>       nw = 2             ! default: 0
  !>       lw = 1.0D0         ! default: 1
  !>     /
  !>     &wave_package_prm
  !>       k = 1, 3           ! default: 1
  !>       a = 1.0D0, 0.5D0   ! default: 1
  !>       s = 0.0D0, 0.2D0   ! default: 0
  !>     /
  !>
  !> If the file is open, it will be read from the current position.
  !> Otherwise, it will be opened for read and closed afterwards.

  type(HarmonicWavePackage) function New_from_File(wave_file) result(this)
    character(len=*), intent(in) :: wave_file !< input file

    integer,   allocatable :: k(:)
    real(RNP), allocatable :: a(:)
    real(RNP), allocatable :: s(:)

    logical   :: opened, exists
    integer   :: io
    integer   :: nw = 0
    real(RNP) :: lw = 1

    namelist /wave_package_dim/ nw, lw
    namelist /wave_package_prm/ k, a, s

    inquire(file=wave_file, opened=opened, exist=exists, number=io)

    if (.not. opened) then
      if (exists) then
        open(newunit=io, file=wave_file, action='READ')
      else
        call Error('New_from_File',                                         &
                   'input file "' // trim(wave_file) // '" does not exist', &
                   'Harmonic_Wave_Package'                                  )
      end if
    end if

    read(io, nml=wave_package_dim)

    allocate(k(nw), source = 1)
    allocate(a(nw), source = ONE)
    allocate(s(nw), source = ZERO)

    read(io, nml=wave_package_prm)

    this = New_from_Args(nw, lw, k, a, s)

    ! close IO unit if file was closed on entry
    if (.not. opened) then
      close(io)
    end if

  end function New_from_File

  !=============================================================================
  ! Get_Amplitude

  !-----------------------------------------------------------------------------
  !> Amplitude at given positions and time -- eXplicit version

  subroutine Get_Amplitude_X(this, v, nu, np, x, t, u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v      !< velocity
    real(RNP), intent(in)  :: nu     !< diffusivity
    integer,   intent(in)  :: np     !< number of points
    real(RNP), intent(in)  :: x(np)  !< points
    real(RNP), intent(in)  :: t      !< time
    real(RNP), intent(out) :: u(np)  !< amplitude u(x,t)

    integer   :: i
    real(RNP) :: c, d, z

    associate(lw => this%lw, k => this%k, a => this%a, s => this%s)

      u = 0
      do i = 1, this % nw
        c = 2 * PI * k(i) / lw
        d = 1 / exp((2 * PI * k(i))**2 * nu * t)
        z = s(i) + v*t
        u = u + a(i) * d * sin(c * (x - z))
      end do

    end associate

  end subroutine Get_Amplitude_X

  !-----------------------------------------------------------------------------
  !> Amplitude at given positions and time -- single point

  subroutine Get_Amplitude_0(this, v, nu, x, t, u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v      !< velocity
    real(RNP), intent(in)  :: nu     !< diffusivity
    real(RNP), intent(in)  :: x      !< point
    real(RNP), intent(in)  :: t      !< time
    real(RNP), intent(out) :: u      !< amplitudes u(x,t)

    real(RNP) :: x_(1), u_(1)

    x_ = x
    call Get_Amplitude_X(this, v, nu, 1, x_, t, u_)
    u = u_(1)

  end subroutine Get_Amplitude_0

  !-----------------------------------------------------------------------------
  !> Amplitude at given positions and time -- 1d array of points

  subroutine Get_Amplitude_1(this, v, nu, x, t, u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v      !< velocity
    real(RNP), intent(in)  :: nu     !< diffusivity
    real(RNP), intent(in)  :: x(:)   !< points
    real(RNP), intent(in)  :: t      !< time
    real(RNP), intent(out) :: u(:)   !< amplitudes u(x,t)

    call Get_Amplitude_X(this, v, nu, size(x), x, t, u)

  end subroutine Get_Amplitude_1

  !-----------------------------------------------------------------------------
  !> Amplitude at given positions and time -- 2d array of points

  subroutine Get_Amplitude_2(this, v, nu, x, t, u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v      !< velocity
    real(RNP), intent(in)  :: nu     !< diffusivity
    real(RNP), intent(in)  :: x(:,:) !< points
    real(RNP), intent(in)  :: t      !< time
    real(RNP), intent(out) :: u(:,:) !< amplitudes u(x,t)

    call Get_Amplitude_X(this, v, nu, size(x), x, t, u)

  end subroutine Get_Amplitude_2

  !=============================================================================
  ! Get_FirstDerivative

  !-----------------------------------------------------------------------------
  !> Derivative at given positions and time -- eXplicit version
  !>
  !>     ∂u/∂x = ∑ᵢ aᵢ 2πkᵢ/lw cos[2πkᵢ/lw (x - sᵢ - vt)] exp[-(2πkᵢ)² nu t]

  subroutine Get_FirstDerivative_X(this, v, nu, np, x, t, dx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    integer,   intent(in)  :: np        !< number of points
    real(RNP), intent(in)  :: x(np)     !< points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dx_u(np)  !< derivatives ∂u/∂x(x,t)

    integer   :: i
    real(RNP) :: c, d, z

    associate(lw => this%lw, k => this%k, a => this%a, s => this%s)

      dx_u = 0
      do i = 1, this % nw
        c = 2 * PI * k(i) / lw
        d = 1 / exp((2 * PI * k(i))**2 * nu * t)
        z = s(i) + v*t
        dx_u = dx_u + a(i) * c * d * cos(c * (x - z))
      end do

    end associate

  end subroutine Get_FirstDerivative_X

  !-----------------------------------------------------------------------------
  !> Derivative at given positions and time -- single point

  subroutine Get_FirstDerivative_0(this, v, nu, x, t, dx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(in)  :: x         !< point
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dx_u      !< derivative ∂u/∂x(x,t)

    real(RNP) :: x_(1), dx_u_(1)

    x_ = x
    call Get_FirstDerivative_X(this, v, nu, 1, x_, t, dx_u_)
    dx_u = dx_u_(1)

  end subroutine Get_FirstDerivative_0

  !-----------------------------------------------------------------------------
  !> Derivative at given positions and time -- 1d array of points

  subroutine Get_FirstDerivative_1(this, v, nu, x, t, dx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(in)  :: x(:)      !< points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dx_u(:)   !< derivatives ∂u/∂x(x,t)

    call Get_FirstDerivative_X(this, v, nu, size(x), x, t, dx_u)

  end subroutine Get_FirstDerivative_1

  !-----------------------------------------------------------------------------
  !> Derivative at given positions and time -- 2d array of points

  subroutine Get_FirstDerivative_2(this, v, nu, x, t, dx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(in)  :: x(:,:)    !< points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dx_u(:,:) !< derivatives ∂u/∂x(x,t)

    call Get_FirstDerivative_X(this, v, nu, size(x), x, t, dx_u)

  end subroutine Get_FirstDerivative_2

  !=============================================================================

end module Harmonic_Wave_Package
