!> summary:  1D harmonic wave packages
!> author:   Joerg Stiller
!> date:     2018/09/26, extended 2023/09/01
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### 1D harmonic wave packages
!>
!> This module provides an exact wave solution for the convection-diffusion
!> equation with constant velocity `v` and diffusivity `nu`. The wave is defined
!> as a package of superimposed sine waves with different wave numbers and phase
!> shifts. For details see the description of the type `HarmonicWavePackage`.
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
  !>     u(x,t) = ∑ᵢ aᵢ sin[κᵢ(x - sᵢ - vt)] exp[-(κᵢ)² nu t]
  !>
  !> with wave numbers
  !>
  !>     κᵢ = 2πkᵢ/lw
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

    generic,   public  :: Get_SpatialDerivative =>   &
                            Get_SpatialDerivative_0, &
                            Get_SpatialDerivative_1, &
                            Get_SpatialDerivative_2
    procedure, private :: Get_SpatialDerivative_0
    procedure, private :: Get_SpatialDerivative_1
    procedure, private :: Get_SpatialDerivative_2

    generic,   public  :: Get_TimeDerivative =>   &
                            Get_TimeDerivative_0, &
                            Get_TimeDerivative_1, &
                            Get_TimeDerivative_2
    procedure, private :: Get_TimeDerivative_0
    procedure, private :: Get_TimeDerivative_1
    procedure, private :: Get_TimeDerivative_2

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
  ! Get_SpatialDerivative

  !-----------------------------------------------------------------------------
  !> Spatial derivatives at given positions and time -- eXplicit version
  !>
  !>     ∂u/∂x   =  ∑ᵢ aᵢ κᵢ    cos[κᵢ(x - sᵢ - vt)] exp[-(κᵢ)² nu t]
  !>     ∂²u/∂x² = -∑ᵢ aᵢ (κᵢ)² sin[κᵢ(x - sᵢ - vt)] exp[-(κᵢ)² nu t]

  subroutine Get_SpatialDerivative_X(this, v, nu, np, x, t, dx_u, dxx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP),           intent(in)  :: v         !< velocity
    real(RNP),           intent(in)  :: nu        !< diffusivity
    integer,             intent(in)  :: np        !< number of points
    real(RNP),           intent(in)  :: x(np)     !< points
    real(RNP),           intent(in)  :: t         !< time
    real(RNP), optional, intent(out) :: dx_u(np)  !< 1st derivative, ∂u/∂x(x,t)
    real(RNP), optional, intent(out) :: dxx_u(np) !< 2nd derivative, ∂²u/∂x²(x,t)

    real(RNP) :: kappa, c, d, z
    integer   :: i

    associate(lw => this%lw, k => this%k, a => this%a, s => this%s)

      if (present(dx_u)) then
        dx_u = 0
        do i = 1, this % nw
          kappa = 2 * PI * k(i) / lw
          c = -kappa ** 2
          d = exp(c * nu * t)
          z = s(i) + v*t
          dx_u = dx_u + a(i) * kappa * d * cos(kappa * (x - z))
        end do
      end if

      if (present(dxx_u)) then
        dxx_u = 0
        do i = 1, this % nw
          kappa = 2 * PI * k(i) / lw
          c = -kappa ** 2
          d = exp(c * nu * t)
          z = s(i) + v*t
          dxx_u = dxx_u + a(i) * c * d * sin(kappa * (x - z))
        end do
      end if

    end associate

  end subroutine Get_SpatialDerivative_X

  !-----------------------------------------------------------------------------
  !> Spatial derivatives at given positions and time -- single point

  subroutine Get_SpatialDerivative_0(this, v, nu, x, t, dx_u, dxx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP),           intent(in)  :: v     !< velocity
    real(RNP),           intent(in)  :: nu    !< diffusivity
    real(RNP),           intent(in)  :: x     !< point
    real(RNP),           intent(in)  :: t     !< time
    real(RNP), optional, intent(out) :: dx_u  !< 1st derivative, ∂u/∂x(x,t)
    real(RNP), optional, intent(out) :: dxx_u !< 2nd derivative, ∂²u/∂x²(x,t)

    real(RNP) :: x_(1), dx_u_(1), dxx_u_(1)

    x_ = x

    call Get_SpatialDerivative_X(this, v, nu, 1, x_, t, dx_u_, dxx_u_)

    if (present( dx_u  ))  dx_u  = dx_u_(1)
    if (present( dxx_u ))  dxx_u = dxx_u_(1)

  end subroutine Get_SpatialDerivative_0

  !-----------------------------------------------------------------------------
  !> Spatial derivatives at given positions and time -- 1d array of points

  subroutine Get_SpatialDerivative_1(this, v, nu, x, t, dx_u, dxx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP),           intent(in)  :: v        !< velocity
    real(RNP),           intent(in)  :: nu       !< diffusivity
    real(RNP),           intent(in)  :: x(:)     !< points
    real(RNP),           intent(in)  :: t        !< time
    real(RNP), optional, intent(out) :: dx_u(:)  !< 1st derivative, ∂u/∂x(x,t)
    real(RNP), optional, intent(out) :: dxx_u(:) !< 2nd derivative, ∂²u/∂x²(x,t)

    call Get_SpatialDerivative_X(this, v, nu, size(x), x, t, dx_u, dxx_u)

  end subroutine Get_SpatialDerivative_1

  !-----------------------------------------------------------------------------
  !> Spatial derivatives at given positions and time -- 2d array of points

  subroutine Get_SpatialDerivative_2(this, v, nu, x, t, dx_u, dxx_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP),           intent(in)  :: v          !< velocity
    real(RNP),           intent(in)  :: nu         !< diffusivity
    real(RNP),           intent(in)  :: x(:,:)     !< points
    real(RNP),           intent(in)  :: t          !< time
    real(RNP), optional, intent(out) :: dx_u(:,:)  !< 1st derivative, ∂u/∂x(x,t)
    real(RNP), optional, intent(out) :: dxx_u(:,:) !< 2nd derivative, ∂²u/∂x²(x,t)

    call Get_SpatialDerivative_X(this, v, nu, size(x), x, t, dx_u, dxx_u)

  end subroutine Get_SpatialDerivative_2

  !=============================================================================
  ! Get_TimeDerivatives

  !-----------------------------------------------------------------------------
  !> Time derivative at given positions and time -- eXplicit version
  !>
  !>     ∂u/∂t   = -∑ᵢ aᵢ ( κᵢ v   cos[κᵢ (x - sᵢ - vt)]
  !>                      + κᵢ² nu sin[κᵢ (x - sᵢ - vt)] ) exp[-κᵢ² nu t]

  subroutine Get_TimeDerivative_X(this, v, nu, np, x, t, dt_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    integer,   intent(in)  :: np        !< number of points
    real(RNP), intent(in)  :: x(np)     !< points
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dt_u(np)  !< time derivative, ∂u/∂t(x,t)

    real(RNP) :: kappa, c, d, z
    integer   :: i

    associate(lw => this%lw, k => this%k, a => this%a, s => this%s)

      dt_u = 0
      do i = 1, this % nw
        kappa = 2 * PI * k(i) / lw
        c = -kappa ** 2
        d = exp(c * nu * t)
        z = s(i) + v*t
        dt_u = dt_u - a(i) * kappa * d * ( cos(kappa * (x - z)) * v          &
                                         + sin(kappa * (x - z)) * kappa * nu )
      end do

    end associate

  end subroutine Get_TimeDerivative_X

  !-----------------------------------------------------------------------------
  !> Time derivative at given positions and time -- single point

  subroutine Get_TimeDerivative_0(this, v, nu, x, t, dt_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v     !< velocity
    real(RNP), intent(in)  :: nu    !< diffusivity
    real(RNP), intent(in)  :: x     !< point
    real(RNP), intent(in)  :: t     !< time
    real(RNP), intent(out) :: dt_u  !< time derivative, ∂u/∂t(x,t)

    real(RNP) :: x_(1), dt_u_(1)

    x_ = x
    call Get_TimeDerivative_X(this, v, nu, 1, x_, t, dt_u_)
    dt_u = dt_u_(1)

  end subroutine Get_TimeDerivative_0

  !-----------------------------------------------------------------------------
  !> Time derivative at given positions and time -- 1d array of points

  subroutine Get_TimeDerivative_1(this, v, nu, x, t, dt_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v       !< velocity
    real(RNP), intent(in)  :: nu      !< diffusivity
    real(RNP), intent(in)  :: x(:)    !< point
    real(RNP), intent(in)  :: t       !< time
    real(RNP), intent(out) :: dt_u(:) !< time derivative, ∂u/∂t(x,t)

    call Get_TimeDerivative_X(this, v, nu, size(x), x, t, dt_u)

  end subroutine Get_TimeDerivative_1

  !-----------------------------------------------------------------------------
  !> Time derivative at given positions and time -- 2d array of points

  subroutine Get_TimeDerivative_2(this, v, nu, x, t, dt_u)
    class(HarmonicWavePackage), intent(in) :: this
    real(RNP), intent(in)  :: v         !< velocity
    real(RNP), intent(in)  :: nu        !< diffusivity
    real(RNP), intent(in)  :: x(:,:)    !< point
    real(RNP), intent(in)  :: t         !< time
    real(RNP), intent(out) :: dt_u(:,:) !< time derivative, ∂u/∂t(x,t)

    call Get_TimeDerivative_X(this, v, nu, size(x), x, t, dt_u)

  end subroutine Get_TimeDerivative_2

  !=============================================================================

end module Harmonic_Wave_Package
