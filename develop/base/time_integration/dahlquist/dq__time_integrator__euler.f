!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Euler method for Dahlquist equation
!> author:   Joerg Stiller
!> date:     2024/07/15
!===============================================================================

module DQ__Time_Integrator__Euler

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use DQ__Time_Integrator

  implicit none
  private

  public :: DQ_TimeIntegrator_Euler
  public :: DQ_TimeIntegrator_Options_Euler

  !-----------------------------------------------------------------------------
  !> IMEX Euler method for Dahlquist equation

  type, extends(DQ_TimeIntegrator) :: DQ_TimeIntegrator_Euler
  contains
    procedure :: Init_DQ_TimeIntegrator_Euler
    procedure :: Show => Show_DQ_TimeIntegrator_Euler
    procedure :: TimeStep
  end type DQ_TimeIntegrator_Euler

  ! overloading the constructor
  interface DQ_TimeIntegrator_Euler
    module procedure New_DQ_TimeIntegrator_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing Euler time-integrator options (none, so far)

  type, extends(DQ_TimeIntegrator_Options) :: DQ_TimeIntegrator_Options_Euler
  end type DQ_TimeIntegrator_Options_Euler

contains

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type DQ_TimeIntegrator_Euler with options

  function New_DQ_TimeIntegrator_Euler(opt) result(this)
    class(DQ_TimeIntegrator_Options_Euler), optional, intent(in) :: opt
    type(DQ_TimeIntegrator_Euler) :: this

    call Init_DQ_TimeIntegrator_Euler(this, opt)

  end function New_DQ_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of a Init_DQ_TimeIntegrator_Euler object

  subroutine Init_DQ_TimeIntegrator_Euler(this, opt)
    class(DQ_TimeIntegrator_Euler),                   intent(inout) :: this
    class(DQ_TimeIntegrator_Options_Euler), optional, intent(in)    :: opt

    ! intialize parent type
    call this % Init_DQ_TimeIntegrator(opt)
    this % name = 'Euler method'

  end subroutine Init_DQ_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Output of DQ_TimeIntegrator_Euler settings

  subroutine Show_DQ_TimeIntegrator_Euler(this, unit)
    class(DQ_TimeIntegrator_Euler), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    ! show parent settings
    call this % Show_DQ_TimeIntegrator(unit)

    write(io,'(2X,A,T15,G0)')  'name:', trim(this % name)
    write(io,'(2X,A,T15,G0)')  'impl:', this % impl

  end subroutine Show_DQ_TimeIntegrator_Euler

  !-----------------------------------------------------------------------------
  !> Performs an IMEX Euler step: u₁ = u₀ + ∆t (iλᵢ u₀ + λᵣ u₁)

  subroutine TimeStep(this, lambda, dt, u)
    class(DQ_TimeIntegrator_Euler), intent(inout) :: this
    real(RNP),    intent(in)    :: dt      !< step size ∆t
    complex(RNP), intent(in)    :: lambda  !< ...
    complex(RNP), intent(inout) :: u       !< u(t) → u(t+ ∆t)

    select case(this%impl)

    case(0)

      ! explicit Euler step: u₁ = u₀ + ∆t λ u₀ .................................

      u = u + dt * lambda * u

    case(1)

      ! implicit Euler step: u₁ = u₀ + ∆t λ u₀ .................................

      u = u / (ONE - dt * lambda)

    case default

      ! IMEX Euler step: u₁ = u₀ + ∆t λ u₀ .....................................

      u = u + dt * (ZERO, ONE) * lambda%im * u  ! explicit convection with i λᵢ
      u = u / (ONE - dt * lambda%re)            ! implicit diffusion with λᵣ

    end select

  end subroutine TimeStep

  !=============================================================================

end module DQ__Time_Integrator__Euler
