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

!> summary:  Picard method for Dahlquist
!> author:   Joerg Stiller
!> date:     2024/07/15
!===============================================================================

module DQ__SDC__Method__Picard

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_Picard
  public :: DQ_SDC_Options_Picard

  !-----------------------------------------------------------------------------
  !> IMEX Picard SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_Picard
  contains
    procedure :: Init_DQ_SDC_Method_Picard
    procedure :: Show => Show_DQ_SDC_Method_Picard
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_Picard

  ! overloading the constructor
  interface DQ_SDC_Method_Picard
    module procedure New_DQ_SDC_Method_Picard
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Picard options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_Picard
  end type DQ_SDC_Options_Picard

contains

  !=============================================================================
  ! SDC_Corrector_Picard: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_Picard

  function New_DQ_SDC_Method_Picard(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Picard),     intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_Picard) :: this

    call Init_DQ_SDC_Method_Picard(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_Picard

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_Picard object

  subroutine Init_DQ_SDC_Method_Picard(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_Picard),      intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_Picard),     intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % corrector_name = 'Picard method'

  end subroutine Init_DQ_SDC_Method_Picard

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Picard settings

  subroutine Show_DQ_SDC_Method_Picard(this, unit)
    class(DQ_SDC_Method_Picard), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_DQ_SDC_Method_Picard

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_Picard), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    select case (this % impl)

    case(0) ! explicit
      F_im = 0
      F_ex = 0 * u

    case(1) ! implicit
      F_im = 0 * u
      F_ex = 0

    case(2) ! IMEX
      F_im =     0 * u
      F_ex = i * 0 * u

    end select

    if (dt > 0) return  ! just to avoid compiler warnings !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, k, t, u , F   &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_Picard), intent(inout) :: this
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

    complex(RNP) :: S
    real(RNP)    :: t0, t1, dt
    real(RNP)    :: delta
    integer      :: j

    ! initialization ..........................................................

    t0 = t(m-1)
    t1 = t(m)
    dt = t1 - t0    ! subinterval

    ! SDC quadrature ..........................................................

    associate(n_sub => this % n_sub, w_nn => this % w_nn)

      delta = t(n_sub) - t(0)
      S = 0
      do j = 0, n_sub
        S = S + delta * F(j) * w_nn(m,j)
      end do

      ! u = u₀ + Sᵏ
      u(m) = u(m-1) + S

    end associate

    ! update RHS
    call this % CorrectorRHS(lambda, dt, u(m), F_ex_new(m), F_im_new(m))

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__Picard
