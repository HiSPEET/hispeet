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

!> summary:  SDC method for Dahlquist using implicit streamline-diffusion
!> author:   Joerg Stiller
!> date:     2024/07/15
!===============================================================================

module DQ__SDC__Method__ISD

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE
  use DQ__Time_Integrator
  use DQ__SDC__Method

  implicit none
  private

  public :: DQ_SDC_Method_ISD
  public :: DQ_SDC_Options_ISD

  !-----------------------------------------------------------------------------
  !> IMEX ISD SDC ...

  type, extends(DQ_SDC_Method) :: DQ_SDC_Method_ISD
    integer   :: n_stage !< number of stages
    integer   :: ad_typ  !< AD type: 0/1 = none/LW
    real(RNP) :: s_base  !< AD baseline scaling factor
    real(RNP) :: s_grow  !< AD per-sweep growth factor
    real(RNP) :: s_max   !< AD maximum scaling factor
  contains
    procedure :: Init_DQ_SDC_Method_ISD
    procedure :: Show => Show_DQ_SDC_Method_ISD
    procedure :: CorrectorRHS
    procedure :: CorrectorStep
  end type DQ_SDC_Method_ISD

  ! overloading the constructor
  interface DQ_SDC_Method_ISD
    module procedure New_DQ_SDC_Method_ISD
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-ISD options

  type, extends(DQ_SDC_Options) :: DQ_SDC_Options_ISD
    integer   :: n_stage = 2  !< number of stages
    integer   :: ad_typ  = 1  !< AD type: 0/1 = none/LW
    real(RNP) :: s_base  = 1  !< AD baseline scaling factor
    real(RNP) :: s_grow  = 0  !< AD per-sweep growth factor
    real(RNP) :: s_max   = 3  !< AD maximum scaling factor
  end type DQ_SDC_Options_ISD

contains

  !=============================================================================
  ! SDC_Corrector_ISD: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for objects of type SDC_Corrector_ISD

  function New_DQ_SDC_Method_ISD(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),        intent(in) :: sdc_opt !< SDC options
    type(DQ_SDC_Method_ISD) :: this

    call Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)

  end function New_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Initialization of a SDC_Corrector_ISD object

  subroutine Init_DQ_SDC_Method_ISD(this, pre_opt, sdc_opt)
    class(DQ_SDC_Method_ISD),         intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_SDC_Options_ISD),        intent(in) :: sdc_opt !< SDC options

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % n_stage = sdc_opt % n_stage
    this % ad_typ  = max(0, min(2, sdc_opt % ad_typ))
    this % s_base  = sdc_opt % s_base
    this % s_grow  = sdc_opt % s_grow
    this % s_max   = sdc_opt % s_max

    write(this%corrector_name,'(A,G0,A)') &
        'ISD method of order 1 with ', this%n_stage, ' stage(s)'

  end subroutine Init_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_DQ_SDC_Method_ISD(this, unit)
    class(DQ_SDC_Method_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name
    select case(this % ad_typ)
    case(1)
      write(io,'(2X,A,T15,G0)') 'AD type:', 'LW'
    case default
      write(io,'(2X,A,T15,G0)') 'AD type:', 'none'
    end select
    write(io,'(2X,A,T15,G0)') 's_base =', this % s_base
    write(io,'(2X,A,T15,G0)') 's_grow =', this % s_grow
    write(io,'(2X,A,T15,G0)') 's_max  =', this % s_max

  end subroutine Show_DQ_SDC_Method_ISD

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector, but without AD

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_SDC_Method_ISD), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    F_im = lambda % re * u
    F_ex = i * lambda % im * u

    if (this % impl == 0 .or. dt > 0) return ! just to avoid compiler warning !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStep( this, lambda, m, k, t, u , F   &
                          , F_ex, F_im, F_ex_new, F_im_new )

    class(DQ_SDC_Method_ISD), intent(inout) :: this
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
    complex(RNP) :: u0, ui, uj, S
    real(RNP)    :: a_inv, dt_step, dt_sub, mu
    integer      :: j

    associate( n_sub  => this % n_sub  &
             , w_nn   => this % w_nn   &
             , ad_typ => this % ad_typ &
             , s_base => this % s_base &
             , s_grow => this % s_grow &
             , s_max  => this % s_max  )

      ! initialization .........................................................

      u0 = u(m)

      dt_step = t(n_sub) - t(0  )
      dt_sub  = t(m    ) - t(m-1)

      ! artificial diffusivity
      select case(ad_typ)
      case(1)
        ! LW (ISD)
        mu = dt_sub * lambda%im**2 / 2
      case default
        mu = 0
      end select

      ! scale
      mu = mu * min(s_max, s_base + k * s_grow)

      ! SDC quadrature .........................................................

      S = 0
      do j = 0, n_sub
        S = S + dt_step * F(j) * w_nn(m,j)
      end do

      ! u' = u₀ + Sᵏ
      ui = u(m-1) + S

      ! correction .............................................................

      a_inv = ONE / (ONE - dt_sub * (lambda%re - mu))

      uj = a_inv * ( ui + dt_sub * ( F_ex_new(m-1)       &
                                   - F_ex(m-1)           &
                                   - F_im(m) + mu * u0 ) )
      do j = 2, this%n_stage
        uj = a_inv * ( ui + dt_sub * ( i * lambda%im * uj  &
                                     - F_ex(m)             &
                                     - F_im(m) + mu * u0 ) )
      end do
      u(m) = uj

      ! update RHS .............................................................

      call this % CorrectorRHS(lambda, dt_sub, u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStep

  !=============================================================================

end module DQ__SDC__Method__ISD
