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

!> summary:  ISD corrector for Dahlquist MLSDC
!> author:   Joerg Stiller
!> date:     2024/07/15
!===============================================================================

module DQ__MLSDC__Corrector__ISD

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, HALF, ONE
  use DQ__Time_Integrator
  use DQ__SDC__Method
  use DQ__SDC__Method__ISD

  implicit none
  private

  public :: DQ_MLSDC_Corrector_ISD
  public :: DQ_MLSDC_Corrector_Options_ISD

  !-----------------------------------------------------------------------------
  !> IMEX ISD SDC ...

  type, extends(DQ_SDC_Method_ISD) :: DQ_MLSDC_Corrector_ISD
    real(RNP) :: c_isd   !< artificial diffusivity factor
    integer   :: scaling !< artificial diffusivity scaling: 0/1/2 none/k/double
  contains
    procedure :: Init_DQ_MLSDC_Corrector_ISD
    procedure :: Show => Show_DQ_MLSDC_Corrector_ISD
    procedure :: CorrectorRHS
    procedure :: CorrectorStepMLSDC
  end type DQ_MLSDC_Corrector_ISD

  ! overloading the constructor
  interface DQ_MLSDC_Corrector_ISD
    module procedure New_DQ_MLSDC_Corrector_ISD
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-ISD options

  type, extends(DQ_SDC_Options_ISD) :: DQ_MLSDC_Corrector_Options_ISD
    real(RNP) :: c_isd   = 1  !< AD amplitude factor
    integer   :: scaling = 0  !< AD scaling: 0/1/2 none/k/double
  end type DQ_MLSDC_Corrector_Options_ISD

contains

  !=============================================================================
  ! MLSDC ISD corrector: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for MLSDC ISD correctors

  function New_DQ_MLSDC_Corrector_ISD(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_MLSDC_Corrector_Options_ISD), intent(in) :: sdc_opt
    type(DQ_MLSDC_Corrector_ISD) :: this

    call Init_DQ_MLSDC_Corrector_ISD(this, pre_opt, sdc_opt)

  end function New_DQ_MLSDC_Corrector_ISD

  !-----------------------------------------------------------------------------
  !> Initialization of an MLSDC ISD corrector

  subroutine Init_DQ_MLSDC_Corrector_ISD(this, pre_opt, sdc_opt)
    class(DQ_MLSDC_Corrector_ISD),         intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_MLSDC_Corrector_Options_ISD), intent(in) :: sdc_opt

    ! intialize parent type
    call this % Init_DQ_SDC_Method_ISD(pre_opt, sdc_opt)

    this % n_stage = sdc_opt % n_stage
    this % c_isd   = sdc_opt % c_isd
    this % scaling = sdc_opt % scaling

    write(this%corrector_name,'(A,G0,A)') &
        'ISD method of order 1 with ', this%n_stage, ' stage(s)'

  end subroutine Init_DQ_MLSDC_Corrector_ISD

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_ISD settings

  subroutine Show_DQ_MLSDC_Corrector_ISD(this, unit)
    class(DQ_MLSDC_Corrector_ISD), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:',  this % corrector_name

  end subroutine Show_DQ_MLSDC_Corrector_ISD

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_MLSDC_Corrector_ISD), intent(in) :: this
    complex(RNP), intent(in)  :: lambda !< λ
    real   (RNP), intent(in)  :: dt     !< step size, used with ISD only
    complex(RNP), intent(in)  :: u      !< u
    complex(RNP), intent(out) :: F_ex   !< explicit RHS for corrector
    complex(RNP), intent(out) :: F_im   !< implicit RHS for corrector

    complex(RNP), parameter :: i = (ZERO, ONE)

    F_im = (lambda % re - this%c_isd * HALF * dt * lambda%im ** 2) * u
    F_ex = i * lambda % im * u

    if (this % impl == 0) return ! just to avoid compiler warning !

  end subroutine CorrectorRHS

  !-----------------------------------------------------------------------------
  !> Execution of a single correction step

  subroutine CorrectorStepMLSDC( this, lambda, m, k, t, u, F, F_ex, F_im &
                          , F_ex_new, F_im_new, G, incremental       )

    class(DQ_MLSDC_Corrector_ISD), intent(inout) :: this
    complex(RNP), intent(in)    :: lambda       !< λ
    integer     , intent(in)    :: m            !< current SDC interval index
    integer     , intent(in)    :: k            !< current corrector sweep
    real   (RNP), intent(in)    :: t(0:)        !< SDC time nodes
    complex(RNP), intent(inout) :: u(0:)        !< SDC solution
    complex(RNP), intent(in)    :: F(0:)        !< Fᵏ
    complex(RNP), intent(in)    :: F_ex(0:)     !< F_exᵏ
    complex(RNP), intent(in)    :: F_im(0:)     !< F_imᵏ
    complex(RNP), intent(inout) :: F_ex_new(0:) !< F_exᵏ⁺¹(0:m-1) → F_exᵏ⁺¹(0:m)
    complex(RNP), intent(inout) :: F_im_new(0:) !< F_imᵏ⁺¹(0:m-1) → F_imᵏ⁺¹(0:m)
    complex(RNP), optional, intent(in) :: G(0:) !< FAS defect correction term
    logical, optional, intent(in) :: incremental !< use incremental form

    ! auxiliary variables .....................................................

    complex(RNP), parameter :: i = (ZERO, ONE)
    complex(RNP)            :: ui, uj, S
    real(RNP), allocatable  :: a_inv(:), c_im(:)
    real(RNP), allocatable  :: dt_isd(:), dt_step(:)
    real(RNP)               :: delta
    integer                 :: j
    logical                 :: use_incremental

    use_incremental = .true.
    if (present(incremental)) use_incremental = incremental

    associate( n_sub => this % n_sub &
             , w_nn  => this % w_nn  &
             , w_0n  => this % w_0n  &
             , c_isd => this % c_isd )

      ! initialization .........................................................

      allocate(a_inv  (1:n_sub))
      allocate(c_im   (1:n_sub))
      allocate(dt_step(1:n_sub))
      allocate(dt_isd (1:n_sub))

      do j = 1, n_sub
        dt_step(j) = t(j) - t(j-1)
        dt_isd(j)  = c_isd * HALF * dt_step(j)
      end do
      delta = t(n_sub) - t(0)

      do j = 1, n_sub
        c_im(j) = 1
        if (this % scaling > 0) then
          c_im(j) = lambda % re - dt_isd(j) * lambda%im**2
          if (abs(c_im(j)) > 1000 * tiny(ONE)) then
            select case(this % scaling)
            case(1)
              c_im(j) = (lambda % re - k * dt_isd(j) * lambda%im**2) / c_im(j)
            case(2)
              c_im(j) = ( lambda % re &
                        - 2**(k-1) * dt_isd(j) * lambda % im**2 ) / c_im(j)
            end select
          end if
        end if

        a_inv(j) = ONE / (ONE - dt_step(j) * c_im(j)     &
                 * (lambda%re - dt_isd(j) * lambda%im**2))

      end do

      ! SDC quadrature .........................................................

      S = 0

      if (use_incremental) then
        do j = 0, n_sub
          S = S + delta * F(j) * w_nn(m,j)
        end do
        ! u' = u₀ + Sᵏ
        ui = u(m-1) + S
      else
        do j = 0, n_sub
          S = S + delta * F(j) * w_0n(m,j)
        end do
        ! u' = u₀ + Sᵏ
        ui = u(0) + S
      end if

      ! correction .............................................................

      if (use_incremental) then

        if (present(G)) then
          uj = ( ui + dt_step(m) * ( F_ex_new(m-1) - F_ex(m-1) &
               - c_im(m) * F_im(m) ) + G(m) ) * a_inv(m)
          do j = 2, this % n_stage
            uj = ( ui + dt_step(m) * ( i * lambda % im * uj - F_ex(m) &
                 - c_im(m) * F_im(m) ) + G(m) ) * a_inv(m)
          end do
        else
          uj = ( ui + dt_step(m) * ( F_ex_new(m-1) - F_ex(m-1) &
               - c_im(m) * F_im(m) ) ) * a_inv(m)
          do j = 2, this % n_stage
            uj = ( ui + dt_step(m) * ( i * lambda % im * uj - F_ex(m) &
                 - c_im(m) * F_im(m) ) ) * a_inv(m)
          end do
        end if

      else

        uj = ui + sum((-ONE/a_inv(1:m-1) + ONE) * u(1:m-1))         &
                + sum(dt_step(1:m) * (F_ex_new(0:m-1) - F_ex(0:m-1) &
                - c_im(1:m) * F_im(1:m)))
        if (present(G)) then
          uj = (uj + G(m)) * a_inv(m)
        else
          uj = uj * a_inv(m)
        end if
        if (present(G)) then
          u(m) = uj
          do j = 2, this % n_stage
            uj = ui + sum((-ONE/a_inv(1:m-1) + ONE) * u(1:m-1))  &
                    + sum(dt_step(1:m) * (i * lambda%im * u(1:m) &
                    - F_ex(1:m) - c_im(1:m) * F_im(1:m))) + G(m)
            uj = uj * a_inv(m)
          end do
        else
          u(m) = uj
          do j = 2, this % n_stage
            uj = ui + sum((-ONE/a_inv(1:m-1) + ONE) * u(1:m-1))  &
                    + sum(dt_step(1:m) * (i * lambda%im * u(1:m) &
                    - F_ex(1:m) - c_im(1:m) * F_im(1:m)))
            uj = uj * a_inv(m)
          end do
        end if

      end if

      u(m) = uj

      ! update RHS .............................................................

      call this % CorrectorRHS( &
          lambda, dt_step(m), u(m), F_ex_new(m), F_im_new(m))

    end associate

  end subroutine CorrectorStepMLSDC

  !=============================================================================

end module DQ__MLSDC__Corrector__ISD
