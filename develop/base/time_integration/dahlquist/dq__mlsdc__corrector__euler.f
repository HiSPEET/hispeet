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

!> summary:  Euler corrector for Dahlquist MLSDC
!> author:   Joerg Stiller
!> date:     2024/07/15
!===============================================================================

module DQ__MLSDC__Corrector__Euler

  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT

  use Kind_Parameters, only: RNP
  use Constants,       only: ONE, ZERO
  use DQ__Time_Integrator
  use DQ__SDC__Method
  use DQ__SDC__Method__Euler

  implicit none
  private

  public :: DQ_MLSDC_Corrector_Euler
  public :: DQ_MLSDC_Corrector_Options_Euler

  !-----------------------------------------------------------------------------
  !> IMEX Euler SDC ...

  type, extends(DQ_SDC_Method_Euler) :: DQ_MLSDC_Corrector_Euler
  contains
    procedure :: Init_DQ_MLSDC_Corrector_Euler
    procedure :: Show => Show_DQ_MLSDC_Corrector_Euler
    procedure :: CorrectorRHS
    procedure :: CorrectorStepMLSDC
  end type DQ_MLSDC_Corrector_Euler

  ! overloading the constructor
  interface DQ_MLSDC_Corrector_Euler
    module procedure New_DQ_MLSDC_Corrector_Euler
  end interface

  !-----------------------------------------------------------------------------
  !> Type for providing SDC-Euler options

  type, extends(DQ_SDC_Options_Euler) :: DQ_MLSDC_Corrector_Options_Euler
  end type DQ_MLSDC_Corrector_Options_Euler

contains

  !=============================================================================
  ! MLSDC Euler corrector: type-bound procedures

  !-----------------------------------------------------------------------------
  !> Constructor for MLSDC Euler correctors

  function New_DQ_MLSDC_Corrector_Euler(pre_opt, sdc_opt) result(this)
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_MLSDC_Corrector_Options_Euler), intent(in) :: sdc_opt
    type(DQ_MLSDC_Corrector_Euler) :: this

    call Init_DQ_MLSDC_Corrector_Euler(this, pre_opt, sdc_opt)

  end function New_DQ_MLSDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Initialization of an MLSDC Euler corrector

  subroutine Init_DQ_MLSDC_Corrector_Euler(this, pre_opt, sdc_opt)
    class(DQ_MLSDC_Corrector_Euler),       intent(inout) :: this
    class(DQ_TimeIntegrator_Options), intent(in) :: pre_opt !< predictor options
    class(DQ_MLSDC_Corrector_Options_Euler), intent(in) :: sdc_opt

    ! intialize parent type
    call this % Init_DQ_SDC_Method(pre_opt, sdc_opt)

    this % corrector_name = 'IMEX Euler method'

  end subroutine Init_DQ_MLSDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Output of SDC_Corrector_Euler settings

  subroutine Show_DQ_MLSDC_Corrector_Euler(this, unit)
    class(DQ_MLSDC_Corrector_Euler), intent(in) :: this
    integer, optional, intent(in) :: unit  !< output unit

    integer :: io

    if (present(unit)) then
      io = unit
    else
      io = OUTPUT_UNIT
    end if

    call this % Show_DQ_SDC_Method(unit)

    write(io,'(2X,A,T15,G0)') 'name:', this % corrector_name

  end subroutine Show_DQ_MLSDC_Corrector_Euler

  !-----------------------------------------------------------------------------
  !> Computes F_ex and F_im as defined in the corrector

  elemental subroutine CorrectorRHS(this, lambda, dt, u, F_ex, F_im)
    class(DQ_MLSDC_Corrector_Euler), intent(in) :: this
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

  subroutine CorrectorStepMLSDC( this, lambda, m, k, t, u, F, F_ex, F_im &
                          , F_ex_new, F_im_new, G, incremental       )

    class(DQ_MLSDC_Corrector_Euler), intent(inout) :: this
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
    logical,      optional, intent(in) :: incremental !< use incremental form

    ! auxiliary variables .....................................................

    complex(RNP), parameter :: i = (ZERO, ONE)

    real(RNP), allocatable  :: dt(:)
    complex(RNP) :: S, u1, u2
    real(RNP)    :: delta
    integer      :: j
    logical      :: use_incremental

    ! initialization ..........................................................

    use_incremental = .true.
    if (present(incremental)) use_incremental = incremental

    associate( n_sub => this % n_sub &
             , w_nn  => this % w_nn  &
             , w_0n  => this % w_0n  )

      allocate(dt(1:n_sub))
      do j = 1, n_sub
        dt(j) = t(j) - t(j-1)
      end do

      ! SDC quadrature ........................................................

      delta = t(n_sub) - t(0)
      S = 0

      if (use_incremental) then
        do j = 0, n_sub
          S = S + delta * F(j) * w_nn(m,j)
        end do
        ! u' = u₀ + Sᵏ
        u1 = u(m-1) + S
      else
        do j = 0, n_sub
          S = S + delta * F(j) * w_0n(m,j)
        end do
        ! u' = u₀ + Sᵏ
        u1 = u(0) + S
      end if


    end associate

    ! correction ..............................................................

    if (use_incremental) then

      if (present(G)) then
        select case (this % impl)
        case(0) ! explicit
          u2 = u1 + dt(m) * (F_ex_new(m-1) - F_ex(m-1)) + G(m)
        case(1) ! implicit
          u2 = (u1 - dt(m) * F_im(m) + G(m)) / (ONE - dt(m) * lambda)
        case(2) ! IMEX
          u1 = u1 + dt(m) * (F_ex_new(m-1) - F_ex(m-1) - F_im(m)) + G(m)
          u2 = u1 / (ONE - dt(m) * lambda%re)
        case(3) ! IMEX partitioned
          u2 = u1 + dt(m) * ( F_ex_new(m-1) - F_ex(m-1) &
                            + F_im_new(m-1) - F_im(m-1) )
          u2 = u1 + dt(m) * (i * lambda%im * u2 - F_ex(m) - F_im(m)) + G(m)
          u2 = u2 / (ONE - dt(m) * lambda%re)
        end select
      else
        select case (this % impl)
        case(0) ! explicit
          u2 = u1 + dt(m) * (F_ex_new(m-1) - F_ex(m-1))
        case(1) ! implicit
          u2 = (u1 - dt(m) * F_im(m)) / (ONE - dt(m) * lambda)
        case(2) ! IMEX
          u1 = u1 + dt(m) * (F_ex_new(m-1) - F_ex(m-1) - F_im(m))
          u2 = u1 / (ONE - dt(m) * lambda%re)
        case(3) ! IMEX partitioned
          u2 = u1 + dt(m) * ( F_ex_new(m-1) - F_ex(m-1) &
                            + F_im_new(m-1) - F_im(m-1) )
          u2 = u1 + dt(m) * (i * lambda%im * u2 - F_ex(m) - F_im(m))
          u2 = u2 / (ONE - dt(m) * lambda%re)
        end select
      end if

    else

      if (present(G)) then
        select case (this % impl)
        case(0) ! explicit
          u2 = u1 + sum(dt(1:m) * (F_ex_new(0:m-1) - F_ex(0:m-1))) + G(m)
        case(1) ! implicit
          u2 = u1 + lambda * sum(dt(1:m-1) * u(1:m-1)) &
                  - sum(dt(1:m) * F_im(1:m)) + G(m)
          u2 = u2 / (ONE - dt(m) * lambda)
        case(2) ! IMEX
          u2 = u1 + lambda%re * sum(dt(1:m-1)*u(1:m-1)) + G(m)                 &
                  + sum(dt(1:m) * (F_ex_new(0:m-1) - F_ex(0:m-1) - F_im(1:m)))
          u2 = u2 / (ONE - dt(m) * lambda%re)
        case(3) ! IMEX partitioned  TODO
          u2 = u1 + dt(m) * (F_ex_new(0) - F_ex(0) + F_im_new(0) - F_im(0))
          u2 = u1 + dt(m) * (i * lambda%im * u2 - F_ex(m) - F_im(m))
          u2 = u2 / (ONE - dt(m) * lambda%re)
        end select
      else
        select case (this % impl)
        case(0) ! explicit
          u2 = u1 + sum(dt(1:m) * (F_ex_new(0:m-1) - F_ex(0:m-1)))
        case(1) ! implicit
          u2 = u1 + lambda * sum(dt(1:m-1)*u(1:m-1)) - sum(dt(1:m)*F_im(1:m))
          u2 = u2 / (ONE - dt(m) * lambda)
        case(2) ! IMEX
          u2 = u1 + lambda%re * sum(dt(1:m-1)*u(1:m-1))                       &
                  + sum(dt(1:m) * (F_ex_new(0:m-1) - F_ex(0:m-1) - F_im(1:m)))
          u2 = u2 / (ONE - dt(m) * lambda%re)
        case(3) ! IMEX partitioned
          u2 = u1 + dt(m) * (F_ex_new(0) - F_ex(0) + F_im_new(0) - F_im(0))
          u2 = u1 + dt(m) * (i * lambda%im * u2 - F_ex(m) - F_im(m))
          u2 = u2 / (ONE - dt(m) * lambda%re)
        end select
      end if

    end if

    u(m) = u2

    ! update RHS
    call this % CorrectorRHS(lambda, dt(m), u(m), F_ex_new(m), F_im_new(m))

  end subroutine CorrectorStepMLSDC

  !=============================================================================

end module DQ__MLSDC__Corrector__Euler
