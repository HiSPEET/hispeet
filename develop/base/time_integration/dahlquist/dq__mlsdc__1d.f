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

!> summary:  Multilevel SDC method for Dahlquist
!> author:   Erik Pfister
!> date:     2024/03/22
!===============================================================================

module DQ__MLSDC__1D
  use Kind_Parameters
  use Standard_Element_Operators__1D
  use HP__Refinement_Operator__1D
  use HP__Coarsening_Operator__1D
  use DQ__Time_Integrator
  use DQ__Time_Integrator__ISD
  use DQ__SDC__Method
  use DQ__MLSDC__Corrector__ISD
  use DQ__MLSDC__Level__1D

  implicit none
  private

  !-----------------------------------------------------------------------------
  !> MLSDC for 1D dahlquist

  type, public :: DQ_MLSDC_1D
    type(DQ_MLSDC_Level_1D), allocatable :: level(:)
  end type DQ_MLSDC_1D

  ! constructor
  interface DQ_MLSDC_1D
    procedure New_DQ_MLSDC_1D
  end interface

  !-----------------------------------------------------------------------------
  !> Options for initializing DQ_MLSDC_1D objects

  type, public :: DQ_MLSDC_Options_1D

    integer :: n_level = -1  !< number of space-time levels

    ! level parameters
    integer, allocatable :: p_time(:)  !< polynomial degree of time step
    integer, allocatable :: n_time(:)  !< number of time steps in one slice
    integer, allocatable :: n_stage_p(:)
    integer, allocatable :: n_stage_c(:)

    character :: projection_method = 'I' !< fine-to-coarse projection method:
                                         !! 'P'  L² projection,
                                         !! 'I'  interpolation

    integer :: projection_smoothing = 0  !< discontinuities in 2:1 projection:
                                         !!  0   no smoothing

  end type 

  ! constructor
  interface DQ_MLSDC_Options_1D
    procedure New_DQ_MLSDC_Options_1D
  end interface

contains

  !=============================================================================
  ! TBP for DQ_MLSDC_1D

  !-----------------------------------------------------------------------------
  !> New DQ_MLSDC_1D object

  function New_DQ_MLSDC_1D(opt, opt_pre, opt_sdc) result(this)
    class(DQ_MLSDC_Options_1D),       intent(in) :: opt
    class(DQ_TimeIntegrator_Options), intent(in) :: opt_pre
    class(DQ_SDC_Options),            intent(in) :: opt_sdc
    type(DQ_MLSDC_1D) :: this

    call Init_DQ_MLSDC_1D(this, opt, opt_pre, opt_sdc)

  end function New_DQ_MLSDC_1D

  !-----------------------------------------------------------------------------
  !> Initialization of a DQ_MLSDC_1D object

  subroutine Init_DQ_MLSDC_1D(this, opt, opt_pre, opt_sdc)
    class(DQ_MLSDC_1D),               intent(inout) :: this
    class(DQ_MLSDC_Options_1D),       intent(in)    :: opt
    class(DQ_TimeIntegrator_Options), intent(in)    :: opt_pre
    class(DQ_SDC_Options),            intent(in)    :: opt_sdc

    type(DQ_MLSDC_Level_Options_1D) :: opt_level
    integer :: l

    ! preliminaries ............................................................

    if (allocated(this % level)) then
      deallocate(this % level)
    end if

    allocate(this % level(opt % n_level))

    ! consistency check ........................................................

    ! TBD

    ! generate levels ..........................................................

    do l = 1, opt % n_level

      ! dimensions . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      opt_level % p_time  = opt % p_time (l)
      opt_level % n_time  = opt % n_time (l)

      ! time options . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      opt_level % pre_opt = opt_pre
      opt_level % sdc_opt = opt_sdc

      ! number of collocation points = polynomial degree + 1
      opt_level % sdc_opt % n_col = opt % p_time(l) + 1

      ! n_stage as ML variable
      select type(pre_opt => opt_level % pre_opt)
      type is(DQ_TimeIntegrator_Options_ISD)
        pre_opt % n_stage = opt % n_stage_p(l)
      end select

      select type(sdc_opt => opt_level % sdc_opt)
      type is(DQ_MLSDC_Corrector_Options_ISD)
        sdc_opt % n_stage = opt % n_stage_c(l)
      end select

      ! transfer options . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      if (l < opt % n_level) then

        opt_level % iop_cf_t =                              &
            HP_RefinementOptions_1D(                        &
                nodes = opt_sdc % nodes,                    &
                po_c  = opt % p_time(l),                    &
                po_f  = opt % p_time(l+1),                  &
                mode  = opt % n_time(l+1) / opt % n_time(l) )

      else

        ! reset interpolation options
        opt_level % iop_cf_t = HP_RefinementOptions_1D()

      end if

      if (l > 1) then

        opt_level % pop_fc_t =                                   &
            HP_CoarseningOptions_1D(                             &
                nodes     = opt_sdc % nodes,                     &
                po_f      = opt % p_time(l),                     &
                po_c      = opt % p_time(l-1),                   &
                mode      = opt % n_time(l) / opt % n_time(l-1), &
                method    = opt % projection_method              )

      else

        ! reset projection options
        opt_level % pop_fc_t = HP_CoarseningOptions_1D()

      end if

      ! create level  . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

      this % level(l) = DQ_MLSDC_Level_1D(opt_level)

    end do

  end subroutine Init_DQ_MLSDC_1D

  !=============================================================================
  ! TBP for DQ_MLSDC_Options_1D

  !-----------------------------------------------------------------------------
  !> New DQ_MLSDC_1D options

  function New_DQ_MLSDC_Options_1D(n_level) result(opt)
    integer, intent(in) :: n_level !< number of space-time levels
    type(DQ_MLSDC_Options_1D) :: opt

    opt % n_level = n_level

    allocate(opt % p_time (n_level))
    allocate(opt % n_time (n_level))
    allocate(opt % n_stage_p (n_level))
    allocate(opt % n_stage_c (n_level))

  end function New_DQ_MLSDC_Options_1D

  !=============================================================================

end module DQ__MLSDC__1D
