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

!> summary:  Evaluation of multilevel flow characteristics
!> author:   Joerg Stiller
!> date:     2025/07/09
!===============================================================================

module ML__INS__Flow_Characteristics__3D
  use Kind_Parameters
  use INS__Flow_Characteristics__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D

  implicit none
  private

  public :: ML_INS_FlowCharacteristics_3D

  !-----------------------------------------------------------------------------
  !> Multilevel incompressible flow statistics

  type ML_INS_FlowCharacteristics_3D
    type(INS_FlowCharacteristics_3D), allocatable :: flow_char(:)
    integer, private :: proc = -1
    logical, private :: has_errors = .false.
    logical, private :: has_dissipation = .false.
  contains
    procedure :: Evaluate
    procedure :: PrintHeader
    procedure :: PrintValues
  end type ML_INS_FlowCharacteristics_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of multilevel incompressible flow characteristics

  subroutine Evaluate(this, ml_ins, t, mu, nu, u, dt, volume, diss, leaf)
    class(ML_INS_FlowCharacteristics_3D), intent(inout) :: this
    class(ML_INS_Operator_3D), intent(in) :: ml_ins
    real(RNP),                 intent(in) :: t
    class(ML_MeshVariable_3D), intent(in) :: mu
    class(ML_MeshVariable_3D), intent(in) :: nu
    class(ML_MeshVariable_3D), intent(in) :: u
    real(RNP), optional,       intent(in) :: dt
    real(RNP), optional,       intent(in) :: volume
    logical,   optional,       intent(in) :: diss !< compute disspation [F]
    logical,   optional,       intent(in) :: leaf !< only leaf elements [F]

    integer :: l, l_top

    ! initialization ...........................................................

    l_top = size(ml_ins % ins_op)

    if (allocated(this % flow_char)) then
      if (size(this % flow_char) /= l_top) then
        deallocate(this % flow_char)
      end if
    end if

    if (.not. allocated(this % flow_char)) then
      allocate(this % flow_char(l_top))
    end if

    this % proc       = ml_ins % ins_op(1) % mesh % proc
    this % has_errors = ml_ins % problem % HasExactSolution()

    if (present(diss)) then
      this % has_dissipation = diss
    else
      this % has_dissipation = .false.
    end if

    ! evaluation ...............................................................

    do l = 1, l_top
      call this % flow_char(l) &
                      % Evaluate( ins_op = ml_ins % ins_op(l)             &
                                , t      = t                              &
                                , mu     = mu % level(l) % val(:,:,:,:,1) &
                                , nu     = nu % level(l) % val(:,:,:,:,1) &
                                , u      = u  % level(l) % val            &
                                , dt     = dt                             &
                                , volume = volume                         &
                                , diss   = diss                           &
                                , leaf   = leaf                           )
    end do

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print header for multilevel incompressible flow characteristics

  subroutine PrintHeader(this, tag)
    class(ML_INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    character(len=*), parameter :: fmt = '(1X,A6,4X)' !! short
  ! character(len=*), parameter :: fmt = '(2X,A6,5X)' !! long

    if (this % proc == 0) then

      !$omp master
      write(*,'(A)',advance='NO') '#'
      write(*,fmt,advance='NO') '  t   '
      if (this % flow_char(1) % dt >= 0) then
        write(*,fmt,advance='NO') '  dt  '
      end if
      write(*,fmt,advance='NO') 'dx_min'
      write(*,fmt,advance='NO') 'dx_max'
      write(*,fmt,advance='NO') 'v_max '
      write(*,fmt,advance='NO') 'e_kin '
      write(*,fmt,advance='NO') 'div_v '

      if (this % has_errors) then
        write(*,fmt,advance='NO') 'err_v '
        write(*,fmt,advance='NO') 'err_p '
      end if

      if (this % has_dissipation) then
        write(*,fmt,advance='NO') 'phi_c '
        write(*,fmt,advance='NO') 'phi_d '
        write(*,fmt,advance='NO') 'phi_ds'
        write(*,fmt,advance='NO') 'phi_d0'
        write(*,fmt,advance='NO') 'phi_s '
      end if

      if (present(tag)) then
        write(*,'(A)') tag
      else
        write(*,*)
      end if
      !$omp end master

    end if

  end subroutine PrintHeader

  !-----------------------------------------------------------------------------
  !> Print incompressible flow characteristics

  subroutine PrintValues(this, tag)
    class(ML_INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    character(len=*), parameter :: fmt = '(ES10.3,1X)' !! short
  ! character(len=*), parameter :: fmt = '(ES12.5,1X)' !! long

    if (this % proc == 0) then
      associate(flow_char => this % flow_char)

        !$omp master
        write(*,fmt,advance='NO') flow_char(1) % t
        if (flow_char(1) % dt >= 0) then
          write(*,fmt,advance='NO') flow_char(1) % dt
        end if
        write(*,fmt,advance='NO') minval(flow_char % dx_min)
        write(*,fmt,advance='NO') maxval(flow_char % dx_max)
        write(*,fmt,advance='NO') maxval(flow_char % v_max)
        write(*,fmt,advance='NO') sum(flow_char % e_kin)
        write(*,fmt,advance='NO') sqrt(sum(flow_char % div_v ** 2))

        if (this % has_errors) then
          write(*,fmt,advance='NO') sqrt(sum(flow_char % err_v ** 2))
          write(*,fmt,advance='NO') sqrt(sum(flow_char % err_p ** 2))
        end if

        if (this % has_dissipation) then
          write(*,fmt,advance='NO') sum(flow_char % phi_c )
          write(*,fmt,advance='NO') sum(flow_char % phi_d )
          write(*,fmt,advance='NO') sum(flow_char % phi_ds)
          write(*,fmt,advance='NO') sum(flow_char % phi_d0)
          write(*,fmt,advance='NO') sum(flow_char % phi_s )
        end if

        if (present(tag)) then
          write(*,'(X,A)') tag
        else
          write(*,*)
        end if
       !$omp end master

      end associate
    end if

  end subroutine PrintValues

  !=============================================================================

end module ML__INS__Flow_Characteristics__3D
