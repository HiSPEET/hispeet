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

!> summary:  Evaluation of multilevel INS boundary fluxes
!> author:   Joerg Stiller
!> date:     2026/08/21
!===============================================================================

module ML__INS__Boundary_Fluxes__3D
  use Kind_Parameters
  use INS__Boundary_Fluxes__3D
  use ML__Mesh_Variable__3D
  use ML__INS__Operator__3D

  implicit none
  private

  public :: ML_INS_BoundaryFluxes_3D

  !-----------------------------------------------------------------------------
  !> Multilevel incompressible flow statistics

  type ML_INS_BoundaryFluxes_3D
    type(INS_BoundaryFluxes_3D), allocatable :: ins_flux(:)
    integer, private :: proc = -1
  contains
    procedure :: Evaluate
    procedure :: PrintValues
  end type ML_INS_BoundaryFluxes_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of multilevel incompressible flow characteristics

  subroutine Evaluate(this, ml_ins, t, mu, nu, u, leaf)
    class(ML_INS_BoundaryFluxes_3D), intent(inout) :: this
    class(ML_INS_Operator_3D), intent(in) :: ml_ins
    real(RNP),                 intent(in) :: t
    class(ML_MeshVariable_3D), intent(in) :: mu
    class(ML_MeshVariable_3D), intent(in) :: nu
    class(ML_MeshVariable_3D), intent(in) :: u
    logical,   optional,       intent(in) :: leaf !< only leaf elements [F]

    integer :: l, l_top

    ! initialization ...........................................................

    l_top = size(ml_ins % ins_op)

    if (allocated(this % ins_flux)) then
      if (size(this % ins_flux) /= l_top) then
        deallocate(this % ins_flux)
      end if
    end if

    if (.not. allocated(this % ins_flux)) then
      allocate(this % ins_flux(l_top))
    end if

    this % proc = ml_ins % ins_op(1) % mesh % proc

    ! evaluation ...............................................................

    do l = 1, l_top
      call this % ins_flux(l) &
                      % Evaluate( ins_op = ml_ins % ins_op(l)             &
                                , t      = t                              &
                                , mu     = mu % level(l) % val(:,:,:,:,1) &
                                , nu     = nu % level(l) % val(:,:,:,:,1) &
                                , u      = u  % level(l) % val            &
                                , leaf   = leaf                           )
    end do

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print multilevel INS boundary fluxes

  subroutine PrintValues(this)
    class(ML_INS_BoundaryFluxes_3D), intent(in) :: this

    integer   :: b, l, l_top, n_bound
    real(RNP) :: m, f_p(3), f_d(3)

    !$omp master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    if (this % proc == 0) then

      associate(ins_flux => this % ins_flux)

        l_top   = size(ins_flux)
        n_bound = size(ins_flux(1) % m)

        write(*,*)
        write(*,'(A,ES12.5)') 'boundary fluxes at t =', ins_flux(1) % t
        write(*,'(A)',advance='NO') '#'
        write(*,'(2X,A1,2X)',advance='NO') 'b'
        write(*,'(1X,A4,8X)',advance='NO') 'm'
        write(*,'(2X,A4,6X)',advance='NO') 'f_p1'
        write(*,'(2X,A4,6X)',advance='NO') 'f_p2'
        write(*,'(2X,A4,6X)',advance='NO') 'f_p3'
        write(*,'(2X,A4,6X)',advance='NO') 'f_d1'
        write(*,'(2X,A4,6X)',advance='NO') 'f_d2'
        write(*,'(2X,A4,6X)'             ) 'f_d3'

        do b = 1, n_bound

          m   = 0
          f_p = 0
          f_d = 0

          do l = 1, l_top
            m   = m   + ins_flux(l) % m(b)
            f_p = f_p + ins_flux(l) % f_p(b,1:3)
            f_d = f_d + ins_flux(l) % f_d(b,1:3)
          end do

          write(*,'(I4,2X,7(ES10.3,2X))') b, m, f_p, f_d

        end do
      end associate
    end if

    !$omp end master !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  end subroutine PrintValues

  !=============================================================================

end module ML__INS__Boundary_Fluxes__3D
