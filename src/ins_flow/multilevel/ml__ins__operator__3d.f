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

!> summary:  3D multilevel DG operator for incompressible Navier-Stokes problems
!> author:   Joerg Stiller
!> date:     2025/03/01
!===============================================================================

module ML__INS__Operator__3D
  use Kind_Parameters
  use XMPI
  use Execution_Control
  use Volume_Integrals__3D
  use INS__Problem__3D
  use INS__Operator__3D
  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D

  implicit none
  private

  public :: ML_INS_Operator_3D
  public :: ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM INS operator options

  type ML_INS_OperatorOptions_3D
    logical :: mixed = .true. !< T/F: use mixed/equal order for (v,p)
    type(INS_OperatorOptions_3D) :: ins_op !< INS operator options
  contains
    procedure :: Bcast => Bcast_ML_INS_OperatorOptions_3D
  end type ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> Multilevel DG-SEM operators for incompressible Navier-Stokes problems

  type ML_INS_Operator_3D
    class(INS_Problem_3D),         pointer     :: problem
    class(ML_Mesh_3D),             pointer     :: ml_mesh
    class(ML_MeshOperators_3D),    pointer     :: ml_op_u
    type (ML_MeshOperators_3D),    allocatable :: ml_op_p
    type (INS_Operator_3D),        allocatable :: ins_op(:)
  contains
    procedure :: Init_ML_INS_Operator_3D
    procedure :: CalibratePressure
  end type ML_INS_Operator_3D

  ! constructor interface
  interface ML_INS_Operator_3D
    procedure New_ML_INS_Operator_3D
  end interface

contains

  !-----------------------------------------------------------------------------
  !> Constructor of ML_INS_Operator_3D

  function New_ML_INS_Operator_3D(opt, problem, ml_mesh, ml_op) result(this)
    class(ML_INS_OperatorOptions_3D), intent(in) :: opt
    class(INS_Problem_3D),            intent(in) :: problem
    class(ML_Mesh_3D),                intent(in) :: ml_mesh
    class(ML_MeshOperators_3D),       intent(in) :: ml_op
    type(ML_INS_Operator_3D) :: this

    call Init_ML_INS_Operator_3D(this, opt, problem, ml_mesh, ml_op)

  end function New_ML_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Initialization of ML_INS_Operator_3D

  subroutine Init_ML_INS_Operator_3D(this, opt, problem, ml_mesh, ml_op)
    class(ML_INS_Operator_3D),          intent(inout) :: this
    class(ML_INS_OperatorOptions_3D),   intent(in)    :: opt
    class(INS_Problem_3D),      target, intent(in)    :: problem
    class(ML_Mesh_3D),          target, intent(in)    :: ml_mesh
    class(ML_MeshOperators_3D), target, intent(in)    :: ml_op

    integer :: l, l_top

    l_top = size(ml_mesh%mesh)

    this % problem => problem
    this % ml_mesh => ml_mesh
    this % ml_op_u => ml_op

    if (opt%mixed) then
      associate(po_p => max(1, ml_op%sem%std_op%po - 1))
        this % ml_op_p = ML_MeshOperators_3D(ml_mesh, po_p)
      end associate
    end if

    allocate(this%ins_op(l_top))
    do l = 1, l_top
      associate(ins_op => this % ins_op(l), sem_u => ml_op % sem(l))
        if(opt%mixed) then
          associate(sem_p => this % ml_op_p % sem(l))
            ins_op = INS_Operator_3D(opt%ins_op, problem, sem_u, sem_p, level=l)
          end associate
        else
          ins_op = INS_Operator_3D(opt%ins_op, problem, sem_u, level=l)
        end if
      end associate
    end do

  end subroutine Init_ML_INS_Operator_3D

  !-----------------------------------------------------------------------------
  !> Remove the mean pressure unless Dirichlet conditions apply
  !>
  !> The multilevel variable `u` must contain the either the flow variables
  !> `[ v(1:3), p, ... ]` or only the pressure `p`.

  subroutine CalibratePressure(this, u)
    class(ML_INS_Operator_3D), intent(in)    :: this !< ML INS operator
    class(ML_MeshVariable_3D), intent(inout) :: u    !< flow variables or p

    real(RNP) :: p_mean, p_mean_l, volume
    integer   :: c, l, l_top

    if (any(this % problem % bc_p == 'D')) return

    ! identify pressure component
    select case(size(u%name))
    case(1)
      c = 1
    case(4:)
      c = 4
    case default
      call Error( 'CalibratePressure'                                   &
                , 'u must contain either the flow variables or p alone' &
                , 'ML__INS__Operator__3D'                               )
    end select

    call this % ml_op_u % Get_Volume(volume)

    l_top = size(u % level)

    p_mean = 0
    do l = 1, l_top
      associate(p_l => u % level(l) % val(:,:,:,:,c))
        call GetVolumeIntegral(this%ins_op(l)%sem_u, p_l, p_mean_l, leaf=.true.)
        p_mean = p_mean + p_mean_l
      end associate
    end do
    p_mean = p_mean / volume

    do l = 1, l_top
      associate(p_l => u % level(l) % val(:,:,:,:,c))
        p_l = p_l - p_mean
      end associate
    end do

  end subroutine CalibratePressure

  !=============================================================================
  ! Type-bound procedures of ML_INS_OperatorOptions_3D

  !-----------------------------------------------------------------------------
  !> MPI_Bcast for objects of type ML_INS_OperatorOptions_3D

   subroutine Bcast_ML_INS_OperatorOptions_3D(this, root, comm)
    class(ML_INS_OperatorOptions_3D), intent(inout) :: this
    integer,        intent(in) :: root !< rank of broadcast root
    type(MPI_Comm), intent(in) :: comm !< MPI communicator

    call XMPI_Bcast(this % mixed, root, comm)
    call this % ins_op % Bcast(root, comm)

  end subroutine Bcast_ML_INS_OperatorOptions_3D

  !=============================================================================

end module ML__INS__Operator__3D
