!> summary:  Evaluation of multilevel flow characteristics
!> author:   Joerg Stiller
!> date:     2025/07/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
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
  contains
    procedure :: Evaluate
    procedure :: PrintHeader
    procedure :: PrintValues
  end type ML_INS_FlowCharacteristics_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of multilevel incompressible flow characteristics

  subroutine Evaluate(this, ml_ins, t, u, dt, volume, leaf)
    class(ML_INS_FlowCharacteristics_3D), intent(inout) :: this
    class(ML_INS_Operator_3D), intent(in) :: ml_ins
    real(RNP),                 intent(in) :: t
    class(ML_MeshVariable_3D), intent(in) :: u
    real(RNP), optional,       intent(in) :: dt
    real(RNP), optional,       intent(in) :: volume
    logical,   optional,       intent(in) :: leaf !< constrain to leaves [F]

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

    ! evaluation ...............................................................

    do l = 1, l_top
      call this % flow_char(l) % Evaluate( ins_op = ml_ins % ins_op(l) &
                                         , t      = t                  &
                                         , u      = u % level(l) % val &
                                         , dt     = dt                 &
                                         , volume = volume             &
                                         , leaf   = leaf               )
    end do

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Print header for multilevel incompressible flow characteristics

  subroutine PrintHeader(this, tag)
    class(ML_INS_FlowCharacteristics_3D), intent(in) :: this
    character(len=*), optional, intent(in) :: tag !< tag placed at end of line

    if (this % proc == 0) then

      !$omp master
      write(*,'(A,2X)'   ,advance='NO') '#'
      write(*,'(2X,A,7X)',advance='NO') 't'
      if (this % flow_char(1) % dt >= 0) then
        write(*,'(2X,A,6X)',advance='NO') 'dt'
      end if
      write(*,'(1X,A,4X)',advance='NO') 'dx_min'
      write(*,'(1X,A,4X)',advance='NO') 'dx_max'
      write(*,'(1X,A,5X)',advance='NO') 'v_max'
      write(*,'(1X,A,5X)',advance='NO') 'e_kin'
      write(*,'(1X,A,5X)',advance='NO') 'div_v'

      if (this % has_errors) then
        write(*,'(1X,A,5X)',advance='NO') 'err_v'
        write(*,'(1X,A,5X)',advance='NO') 'err_p'
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

    if (this % proc == 0) then
      associate(flow_char => this % flow_char)

        !$omp master
        write(*,'(ES10.3,1X)',advance='NO') flow_char(1) % t
        if (flow_char(1) % dt >= 0) then
          write(*,'(ES10.3,1X)',advance='NO') flow_char(1) % dt
        end if
        write(*,'(ES10.3,1X)',advance='NO') minval(flow_char % dx_min)
        write(*,'(ES10.3,1X)',advance='NO') maxval(flow_char % dx_max)
        write(*,'(ES10.3,1X)',advance='NO') maxval(flow_char % v_max)
        write(*,'(ES10.3,1X)',advance='NO') sum(flow_char % e_kin)
        write(*,'(ES10.3,1X)',advance='NO') sqrt(sum(flow_char % div_v ** 2))

        if (this % has_errors) then
          write(*,'(ES10.3,1X)',advance='NO') sqrt(sum(flow_char % err_v ** 2))
          write(*,'(ES10.3,1X)',advance='NO') sqrt(sum(flow_char % err_p ** 2))
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
