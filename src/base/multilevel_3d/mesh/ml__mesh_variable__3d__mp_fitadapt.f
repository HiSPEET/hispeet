!> summary:  Fit 3D multilevel mesh variable to adapted mesh
!> author:   Joerg Stiller
!> date:     2025/01/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(ML__Mesh_Variable__3D) MP_FitAdapt
  use Parent_To_Child_Interpolation__3D
!### CHECK
use XMPI
!### CHECK END
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fit to adapted mesh

  module subroutine FitAdapt(this, ml_op, x_plan)
    class(ML_MeshVariable_3D),  intent(inout) :: this
    class(ML_MeshOperators_3D), intent(in)    :: ml_op     !< adapted operators
    class(DataExchangePlan_3D), intent(in)    :: x_plan(:) !< reassignment plan

    type(MeshVariable_3D), allocatable, save :: old_level(:)

    integer :: l

    ! preliminaries ............................................................

    ! save old mesh data
    call move_alloc(this % level, old_level)

    ! re-initialize ML mesh variable
    call this % Init(ml_op, size(this%name), this%name)

    ! redistribute root level ..................................................

    call CopyRetainedData(old_level(1), this%level(1), x_plan(1))

    ! interpolate/reassign higher levels .......................................

!!     do l = 2, size(this%level)
!!
!!       call ParentToChildInterpolation_3D( parent = this % level(l-1) % mesh &
!!                                         , child  = this % level(l)   % mesh &
!!                                         , iop    = ml_op % iop_cf_x(l-1)    &
!!                                         , v_p    = this % level(l-1) % val  &
!!                                         , v_c    = this % level(l  ) % val  )
!!
!!       if (l <= size(x_plan)) then
!!         call CopyRetainedData(old_level(l), this%level(l), x_plan(l))
!!       end if
!!
!!     end do

    ! finalization .............................................................

    deallocate(old_level)

  end subroutine FitAdapt

  !-----------------------------------------------------------------------------
  !> Copy retained element attributes

  subroutine CopyRetainedData(old_var, new_var, x_plan)
    class(MeshVariable_3D),     intent(in)    :: old_var
    class(MeshVariable_3D),     intent(inout) :: new_var
    class(DataExchangePlan_3D), intent(in)    :: x_plan

    type(DataExchangeSendBuf_3D), allocatable, save :: send_buf(:)
    type(DataExchangeRecvBuf_3D), allocatable, save :: recv_buf(:)

    integer :: n_recv, n_send
    integer :: i

    ! initialization ...........................................................

    if (allocated(x_plan % send_map)) then
      n_send = size(x_plan % send_map)
    else
      n_send = 0
    end if

    if (allocated(x_plan % recv_map)) then
      n_recv = size(x_plan % recv_map)
    else
      n_recv = 0
    end if

    if (max(n_send, n_recv) > 0) then
      allocate(send_buf(n_send), recv_buf(n_recv))
    else
      return
    end if
!### CHECK
block
  integer :: rank
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)
  print '(99(G0,X))', '### rank',rank,'n_send =',n_send,'n_recv =',n_recv
  do i = 1, n_send
    print '(99(G0,X))', '### rank',rank,'x_plan%send_map(',i,')%id_elem =',x_plan%send_map(i)%id_elem
  end do
  do i = 1, n_recv
    print '(99(G0,X))', '### rank',rank,'x_plan%recv_map(',i,')%id_elem =',x_plan%recv_map(i)%id_elem
  end do
end block
!### CHECK END

    ! transfer retained data ...................................................

    do i = 1, n_send
      call send_buf(i) % Extract_Data(x_plan % send_map(i), v = old_var % val)
      call send_buf(i) % Send_Start()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Init(x_plan % recv_map(i), v = new_var % val)
      call recv_buf(i) % Recv_Start()
    end do

    do i = 1, n_send
      call send_buf(i) % Send_Finish()
    end do

    do i = 1, n_recv
      call recv_buf(i) % Recv_Finish()
      call recv_buf(i) % Assign_Data(v = new_var % val)
    end do

    ! finalization .............................................................

    deallocate(send_buf, recv_buf)

  end subroutine CopyRetainedData

  !=============================================================================

end submodule MP_FitAdapt
