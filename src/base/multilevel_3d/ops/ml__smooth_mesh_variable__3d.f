!> summary:  Utility for removing jumps and smoothing multilevel mesh variables
!> author:   Joerg Stiller
!> date:     2026/08/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module ML__Smooth_Mesh_Variable__3D
  use Smooth_Mesh_Data__3D
  use Child_To_Parent_Projection__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  implicit none
  private

  public :: ML_SmoothMeshVariable_3D

contains

  !-----------------------------------------------------------------------------
  !> Remove jumps and smooth multilevel mesh data

  subroutine ML_SmoothMeshVariable_3D(ml_op, u, filter, order)
    class(ML_MeshOperators_3D), intent(in)    :: ml_op !< ML mesh operators
    class(ML_MeshVariable_3D),  intent(inout) :: u     !< ML mesh variable
    integer, intent(in) :: filter !< 0/1/2/3: none/cut-off/erfc-log/exponential
    integer, intent(in) :: order  !< filter order (0: auto)

    integer :: l

    do l = 1, size(u%level)

      ! smooth current mesh level
      call SmoothMeshData_3D( mesh   = ml_op % sem(l) % mesh   &
                            , eop    = ml_op % sem(l) % std_op &
                            , u      = u % level(l) % val      &
                            , filter = filter                  &
                            , order  = order                   )

      if (l == 1) exit

      ! project to parent level (currently using embedded interpolation)
      call ChildToParentProjection_3D( child  = ml_op % sem(l  ) % mesh &
                                     , parent = ml_op % sem(l-1) % mesh &
                                     , pop    = ml_op % iop_fc_x(l)     &
                                     , v_c    = u % level(l  ) % val    &
                                     , v_p    = u % level(l-1) % val    )

    end do

  end subroutine ML_SmoothMeshVariable_3D

  !=============================================================================

end module ML__Smooth_Mesh_Variable__3D
