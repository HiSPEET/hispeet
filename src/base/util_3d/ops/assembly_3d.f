!> summary:  Assembly of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2021/03/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Assembly_3d
  use Kind_Parameters, only: RNP
  use Constants      , only: ONE
  use Execution_Control
  use Mesh_3d__Partition
  use Element_Transfer_Buffer_3d
  implicit none
  private

  public :: Assembly3d

  interface Assembly3d
    module procedure :: Assembly3d_0
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Unstructured assembly and, optionally, averaging of a scalar variable

    module subroutine Assembly3d_U0(mesh, u, u_buf, avg)
      !> mesh partition
      class(Mesh3d_Partition), intent(in)    :: mesh
      !> mesh variable, including ghost entries
      real(RNP), intent(inout) :: u(:,:,:,:)
      !> MPI transfer buffer
      class(ElementTransferBuffer3d), asynchronous, intent(inout) :: u_buf
      !> switch for averaging over element boundaries [F]
      logical, optional, intent(in) :: avg
    end subroutine Assembly3d_U0

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Assembly and, optionally, averaging of a scalar variable

  subroutine Assembly3d_0(mesh, u, u_buf, avg)
    !> mesh partition
    class(Mesh3d_Partition), intent(in)    :: mesh
    !> mesh variable, including ghost entries
    real(RNP), intent(inout) :: u(:,:,:,:)
    !> MPI transfer buffer
    type(ElementTransferBuffer3d), asynchronous, intent(inout) :: u_buf
    !> switch for averaging over element boundaries [F]
    logical, optional, intent(in) :: avg

    !$omp master
    if (size(u,4) /= mesh%n_elem + mesh%n_ghost) then
      call Error('Assembly3d_0', 'Dimension 4 of u does not match')
    end if
    !$omp end master

    ! insert special treatment of structured case here
    call Assembly3d_U0(mesh, u, u_buf, avg)

  end subroutine Assembly3d_0

  !=============================================================================

end module Assembly_3d
