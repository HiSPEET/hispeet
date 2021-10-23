!> summary:  Assembly of 3D mesh variables
!> author:   Joerg Stiller
!> date:     2021/03/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Assembly__3D
  use Kind_Parameters, only: RNP
  use Constants      , only: ONE
  use Execution_Control
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  implicit none
  private

  public :: Assembly_3D

  interface Assembly_3D
    module procedure :: Assembly_3D_0
  end interface

  interface

    !---------------------------------------------------------------------------
    !> Unstructured assembly and, optionally, averaging of a scalar variable

    module subroutine Assembly_3D_0U(mesh, u, u_buf, avg)
      !> mesh partition
      class(Mesh_3D), intent(in) :: mesh
      !> mesh variable, including ghost entries
      real(RNP), intent(inout) :: u(:,:,:,:)
      !> MPI transfer buffer
      class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: u_buf
      !> switch for averaging over element boundaries [F]
      logical, optional, intent(in) :: avg
    end subroutine Assembly_3D_0U

  end interface

contains

  !-----------------------------------------------------------------------------
  !> Assembly and, optionally, averaging of a scalar variable

  subroutine Assembly_3D_0(mesh, u, u_buf, avg)
    !> mesh partition
    class(Mesh_3D), intent(in) :: mesh
    !> mesh variable, including ghost entries
    real(RNP), intent(inout) :: u(:,:,:,:)
    !> MPI transfer buffer
    type(ElementTransferBuffer_3D), asynchronous, intent(inout) :: u_buf
    !> switch for averaging over element boundaries [F]
    logical, optional, intent(in) :: avg

    !$omp master
    if (size(u,4) /= mesh%n_elem + mesh%n_ghost) then
      call Error('Assembly_3D_0', 'Dimension 4 of u does not match')
    end if
    !$omp end master

    ! insert special treatment of structured case here
    call Assembly_3D_0U(mesh, u, u_buf, avg)

  end subroutine Assembly_3D_0

  !=============================================================================

end module Assembly__3D
