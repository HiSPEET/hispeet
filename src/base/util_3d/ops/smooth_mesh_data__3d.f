!> summary:  Utility for removing jumps and smoothing mesh data
!> author:   Joerg Stiller
!> date:     2026/08/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Smooth_Mesh_Data__3D
  use Kind_Parameters
  use Constants, only: ZERO
  use Array_Assignments
  use Standard_Element_Operators__1D
  use TPO__AAA__3D
  use Assembly__3D
  use Mesh__3D
  use Element_Transfer_Buffer__3D
  implicit none
  private

  public :: SmoothMeshData_3D

contains

  !-----------------------------------------------------------------------------
  !> Remove jumps and smooth mesh data

  subroutine SmoothMeshData_3D(mesh, eop, u, filter, order)

    ! arguments ................................................................

    class(Mesh_3D),                     intent(in) :: mesh  !< mesh partition
    class(StandardElementOperators_1D), intent(in) :: eop   !< element operators
    real(RNP), contiguous, intent(inout) :: u(0:,0:,0:,:,:) !< mesh data

    integer, intent(in) :: filter !< 0/1/2/3: none/cut-off/erfc-log/exponential
    integer, intent(in) :: order  !< filter order (0: auto)


    ! internal variables .......................................................

    type(ElementTransferBuffer_3D), allocatable, asynchronous, save :: buf_s
    real(RNP), allocatable, save :: s(:,:,:,:)

    real(RNP), allocatable :: A(:,:)
    integer :: po, pf
    integer :: c

    ! initialization ...........................................................

    po = ubound(u,1)

    !$omp master
    ! data structures for removal of discontinuities
    allocate(s(0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost), source = ZERO)
    buf_s = ElementTransferBuffer_3D(mesh, s)
    !$omp end master
    !$omp barrier

    ! set filter order
    if (order > 0) then
      pf = order
    else
      select case(filter)
      case(1)
        ! default cut-off degree
        pf = max(po-2, 2)
      case(2)
        ! default order of erfc-log filter
        pf = 4
      case(3)
        ! default order of exponential filter
        pf = 5
      end select
    end if

    if (po > 1) then
      ! get filter matrix
      allocate(A(0:po,0:po))
      select case(filter)
      case(1)
        call eop % Get_Bubble_CutoffFilter(pf, A)
      case(2)
        call eop % Get_ErfcLogFilter(real(pf,RNP), A, modes='B')
      case(3)
        call eop % Get_ExponentialFilter(real(pf,RNP), A, modes='B')
      end select
    end if

    do c = 1, size(u,5)
      ! step 1: remove jumps
      call SetArray(s(:,:,:,1:mesh%n_elem), u(:,:,:,:,c))
      call Assembly_3D(mesh, s, buf_s, avg=.true.)
      ! step 2: filter
      if (allocated(A)) then
        call TPO_AAA(A, s(:,:,:,1::mesh%n_elem), u(:,:,:,:,c))
      else
        call SetArray(u(:,:,:,:,c), s(:,:,:,1:mesh%n_elem))
      end if
    end do

    !$omp master
    deallocate(s, buf_s)
    !$omp end master

  end subroutine SmoothMeshData_3D

  !=============================================================================

end module Smooth_Mesh_Data__3D
