!> summary:  Conversion of interior traces to exterior traces
!> author:   Joerg Stiller
!> date:     2022/01/11
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module Trace_Operators__3D
  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO
  use Mesh__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none
  private

  public :: GetInnerTraces_3D
  public :: GetOuterTraces_3D
  public :: ConvertInnerToOuterTraces_3D

  interface GetInnerTraces_3D
    module procedure GetInnerTraces_S
    module procedure GetInnerTraces_A
  end interface

  interface GetOuterTraces_3D
    module procedure GetOuterTraces_S
    module procedure GetOuterTraces_A
  end interface

  interface ConvertInnerToOuterTraces_3D
    module procedure ConvertInnerToOuterTraces_S
    module procedure ConvertInnerToOuterTraces_A
  end interface

contains

  !=============================================================================
  ! GetInnerTraces_3D variants

  !-----------------------------------------------------------------------------
  !> Extract the inner traces `u⁻` of a scalar mesh variable `u`
  !>
  !> The traces are stored as element-face variables. If `align` is passed `T`
  !> they are aligned with the mesh faces, otherwise with the element faces.
  !> Traces of remote elements are stored in their ghost entries. Therefore,
  !> the trace values must be dimensioned as `um(np,np,6,ne+ng,nc)`, where
  !>
  !>   - `np` is the number of element values per direction, i.e. `size(u,1)`
  !>   - `ne` is the number of local elements, i.e. `mesh % n_elem`
  !>   - `ng` is the number of ghost elements, i.e. `mesh % n_ghost`

  module subroutine GetInnerTraces_S(mesh, u, um, align)
    class(Mesh_3D),        intent(in)    :: mesh
    real(RNP), contiguous, intent(in)    :: u (:,:,:,:) !< u
    real(RNP), contiguous, intent(inout) :: um(:,:,:,:) !< u⁻ with ghosts
    logical, optional,     intent(in)    :: align !< align u⁻ with mesh face [F]

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_um
    integer :: e, o, p
    logical :: as_is

    if (present(align)) then
      as_is = .not. align
    else
      as_is = .true.
    end if

    ! initialization ...........................................................

    !$omp master
    buf_um = ElementFaceTransferBuffer_3D(mesh, um)
    !$omp end master
    !$omp barrier

    ! local traces .............................................................

    o = lbound(u,1)
    p = ubound(u,1)
    !$omp do
    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)

        if (as_is .or. mesh % structured) then
          um(:,:,1,e) = u(o,:,:,e)
          um(:,:,2,e) = u(p,:,:,e)
          um(:,:,3,e) = u(:,o,:,e)
          um(:,:,4,e) = u(:,p,:,e)
          um(:,:,5,e) = u(:,:,o,e)
          um(:,:,6,e) = u(:,:,p,e)
        else
          call face(1) % AlignWithMesh(u(o,:,:,e), um(:,:,1,e))
          call face(2) % AlignWithMesh(u(p,:,:,e), um(:,:,2,e))
          call face(3) % AlignWithMesh(u(:,o,:,e), um(:,:,3,e))
          call face(4) % AlignWithMesh(u(:,p,:,e), um(:,:,4,e))
          call face(5) % AlignWithMesh(u(:,:,o,e), um(:,:,5,e))
          call face(6) % AlignWithMesh(u(:,:,p,e), um(:,:,6,e))
        end if

      end associate
    end do

    ! transfer to/from adjoining partitions ....................................

    call buf_um % Transfer(mesh, um, tag=100)
    call buf_um % Merge(um)

    !$omp master
    deallocate(buf_um)
    !$omp end master

   end subroutine GetInnerTraces_S

  !-----------------------------------------------------------------------------
  !> Extract the inner traces `u⁻` of an array-valued mesh variable `u`
  !>
  !> The traces are stored as element-face variables. If `align` is passed `T`
  !> they are aligned with the mesh faces, otherwise with the element faces.
  !> Traces of remote elements are stored in their ghost entries. Therefore,
  !> the trace values must be dimensioned as `um(np,np,6,ne+ng,nc)`, where
  !>
  !>   - `np` is the number of element values per direction, i.e. `size(u,1)`
  !>   - `ne` is the number of local elements, i.e. `mesh % n_elem`
  !>   - `ng` is the number of ghost elements, i.e. `mesh % n_ghost`
  !>   - `nc` is the number of components

  subroutine GetInnerTraces_A(mesh, u, um, align)
    class(Mesh_3D),        intent(in)    :: mesh
    real(RNP), contiguous, intent(in)    :: u (:,:,:,:,:) !< u
    real(RNP), contiguous, intent(inout) :: um(:,:,:,:,:) !< u⁻ with ghosts
    logical, optional,     intent(in)    :: align !< align u⁻ with mesh face [F]

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: buf_um
    integer :: c, e, nc, o, p
    logical :: as_is

    if (present(align)) then
      as_is = .not. align
    else
      as_is = .true.
    end if

    ! initialization ...........................................................

    !$omp master
    buf_um = ElementFaceTransferBuffer_3D(mesh, um)
    !$omp end master
    !$omp barrier

    nc = size(um,5)

    ! local traces .............................................................

    o = lbound(u,1)
    p = ubound(u,1)
    !$omp do
    do e = 1, mesh % n_elem
      associate(face => mesh % element(e) % face)

        if (as_is .or. mesh % structured) then
          do c = 1, nc
            um(:,:,1,e,c) = u(o,:,:,e,c)
            um(:,:,2,e,c) = u(p,:,:,e,c)
            um(:,:,3,e,c) = u(:,o,:,e,c)
            um(:,:,4,e,c) = u(:,p,:,e,c)
            um(:,:,5,e,c) = u(:,:,o,e,c)
            um(:,:,6,e,c) = u(:,:,p,e,c)
          end do
        else
          do c = 1, nc
            call face(1) % AlignWithMesh(u(o,:,:,e,c), um(:,:,1,e,c))
            call face(2) % AlignWithMesh(u(p,:,:,e,c), um(:,:,2,e,c))
            call face(3) % AlignWithMesh(u(:,o,:,e,c), um(:,:,3,e,c))
            call face(4) % AlignWithMesh(u(:,p,:,e,c), um(:,:,4,e,c))
            call face(5) % AlignWithMesh(u(:,:,o,e,c), um(:,:,5,e,c))
            call face(6) % AlignWithMesh(u(:,:,p,e,c), um(:,:,6,e,c))
          end do
        end if

      end associate
    end do

    ! transfer to/from adjoining partitions ....................................

    call buf_um % Transfer(mesh, um, tag=100)
    call buf_um % Merge(um)

    !$omp master
    deallocate(buf_um)
    !$omp end master

   end subroutine GetInnerTraces_A

  !=============================================================================
  ! GetOuterTraces_3D variants

  !-----------------------------------------------------------------------------
  !> Build the outer traces `u⁺` of an array-valued mesh variable `u`
  !>
  !> The traces are stored as element-face variables. They are aligned with the
  !> faces of the corresponding elements. Ghost entries are assumed to be absent
  !> and hence not touched.

  subroutine GetOuterTraces_A(mesh, u, up)
    class(Mesh_3D),        intent(in)  :: mesh
    real(RNP), contiguous, intent(in)  :: u (:,:,:,:,:) !< u
    real(RNP), contiguous, intent(out) :: up(:,:,:,:,:) !< u⁺

    real(RNP), allocatable, save :: um(:,:,:,:,:)
    integer :: nc, np

    !$omp master
    np = size(u,1)
    nc = size(u,5)
    allocate(um(np, np, np, mesh%n_elem + mesh%n_ghost, nc), source = ZERO)
    !$omp end master
    ! no barrier required ;)

    call GetInnerTraces_A(mesh, u, um)
    call ConvertInnerToOuterTraces_A(mesh, um, up)

    !$omp master
    deallocate(um)
    !$omp end master

  end subroutine GetOuterTraces_A

  !-----------------------------------------------------------------------------
  !> Build the outer traces `u⁺` of a scalar mesh variable `u`
  !>
  !> The traces are stored as element-face variables. They are aligned with the
  !> faces of the corresponding elements. Ghost entries are assumed to be absent
  !> and hence not touched.

  subroutine GetOuterTraces_S(mesh, u, up)
    class(Mesh_3D),        intent(in)  :: mesh
    real(RNP), contiguous, intent(in)  :: u (:,:,:,:) !< u
    real(RNP), contiguous, intent(out) :: up(:,:,:,:) !< u⁺

    real(RNP), allocatable, save :: um(:,:,:,:)
    integer :: np

    !$omp master
    np = size(u,1)
    allocate(um(np, np, np, mesh%n_elem + mesh%n_ghost), source = ZERO)
    !$omp end master
    ! no barrier required ;)

    call GetInnerTraces_S(mesh, u, um)
    call ConvertInnerToOuterTraces_S(mesh, um, up)

    !$omp master
    deallocate(um)
    !$omp end master

  end subroutine GetOuterTraces_S

  !=============================================================================
  ! ConvertInnerToOuterTraces_3D variants

  !-----------------------------------------------------------------------------
  !> Converts interior traces u⁻ to exterior traces u⁺ -- scalar version
  !>
  !> This routine converts the interior traces `u⁻ = um` into exterior traces
  !> `u⁺ = up`. Both are aligned with the element faces. In order to transfer
  !> the traces of remote elements, they must be present in the corresponding
  !> ghost entries of `um`.

  subroutine ConvertInnerToOuterTraces_S(mesh, um, up)
    class(Mesh_3D),        intent(in)  :: mesh
    real(RNP), contiguous, intent(in)  :: um(:,:,:,:) !< u⁻
    real(RNP), contiguous, intent(out) :: up(:,:,:,:) !< u⁺

    integer :: e, i, f, l, m

    !$omp do
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))
        do f = 1, 6
          i = element % face(f) % i_neighbor
          if (i > 0) then ! adopt trace from neighbor
            l = element % neighbor(i) % id
            m = element % neighbor(i) % component
            call element % AlignFromNeighborFace(f, i, um(:,:,m,l), up(:,:,f,e))
          else ! no neighbor
            up(:,:,f,e) = um(:,:,f,e)
          end if
        end do
      end associate
    end do

  end subroutine ConvertInnerToOuterTraces_S

  !-----------------------------------------------------------------------------
  !> Converts interior traces u⁻ to exterior traces u⁺ -- array version

  subroutine ConvertInnerToOuterTraces_A(mesh, um, up)
    class(Mesh_3D),        intent(in)  :: mesh
    real(RNP), contiguous, intent(in)  :: um(:,:,:,:,:) !< u⁻
    real(RNP), contiguous, intent(out) :: up(:,:,:,:,:) !< u⁺

    integer :: c, e, i, f, l, m, nc

    nc = size(um,5)

    !$omp do
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))
        do f = 1, 6
          i = element % face(f) % i_neighbor
          if (i > 0) then ! adopt trace from neighbor
            l = element % neighbor(i) % id
            m = element % neighbor(i) % component
            do c = 1, nc
              call element % AlignFromNeighborFace &
                                 (f, i, um(:,:,m,l,c), up(:,:,f,e,c))
            end do
          else ! no neighbor
            do c = 1, nc
              up(:,:,f,e,c) = um(:,:,f,e,c)
            end do
          end if
        end do
      end associate
    end do

  end subroutine ConvertInnerToOuterTraces_A

  !=============================================================================

end module Trace_Operators__3D
