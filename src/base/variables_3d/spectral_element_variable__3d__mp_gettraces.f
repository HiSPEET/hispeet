!> summary:  3D spectral element variable: TBP for extracting traces
!> author:   Joerg Stiller
!> date:     2021/08/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Spectral_Element_Variable__3D) MP_GetTraces
  use Mesh__3D
  use Element_Face_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Extract the traces of SEV components as an element-face variable aligned
  !> with mesh faces
  !>
  !> The traces are stored as element-face variables. If `aligned`is passed `T`
  !> they are aligned with the mesh faces, otherwise with the element faces.
  !> Traces of remote elements are stored in their ghost entries. Therefore,
  !> the trace values must be dimensioned as `tr_val(np,np,6,nl+ng,nc)`, where
  !>
  !>   - `np` is the number of element values per direction
  !>   - `nl` is the number of local elements, i.e. `this%sem%mesh%n_elem`
  !>   - `ng` is the number of ghost elements, i.e. `this%sem%mesh%n_ghost`
  !>   - `nc` is the number of components

  module subroutine GetTraces(this, tr_val, align)
    class(SpectralElementVariable_3D), intent(in) :: this
    real(RNP), contiguous, intent(inout) :: tr_val(:,:,:,:,:) !< trace of val
    logical, optional, intent(in) :: align !< align traces with mesh face [F]

    type(ElementFaceTransferBuffer_3D), asynchronous, allocatable, save :: tr_buf
    integer :: c, e, nc, o, p
    logical :: as_is

    if (present(align)) then
      as_is = .not. align
    else
      as_is = .true.
    end if

    associate(mesh => this % sem % mesh, val => this % val)

      ! initialization .........................................................

      !$omp master
      tr_buf = ElementFaceTransferBuffer_3D(mesh, tr_val)
      !$omp end master
      !$omp barrier

      nc = size(tr_val,5)

      ! local traces ...........................................................

      o = lbound(val,1)
      p = ubound(val,1)
      !$omp do
      do e = 1, mesh % n_elem
        associate(face => mesh % element(e) % face)

          if (as_is .or. mesh % structured) then
            do c = 1, nc
              tr_val(:,:,1,e,c) = val(o,:,:,e,c)
              tr_val(:,:,2,e,c) = val(p,:,:,e,c)
              tr_val(:,:,3,e,c) = val(:,o,:,e,c)
              tr_val(:,:,4,e,c) = val(:,p,:,e,c)
              tr_val(:,:,5,e,c) = val(:,:,o,e,c)
              tr_val(:,:,6,e,c) = val(:,:,p,e,c)
            end do
          else
            do c = 1, nc
              call face(1) % AlignWithMesh(val(o,:,:,e,c), tr_val(:,:,1,e,c))
              call face(2) % AlignWithMesh(val(p,:,:,e,c), tr_val(:,:,2,e,c))
              call face(3) % AlignWithMesh(val(:,o,:,e,c), tr_val(:,:,3,e,c))
              call face(4) % AlignWithMesh(val(:,p,:,e,c), tr_val(:,:,4,e,c))
              call face(5) % AlignWithMesh(val(:,:,o,e,c), tr_val(:,:,5,e,c))
              call face(6) % AlignWithMesh(val(:,:,p,e,c), tr_val(:,:,6,e,c))
            end do
          end if

        end associate
      end do

      ! transfer to/from adjoining partitions ..................................

      call tr_buf % Transfer(mesh, tr_val, tag=100)
      call tr_buf % Merge(tr_val)

      !$omp barrier
      !$omp master
      deallocate(tr_buf)
      !$omp end master

    end associate

   end subroutine GetTraces

  !=============================================================================

end submodule MP_GetTraces
