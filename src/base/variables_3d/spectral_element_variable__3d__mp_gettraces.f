!> summary:  3D spectral element variable: TBP for extracting traces
!> author:   Joerg Stiller
!> date:     2021/8/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - check if OpenMP parallelization makes sense
!===============================================================================

submodule(Spectral_Element_Variable__3D) MP_GetTraces
  use Mesh__3D
  use Trace_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Extract the traces of SEV components on mesh faces
  !>
  !> The traces are shaped as `tr_val(np,np,2,nf,nc)` where
  !>    - `np` is the number of element points per direction, i.e. `po+1`
  !>    - `nf` is the number of local mesh faces
  !>    - `nc` is the number of components, i.e. `size(this%val,5)`

  module subroutine GetTraces(this, tr_val)
    class(SpectralElementVariable_3D), intent(in) :: this
    real(RNP), intent(out) :: tr_val(:,:,:,:,:) !< trace of val

    type(TraceTransferBuffer_3D), asynchronous, allocatable, save :: trace_buf
    integer :: c, e, i, f(6), s(6), nc, o, p

    associate(mesh => this % sem % mesh, val => this % val)

      ! initialization .........................................................

      !$omp master
      trace_buf = TraceTransferBuffer_3D(mesh, tr_val)
      !$omp master
      !$omp barrier

      nc = size(tr_val,5)

      ! local traces ...........................................................

      o = lbound(val,1)
      p = ubound(val,1)
      !$omp do
      do e = 1, mesh % n_elem
        associate(face => mesh % element(e) % face)

          do i = 1, 6
            f(i) = face(i) % id      ! mesh face adjacent to element face i
            s(i) = face(i) % Side()  ! mesh face side that is touched {1,2}
          end do

          if (mesh % structured) then
            do c = 1, nc
              tr_val(:,:,s(1),f(1),c) = val(o,:,:,e,c)
              tr_val(:,:,s(2),f(2),c) = val(p,:,:,e,c)
              tr_val(:,:,s(3),f(3),c) = val(:,p,:,e,c)
              tr_val(:,:,s(4),f(4),c) = val(:,p,:,e,c)
              tr_val(:,:,s(5),f(5),c) = val(:,:,o,e,c)
              tr_val(:,:,s(6),f(6),c) = val(:,:,p,e,c)
            end do
          else
            do c = 1, nc
              call face(1) % AlignWithMesh(val(o,:,:,e,c), tr_val(:,:,s(1),f(1),c))
              call face(2) % AlignWithMesh(val(p,:,:,e,c), tr_val(:,:,s(2),f(2),c))
              call face(3) % AlignWithMesh(val(:,p,:,e,c), tr_val(:,:,s(3),f(3),c))
              call face(4) % AlignWithMesh(val(:,p,:,e,c), tr_val(:,:,s(4),f(4),c))
              call face(5) % AlignWithMesh(val(:,:,o,e,c), tr_val(:,:,s(5),f(5),c))
              call face(6) % AlignWithMesh(val(:,:,p,e,c), tr_val(:,:,s(6),f(6),c))
            end do
          end if

        end associate
      end do

      ! transfer to/from adjoining partitions ..................................

      call trace_buf % Transfer(mesh, tr_val, tag=100)
      call trace_buf % Merge(mesh, tr_val, alpha=ZERO, beta=ONE)

      !$omp barrier
      !$omp master
      deallocate(trace_buf)
      !$omp end master

    end associate

  end subroutine GetTraces

  !=============================================================================

end submodule MP_GetTraces
