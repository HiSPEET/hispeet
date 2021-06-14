!> summary:  Evaluation of the valency of mesh points
!> author:   Joerg Stiller
!> date:     2020/12/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(Mesh_3d__Partition) MP_GetPointValency
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Computes the valency of mesh points
  !>
  !> The valency of a point equals the number of adjoining mesh elements.
  !> This procedure assumes that the underlying point set includes the element
  !> boundaries, whichis the case, e.g., for Lobatto points.
  !>
  !> With interior sets like Gauss the valency is 1 for all points and, hence,
  !> does not need to be computed.

  module subroutine GetPointValency(mesh, v)
    class(Mesh3d_Partition), intent(in)  :: mesh          !< mesh parition
    integer,                 intent(out) :: v(0:,0:,0:,:) !< point valency

    integer :: e, po

    po = ubound(v, 1)

    !$omp do
    do e = 1, mesh % n_elem

      v(:,:,:,e) = 1

      associate(face => mesh % element(e)%face)
        v(  0,  :,  :, e)  =  v(  0,  :,  :, e)  +  face(1) % n_neighbor
        v( po,  :,  :, e)  =  v( po,  :,  :, e)  +  face(2) % n_neighbor
        v(  :,  0,  :, e)  =  v(  :,  0,  :, e)  +  face(3) % n_neighbor
        v(  :, po,  :, e)  =  v(  :, po,  :, e)  +  face(4) % n_neighbor
        v(  :,  :,  0, e)  =  v(  :,  :,  0, e)  +  face(5) % n_neighbor
        v(  :,  :, po, e)  =  v(  :,  :, po, e)  +  face(6) % n_neighbor
      end associate

      associate(edge => mesh % element(e)%edge)
        v(  :,  0,  0, e)  =  v(  :,  0,  0, e)  +  edge( 1) % n_neighbor
        v(  :, po,  0, e)  =  v(  :, po,  0, e)  +  edge( 2) % n_neighbor
        v(  :,  0, po, e)  =  v(  :,  0, po, e)  +  edge( 3) % n_neighbor
        v(  :, po, po, e)  =  v(  :, po, po, e)  +  edge( 4) % n_neighbor
        v(  0,  :,  0, e)  =  v(  0,  :,  0, e)  +  edge( 5) % n_neighbor
        v( po,  :,  0, e)  =  v( po,  :,  0, e)  +  edge( 6) % n_neighbor
        v(  0,  :, po, e)  =  v(  0,  :, po, e)  +  edge( 7) % n_neighbor
        v( po,  :, po, e)  =  v( po,  :, po, e)  +  edge( 8) % n_neighbor
        v(  0,  0,  :, e)  =  v(  0,  0,  :, e)  +  edge( 9) % n_neighbor
        v( po,  0,  :, e)  =  v( po,  0,  :, e)  +  edge(10) % n_neighbor
        v(  0, po,  :, e)  =  v(  0, po,  :, e)  +  edge(11) % n_neighbor
        v( po, po,  :, e)  =  v( po, po,  :, e)  +  edge(12) % n_neighbor
      end associate

      associate(vertex => mesh % element(e)%vertex)
        v(  0,  0,  0, e)  =  v(  0,  0,  0, e)  +  vertex(1) % n_neighbor
        v( po,  0,  0, e)  =  v( po,  0,  0, e)  +  vertex(2) % n_neighbor
        v(  0, po,  0, e)  =  v(  0, po,  0, e)  +  vertex(3) % n_neighbor
        v( po, po,  0, e)  =  v( po, po,  0, e)  +  vertex(4) % n_neighbor
        v(  0,  0, po, e)  =  v(  0,  0, po, e)  +  vertex(5) % n_neighbor
        v( po,  0, po, e)  =  v( po,  0, po, e)  +  vertex(6) % n_neighbor
        v(  0, po, po, e)  =  v(  0, po, po, e)  +  vertex(7) % n_neighbor
        v( po, po, po, e)  =  v( po, po, po, e)  +  vertex(8) % n_neighbor
      end associate

    end do

  end subroutine GetPointValency

  !=============================================================================

end submodule MP_GetPointValency
