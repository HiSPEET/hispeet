!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  SFC based ordering of mesh elements
!> author:   Joerg Stiller
!> date:     2026/01/16
!>
!> The mesh elements are ordered corresponding to their coordinates on a Hilbert
!> curve defined in a bounding box that contains the element centers. Connecting
!> the centers accordingly yields the approximation of a space filling curve
!> through the computational domain. Although the curve may be disrupted in
!> complex configurations, it can still be used for partitioning. Furthermore,
!> the path of the SFC is identified within each element, allowing it to be
!> subdivided without recursion when the mesh is refined.
!>
!> @todo
!>   Generalize SFC generation to partitioned meshes
!===============================================================================

submodule(Mesh__3D) MP_BuildSFC
  use Execution_Control
  use Hilbert_Curve
  use Quick_Sort
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Arranges the elements along a Hilbert curve defined in the bounding box
  !>
  !> Requires
  !>   - mesh % n_elem
  !>   - mesh % element % vertex
  !>   - mesh % element % edge
  !>   - mesh % element % face
  !>   - mesh % element % neighbor
  !>   - mesh % element % geometry % x
  !>
  !> Generates
  !>   - mesh % element % sfc_rank
  !>   - mesh % element % sfc_path
  !>
  !> This routine is usually applied only to the root mesh.
  !>
  !> The current version requires the mesh to be unpartitioned.
  !> If the mesh is partitioned, no SFC will be generated.

  module subroutine BuildSFC(mesh)
    class(Mesh_3D), intent(inout)  :: mesh !< mesh partition

    real(RNP), allocatable, save :: x(:,:) ! element midpoint coordinates
    real(RNP), allocatable, save :: c(:)   ! element curve coordinates
    integer,   allocatable, save :: map(:) ! map from element index to rank

    real(RNP), save :: dx_min(3) ! min mesh spacing in directions 1:3
    real(RNP), save :: x_min(3)  ! corner of bounding box closest to -∞
    real(RNP), save :: x_max(3)  ! corner of bounding box closest to +∞

    logical :: uniform = .true.  ! T/F for uniform/direction-dependent scaling

    integer      :: lc    ! Hilbert curve level
    integer(IXL) :: ix(3) ! 3D index in scaled bounding box
    integer(IXL) :: ic    ! 1D index on Hilbert curve

    real(RNP) :: dx(3), lx(3), sc(3)
    integer   :: e_pre, e_suc, v_pre, v_start, v_exit
    integer   :: e, i

    associate(ne => mesh % n_elem)

      ! initial checks .........................................................

      if (mesh % n_parts > 1) then
        call Warning('BuildSFC', 'SFC not built for distributed mesh')
        return
      end if

      if (ne < 1) then
        mesh % has_sfc = .true.
        return
      end if

      ! initialization .........................................................

      !$omp master
      dx_min =  huge(dx_min)
      x_min  =  huge(x_min)
      x_max  = -huge(x_max)
      allocate(x(ne,3), c(ne), map(ne))
      !$omp end master
      !$omp barrier

      ! bounding box and mesh spacing ..........................................

      !$omp do reduction(min:dx_min,x_min) reduction(max:x_max)
      do e = 1, ne

        call mesh % element(e) % GetCuboidDimensions(dx)
        dx_min = min(dx_min, dx)

        do i = 1, 3
          x(e,i)   = mesh % element(e) % geometry % x_c(0,i)
          x_min(i) = min(x_min(i), x(e,i))
          x_max(i) = max(x_max(i), x(e,i))
        end do

      end do

      ! curve level and scaling factor .........................................

      if (uniform) then
        lx = maxval(x_max - x_min)
        dx = minval(dx_min)
        x_min = (x_min + x_max - lx) / 2
      else
        lx = x_max - x_min
        dx = dx_min
      end if
      lc = int( log(maxval(lx/dx)) / log(real(2,RNP)) ) + 1
      sc = (2**lc - 1) / lx

      ! midpoint curve positions ...............................................

      !$omp do
      do e = 1, ne
        ix = nint(sc * (x(e,1:3) - x_min))
        call Hilbert_C2I(lc, ix, ic)
        c(e) = real(ic, RNP)
      end do

      ! sfc ranking ............................................................

      ! sort elements WRT to curve positions
      !$omp master
      call SortIndex(c, map)
      !$omp end master

      !$omp do
      do i = 1, ne
        mesh % element(map(i)) % sfc_rank = i
      end do

      ! sfc path inside elements ...............................................

      !$omp master
      do i = 1, ne
        e = map(i)

        ! predecessor ID and connected vertex
        if (i > 1) then
          e_pre = map(i-1)
          v_pre = mod(mesh%element(e_pre)%sfc_path, 10)
        else
          e_pre = 0
          v_pre = 0
        end if
        v_start = StartPosition(mesh%element(e), e_pre, v_pre)

        ! successor
        if (i < ne) then
          e_suc = map(i+1)
        else
          e_suc = 0
        end if
        v_exit = ExitPosition(mesh%element(e), v_start, e_suc)

        mesh % element(e) % sfc_path = 10 * v_start + v_exit

      end do

      mesh % has_sfc = .true.
      !$omp end master

      ! finalization ...........................................................

      !$omp master
      deallocate(x, c, map)
      !$omp end master

    end associate

  end subroutine BuildSFC

  !-----------------------------------------------------------------------------
  !> Defines the start position of the SFC inside the given element

  function StartPosition(element, e_pre, v_pre) result(v_start)
    class(MeshElement_3D), intent(in) :: element !< element
    integer, intent(in) :: v_pre   !< connected predecessor vertex, 0 if none
    integer, intent(in) :: e_pre   !< ID of predecessor, 0 if none
    integer             :: v_start !< vertex from which the SFC starts

    integer :: m_pre(8)     ! marker for predecessor element vertices
    integer :: mf_pre(2,2)  ! marker for predecessor element face vertices
    integer :: mf(2,2)      ! marker for current element face vertices
    integer :: me_pre(2)    ! marker for predecessor element edge vertices
    integer :: me(2)        ! marker for current element edge vertices

    integer :: i, j, k, l, n, n0, nn

    associate(neighbor => element % neighbor)

      ! initialization .........................................................

      v_start = 1

      PREDECESSOR: if (e_pre > 0) then

        ! mark connected predecessor element vertex .............................

        m_pre = 0
        select case(v_pre)
        case(1:8)
          m_pre(v_pre) = 1
        case default
          exit PREDECESSOR
        end select

        ! check for predecessor on element face ................................

        FACES: do k = 1, 6
          n0 = element % face(k) % i_neighbor
          nn = element % face(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_pre) then
              l = ElementFaceID(neighbor(n)%component)
              exit FACES
            end if
          end do
        end do FACES

        if (k <= 6) then
          mf_pre = reshape(m_pre(V_FACE(:,l)), [2,2])
          call element % AlignFromNeighborFace(k, n, mf_pre, mf)
          do j = 1, 2
          do i = 1, 2
            if (mf(i,j) == 1) then
              v_start = V_FACE(i + 2*(j-1), k)
              exit PREDECESSOR
            end if
          end do
          end do
        end if

        ! check for predecessor on element edge ................................

        EDGES: do k = 1, 12
          n0 = element % edge(k) % i_neighbor
          nn = element % edge(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_pre) then
              l = ElementEdgeID(neighbor(n)%component)
              exit EDGES
            end if
          end do
        end do EDGES

        if (k <= 12) then
          me_pre = m_pre(V_EDGE(:,l))
          call element % AlignFromNeighborEdge(k, n, me_pre, me)
          do i = 1, 2
            if (me(i) == 1) then
              v_start = V_EDGE(i, k)
              exit PREDECESSOR
            end if
          end do
        end if

        ! check for predecessor on element vertex ..............................

        VERTICES: do k = 1, 8
          n0 = element % vertex(k) % i_neighbor
          nn = element % vertex(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_pre) exit VERTICES
          end do
        end do VERTICES

        if (k <= 8) then
          v_start = k
        end if

      end if PREDECESSOR

    end associate

  end function StartPosition

  !-----------------------------------------------------------------------------
  !> Defines the end position of the SFC inside the given element

  function ExitPosition(element,  v_start, e_suc) result(v_exit)
    class(MeshElement_3D), intent(in) :: element !< element
    integer, intent(in) :: v_start !< vertex from which the SFC starts
    integer, intent(in) :: e_suc   !< ID of successor, 0 if none
    integer             :: v_exit  !< vertex in which the SFC exits

    integer :: v1, v2, v3
    integer :: k, n, n0, nn

    associate(neighbor => element % neighbor)

      ! identify opposite element vertices .....................................

      select case(v_start)
      case(1)
        v1 = 2
        v2 = 3
        v3 = 5
      case(2)
        v1 = 1
        v2 = 4
        v3 = 6
      case(3)
        v1 = 4
        v2 = 1
        v3 = 7
      case(4)
        v1 = 3
        v2 = 2
        v3 = 8
      case(5)
        v1 = 6
        v2 = 7
        v3 = 1
      case(6)
        v1 = 5
        v2 = 8
        v3 = 2
      case(7)
        v1 = 8
        v2 = 5
        v3 = 3
      case(8)
        v1 = 7
        v2 = 6
        v3 = 4
      end select

      ! default ................................................................

      v_exit = v1

      SUCCESSOR: if (e_suc > 0) then

        ! check for successor on element face ..................................

        FACES: do k = 1, 6
          n0 = element % face(k) % i_neighbor
          nn = element % face(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_suc) exit FACES
          end do
        end do FACES

        if (k <= 6) then
          if (any(V_FACE(:,k) == v1)) then
            v_exit = v1
          else if (any(V_FACE(:,k) == v2)) then
            v_exit = v2
          else
            v_exit = v3
          end if
          exit SUCCESSOR
        end if

        ! check for successor on element edge ..................................

        EDGES: do k = 1, 12
          n0 = element % edge(k) % i_neighbor
          nn = element % edge(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_suc) exit EDGES
          end do
        end do EDGES

        if (k <= 12) then
          if (any(V_EDGE(:,k) == v1)) then
            v_exit = v1
          else if (any(V_EDGE(:,k) == v2)) then
            v_exit = v2
          else
            v_exit = v3
          end if
          exit SUCCESSOR
        end if

        ! check for successor on element vertex ................................

        VERTICES: do k = 1, 8
          n0 = element % vertex(k) % i_neighbor
          nn = element % vertex(k) % n_neighbor
          if (nn < 1) cycle
          do n = n0, n0+nn-1
            if (neighbor(n) % id == e_suc) exit VERTICES
          end do
        end do VERTICES

        if (k == v1 .or. k == v2 .or. k == v3) then
          v_exit = k
        end if

      end if SUCCESSOR

    end associate

  end function ExitPosition

  !=============================================================================

end submodule MP_BuildSFC
