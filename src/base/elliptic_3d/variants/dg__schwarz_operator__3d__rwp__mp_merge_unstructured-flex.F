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

  !-----------------------------------------------------------------------------
  !> Merge subdomain contributions into unstructured mesh variable
  !>
  !> The strategy is as a follows:
  !>
  !>   1. Send `us` from linked local elements to ghosts in remote partitions
  !>   2. Add `us` from subdomain core regions to `u`
  !>   3. Merge received remote data into ghosts
  !>   4. Merge `us` from local and ghost overlap regions into `u`
  !>
  !> The subdomain variable must be dimensioned as `u(ns,ns,ns,ne+ng)`, where
  !> `ns` is the number of subdomain points per direction, `ne` the number of
  !> local elements and `ng` the number of ghosts.

  module subroutine Merge_Unstructured(this, mesh, cfg, buf_us, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    integer,        intent(in)    :: cfg(:,:)      !< subdomain configurations
    real(RWP),      intent(inout) :: us(:,:,:,:)   !< subdomain solutions
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< mesh variable

    real(RNP) :: du(0:ubound(u,1), 0:ubound(u,1), 0:ubound(u,1))
    real(RNP) :: ue(0:ubound(u,1)), uv

    real(RNP) :: wa(0:ubound(u,1), 0:ubound(u,1), 0:ubound(u,1))
    real(RNP) :: ws(size(us,1), size(us,1), size(us,1))

    integer(IXS) :: orientation
    integer :: e, i, j, k, l, m, n
    integer :: ii, jj, kk, nn, np, ns
    integer :: c0, c1, lo, ro
    integer :: o(3), s(3,3), t(3)
    logical :: merge_edge(12), merge_vertex(8)

    associate(po => this%po, no => this%no)

      if (no <= 0) then
        call Warning('Merge_Unstructured', 'Expecting no > 0')
      end if

      ! prerequisites ..........................................................

      np = po + 1       ! number of element   points per direction
      ns = np + 2 * no  ! number of subdomain points per direction

      lo = no - 1       ! left overlap end
      ro = po - lo      ! right overlap start

      c0 = no + 1       ! first subdomain point inside core
      c1 = c0 + np      ! first subdomain point behind core

      ! get ghost corrections ..................................................

      call buf_us % Transfer(mesh, us, tag=2000)
      call buf_us % Merge(us)

      ! merge corrections ......................................................

      !$omp do
      do e = 1, mesh % n_elem_active
        associate(element => mesh % element)

          merge_edge   = this % edges
          merge_vertex = this % vertices

          ! subdomain weights
          l = cfg(1,e)
          m = cfg(2,e)
          n = cfg(3,e)
          select case(this % wp)
          case(RSP)
            associate(W => this % ops_sp % W)
              do k = 1, ns
              do j = 1, ns
              do i = 1, ns
                ws(i,j,k) = W(i,l) * W(j,m) * W(k,n)
              end do
              end do
              end do
            end associate
          case(RDP)
            associate(W => this % ops_dp % W)
              do k = 1, ns
              do j = 1, ns
              do i = 1, ns
                ws(i,j,k) = W(i,l) * W(j,m) * W(k,n)
              end do
              end do
              end do
            end associate
          end select

          ! core
          do k = 0, po
          do j = 0, po
          do i = 0, po
            du(i,j,k) = us(c0 + i, c0 + j, c0 + k, e)
            wa(i,j,k) = ws(c0 + i, c0 + j, c0 + k)
          end do
          end do
          end do

          ! face 1: ξ = -1 <--> west
          m = element(e) % face(1) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [ c1, c0, c0 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, po
            do j = 0, po
            do i = 0, lo
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
              wa(i, j, k) = wa(i, j, k) + ws(c0-1-i, c0+j, c0+k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,1)) = .false.
            merge_vertex (V_FACE(:,1)) = .false.
          end if

          ! face 2:  ξ = 1 <--> east
          m = element(e) % face(2) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [  1 , c0, c0 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, po
            do j = 0, po
            do i = 0, lo
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(ro+i, j, k) = du(ro+i, j, k) + us(ii, jj, kk, l)
              wa(ro+i, j, k) = wa(ro+i, j, k) + ws(ns-i, c0+j, c0+k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,2)) = .false.
            merge_vertex (V_FACE(:,2)) = .false.
          end if

          ! face 3: η = -1 <--> south
          m = element(e) % face(3) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [ c0 , c1, c0 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, po
            do j = 0, lo
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
              wa(i, j, k) = wa(i, j, k) + ws(c0+i, c0-1-j, c0+k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,3)) = .false.
            merge_vertex (V_FACE(:,3)) = .false.
          end if

          ! face 4: η = 1 <--> north
          m = element(e) % face(4) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [ c0,  1, c0 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, po
            do j = 0, lo
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(i, ro+j, k) = du(i, ro+j, k) + us(ii, jj, kk, l)
              wa(i, ro+j, k) = wa(i, ro+j, k) + ws(c0+i, ns-j, c0+k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,4)) = .false.
            merge_vertex (V_FACE(:,4)) = .false.
          end if

          ! face 5: ζ = -1 <--> bottom
          m = element(e) % face(5) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [ c0, c0, c1 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, lo
            do j = 0, po
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
              wa(i, j, k) = wa(i, j, k) + ws(c0+i, c0+j, c0-1-k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,5)) = .false.
            merge_vertex (V_FACE(:,5)) = .false.
          end if

          ! face 6: ζ = 1 <--> top
          m = element(e) % face(6) % i_neighbor
          if (m > 0) then

            orientation = element(e) % neighbor(m) % orientation
            l = element(e) % neighbor(m) % id
            o = [ c0, c0,  1 ]
            call TransformIndex(orientation, ns, 1, o, 1, t, s)
            do k = 0, lo
            do j = 0, po
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              du(i, j, ro+k) = du(i, j, ro+k) + us(ii, jj, kk, l)
              wa(i, j, ro+k) = wa(i, j, ro+k) + ws(c0+i, c0+j, ns-k)
            end do
            end do
            end do

          else
            merge_edge   (E_FACE(:,6)) = .false.
            merge_vertex (V_FACE(:,6)) = .false.
          end if

          ! edge  1: (x1  , south, bottom)
          nn = element(e) % edge( 1) % n_neighbor
          if (nn > 0 .and. merge_edge( 1)) then
            m = element(e) % edge( 1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c0, c1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
                wa(i, j, k) = wa(i, j, k) + ws(c0+i, c0-1-j, c0-1-k)
              end do
              end do
              end do

            else

              j = 0
              k = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 1, m, n, no, us, ue)
              do i = 0, po
                du(i, j, k) = du(i, j, k) + ue(i)
                wa(i, j, k) = wa(i, j, k) + ws(c0+i, c0-1, c0-1)
              end do

            end if
          end if

          ! edge  2: (x1  , north, bottom)
          nn = element(e) % edge( 2) % n_neighbor
          if (nn > 0 .and. merge_edge( 2)) then
            m = element(e) % edge( 2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c0,  1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, ro+j, k) = du(i, ro+j, k) + us(ii, jj, kk, l)
                wa(i, ro+j, k) = wa(i, ro+j, k) + ws(c0+i, ns-j, c0-1-k)
              end do
              end do
              end do

            else

              j = po
              k = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 2, m, n, no, us, ue)
              do i = 0, po
                du(i, j, k) = du(i, j, k) + ue(i)
                wa(i, j, k) = wa(i, j, k) + ws(c0+i, c1, c0-1)
              end do

            end if
          end if

          ! edge  3: (x1  , south, top   )
          nn = element(e) % edge( 3) % n_neighbor
          if (nn > 0 .and. merge_edge( 3)) then
            m = element(e) % edge( 3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c0, c1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, ro+k) = du(i, j, ro+k) + us(ii, jj, kk, l)
                wa(i, j, ro+k) = wa(i, j, ro+k) + ws(c0+i, c0-1-j, ns-k)
              end do
              end do
              end do

            else

              j = 0
              k = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 3, m, n, no, us, ue)
              do i = 0, po
                du(i, j, k) = du(i, j, k) + ue(i)
                wa(i, j, k) = wa(i, j, k) + ws(c0+i, c0-1, c1)
              end do

            end if
          end if

          ! edge  4: (x1  , north, top   )
          nn = element(e) % edge( 4) % n_neighbor
          if (nn > 0 .and. merge_edge( 4)) then
            m = element(e) % edge( 4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c0,  1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, ro+j, ro+k) = du(i, ro+j, ro+k) + us(ii, jj, kk, l)
                wa(i, ro+j, ro+k) = wa(i, ro+j, ro+k) + ws(c0+i, ns-j, ns-k)
              end do
              end do
              end do

            else

              j = po
              k = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 4, m, n, no, us, ue)
              do i = 0, po
                du(i, j, k) = du(i, j, k) + ue(i)
                wa(i, j, k) = wa(i, j, k) + ws(c0+i, c1, c1)
              end do

            end if
          end if

          ! edge  5: (west, x2   , bottom)
          nn = element(e) % edge( 5) % n_neighbor
          if (nn > 0 .and. merge_edge( 5)) then
            m = element(e) % edge( 5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1, c0, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1-i, c0+j, c0-1-k)
              end do
              end do
              end do

            else

              i = 0
              k = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 5, m, n, no, us, ue)
              do j = 0, po
                du(i, j, k) = du(i, j, k) + ue(j)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1, c0+j, c0-1)
              end do

            end if
          end if

          ! edge  6: (east, x2   , bottom)
          nn = element(e) % edge( 6) % n_neighbor
          if (nn > 0 .and. merge_edge( 6)) then
            m = element(e) % edge( 6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1, c0, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, j, k) = du(ro+i, j, k) + us(ii, jj, kk, l)
                wa(ro+i, j, k) = wa(ro+i, j, k) + ws(ns-i, c0+j, c0-1-k)
              end do
              end do
              end do

            else

              i = po
              k = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 6, m, n, no, us, ue)
              do j = 0, po
                du(i, j, k) = du(i, j, k) + ue(j)
                wa(i, j, k) = wa(i, j, k) + ws(c1, c0+j, c0-1)
              end do

            end if
          end if

          ! edge  7: (west, x2   , top   )
          nn = element(e) % edge( 7) % n_neighbor
          if (nn > 0 .and. merge_edge( 7)) then
            m = element(e) % edge( 7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1, c0,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, ro+k) = du(i, j, ro+k) + us(ii, jj, kk, l)
                wa(i, j, ro+k) = wa(i, j, ro+k) + ws(c0-1-i, c0+j, ns-k)
              end do
              end do
              end do

            else

              i = 0
              k = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 7, m, n, no, us, ue)
              do j = 0, po
                du(i, j, k) = du(i, j, k) + ue(j)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1, c0+j, c1)
              end do

            end if
          end if

          ! edge  8: (east, x2   , top   )
          nn = element(e) % edge( 8) % n_neighbor
          if (nn > 0 .and. merge_edge( 8)) then
            m = element(e) % edge( 8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1, c0,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, j, ro+k) = du(ro+i, j, ro+k) + us(ii, jj, kk, l)
                wa(ro+i, j, ro+k) = wa(ro+i, j, ro+k) + ws(ns-i, c0+j, ns-k)
              end do
              end do
              end do

            else

              i = po
              k = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 8, m, n, no, us, ue)
              do j = 0, po
                du(i, j, k) = du(i, j, k) + ue(j)
                wa(i, j, k) = wa(i, j, k) + ws(c1, c0+j, c1)
              end do

            end if
          end if

          ! edge  9: (west, south, x3    )
          nn = element(e) % edge( 9) % n_neighbor
          if (nn > 0 .and. merge_edge( 9)) then
            m = element(e) % edge( 9) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1, c1, c0 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1-i, c0-1-j, c0+k)
              end do
              end do
              end do

            else

              i = 0
              j = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 9, m, n, no, us, ue)
              do k = 0, po
                du(i, j, k) = du(i, j, k) + ue(k)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1, c0-1, c0+k)
              end do

            end if
          end if

          ! edge 10: (east, south, x3    )
          nn = element(e) % edge(10) % n_neighbor
          if (nn > 0 .and. merge_edge(10)) then
            m = element(e) % edge(10) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1, c1, c0 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, j, k) = du(ro+i, j, k) + us(ii, jj, kk, l)
                wa(ro+i, j, k) = wa(ro+i, j, k) + ws(ns-i, c0-1-j, c0+k)
              end do
              end do
              end do

            else

              i = po
              j = 0
              n = m + nn - 1
              call CollectEdgeData(element(e), 10, m, n, no, us, ue)
              do k = 0, po
                du(i, j, k) = du(i, j, k) + ue(k)
                wa(i, j, k) = wa(i, j, k) + ws(c1, c0-1, c0+k)
              end do

            end if
          end if

          ! edge 11: (west, north, x3    )
          nn = element(e) % edge(11) % n_neighbor
          if (nn > 0 .and. merge_edge(11)) then
            m = element(e) % edge(11) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1,  1, c0 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, ro+j, k) = du(i, ro+j, k) + us(ii, jj, kk, l)
                wa(i, ro+j, k) = wa(i, ro+j, k) + ws(c0-1-i, ns-j, c0+k)
              end do
              end do
              end do

            else

              i = 0
              j = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 11, m, n, no, us, ue)
              do k = 0, po
                du(i, j, k) = du(i, j, k) + ue(k)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1, c1, c0+k)
              end do

            end if
          end if

          ! edge 12: (east, north, x3    )
          nn = element(e) % edge(12) % n_neighbor
          if (nn > 0 .and. merge_edge(12)) then
            m = element(e) % edge(12) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1,  1, c0 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, ro+j, k) = du(ro+i, ro+j, k) + us(ii, jj, kk, l)
                wa(ro+i, ro+j, k) = wa(ro+i, ro+j, k) + ws(ns-i, ns-j, c0+k)
              end do
              end do
              end do

            else

              i = po
              j = po
              n = m + nn - 1
              call CollectEdgeData(element(e), 12, m, n, no, us, ue)
              do k = 0, po
                du(i, j, k) = du(i, j, k) + ue(k)
                wa(i, j, k) = wa(i, j, k) + ws(c1, c1, c0+k)
              end do

            end if
          end if

          ! vertex 1: (west, south, bottom)
          nn = element(e) % vertex(1) % n_neighbor
          if (nn > 0 .and. merge_vertex(1)) then
            m = element(e) % vertex(1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1, c1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, k) = du(i, j, k) + us(ii, jj, kk, l)
                wa(i, j, k) = wa(i, j, k) + ws(c0-1-i, c0-1-j, c0-1-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du(  0,  0,  0) = du(  0,  0,  0) + uv
              wa(  0,  0,  0) = wa(  0,  0,  0) + ws(c0-1, c0-1, c0-1)

            end if
          end if

          ! vertex 2: (east, south, bottom)
          nn = element(e) % vertex(2) % n_neighbor
          if (nn > 0 .and. merge_vertex(2)) then
            m = element(e) % vertex(2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1, c1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, j, k) = du(ro+i, j, k) + us(ii, jj, kk, l)
                wa(ro+i, j, k) = wa(ro+i, j, k) + ws(ns-i, c0-1-j, c0-1-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du( po,  0,  0) = du( po,  0,  0) + uv
              wa( po,  0,  0) = wa( po,  0,  0) + ws(  c1, c0-1, c0-1)

            end if
          end if

          ! vertex 3: (west, north, bottom)
          nn = element(e) % vertex(3) % n_neighbor
          if (nn > 0 .and. merge_vertex(3)) then
            m = element(e) % vertex(3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1,  1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, ro+j, k) = du(i, ro+j, k) + us(ii, jj, kk, l)
                wa(i, ro+j, k) = wa(i, ro+j, k) + ws(c0-1-i, ns-j, c0-1-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du(  0, po,  0) = du(  0, po,  0) + uv
              wa(  0, po,  0) = wa(  0, po,  0) + ws(c0-1,   c1, c0-1)

            end if
          end if

          ! vertex 4: (east, north, bottom)
          nn = element(e) % vertex(4) % n_neighbor
          if (nn > 0 .and. merge_vertex(4)) then
            m = element(e) % vertex(4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1,  1, c1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, ro+j, k) = du(ro+i, ro+j, k) + us(ii, jj, kk, l)
                wa(ro+i, ro+j, k) = wa(ro+i, ro+j, k) + ws(ns-i, ns-j, c0-1-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du( po, po,  0) = du( po, po,  0) + uv
              wa( po, po,  0) = wa( po, po,  0) + ws(  c1,   c1, c0-1)

            end if
          end if

          ! vertex 5: (west, south, top   )
          nn = element(e) % vertex(5) % n_neighbor
          if (nn > 0 .and. merge_vertex(5)) then
            m = element(e) % vertex(5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1, c1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, j, ro+k) = du(i, j, ro+k) + us(ii, jj, kk, l)
                wa(i, j, ro+k) = wa(i, j, ro+k) + ws(c0-1-i, c0-1-j, ns-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du(  0,  0, po) = du(  0,  0, po) + uv
              wa(  0,  0, po) = wa(  0,  0, po) + ws(c0-1, c0-1,   c1)

            end if
          end if

          ! vertex 6: (east, south, top   )
          nn = element(e) % vertex(6) % n_neighbor
          if (nn > 0 .and. merge_vertex(6)) then
            m = element(e) % vertex(6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1, c1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, j, ro+k) = du(ro+i, j, ro+k) + us(ii, jj, kk, l)
                wa(ro+i, j, ro+k) = wa(ro+i, j, ro+k) + ws(ns-i, c0-1-j, ns-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du( po,  0, po) = du( po,  0, po) + uv
              wa( po,  0, po) = wa( po,  0, po) + ws(  c1, c0-1,   c1)

            end if
          end if

          ! vertex 7: (west, north, top   )
          nn = element(e) % vertex(7) % n_neighbor
          if (nn > 0 .and. merge_vertex(7)) then
            m = element(e) % vertex(7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ c1,  1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(i, ro+j, ro+k) = du(i, ro+j, ro+k) + us(ii, jj, kk, l)
                wa(i, ro+j, ro+k) = wa(i, ro+j, ro+k) + ws(c0-1-i, ns-j, ns-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du(  0, po, po) = du(  0, po, po) + uv
              wa(  0, po, po) = wa(  0, po, po) + ws(c0-1,   c1,   c1)

            end if
          end if

          ! vertex 8: (east, north, top   )
          nn = element(e) % vertex(8) % n_neighbor
          if (nn > 0 .and. merge_vertex(8)) then
            m = element(e) % vertex(8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  1,  1,  1 ]
              call TransformIndex(orientation, ns, 1, o, 1, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                du(ro+i, ro+j, ro+k) = du(ro+i, ro+j, ro+k) + us(ii, jj, kk, l)
                wa(ro+i, ro+j, ro+k) = wa(ro+i, ro+j, ro+k) + ws(ns-i, ns-j, ns-k)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              du( po, po, po) = du( po, po, po) + uv
              wa( po, po, po) = wa( po, po, po) + ws(  c1,   c1,   c1)

            end if
          end if

          ! add correction
          if (this % renormalize) then
            u(:,:,:,e) = u(:,:,:,e) + du / wa
          else
            u(:,:,:,e) = u(:,:,:,e) + du
          end if

        end associate
      end do

    end associate

  end subroutine Merge_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(element, e, m, n, no, us, ue)
    class(MeshElement_3D), intent(in)  :: element     !< current element
    integer,               intent(in)  :: e           !< element edge ID
    integer,               intent(in)  :: m           !< first neighbor ID
    integer,               intent(in)  :: n           !< last neighbor ID
    integer,               intent(in)  :: no          !< overlapped layers
    real(RWP),             intent(in)  :: us(:,:,:,:) !< element data
    real(RNP),             intent(out) :: ue(:)       !< avg neighbor edge data

    real(RNP) :: uk(size(ue))
    integer   :: i0, ip, np
    integer   :: k, l

    np = size(ue)
    i0 = no + 1
    ip = no + np
    ue = 0

    do k = m, n

      l = element % neighbor(k) % id

      select case(ElementEdgeID(element % neighbor(k) % component))
      case( 1)
        uk = us( i0:ip, i0, i0, l)
      case( 2)
        uk = us( i0:ip, ip, i0, l)
      case( 3)
        uk = us( i0:ip, i0, ip, l)
      case( 4)
        uk = us( i0:ip, ip, ip, l)
      case( 5)
        uk = us( i0, i0:ip, i0, l)
      case( 6)
        uk = us( ip, i0:ip, i0, l)
      case( 7)
        uk = us( i0, i0:ip, ip, l)
      case( 8)
        uk = us( ip, i0:ip, ip, l)
      case( 9)
        uk = us( i0, i0, i0:ip, l)
      case(10)
        uk = us( ip, i0, i0:ip, l)
      case(11)
        uk = us( i0, ip, i0:ip, l)
      case(12)
        uk = us( ip, ip, i0:ip, l)
      end select

      if (element % NeigborEdgeIsAligned(e, k)) then
        ue = ue + uk
      else
        ue = ue + uk(np:1:-1)
      end if

    end do

    ue = ue / (n - m + 1)

  end subroutine CollectEdgeData

  !-----------------------------------------------------------------------------
  !> Collect and average neighbor vertex data

  pure subroutine CollectVertexData(neighbor, no, us, uv)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    integer,                       intent(in)  :: no
    real(RWP),                     intent(in)  :: us(:,:,:,:)
    real(RNP),                     intent(out) :: uv

    integer :: i0, ip, nn
    integer :: k, l, n

    nn = size(neighbor)
    i0 = no + 1
    ip = size(us,1) - no

    uv = 0

    do n = 1, nn

      l = neighbor(n) % id
      k = neighbor(n) % component - 18

      select case(k)
      case( 1)
        uv = uv + us( i0, i0, i0, l)
      case( 2)
        uv = uv + us( ip, i0, i0, l)
      case( 3)
        uv = uv + us( i0, ip, i0, l)
      case( 4)
        uv = uv + us( ip, ip, i0, l)
      case( 5)
        uv = uv + us( i0, i0, ip, l)
      case( 6)
        uv = uv + us( ip, i0, ip, l)
      case( 7)
        uv = uv + us( i0, ip, ip, l)
      case( 8)
        uv = uv + us( ip, ip, ip, l)
      end select

    end do

    uv = ONE/nn * uv

  end subroutine CollectVertexData
