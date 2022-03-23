  !-----------------------------------------------------------------------------
  !> Restrict an unstructured mesh variable `u` to subdomains
  !>
  !> The strategy is as a follows:
  !>
  !>   1. Send `u` from linked local elements to ghosts in remote partitions
  !>   2. Assign `u` to subdomain core regions, i.e. 1:1 copy to `us`
  !>   3. Merge received remote data into ghosts
  !>   4. Copy local and ghost data to overlap zones of subdomains
  !>
  !> The mesh variable must be dimensioned as `u(0:po,0:po,0:po,ne+ng)`, where
  !> `ne` is the number of local elements and `ng` the number of ghosts.

  module subroutine Restrict_Unstructured(this, mesh, buf_u, u, us)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_u
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: us(:,:,:,:)   !< restricted variable

    real(RNP) :: ue(size(u,1)), uv

    integer(IXS) :: orientation
    integer :: e, i, j, k, l, m, n
    integer :: lo, ro, nn, np
    integer :: ii, jj, kk
    integer :: c0, c1, o(3), s(3,3), t(3)
    integer :: fill_edge(12), fill_vertex(8)

    associate(po => this%po, no => this%no)

      ! start transfer .........................................................

      np = po + 1       ! number of element   points per direction

      lo = no - 1       ! last  point of left  element overlap zone
      ro = po - lo      ! first point of right element overlap zone

      c0 = no + 1       ! first subdomain point inside core
      c1 = c0 + np      ! first subdomain point behind core

      if (no > 0) then
        call buf_u % Transfer(mesh, u, tag=1000)
      end if

      ! assign core regions ....................................................

      !$omp do
      do e = 1, mesh % n_elem

        ! initialization and removal of NANs
        us(:,:,:,e) = 0

        do k = 0, po
        do j = 0, po
        do i = 0, po
          us(c0 + i, c0 + j, c0 + k, e) = real(u(i,j,k,e), RWP)
        end do
        end do
        end do

      end do

      ! merge remote data to ghosts ............................................

      if (no > 0) then
        call buf_u % Merge(u)
      else
        return
      end if

      ! fill overlap regions ...................................................

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element)

          fill_edge   = 2 ! preset edge   overlap zones to complete autofill
          fill_vertex = 2 ! preset vertex overlap zones to complete autofill

          ! face 1: ξ = -1 <--> west
          m = element(e) % face(1) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [ ro,  0,  0 ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, po
            do j = 0, po
            do i = 0, lo
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(1 + i, c0 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,1)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,1)) = 0 ! behind boundaries
          end if

          ! face 2:  ξ = 1 <--> east
          m = element(e) % face(2) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [  0,  0,  0 ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, po
            do j = 0, po
            do i = 0, lo
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(c1 + i, c0 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,2)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,2)) = 0 ! behind boundaries
          end if

          ! face 3: η = -1 <--> south
          m = element(e) % face(3) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [  0, ro,  0 ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, po
            do j = 0, lo
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(c0 + i, 1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,3)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,3)) = 0 ! behind boundaries
          end if

          ! face 4: η = 1 <--> north
          m = element(e) % face(4) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [  0,  0,  0 ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, po
            do j = 0, lo
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(c0 + i, c1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,4)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,4)) = 0 ! behind boundaries
          end if

          ! face 5: ζ = -1 <--> bottom
          m = element(e) % face(5) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [  0,  0, ro ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, lo
            do j = 0, po
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(c0 + i, c0 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,5)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,5)) = 0 ! behind boundaries
          end if

          ! face 6: ζ = 1 <--> top
          m = element(e) % face(6) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            o = [  0,  0,  0 ]
            call TransformIndex(orientation, np, 0, o, 0, t, s)
            do k = 0, lo
            do j = 0, po
            do i = 0, po
              ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
              jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
              kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
              us(c0 + i, c0 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
            end do
            end do
            end do

          else
            fill_edge   (E_FACE(:,6)) = 0 ! disable autofill
            fill_vertex (V_FACE(:,6)) = 0 ! behind boundaries
         end if

          ! edge  1: (x1  , south, bottom)
          nn = element(e) % edge( 1) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0, ro, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c0 + i, 1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(1) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 1) % orientation, &
                                    u, ue                                )
              us(c0 : c0 + po, no, no, e) = real(ue, RWP)
              fill_edge(1) = min(fill_edge(1), 1)

            end if
          end if

          ! edge  2: (x1  , north, bottom)
          nn = element(e) % edge( 2) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c0 + i, c1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(2) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 2) % orientation, &
                                    u, ue                                )
              us(c0 : c0 + po, c1 + 1, no, e) = real(ue, RWP)
              fill_edge(2) = min(fill_edge(2), 1)

            end if
          end if

          ! edge  3: (x1  , south, top   )
          nn = element(e) % edge( 3) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0, ro,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c0 + i, 1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(3) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 3) % orientation, &
                                    u, ue                                )
              us(c0 : c0 + po, no, c1 + 1, e) = real(ue, RWP)
              fill_edge(3) = min(fill_edge(3), 1)

            end if
          end if

          ! edge  4: (x1  , north, top   )
          nn = element(e) % edge( 4) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, po
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c0 + i, c1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(4) = 0

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 4) % orientation, &
                                    u, ue                                )
              us(c0 : c0 + po, c1 + 1, c1 + 1, e) = real(ue, RWP)
              fill_edge(4) = min(fill_edge(4), 1)

            end if
          end if

          ! edge  5: (west, x2   , bottom)
          nn = element(e) % edge( 5) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro,  0, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, c0 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(5) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 5) % orientation, &
                                    u, ue                                )
              us(no, c0 : c0 + po, no, e) = real(ue, RWP)
              fill_edge(5) = min(fill_edge(5), 1)

            end if
          end if

          ! edge  6: (east, x2   , bottom)
          nn = element(e) % edge( 6) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, c0 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(6) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 6) % orientation, &
                                    u, ue                                )
              us(c1 + 1, c0 : c0 + po, no, e) = real(ue, RWP)
              fill_edge(6) = min(fill_edge(6), 1)

            end if
          end if

          ! edge  7: (west, x2   , top   )
          nn = element(e) % edge( 7) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, c0 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(7) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 7) % orientation, &
                                    u, ue                                )
              us(no, c0 : c0 + po, c1 + 1, e) = real(ue, RWP)
              fill_edge(7) = min(fill_edge(7), 1)

            end if
          end if

          ! edge  8: (east, x2   , top   )
          nn = element(e) % edge( 8) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, po
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, c0 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(8) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 8) % orientation, &
                                    u, ue                                )
              us(c1 + 1, c0 : c0 + po, c1 + 1, e) = real(ue, RWP)
              fill_edge(8) = min(fill_edge(8), 1)

            end if
          end if

          ! edge  9: (west, south, x3    )
          nn = element(e) % edge( 9) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 9) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro, ro,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, 1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(9) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 9) % orientation, &
                                    u, ue                                )
              us(no, no, c0 : c0 + po, e) = real(ue, RWP)
              fill_edge(9) = min(fill_edge(9), 1)

            end if
          end if

          ! edge 10: (east, south, x3    )
          nn = element(e) % edge(10) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(10) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0, ro,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, 1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(10) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(10) % orientation, &
                                    u, ue                                )
              us(c1 + 1, no, c0 : c0 + po, e) = real(ue, RWP)
              fill_edge(10) = min(fill_edge(10), 1)

            end if
          end if

          ! edge 11: (west, north, x3    )
          nn = element(e) % edge(11) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(11) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, c1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(11) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(11) % orientation, &
                                    u, ue                                )
              us(no, c1 + 1, c0 : c0 + po, e) = real(ue, RWP)
              fill_edge(11) = min(fill_edge(11), 1)

            end if
          end if

          ! edge 12: (east, north, x3    )
          nn = element(e) % edge(12) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(12) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, po
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, c1 + j, c0 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(12) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(12) % orientation, &
                                    u, ue                                )
              us(c1 + 1, c1 + 1, c0 : c0 + po, e) = real(ue, RWP)
              fill_edge(12) = min(fill_edge(12), 1)

            end if
          end if

          ! vertex 1: (west, south, bottom)
          nn = element(e) % vertex(1) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro, ro, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, 1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(1) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(no, no, no, e) = real(uv, RWP)
              fill_vertex(1) = min(fill_vertex(1), 1)

            end if
          end if

          ! vertex 2: (east, south, bottom)
          nn = element(e) % vertex(2) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0, ro, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, 1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(2) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(c1 + 1, no, no, e) = real(uv, RWP)
              fill_vertex(2) = min(fill_vertex(2), 1)

            end if
          end if

          ! vertex 3: (west, north, bottom)
          nn = element(e) % vertex(3) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro,  0, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, c1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(3) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(no, c1 + 1, no, e) = real(uv, RWP)
              fill_vertex(3) = min(fill_vertex(3), 1)

            end if
          end if

          ! vertex 4: (east, north, bottom)
          nn = element(e) % vertex(4) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0, ro ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, c1 + j, 1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(4) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(c1 + 1, c1 + 1, no, e) = real(uv, RWP)
              fill_vertex(4) = min(fill_vertex(4), 1)

            end if
          end if

          ! vertex 5: (west, south, top   )
          nn = element(e) % vertex(5) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro, ro,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, 1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(5) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(no, no, c1 + 1, e) = real(uv, RWP)
              fill_vertex(5) = min(fill_vertex(5), 1)

            end if
          end if

          ! vertex 6: (east, south, top   )
          nn = element(e) % vertex(6) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0, ro,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, 1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(6) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(c1 + 1, no, c1 + 1, e) = real(uv, RWP)
              fill_vertex(6) = min(fill_vertex(6), 1)

            end if
          end if

          ! vertex 7: (west, north, top   )
          nn = element(e) % vertex(7) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [ ro,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(1 + i, c1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(7) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(no, c1 + 1, c1 + 1, e) = real(uv, RWP)
              fill_vertex(7) = min(fill_vertex(7), 1)

            end if
          end if

          ! vertex 8: (east, north, top   )
          nn = element(e) % vertex(8) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              l = element(e) % neighbor(m) % id
              o = [  0,  0,  0 ]
              call TransformIndex(orientation, np, 0, o, 0, t, s)
              do k = 0, lo
              do j = 0, lo
              do i = 0, lo
                ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                us(c1 + i, c1 + j, c1 + k, e) = real(u(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(8) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), u, uv)
              us(c1 + 1, c1 + 1, c1 + 1, e) = real(uv, RWP)
              fill_vertex(8) = min(fill_vertex(8), 1)

            end if
          end if

          ! fill irregular overlap zones if requested
          if (this % fill) then
            call FillIrregularZones(no, fill_edge, fill_vertex, us(:,:,:,e))
          end if

        end associate
      end do

    end associate

  end subroutine Restrict_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(neighbor, element, o, u, ue)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    class(MeshElement_3D),         intent(in)  :: element(:)
    integer(IXS),                  intent(in)  :: o
    real(RNP),                     intent(in)  :: u(:,:,:,:)
    real(RNP),                     intent(out) :: ue(:)

    real(RNP) :: un(size(ue))
    integer :: nn, np
    integer :: e, l, n

    nn = size(neighbor)
    np = size(ue)

    ue = 0

    do n = 1, nn

      l = neighbor(n) % id
      e = neighbor(n) % component - 6

      select case(e)
      case( 1)
        un = u( :, 1, 1,l)
      case( 2)
        un = u( :,np, 1,l)
      case( 3)
        un = u( :, 1,np,l)
      case( 4)
        un = u( :,np,np,l)
      case( 5)
        un = u( 1, :, 1,l)
      case( 6)
        un = u(np, :, 1,l)
      case( 7)
        un = u( 1, :,np,l)
      case( 8)
        un = u(np, :,np,l)
      case( 9)
        un = u( 1, 1, :,l)
      case(10)
        un = u(np, 1, :,l)
      case(11)
        un = u( 1,np, :,l)
      case(12)
        un = u(np,np, :,l)
      end select

      if (element(l) % edge(e) % orientation == o) then
        ue = ue + un
      else
        ue = ue + un(np:1:-1)
      end if

    end do

    ue = ONE/nn * ue

  end subroutine CollectEdgeData

  !-----------------------------------------------------------------------------
  !> Collect and average neighbor vertex data

  pure subroutine CollectVertexData(neighbor, u, uv)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    real(RNP),                     intent(in)  :: u(:,:,:,:)
    real(RNP),                     intent(out) :: uv

    integer :: nn, np
    integer :: k, l, n

    nn = size(neighbor)
    np = size(u,1)

    uv = 0

    do n = 1, nn

      l = neighbor(n) % id
      k = neighbor(n) % component - 18

      select case(k)
      case( 1)
        uv = uv + u( 1, 1, 1,l)
      case( 2)
        uv = uv + u(np, 1, 1,l)
      case( 3)
        uv = uv + u( 1,np, 1,l)
      case( 4)
        uv = uv + u(np,np, 1,l)
      case( 5)
        uv = uv + u( 1, 1,np,l)
      case( 6)
        uv = uv + u(np, 1,np,l)
      case( 7)
        uv = uv + u( 1,np,np,l)
      case( 8)
        uv = uv + u(np,np,np,l)
      end select

    end do

    uv = ONE/nn * uv

  end subroutine CollectVertexData

  !-----------------------------------------------------------------------------
  !> Fill irregular overlap zones by interpolation

  pure subroutine FillIrregularZones(no, fill_edge, fill_vertex, us)
    integer, intent(in) :: no             !< overlap
    integer, intent(in) :: fill_edge(12)  !< edge   type  {0,1,2}
    integer, intent(in) :: fill_vertex(8) !< vertex type  {0,1,2}
    real(RWP), contiguous, intent(inout) :: us(:,:,:) !< subdomain variable

    integer   :: lc, lo, rc, ro, ns
    integer   :: i, i0, j, j0, k, k0
    real(RWP) :: c1, c2, c3

    if (no < 1) return

    ns = size(us,1)
    lo = no         ! end of left overlap zone
    lc = no + 1     ! leftmost  core point
    rc = ns - no    ! rightmost core point
    ro = rc + 1     ! start of right overlap zone

    EDGE_1: if (fill_edge( 1) > 0) then
      ! adopt from core
      if (fill_edge( 1) > 1) us(lc:rc,lo,lo) = us(lc:rc,lc,lc)
      if (no == 1) exit EDGE_1
      ! adopt from face overlap zones
      us(lc:rc,  1: lo-1,     lo  ) = us(lc:rc,  1: lo-1,     lc  )
      us(lc:rc,     lo  ,  1: lo-1) = us(lc:rc,     lc  ,  1: lo-1)
      ! interpolate remaining points
      j0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        us(i,j,k) = (c2 * us(i,j,k0) + c3 * us(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_1

    EDGE_2: if (fill_edge( 2) > 0) then
      ! adopt from core
      if (fill_edge( 2) > 1) us(lc:rc,ro,lo) = us(lc:rc,rc,lc)
      if (no == 1) exit EDGE_2
      ! adopt from face overlap zones
      us(lc:rc, ro+1: ns,     lo  ) = us(lc:rc, ro+1: ns,     lc  )
      us(lc:rc, ro      ,  1: lo-1) = us(lc:rc, rc      ,  1: lo-1)
      ! interpolate remaining points
      j0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        us(i,j,k) = (c2 * us(i,j,k0) + c3 * us(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_2

    EDGE_3: if (fill_edge( 3) > 0) then
      ! adopt from core
      if (fill_edge( 3) > 1) us(lc:rc,lo,ro) = us(lc:rc,lc,rc)
      if (no == 1) exit EDGE_3
      ! adopt from face overlap zones
      us(lc:rc,  1: lo-1, ro      ) = us(lc:rc,  1: lo-1, rc      )
      us(lc:rc,     lo  , ro+1: ns) = us(lc:rc,     lc  , ro+1: ns)
      ! interpolate remaining points
      j0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        us(i,j,k) = (c2 * us(i,j,k0) + c3 * us(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_3

    EDGE_4: if (fill_edge( 4) > 0) then
      ! adopt from core
      if (fill_edge( 4) > 1) us(lc:rc,ro,ro) = us(lc:rc,rc,rc)
      if (no == 1) exit EDGE_4
      ! adopt from face overlap zones
      us(lc:rc, ro+1: ns, ro      ) = us(lc:rc, ro+1: ns, rc      )
      us(lc:rc, ro      , ro+1: ns) = us(lc:rc, rc      , ro+1: ns)
      ! interpolate remaining points
      j0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        us(i,j,k) = (c2 * us(i,j,k0) + c3 * us(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_4

    EDGE_5: if (fill_edge( 5) > 0) then
      ! adopt from core
      if (fill_edge( 5) > 1) us(lo,lc:rc,lo) = us(lc,lc:rc,lc)
      if (no == 1) exit EDGE_5
      ! adopt from face overlap zones
      us( 1: lo-1, lc:rc,     lo  ) = us( 1: lo-1, lc:rc,     lc  )
      us(    lo  , lc:rc,  1: lo-1) = us(    lc  , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j,k0) + c3 * us(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_5

    EDGE_6: if (fill_edge( 6) > 0) then
      ! adopt from core
      if (fill_edge( 6) > 1) us(ro,lc:rc,lo) = us(rc,lc:rc,lc)
      if (no == 1) exit EDGE_6
      ! adopt from face overlap zones
      us(ro+1: ns, lc:rc,     lo  ) = us(ro+1: ns, lc:rc,     lc  )
      us(ro      , lc:rc,  1: lo-1) = us(rc      , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j,k0) + c3 * us(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_6

    EDGE_7: if (fill_edge( 7) > 0) then
      ! adopt from core
      if (fill_edge( 7) > 1) us(lo,lc:rc,ro) = us(lc,lc:rc,rc)
      if (no == 1) exit EDGE_7
      ! adopt from face overlap zones
      us( 1: lo-1, lc:rc, ro      ) = us( 1: lo-1, lc:rc, rc      )
      us(    lo  , lc:rc, ro+1: ns) = us(    lc  , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j,k0) + c3 * us(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_7

    EDGE_8: if (fill_edge( 8) > 0) then
      ! adopt from core
      if (fill_edge( 8) > 1) us(ro,lc:rc,ro) = us(rc,lc:rc,rc)
      if (no == 1) exit EDGE_8
      ! adopt from face overlap zones
      us(ro+1: ns, lc:rc, ro      ) = us(ro+1: ns, lc:rc, rc      )
      us(ro      , lc:rc, ro+1: ns) = us(rc      , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j,k0) + c3 * us(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_8

    EDGE_9: if (fill_edge( 9) > 0) then
      ! adopt from core
      if (fill_edge( 9) > 1) us(lo,lo,lc:rc) = us(lc,lc,lc:rc)
      if (no == 1) exit EDGE_9
      ! adopt from face overlap zones
      us( 1: lo-1,     lo  , lc:rc) = us( 1: lo-1,     lc  , lc:rc)
      us(    lo  ,  1: lo-1, lc:rc) = us(    lc  ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j0,k) + c2 * us(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_9

    EDGE_10: if (fill_edge(10) > 0) then
      ! adopt from core
      if (fill_edge(10) > 1) us(ro,lo,lc:rc) = us(rc,lc,lc:rc)
      if (no == 1) exit EDGE_10
      ! adopt from face overlap zones
      us(ro+1: ns,     lo  , lc:rc) = us(ro+1: ns,     lc  , lc:rc)
      us(ro      ,  1: lo-1, lc:rc) = us(rc      ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j0,k) + c2 * us(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_10

    EDGE_11: if (fill_edge(11) > 0) then
      ! adopt from core
      if (fill_edge(11) > 1) us(lo,ro,lc:rc) = us(lc,rc,lc:rc)
      if (no == 1) exit EDGE_11
      ! adopt from face overlap zones
      us( 1: lo-1, ro      , lc:rc) = us( 1: lo-1, rc      , lc:rc)
      us(    lo  , ro+1: ns, lc:rc) = us(    lc  , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j0,k) + c2 * us(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_11

    EDGE_12: if (fill_edge(12) > 0) then
      ! adopt from core
      if (fill_edge(12) > 1) us(ro,ro,lc:rc) = us(rc,rc,lc:rc)
      if (no == 1) exit EDGE_12
      ! adopt from face overlap zones
      us(ro+1: ns, ro      , lc:rc) = us(ro+1: ns, rc      , lc:rc)
      us(ro      , ro+1: ns, lc:rc) = us(rc      , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        us(i,j,k) = (c1 * us(i,j0,k) + c2 * us(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_12

    VERTEX_1: if (fill_vertex(1) > 0) then
      ! adopt from core
      if (fill_vertex(1) > 1) us(lo,lo,lo) = us(lc,lc,lc)
      if (no == 1) exit VERTEX_1
      ! adopt from face overlap zones
      us(  1:lo-1,    lo  ,    lo  ) = us(  1:lo-1,    lc  ,    lc  )
      us(    lo  ,  1:lo-1,    lo  ) = us(    lc  ,  1:lo-1,    lc  )
      us(    lo  ,    lo  ,  1:lo-1) = us(    lc  ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      us(  1:lo-1,  1:lo-1,    lo  ) = us(  1:lo-1,  1:lo-1,    lc  )
      us(  1:lo-1,    lo  ,  1:lo-1) = us(  1:lo-1,    lc  ,  1:lo-1)
      us(    lo  ,  1:lo-1,  1:lo-1) = us(    lc  ,  1:lo-1,  1:lo-1)
      ! interpolate remaining points
      i0 = lo
      j0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lo-1,  1, -1
      do i = lo-1,  1, -1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_1

    VERTEX_2: if (fill_vertex(2) > 0) then
      ! adopt from core
      if (fill_vertex(2) > 1) us(ro,lo,lo) = us(rc,lc,lc)
      if (no == 1) exit VERTEX_2
      ! adopt from face overlap zones
      us( ro+1:ns,    lo  ,    lo  ) = us( ro+1:ns,    lc  ,    lc  )
      us( ro     ,  1:lo-1,    lo  ) = us( rc     ,  1:lo-1,    lc  )
      us( ro     ,    lo  ,  1:lo-1) = us( rc     ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      us( ro+1:ns,  1:lo-1,    lo  ) = us( ro+1:ns,  1:lo-1,    lc  )
      us( ro+1:ns,    lo  ,  1:lo-1) = us( ro+1:ns,    lc  ,  1:lo-1)
      us( ro     ,  1:lo-1,  1:lo-1) = us( rc     ,  1:lo-1,  1:lo-1)
      ! interpolate remaining points
      i0 = ro
      j0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lo-1,  1, -1
      do i = ro+1, ns,  1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_2

    VERTEX_3: if (fill_vertex(3) > 0) then
      ! adopt from core
      if (fill_vertex(3) > 1) us(lo,ro,lo) = us(lc,rc,lc)
      if (no == 1) exit VERTEX_3
      ! adopt from face overlap zones
      us(    lo  , ro+1:ns,    lo  ) = us(    lc  , ro+1:ns,    lc  )
      us(  1:lo-1, ro     ,    lo  ) = us(  1:lo-1, rc     ,    lc  )
      us(    lo  , ro     ,  1:lo-1) = us(    lc  , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      us(  1:lo-1, ro+1:ns,    lo  ) = us(  1:lo-1, ro+1:ns,    lc  )
      us(  1:lo-1, ro     ,  1:lo-1) = us(  1:lo-1, rc     ,  1:lo-1)
      us(    lo  , ro+1:ns,  1:lo-1) = us(    lc  , ro+1:ns,  1:lo-1)
      ! interpolate remaining points
      i0 = lo
      j0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = ro+1, ns,  1
      do i = lo-1,  1, -1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_3

    VERTEX_4: if (fill_vertex(4) > 0) then
      ! adopt from core
      if (fill_vertex(4) > 1) us(ro,ro,lo) = us(rc,rc,lc)
      if (no == 1) exit VERTEX_4
      ! adopt from face overlap zones
      us( ro     , ro+1:ns,    lo  ) = us( rc     , ro+1:ns,    lc  )
      us( ro+1:ns, ro     ,    lo  ) = us( ro+1:ns, rc     ,    lc  )
      us( ro     , ro     ,  1:lo-1) = us( rc     , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      us( ro+1:ns, ro+1:ns,    lo  ) = us( ro+1:ns, ro+1:ns,    lc  )
      us( ro+1:ns, ro     ,  1:lo-1) = us( ro+1:ns, rc     ,  1:lo-1)
      us( ro     , ro+1:ns,  1:lo-1) = us( rc     , ro+1:ns,  1:lo-1)
      ! interpolate remaining points
      i0 = ro
      j0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = ro+1, ns,  1
      do i = ro+1, ns,  1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_4

    VERTEX_5: if (fill_vertex(5) > 0) then
      ! adopt from core
      if (fill_vertex(5) > 1) us(lo,lo,ro) = us(lc,lc,rc)
      if (no == 1) exit VERTEX_5
      ! adopt from face overlap zones
      us(  1:lo-1,    lo  , ro     ) = us(  1:lo-1,    lc  , rc     )
      us(    lo  ,  1:lo-1, ro     ) = us(    lc  ,  1:lo-1, rc     )
      us(    lo  ,    lo  , ro+1:ns) = us(    lc  ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      us(  1:lo-1,  1:lo-1, ro     ) = us(  1:lo-1,  1:lo-1, rc     )
      us(  1:lo-1,    lo  , ro+1:ns) = us(  1:lo-1,    lc  , ro+1:ns)
      us(    lo  ,  1:lo-1, ro+1:ns) = us(    lc  ,  1:lo-1, ro+1:ns)
      ! interpolate remaining points
      i0 = lo
      j0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lo-1,  1, -1
      do i = lo-1,  1, -1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_5

    VERTEX_6: if (fill_vertex(6) > 0) then
      ! adopt from core
      if (fill_vertex(6) > 1) us(ro,lo,ro) = us(rc,lc,rc)
      if (no == 1) exit VERTEX_6
      ! adopt from face overlap zones
      us( ro+1:ns,    lo  , ro     ) = us( ro+1:ns,    lc  , rc     )
      us( ro     ,  1:lo-1, ro     ) = us( rc     ,  1:lo-1, rc     )
      us( ro     ,    lo  , ro+1:ns) = us( rc     ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      us( ro+1:ns,  1:lo-1, ro     ) = us( ro+1:ns,  1:lo-1, rc     )
      us( ro+1:ns,    lo  , ro+1:ns) = us( ro+1:ns,    lc  , ro+1:ns)
      us( ro     ,  1:lo-1, ro+1:ns) = us( rc     ,  1:lo-1, ro+1:ns)
      ! interpolate remaining points
      i0 = ro
      j0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lo-1,  1, -1
      do i = ro+1, ns,  1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_6

    VERTEX_7: if (fill_vertex(7) > 0) then
      ! adopt from core
      if (fill_vertex(7) > 1) us(lo,ro,ro) = us(lc,rc,rc)
      if (no == 1) exit VERTEX_7
      ! adopt from face overlap zones
      us(    lo  , ro+1:ns, ro     ) = us(    lc  , ro+1:ns, rc     )
      us(  1:lo-1, ro     , ro     ) = us(  1:lo-1, rc     , rc     )
      us(    lo  , ro     , ro+1:ns) = us(    lc  , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      us(  1:lo-1, ro+1:ns, ro     ) = us(  1:lo-1, ro+1:ns, rc     )
      us(  1:lo-1, ro     , ro+1:ns) = us(  1:lo-1, rc     , ro+1:ns)
      us(    lo  , ro+1:ns, ro+1:ns) = us(    lc  , ro+1:ns, ro+1:ns)
      ! interpolate remaining points
      i0 = lo
      j0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = ro+1, ns,  1
      do i = lo-1,  1, -1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_7

    VERTEX_8: if (fill_vertex(8) > 0) then
      ! adopt from core
      if (fill_vertex(8) > 1) us(ro,ro,ro) = us(rc,rc,rc)
      if (no == 1) exit VERTEX_8
      ! adopt from face overlap zones
      us( ro     , ro+1:ns, ro     ) = us( rc     , ro+1:ns, rc     )
      us( ro+1:ns, ro     , ro     ) = us( ro+1:ns, rc     , rc     )
      us( ro     , ro     , ro+1:ns) = us( rc     , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      us( ro+1:ns, ro+1:ns, ro     ) = us( ro+1:ns, ro+1:ns, rc     )
      us( ro+1:ns, ro     , ro+1:ns) = us( ro+1:ns, rc     , ro+1:ns)
      us( ro     , ro+1:ns, ro+1:ns) = us( rc     , ro+1:ns, ro+1:ns)
      ! interpolate remaining points
      i0 = ro
      j0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = ro+1, ns,  1
      do i = ro+1, ns,  1
        c1 = 1 / sqrt(real((j - j0)**2 + (k - k0)**2, RWP))
        c2 = 1 / sqrt(real((i - i0)**2 + (k - k0)**2, RWP))
        c3 = 1 / sqrt(real((i - i0)**2 + (j - j0)**2, RWP))
        us(i,j,k) = (c1* us(i,j0,k0) + c2 * us(i0,j,k0) + c3 * us(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_8

  end subroutine FillIrregularZones
