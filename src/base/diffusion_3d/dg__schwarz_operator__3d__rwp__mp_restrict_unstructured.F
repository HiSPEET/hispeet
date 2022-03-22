  !-----------------------------------------------------------------------------
  !> Restrict an unstructured mesh variable `v` to subdomains
  !>
  !> The strategy is as a follows:
  !>
  !>   1. Send `v` from linked local elements to ghosts in remote partitions
  !>   2. Assign `v` to subdomain core regions, i.e. 1:1 copy to `vs`
  !>   3. Merge received remote data into ghosts
  !>   4. Copy local and ghost data to overlap zones of subdomains
  !>
  !> The mesh variable must be dimensioned as `v(0:po,0:po,0:po,ne+ng)`, where
  !> `ne` is the number of local elements and `ng` the number of ghosts.

  module subroutine Restrict_Unstructured(this, mesh, buf_v, v, vs)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_v
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: vs(:,:,:,:)   !< restricted variable

    real(RNP) :: ve(size(v,1)), vv

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
        call buf_v % Transfer(mesh, v, tag=1000)
      end if

      ! assign core regions ....................................................

      !$omp do
      do e = 1, mesh % n_elem

        ! initialization and removal of NANs
        vs(:,:,:,e) = 0

        do k = 0, po
        do j = 0, po
        do i = 0, po
          vs(c0 + i, c0 + j, c0 + k, e) = real(v(i,j,k,e), RWP)
        end do
        end do
        end do

      end do

      ! merge remote data to ghosts ............................................

      if (no > 0) then
        call buf_v % Merge(v)
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
              vs(1 + i, c0 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
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
              vs(c1 + i, c0 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
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
              vs(c0 + i, 1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
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
              vs(c0 + i, c1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
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
              vs(c0 + i, c0 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
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
              vs(c0 + i, c0 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
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
                vs(c0 + i, 1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(1) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 1) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, no, no, e) = real(ve, RWP)
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
                vs(c0 + i, c1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(2) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 2) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, c1 + 1, no, e) = real(ve, RWP)
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
                vs(c0 + i, 1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(3) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 3) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, no, c1 + 1, e) = real(ve, RWP)
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
                vs(c0 + i, c1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(4) = 0

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 4) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, c1 + 1, c1 + 1, e) = real(ve, RWP)
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
                vs(1 + i, c0 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(5) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 5) % orientation, &
                                    v, ve                                )
              vs(no, c0 : c0 + po, no, e) = real(ve, RWP)
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
                vs(c1 + i, c0 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(6) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 6) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c0 : c0 + po, no, e) = real(ve, RWP)
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
                vs(1 + i, c0 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(7) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 7) % orientation, &
                                    v, ve                                )
              vs(no, c0 : c0 + po, c1 + 1, e) = real(ve, RWP)
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
                vs(c1 + i, c0 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(8) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 8) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c0 : c0 + po, c1 + 1, e) = real(ve, RWP)
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
                vs(1 + i, 1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(9) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 9) % orientation, &
                                    v, ve                                )
              vs(no, no, c0 : c0 + po, e) = real(ve, RWP)
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
                vs(c1 + i, 1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(10) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(10) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, no, c0 : c0 + po, e) = real(ve, RWP)
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
                vs(1 + i, c1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(11) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(11) % orientation, &
                                    v, ve                                )
              vs(no, c1 + 1, c0 : c0 + po, e) = real(ve, RWP)
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
                vs(c1 + i, c1 + j, c0 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(12) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(12) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c1 + 1, c0 : c0 + po, e) = real(ve, RWP)
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
                vs(1 + i, 1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(1) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, no, no, e) = real(vv, RWP)
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
                vs(c1 + i, 1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(2) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, no, no, e) = real(vv, RWP)
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
                vs(1 + i, c1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(3) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, c1 + 1, no, e) = real(vv, RWP)
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
                vs(c1 + i, c1 + j, 1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(4) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, c1 + 1, no, e) = real(vv, RWP)
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
                vs(1 + i, 1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(5) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, no, c1 + 1, e) = real(vv, RWP)
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
                vs(c1 + i, 1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(6) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, no, c1 + 1, e) = real(vv, RWP)
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
                vs(1 + i, c1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(7) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, c1 + 1, c1 + 1, e) = real(vv, RWP)
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
                vs(c1 + i, c1 + j, c1 + k, e) = real(v(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(8) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, c1 + 1, c1 + 1, e) = real(vv, RWP)
              fill_vertex(8) = min(fill_vertex(8), 1)

            end if
          end if

          ! fill irregular overlap zones if requested
          if (this % fill) then
            call FillIrregularZones(no, fill_edge, fill_vertex, vs(:,:,:,e))
          end if

        end associate
      end do

    end associate

  end subroutine Restrict_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(neighbor, element, o, v, ve)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    class(MeshElement_3D),         intent(in)  :: element(:)
    integer(IXS),                  intent(in)  :: o
    real(RNP),                     intent(in)  :: v(:,:,:,:)
    real(RNP),                     intent(out) :: ve(:)

    real(RNP) :: vn(size(ve))
    integer :: nn, np
    integer :: e, l, n

    nn = size(neighbor)
    np = size(ve)

    ve = 0

    do n = 1, nn

      l = neighbor(n) % id
      e = neighbor(n) % component - 6

      select case(e)
      case( 1)
        vn = v( :, 1, 1,l)
      case( 2)
        vn = v( :,np, 1,l)
      case( 3)
        vn = v( :, 1,np,l)
      case( 4)
        vn = v( :,np,np,l)
      case( 5)
        vn = v( 1, :, 1,l)
      case( 6)
        vn = v(np, :, 1,l)
      case( 7)
        vn = v( 1, :,np,l)
      case( 8)
        vn = v(np, :,np,l)
      case( 9)
        vn = v( 1, 1, :,l)
      case(10)
        vn = v(np, 1, :,l)
      case(11)
        vn = v( 1,np, :,l)
      case(12)
        vn = v(np,np, :,l)
      end select

      if (element(l) % edge(e) % orientation == o) then
        ve = ve + vn
      else
        ve = ve + vn(np:1:-1)
      end if

    end do

    ve = ONE/nn * ve

  end subroutine CollectEdgeData

  !-----------------------------------------------------------------------------
  !> Collect and average neighbor vertex data

  pure subroutine CollectVertexData(neighbor, v, vv)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    real(RNP),                     intent(in)  :: v(:,:,:,:)
    real(RNP),                     intent(out) :: vv

    integer :: nn, np
    integer :: k, l, n

    nn = size(neighbor)
    np = size(v,1)

    vv = 0

    do n = 1, nn

      l = neighbor(n) % id
      k = neighbor(n) % component - 18

      select case(k)
      case( 1)
        vv = vv + v( 1, 1, 1,l)
      case( 2)
        vv = vv + v(np, 1, 1,l)
      case( 3)
        vv = vv + v( 1,np, 1,l)
      case( 4)
        vv = vv + v(np,np, 1,l)
      case( 5)
        vv = vv + v( 1, 1,np,l)
      case( 6)
        vv = vv + v(np, 1,np,l)
      case( 7)
        vv = vv + v( 1,np,np,l)
      case( 8)
        vv = vv + v(np,np,np,l)
      end select

    end do

    vv = ONE/nn * vv

  end subroutine CollectVertexData

  !-----------------------------------------------------------------------------
  !> Fill irregular overlap zones by interpolation

  pure subroutine FillIrregularZones(no, fill_edge, fill_vertex, vs)
    integer, intent(in) :: no             !< overlap
    integer, intent(in) :: fill_edge(12)  !< edge   type  {0,1,2}
    integer, intent(in) :: fill_vertex(8) !< vertex type  {0,1,2}
    real(RWP), contiguous, intent(inout) :: vs(:,:,:) !< subdomain variable

    integer   :: lc, lo, rc, ro, ns
    integer   :: i, i0, j, j0, k, k0
    real(RWP) :: c1, c2, c3

    if (no < 2) return

    ns = size(vs,1)
    lo = no         ! end of left overlap zone
    lc = no + 1     ! leftmost  core point
    rc = ns - no    ! rightmost core point
    ro = rc + 1     ! start of right overlap zone

    if (fill_edge( 1) > 0) then
      ! adopt from core
      if (fill_edge( 1) > 1) vs(lc:rc,lo,lo) = vs(lc:rc,lc,lc)
      ! adopt from face overlap zones
      vs(lc:rc,  1: lo-1,     lo  ) = vs(lc:rc,  1: lo-1,     lc  )
      vs(lc:rc,     lo  ,  1: lo-1) = vs(lc:rc,     lc  ,  1: lo-1)
      ! interpolate remaining points
      j0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        vs(i,j,k) = (c2 * vs(i,j,k0) + c3 * vs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 2) > 0) then
      ! adopt from core
      if (fill_edge( 2) > 1) vs(lc:rc,ro,lo) = vs(lc:rc,rc,lc)
      ! adopt from face overlap zones
      vs(lc:rc, ro+1: ns,     lo  ) = vs(lc:rc, ro+1: ns,     lc  )
      vs(lc:rc, ro      ,  1: lo-1) = vs(lc:rc, rc      ,  1: lo-1)
      ! interpolate remaining points
      j0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        vs(i,j,k) = (c2 * vs(i,j,k0) + c3 * vs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 3) > 0) then
      ! adopt from core
      if (fill_edge( 3) > 1) vs(lc:rc,lo,ro) = vs(lc:rc,lc,rc)
      ! adopt from face overlap zones
      vs(lc:rc,  1: lo-1, ro      ) = vs(lc:rc,  1: lo-1, rc      )
      vs(lc:rc,     lo  , ro+1: ns) = vs(lc:rc,     lc  , ro+1: ns)
      ! interpolate remaining points
      j0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        vs(i,j,k) = (c2 * vs(i,j,k0) + c3 * vs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 4) > 0) then
      ! adopt from core
      if (fill_edge( 4) > 1) vs(lc:rc,ro,ro) = vs(lc:rc,rc,rc)
      ! adopt from face overlap zones
      vs(lc:rc, ro+1: ns, ro      ) = vs(lc:rc, ro+1: ns, rc      )
      vs(lc:rc, ro      , ro+1: ns) = vs(lc:rc, rc      , ro+1: ns)
      ! interpolate remaining points
      j0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        vs(i,j,k) = (c2 * vs(i,j,k0) + c3 * vs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 5) > 0) then
      ! adopt from core
      if (fill_edge( 5) > 1) vs(lo,lc:rc,lo) = vs(lc,lc:rc,lc)
      ! adopt from face overlap zones
      vs( 1: lo-1, lc:rc,     lo  ) = vs( 1: lo-1, lc:rc,     lc  )
      vs(    lo  , lc:rc,  1: lo-1) = vs(    lc  , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j,k0) + c3 * vs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 6) > 0) then
      ! adopt from core
      if (fill_edge( 6) > 1) vs(ro,lc:rc,lo) = vs(rc,lc:rc,lc)
      ! adopt from face overlap zones
      vs(ro+1: ns, lc:rc,     lo  ) = vs(ro+1: ns, lc:rc,     lc  )
      vs(ro      , lc:rc,  1: lo-1) = vs(rc      , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j,k0) + c3 * vs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 7) > 0) then
      ! adopt from core
      if (fill_edge( 7) > 1) vs(lo,lc:rc,ro) = vs(lc,lc:rc,rc)
      ! adopt from face overlap zones
      vs( 1: lo-1, lc:rc, ro      ) = vs( 1: lo-1, lc:rc, rc      )
      vs(    lo  , lc:rc, ro+1: ns) = vs(    lc  , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j,k0) + c3 * vs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 8) > 0) then
      ! adopt from core
      if (fill_edge( 8) > 1) vs(ro,lc:rc,ro) = vs(rc,lc:rc,rc)
      ! adopt from face overlap zones
      vs(ro+1: ns, lc:rc, ro      ) = vs(ro+1: ns, lc:rc, rc      )
      vs(ro      , lc:rc, ro+1: ns) = vs(rc      , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j,k0) + c3 * vs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if

    if (fill_edge( 9) > 0) then
      ! adopt from core
      if (fill_edge( 9) > 1) vs(lo,lo,lc:rc) = vs(lc,lc,lc:rc)
      ! adopt from face overlap zones
      vs( 1: lo-1,     lo  , lc:rc) = vs( 1: lo-1,     lc  , lc:rc)
      vs(    lo  ,  1: lo-1, lc:rc) = vs(    lc  ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j0,k) + c2 * vs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if

    if (fill_edge(10) > 0) then
      ! adopt from core
      if (fill_edge(10) > 1) vs(ro,lo,lc:rc) = vs(rc,lc,lc:rc)
      ! adopt from face overlap zones
      vs(ro+1: ns,     lo  , lc:rc) = vs(ro+1: ns,     lc  , lc:rc)
      vs(ro      ,  1: lo-1, lc:rc) = vs(rc      ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j0,k) + c2 * vs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if

    if (fill_edge(11) > 0) then
      ! adopt from core
      if (fill_edge(11) > 1) vs(lo,ro,lc:rc) = vs(lc,rc,lc:rc)
      ! adopt from face overlap zones
      vs( 1: lo-1, ro      , lc:rc) = vs( 1: lo-1, rc      , lc:rc)
      vs(    lo  , ro+1: ns, lc:rc) = vs(    lc  , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j0,k) + c2 * vs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if

    if (fill_edge(12) > 0) then
      ! adopt from core
      if (fill_edge(12) > 1) vs(ro,ro,lc:rc) = vs(rc,rc,lc:rc)
      ! adopt from face overlap zones
      vs(ro+1: ns, ro      , lc:rc) = vs(ro+1: ns, rc      , lc:rc)
      vs(ro      , ro+1: ns, lc:rc) = vs(rc      , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        vs(i,j,k) = (c1 * vs(i,j0,k) + c2 * vs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if

    if (fill_vertex(1) > 0) then
      ! adopt from core
      if (fill_vertex(1) > 1) vs(lo,lo,lo) = vs(lc,lc,lc)
      ! adopt from face overlap zones
      vs(  1:lo-1,    lo  ,    lo  ) = vs(  1:lo-1,    lc  ,    lc  )
      vs(    lo  ,  1:lo-1,    lo  ) = vs(    lc  ,  1:lo-1,    lc  )
      vs(    lo  ,    lo  ,  1:lo-1) = vs(    lc  ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      vs(  1:lo-1,  1:lo-1,    lo  ) = vs(  1:lo-1,  1:lo-1,    lc  )
      vs(  1:lo-1,    lo  ,  1:lo-1) = vs(  1:lo-1,    lc  ,  1:lo-1)
      vs(    lo  ,  1:lo-1,  1:lo-1) = vs(    lc  ,  1:lo-1,  1:lo-1)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(2) > 0) then
      ! adopt from core
      if (fill_vertex(2) > 1) vs(ro,lo,lo) = vs(rc,lc,lc)
      ! adopt from face overlap zones
      vs( ro+1:ns,    lo  ,    lo  ) = vs( ro+1:ns,    lc  ,    lc  )
      vs( ro     ,  1:lo-1,    lo  ) = vs( rc     ,  1:lo-1,    lc  )
      vs( ro     ,    lo  ,  1:lo-1) = vs( rc     ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      vs( ro+1:ns,  1:lo-1,    lo  ) = vs( ro+1:ns,  1:lo-1,    lc  )
      vs( ro+1:ns,    lo  ,  1:lo-1) = vs( ro+1:ns,    lc  ,  1:lo-1)
      vs( ro     ,  1:lo-1,  1:lo-1) = vs( rc     ,  1:lo-1,  1:lo-1)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(3) > 0) then
      ! adopt from core
      if (fill_vertex(3) > 1) vs(lo,ro,lo) = vs(lc,rc,lc)
      ! adopt from face overlap zones
      vs(    lo  , ro+1:ns,    lo  ) = vs(    lc  , ro+1:ns,    lc  )
      vs(  1:lo-1, ro     ,    lo  ) = vs(  1:lo-1, rc     ,    lc  )
      vs(    lo  , ro     ,  1:lo-1) = vs(    lc  , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      vs(  1:lo-1, ro+1:ns,    lo  ) = vs(  1:lo-1, ro+1:ns,    lc  )
      vs(  1:lo-1, ro     ,  1:lo-1) = vs(  1:lo-1, rc     ,  1:lo-1)
      vs(    lo  , ro+1:ns,  1:lo-1) = vs(    lc  , ro+1:ns,  1:lo-1)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(4) > 0) then
      ! adopt from core
      if (fill_vertex(4) > 1) vs(ro,ro,lo) = vs(rc,rc,lc)
      ! adopt from face overlap zones
      vs( ro     , ro+1:ns,    lo  ) = vs( rc     , ro+1:ns,    lc  )
      vs( ro+1:ns, ro     ,    lo  ) = vs( ro+1:ns, rc     ,    lc  )
      vs( ro     , ro     ,  1:lo-1) = vs( rc     , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      vs( ro+1:ns, ro+1:ns,    lo  ) = vs( ro+1:ns, ro+1:ns,    lc  )
      vs( ro+1:ns, ro     ,  1:lo-1) = vs( ro+1:ns, rc     ,  1:lo-1)
      vs( ro     , ro+1:ns,  1:lo-1) = vs( rc     , ro+1:ns,  1:lo-1)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(5) > 0) then
      ! adopt from core
      if (fill_vertex(5) > 1) vs(lo,lo,ro) = vs(lc,lc,rc)
      ! adopt from face overlap zones
      vs(  1:lo-1,    lo  , ro     ) = vs(  1:lo-1,    lc  , rc     )
      vs(    lo  ,  1:lo-1, ro     ) = vs(    lc  ,  1:lo-1, rc     )
      vs(    lo  ,    lo  , ro+1:ns) = vs(    lc  ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      vs(  1:lo-1,  1:lo-1, ro     ) = vs(  1:lo-1,  1:lo-1, rc     )
      vs(  1:lo-1,    lo  , ro+1:ns) = vs(  1:lo-1,    lc  , ro+1:ns)
      vs(    lo  ,  1:lo-1, ro+1:ns) = vs(    lc  ,  1:lo-1, ro+1:ns)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(6) > 0) then
      ! adopt from core
      if (fill_vertex(6) > 1) vs(ro,lo,ro) = vs(rc,lc,rc)
      ! adopt from face overlap zones
      vs( ro+1:ns,    lo  , ro     ) = vs( ro+1:ns,    lc  , rc     )
      vs( ro     ,  1:lo-1, ro     ) = vs( rc     ,  1:lo-1, rc     )
      vs( ro     ,    lo  , ro+1:ns) = vs( rc     ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      vs( ro+1:ns,  1:lo-1, ro     ) = vs( ro+1:ns,  1:lo-1, rc     )
      vs( ro+1:ns,    lo  , ro+1:ns) = vs( ro+1:ns,    lc  , ro+1:ns)
      vs( ro     ,  1:lo-1, ro+1:ns) = vs( rc     ,  1:lo-1, ro+1:ns)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(7) > 0) then
      ! adopt from core
      if (fill_vertex(7) > 1) vs(lo,ro,ro) = vs(lc,rc,rc)
      ! adopt from face overlap zones
      vs(    lo  , ro+1:ns, ro     ) = vs(    lc  , ro+1:ns, rc     )
      vs(  1:lo-1, ro     , ro     ) = vs(  1:lo-1, rc     , rc     )
      vs(    lo  , ro     , ro+1:ns) = vs(    lc  , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      vs(  1:lo-1, ro+1:ns, ro     ) = vs(  1:lo-1, ro+1:ns, rc     )
      vs(  1:lo-1, ro     , ro+1:ns) = vs(  1:lo-1, rc     , ro+1:ns)
      vs(    lo  , ro+1:ns, ro+1:ns) = vs(    lc  , ro+1:ns, ro+1:ns)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

    if (fill_vertex(8) > 0) then
      ! adopt from core
      if (fill_vertex(8) > 1) vs(ro,ro,ro) = vs(rc,rc,rc)
      ! adopt from face overlap zones
      vs( ro     , ro+1:ns, ro     ) = vs( rc     , ro+1:ns, rc     )
      vs( ro+1:ns, ro     , ro     ) = vs( ro+1:ns, rc     , rc     )
      vs( ro     , ro     , ro+1:ns) = vs( rc     , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      vs( ro+1:ns, ro+1:ns, ro     ) = vs( ro+1:ns, ro+1:ns, rc     )
      vs( ro+1:ns, ro     , ro+1:ns) = vs( ro+1:ns, rc     , ro+1:ns)
      vs( ro     , ro+1:ns, ro+1:ns) = vs( rc     , ro+1:ns, ro+1:ns)
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
        vs(i,j,k) = (c1* vs(i,j0,k0) + c2 * vs(i0,j,k0) + c3 * vs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if

  end subroutine FillIrregularZones
