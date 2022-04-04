  !-----------------------------------------------------------------------------
  !> Restrict an unstructured mesh variable `r` to subdomains
  !>
  !> The strategy is as a follows:
  !>
  !>   1. Send `r` from linked local elements to ghosts in remote partitions
  !>   2. Assign `r` to subdomain core regions, i.e. 1:1 copy to `rs`
  !>   3. Merge received remote data into ghosts
  !>   4. Copy local and ghost data to overlap zones of subdomains
  !>
  !> The mesh variable must be dimensioned as `r(0:po,0:po,0:po,ne+ng)`, where
  !> `ne` is the number of local elements and `ng` the number of ghosts.

  module subroutine Restrict_Unstructured(this, mesh, buf_r, r, rs, sgn)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
    integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]

    real(RNP) :: re(size(r,1)), rv
    real(RWP) :: c_sgn

    integer(IXS) :: orientation
    integer :: e, i, j, k, l, m, n
    integer :: lo, ro, nn, np
    integer :: ii, jj, kk
    integer :: c0, c1, o(3), s(3,3), t(3)
    integer :: fill_edge(12), fill_vertex(8)

    if (present(sgn)) then
      c_sgn = sgn
    else
      c_sgn = 1
    end if

    associate(po => this%po, no => this%no)

      ! start transfer .........................................................

      np = po + 1       ! number of element   points per direction

      lo = no - 1       ! last  point of left  element overlap zone
      ro = po - lo      ! first point of right element overlap zone

      c0 = no + 1       ! first subdomain point inside core
      c1 = c0 + np      ! first subdomain point behind core

      if (no > 0) then
        call buf_r % Transfer(mesh, r, tag=1000)
      end if

      ! assign core regions ....................................................

      !$omp do
      do e = 1, mesh % n_elem

        ! initialization and removal of NANs
        rs(:,:,:,e) = 0

        do k = 0, po
        do j = 0, po
        do i = 0, po
          rs(c0 + i, c0 + j, c0 + k, e) = c_sgn * real(r(i,j,k,e), RWP)
        end do
        end do
        end do

      end do

      ! merge remote data to ghosts ............................................

      if (no > 0) then
        call buf_r % Merge(r)
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
              rs(1 + i, c0 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
              rs(c1 + i, c0 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
              rs(c0 + i, 1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
              rs(c0 + i, c1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
              rs(c0 + i, c0 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
              rs(c0 + i, c0 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
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
                rs(c0 + i, 1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(1) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 1) % orientation, &
                                    r, re                                )
              rs(c0 : c0 + po, no, no, e) = c_sgn * real(re, RWP)
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
                rs(c0 + i, c1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(2) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 2) % orientation, &
                                    r, re                                )
              rs(c0 : c0 + po, c1, no, e) = c_sgn * real(re, RWP)
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
                rs(c0 + i, 1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(3) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 3) % orientation, &
                                    r, re                                )
              rs(c0 : c0 + po, no, c1, e) = c_sgn * real(re, RWP)
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
                rs(c0 + i, c1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(4) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 4) % orientation, &
                                    r, re                                )
              rs(c0 : c0 + po, c1, c1, e) = c_sgn * real(re, RWP)
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
                rs(1 + i, c0 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(5) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 5) % orientation, &
                                    r, re                                )
              rs(no, c0 : c0 + po, no, e) = c_sgn * real(re, RWP)
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
                rs(c1 + i, c0 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(6) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 6) % orientation, &
                                    r, re                                )
              rs(c1, c0 : c0 + po, no, e) = c_sgn * real(re, RWP)
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
                rs(1 + i, c0 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(7) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 7) % orientation, &
                                    r, re                                )
              rs(no, c0 : c0 + po, c1, e) = c_sgn * real(re, RWP)
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
                rs(c1 + i, c0 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(8) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 8) % orientation, &
                                    r, re                                )
              rs(c1, c0 : c0 + po, c1, e) = c_sgn * real(re, RWP)
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
                rs(1 + i, 1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(9) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 9) % orientation, &
                                    r, re                                )
              rs(no, no, c0 : c0 + po, e) = c_sgn * real(re, RWP)
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
                rs(c1 + i, 1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(10) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(10) % orientation, &
                                    r, re                                )
              rs(c1, no, c0 : c0 + po, e) = c_sgn * real(re, RWP)
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
                rs(1 + i, c1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(11) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(11) % orientation, &
                                    r, re                                )
              rs(no, c1, c0 : c0 + po, e) = c_sgn * real(re, RWP)
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
                rs(c1 + i, c1 + j, c0 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_edge(12) = 0

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(12) % orientation, &
                                    r, re                                )
              rs(c1, c1, c0 : c0 + po, e) = c_sgn * real(re, RWP)
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
                rs(1 + i, 1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(1) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, no, no, e) = c_sgn * real(rv, RWP)
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
                rs(c1 + i, 1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(2) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, no, no, e) = c_sgn * real(rv, RWP)
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
                rs(1 + i, c1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(3) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, c1, no, e) = c_sgn * real(rv, RWP)
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
                rs(c1 + i, c1 + j, 1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(4) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, c1, no, e) = c_sgn * real(rv, RWP)
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
                rs(1 + i, 1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(5) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, no, c1, e) = c_sgn * real(rv, RWP)
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
                rs(c1 + i, 1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(6) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, no, c1, e) = c_sgn * real(rv, RWP)
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
                rs(1 + i, c1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(7) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, c1, c1, e) = c_sgn * real(rv, RWP)
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
                rs(c1 + i, c1 + j, c1 + k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do
              fill_vertex(8) = 0

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, c1, c1, e) = c_sgn * real(rv, RWP)
              fill_vertex(8) = min(fill_vertex(8), 1)

            end if
          end if

          ! fill irregular overlap zones if requested
          if (this % fill) then
            call FillIrregularZones(no, fill_edge, fill_vertex, rs(:,:,:,e))
          end if

        end associate
      end do

    end associate

  end subroutine Restrict_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(neighbor, element, o, r, re)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    class(MeshElement_3D),         intent(in)  :: element(:)
    integer(IXS),                  intent(in)  :: o
    real(RNP),                     intent(in)  :: r(:,:,:,:)
    real(RNP),                     intent(out) :: re(:)

    real(RNP) :: rn(size(re))
    integer :: nn, np
    integer :: e, l, n

    nn = size(neighbor)
    np = size(re)

    re = 0

    do n = 1, nn

      l = neighbor(n) % id
      e = neighbor(n) % component - 6

      select case(e)
      case( 1)
        rn = r( :, 1, 1,l)
      case( 2)
        rn = r( :,np, 1,l)
      case( 3)
        rn = r( :, 1,np,l)
      case( 4)
        rn = r( :,np,np,l)
      case( 5)
        rn = r( 1, :, 1,l)
      case( 6)
        rn = r(np, :, 1,l)
      case( 7)
        rn = r( 1, :,np,l)
      case( 8)
        rn = r(np, :,np,l)
      case( 9)
        rn = r( 1, 1, :,l)
      case(10)
        rn = r(np, 1, :,l)
      case(11)
        rn = r( 1,np, :,l)
      case(12)
        rn = r(np,np, :,l)
      end select

      if (element(l) % edge(e) % orientation == o) then
        re = re + rn
      else
        re = re + rn(np:1:-1)
      end if

    end do

    re = ONE/nn * re

  end subroutine CollectEdgeData

  !-----------------------------------------------------------------------------
  !> Collect and average neighbor vertex data

  pure subroutine CollectVertexData(neighbor, r, rv)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    real(RNP),                     intent(in)  :: r(:,:,:,:)
    real(RNP),                     intent(out) :: rv

    integer :: nn, np
    integer :: k, l, n

    nn = size(neighbor)
    np = size(r,1)

    rv = 0

    do n = 1, nn

      l = neighbor(n) % id
      k = neighbor(n) % component - 18

      select case(k)
      case( 1)
        rv = rv + r( 1, 1, 1,l)
      case( 2)
        rv = rv + r(np, 1, 1,l)
      case( 3)
        rv = rv + r( 1,np, 1,l)
      case( 4)
        rv = rv + r(np,np, 1,l)
      case( 5)
        rv = rv + r( 1, 1,np,l)
      case( 6)
        rv = rv + r(np, 1,np,l)
      case( 7)
        rv = rv + r( 1,np,np,l)
      case( 8)
        rv = rv + r(np,np,np,l)
      end select

    end do

    rv = ONE/nn * rv

  end subroutine CollectVertexData

  !-----------------------------------------------------------------------------
  !> Fill irregular overlap zones by interpolation

  pure subroutine FillIrregularZones(no, fill_edge, fill_vertex, rs)
    integer, intent(in) :: no             !< overlap
    integer, intent(in) :: fill_edge(12)  !< edge   type  {0,1,2}
    integer, intent(in) :: fill_vertex(8) !< vertex type  {0,1,2}
    real(RWP), contiguous, intent(inout) :: rs(:,:,:) !< subdomain variable

    integer   :: lc, lo, rc, ro, ns
    integer   :: i, i0, j, j0, k, k0
    real(RWP) :: c1, c2, c3

    if (no < 1) return

    ns = size(rs,1)
    lo = no         ! end of left overlap zone
    lc = no + 1     ! leftmost  core point
    rc = ns - no    ! rightmost core point
    ro = rc + 1     ! start of right overlap zone

    EDGE_1: if (fill_edge( 1) > 0) then
      ! adopt from core
      if (fill_edge( 1) > 1) rs(lc:rc,lo,lo) = rs(lc:rc,lc,lc)
      if (no == 1) exit EDGE_1
      ! adopt from face overlap zones
      rs(lc:rc,  1: lo-1,     lo  ) = rs(lc:rc,  1: lo-1,     lc  )
      rs(lc:rc,     lo  ,  1: lo-1) = rs(lc:rc,     lc  ,  1: lo-1)
      ! interpolate remaining points
      j0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        rs(i,j,k) = (c2 * rs(i,j,k0) + c3 * rs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_1

    EDGE_2: if (fill_edge( 2) > 0) then
      ! adopt from core
      if (fill_edge( 2) > 1) rs(lc:rc,ro,lo) = rs(lc:rc,rc,lc)
      if (no == 1) exit EDGE_2
      ! adopt from face overlap zones
      rs(lc:rc, ro+1: ns,     lo  ) = rs(lc:rc, ro+1: ns,     lc  )
      rs(lc:rc, ro      ,  1: lo-1) = rs(lc:rc, rc      ,  1: lo-1)
      ! interpolate remaining points
      j0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        rs(i,j,k) = (c2 * rs(i,j,k0) + c3 * rs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_2

    EDGE_3: if (fill_edge( 3) > 0) then
      ! adopt from core
      if (fill_edge( 3) > 1) rs(lc:rc,lo,ro) = rs(lc:rc,lc,rc)
      if (no == 1) exit EDGE_3
      ! adopt from face overlap zones
      rs(lc:rc,  1: lo-1, ro      ) = rs(lc:rc,  1: lo-1, rc      )
      rs(lc:rc,     lo  , ro+1: ns) = rs(lc:rc,     lc  , ro+1: ns)
      ! interpolate remaining points
      j0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lo-1,  1, -1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        rs(i,j,k) = (c2 * rs(i,j,k0) + c3 * rs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_3

    EDGE_4: if (fill_edge( 4) > 0) then
      ! adopt from core
      if (fill_edge( 4) > 1) rs(lc:rc,ro,ro) = rs(lc:rc,rc,rc)
      if (no == 1) exit EDGE_4
      ! adopt from face overlap zones
      rs(lc:rc, ro+1: ns, ro      ) = rs(lc:rc, ro+1: ns, rc      )
      rs(lc:rc, ro      , ro+1: ns) = rs(lc:rc, rc      , ro+1: ns)
      ! interpolate remaining points
      j0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = ro+1, ns,  1
      do i = lc  , rc
        c2 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(j - j0), RWP)
        rs(i,j,k) = (c2 * rs(i,j,k0) + c3 * rs(i,j0,k)) / (c2 + c3)
      end do
      end do
      end do
    end if EDGE_4

    EDGE_5: if (fill_edge( 5) > 0) then
      ! adopt from core
      if (fill_edge( 5) > 1) rs(lo,lc:rc,lo) = rs(lc,lc:rc,lc)
      if (no == 1) exit EDGE_5
      ! adopt from face overlap zones
      rs( 1: lo-1, lc:rc,     lo  ) = rs( 1: lo-1, lc:rc,     lc  )
      rs(    lo  , lc:rc,  1: lo-1) = rs(    lc  , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = lo
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j,k0) + c3 * rs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_5

    EDGE_6: if (fill_edge( 6) > 0) then
      ! adopt from core
      if (fill_edge( 6) > 1) rs(ro,lc:rc,lo) = rs(rc,lc:rc,lc)
      if (no == 1) exit EDGE_6
      ! adopt from face overlap zones
      rs(ro+1: ns, lc:rc,     lo  ) = rs(ro+1: ns, lc:rc,     lc  )
      rs(ro      , lc:rc,  1: lo-1) = rs(rc      , lc:rc,  1: lo-1)
      ! interpolate remaining points
      i0 = ro
      k0 = lo
      do k = lo-1,  1, -1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j,k0) + c3 * rs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_6

    EDGE_7: if (fill_edge( 7) > 0) then
      ! adopt from core
      if (fill_edge( 7) > 1) rs(lo,lc:rc,ro) = rs(lc,lc:rc,rc)
      if (no == 1) exit EDGE_7
      ! adopt from face overlap zones
      rs( 1: lo-1, lc:rc, ro      ) = rs( 1: lo-1, lc:rc, rc      )
      rs(    lo  , lc:rc, ro+1: ns) = rs(    lc  , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = lo
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j,k0) + c3 * rs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_7

    EDGE_8: if (fill_edge( 8) > 0) then
      ! adopt from core
      if (fill_edge( 8) > 1) rs(ro,lc:rc,ro) = rs(rc,lc:rc,rc)
      if (no == 1) exit EDGE_8
      ! adopt from face overlap zones
      rs(ro+1: ns, lc:rc, ro      ) = rs(ro+1: ns, lc:rc, rc      )
      rs(ro      , lc:rc, ro+1: ns) = rs(rc      , lc:rc, ro+1: ns)
      ! interpolate remaining points
      i0 = ro
      k0 = ro
      do k = ro+1, ns,  1
      do j = lc  , rc
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(k - k0), RWP)
        c3 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j,k0) + c3 * rs(i0,j,k)) / (c1 + c3)
      end do
      end do
      end do
    end if EDGE_8

    EDGE_9: if (fill_edge( 9) > 0) then
      ! adopt from core
      if (fill_edge( 9) > 1) rs(lo,lo,lc:rc) = rs(lc,lc,lc:rc)
      if (no == 1) exit EDGE_9
      ! adopt from face overlap zones
      rs( 1: lo-1,     lo  , lc:rc) = rs( 1: lo-1,     lc  , lc:rc)
      rs(    lo  ,  1: lo-1, lc:rc) = rs(    lc  ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j0,k) + c2 * rs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_9

    EDGE_10: if (fill_edge(10) > 0) then
      ! adopt from core
      if (fill_edge(10) > 1) rs(ro,lo,lc:rc) = rs(rc,lc,lc:rc)
      if (no == 1) exit EDGE_10
      ! adopt from face overlap zones
      rs(ro+1: ns,     lo  , lc:rc) = rs(ro+1: ns,     lc  , lc:rc)
      rs(ro      ,  1: lo-1, lc:rc) = rs(rc      ,  1: lo-1, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = lo
      do k = lc  , rc
      do j = lo-1,  1, -1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j0,k) + c2 * rs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_10

    EDGE_11: if (fill_edge(11) > 0) then
      ! adopt from core
      if (fill_edge(11) > 1) rs(lo,ro,lc:rc) = rs(lc,rc,lc:rc)
      if (no == 1) exit EDGE_11
      ! adopt from face overlap zones
      rs( 1: lo-1, ro      , lc:rc) = rs( 1: lo-1, rc      , lc:rc)
      rs(    lo  , ro+1: ns, lc:rc) = rs(    lc  , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = lo
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = lo-1,  1, -1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j0,k) + c2 * rs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_11

    EDGE_12: if (fill_edge(12) > 0) then
      ! adopt from core
      if (fill_edge(12) > 1) rs(ro,ro,lc:rc) = rs(rc,rc,lc:rc)
      if (no == 1) exit EDGE_12
      ! adopt from face overlap zones
      rs(ro+1: ns, ro      , lc:rc) = rs(ro+1: ns, rc      , lc:rc)
      rs(ro      , ro+1: ns, lc:rc) = rs(rc      , ro+1: ns, lc:rc)
      ! interpolate remaining points
      i0 = ro
      j0 = ro
      do k = lc  , rc
      do j = ro+1, ns,  1
      do i = ro+1, ns,  1
        c1 = 1 / real(abs(j - j0), RWP)
        c2 = 1 / real(abs(i - i0), RWP)
        rs(i,j,k) = (c1 * rs(i,j0,k) + c2 * rs(i0,j,k)) / (c1 + c2)
      end do
      end do
      end do
    end if EDGE_12

    VERTEX_1: if (fill_vertex(1) > 0) then
      ! adopt from core
      if (fill_vertex(1) > 1) rs(lo,lo,lo) = rs(lc,lc,lc)
      if (no == 1) exit VERTEX_1
      ! adopt from face overlap zones
      rs(  1:lo-1,    lo  ,    lo  ) = rs(  1:lo-1,    lc  ,    lc  )
      rs(    lo  ,  1:lo-1,    lo  ) = rs(    lc  ,  1:lo-1,    lc  )
      rs(    lo  ,    lo  ,  1:lo-1) = rs(    lc  ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      rs(  1:lo-1,  1:lo-1,    lo  ) = rs(  1:lo-1,  1:lo-1,    lc  )
      rs(  1:lo-1,    lo  ,  1:lo-1) = rs(  1:lo-1,    lc  ,  1:lo-1)
      rs(    lo  ,  1:lo-1,  1:lo-1) = rs(    lc  ,  1:lo-1,  1:lo-1)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_1

    VERTEX_2: if (fill_vertex(2) > 0) then
      ! adopt from core
      if (fill_vertex(2) > 1) rs(ro,lo,lo) = rs(rc,lc,lc)
      if (no == 1) exit VERTEX_2
      ! adopt from face overlap zones
      rs( ro+1:ns,    lo  ,    lo  ) = rs( ro+1:ns,    lc  ,    lc  )
      rs( ro     ,  1:lo-1,    lo  ) = rs( rc     ,  1:lo-1,    lc  )
      rs( ro     ,    lo  ,  1:lo-1) = rs( rc     ,    lc  ,  1:lo-1)
      ! adopt from edge overlap zones
      rs( ro+1:ns,  1:lo-1,    lo  ) = rs( ro+1:ns,  1:lo-1,    lc  )
      rs( ro+1:ns,    lo  ,  1:lo-1) = rs( ro+1:ns,    lc  ,  1:lo-1)
      rs( ro     ,  1:lo-1,  1:lo-1) = rs( rc     ,  1:lo-1,  1:lo-1)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_2

    VERTEX_3: if (fill_vertex(3) > 0) then
      ! adopt from core
      if (fill_vertex(3) > 1) rs(lo,ro,lo) = rs(lc,rc,lc)
      if (no == 1) exit VERTEX_3
      ! adopt from face overlap zones
      rs(    lo  , ro+1:ns,    lo  ) = rs(    lc  , ro+1:ns,    lc  )
      rs(  1:lo-1, ro     ,    lo  ) = rs(  1:lo-1, rc     ,    lc  )
      rs(    lo  , ro     ,  1:lo-1) = rs(    lc  , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      rs(  1:lo-1, ro+1:ns,    lo  ) = rs(  1:lo-1, ro+1:ns,    lc  )
      rs(  1:lo-1, ro     ,  1:lo-1) = rs(  1:lo-1, rc     ,  1:lo-1)
      rs(    lo  , ro+1:ns,  1:lo-1) = rs(    lc  , ro+1:ns,  1:lo-1)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_3

    VERTEX_4: if (fill_vertex(4) > 0) then
      ! adopt from core
      if (fill_vertex(4) > 1) rs(ro,ro,lo) = rs(rc,rc,lc)
      if (no == 1) exit VERTEX_4
      ! adopt from face overlap zones
      rs( ro     , ro+1:ns,    lo  ) = rs( rc     , ro+1:ns,    lc  )
      rs( ro+1:ns, ro     ,    lo  ) = rs( ro+1:ns, rc     ,    lc  )
      rs( ro     , ro     ,  1:lo-1) = rs( rc     , rc     ,  1:lo-1)
      ! adopt from edge overlap zones
      rs( ro+1:ns, ro+1:ns,    lo  ) = rs( ro+1:ns, ro+1:ns,    lc  )
      rs( ro+1:ns, ro     ,  1:lo-1) = rs( ro+1:ns, rc     ,  1:lo-1)
      rs( ro     , ro+1:ns,  1:lo-1) = rs( rc     , ro+1:ns,  1:lo-1)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_4

    VERTEX_5: if (fill_vertex(5) > 0) then
      ! adopt from core
      if (fill_vertex(5) > 1) rs(lo,lo,ro) = rs(lc,lc,rc)
      if (no == 1) exit VERTEX_5
      ! adopt from face overlap zones
      rs(  1:lo-1,    lo  , ro     ) = rs(  1:lo-1,    lc  , rc     )
      rs(    lo  ,  1:lo-1, ro     ) = rs(    lc  ,  1:lo-1, rc     )
      rs(    lo  ,    lo  , ro+1:ns) = rs(    lc  ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      rs(  1:lo-1,  1:lo-1, ro     ) = rs(  1:lo-1,  1:lo-1, rc     )
      rs(  1:lo-1,    lo  , ro+1:ns) = rs(  1:lo-1,    lc  , ro+1:ns)
      rs(    lo  ,  1:lo-1, ro+1:ns) = rs(    lc  ,  1:lo-1, ro+1:ns)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_5

    VERTEX_6: if (fill_vertex(6) > 0) then
      ! adopt from core
      if (fill_vertex(6) > 1) rs(ro,lo,ro) = rs(rc,lc,rc)
      if (no == 1) exit VERTEX_6
      ! adopt from face overlap zones
      rs( ro+1:ns,    lo  , ro     ) = rs( ro+1:ns,    lc  , rc     )
      rs( ro     ,  1:lo-1, ro     ) = rs( rc     ,  1:lo-1, rc     )
      rs( ro     ,    lo  , ro+1:ns) = rs( rc     ,    lc  , ro+1:ns)
      ! adopt from edge overlap zones
      rs( ro+1:ns,  1:lo-1, ro     ) = rs( ro+1:ns,  1:lo-1, rc     )
      rs( ro+1:ns,    lo  , ro+1:ns) = rs( ro+1:ns,    lc  , ro+1:ns)
      rs( ro     ,  1:lo-1, ro+1:ns) = rs( rc     ,  1:lo-1, ro+1:ns)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_6

    VERTEX_7: if (fill_vertex(7) > 0) then
      ! adopt from core
      if (fill_vertex(7) > 1) rs(lo,ro,ro) = rs(lc,rc,rc)
      if (no == 1) exit VERTEX_7
      ! adopt from face overlap zones
      rs(    lo  , ro+1:ns, ro     ) = rs(    lc  , ro+1:ns, rc     )
      rs(  1:lo-1, ro     , ro     ) = rs(  1:lo-1, rc     , rc     )
      rs(    lo  , ro     , ro+1:ns) = rs(    lc  , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      rs(  1:lo-1, ro+1:ns, ro     ) = rs(  1:lo-1, ro+1:ns, rc     )
      rs(  1:lo-1, ro     , ro+1:ns) = rs(  1:lo-1, rc     , ro+1:ns)
      rs(    lo  , ro+1:ns, ro+1:ns) = rs(    lc  , ro+1:ns, ro+1:ns)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_7

    VERTEX_8: if (fill_vertex(8) > 0) then
      ! adopt from core
      if (fill_vertex(8) > 1) rs(ro,ro,ro) = rs(rc,rc,rc)
      if (no == 1) exit VERTEX_8
      ! adopt from face overlap zones
      rs( ro     , ro+1:ns, ro     ) = rs( rc     , ro+1:ns, rc     )
      rs( ro+1:ns, ro     , ro     ) = rs( ro+1:ns, rc     , rc     )
      rs( ro     , ro     , ro+1:ns) = rs( rc     , rc     , ro+1:ns)
      ! adopt from edge overlap zones
      rs( ro+1:ns, ro+1:ns, ro     ) = rs( ro+1:ns, ro+1:ns, rc     )
      rs( ro+1:ns, ro     , ro+1:ns) = rs( ro+1:ns, rc     , ro+1:ns)
      rs( ro     , ro+1:ns, ro+1:ns) = rs( rc     , ro+1:ns, ro+1:ns)
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
        rs(i,j,k) = (c1* rs(i,j0,k0) + c2 * rs(i0,j,k0) + c3 * rs(i0,j0,k)) / (c1+c2+c3)
      end do
      end do
      end do
    end if VERTEX_8

  end subroutine FillIrregularZones
