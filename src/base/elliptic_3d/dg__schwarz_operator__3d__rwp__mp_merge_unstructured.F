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

  module subroutine Merge_Unstructured(this, mesh, buf_us, us, u)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_us
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: u(0:,0:,0:,:) !< subdomain solutions
    real(RWP),      intent(inout) :: us(:,:,:,:)   !< mesh variable

    real(RNP) :: ue(size(u,1)), uv

    integer(IXS) :: orientation
    integer :: e, i, j, k, l, m, n
    integer :: ii, jj, kk, nn, np, ns
    integer :: c0, c1, lo, ro
    integer :: o(3), s(3,3), t(3)
    logical :: merge_edge(12), merge_vertex(8)

    associate(po => this%po, no => this%no)

      np = po + 1       ! number of element   points per direction
      ns = np + 2 * no  ! number of subdomain points per direction

      lo = no - 1       ! last  point of left  element overlap zone
      ro = po - lo      ! first point of right element overlap zone

      c0 = no + 1       ! first subdomain point inside core
      c1 = c0 + np      ! first subdomain point behind core

      ! start transfer .........................................................

      if (no > 0) then
        call buf_us % Transfer(mesh, us, tag=2000)
      end if

      ! core regions ...........................................................

      ! indices of first and last core points within subdomain

      !$omp do
      do e = 1, mesh % n_elem
        do k = 0, po
        do j = 0, po
        do i = 0, po
          u(i,j,k,e) = u(i,j,k,e) + us(c0 + i, c0 + j, c0 + k, e)
        end do
        end do
        end do
      end do

      ! merge remote data to ghosts ............................................

      if (no > 0) then
        call buf_us % Merge(us)
      else
        return
      end if

      ! overlap regions ........................................................

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element)

          merge_edge   = .true.
          merge_vertex = .true.

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
              u(i, j, k, e) = u(i, j, k, e) + us(ii, jj, kk, l)
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
              u(ro + i, j, k, e) = u(ro + i, j, k, e) + us(ii, jj, kk, l)
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
              u(i, j, k, e) = u(i, j, k, e) + us(ii, jj, kk, l)
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
              u(i, ro + j, k, e) = u(i, ro + j, k, e) + us(ii, jj, kk, l)
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
              u(i, j, k, e) = u(i, j, k, e) + us(ii, jj, kk, l)
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
              u(i, j, ro + k, e) = u(i, j, ro + k, e) + us(ii, jj, kk, l)
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
                u(i, j, k, e) = u(i, j, k, e) + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 1) % orientation, &
                                    no, us, ue                           )
              u( :,  0,  0, e) = u( :,  0,  0, e) + ue

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
                u(i, ro + j, k, e) = u(i, ro + j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 2) % orientation, &
                                    no, us, ue                           )
              u( :, po,  0, e) = u( :, po,  0, e) + ue

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
                u(i, j, ro + k, e) = u(i, j, ro + k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 3) % orientation, &
                                    no, us, ue                           )
              u( :,  0, po, e) = u( :,  0, po, e) + ue

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
                u(i, ro + j, ro + k, e) = u(i, ro + j, ro + k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 4) % orientation, &
                                    no, us, ue                           )
              u( :, po, po, e) = u( :, po, po, e) + ue

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
                u(i, j, k, e) = u(i, j, k, e) &
                              + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 5) % orientation, &
                                    no, us, ue                           )
              u(  0, :,  0, e) = u(  0, :,  0, e) + ue

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
                u(ro + i, j, k, e) = u(ro + i, j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 6) % orientation, &
                                    no, us, ue                           )
              u( po, :,  0, e) = u( po, :,  0, e) + ue

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
                u(i, j, ro + k, e) = u(i, j, ro + k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 7) % orientation, &
                                    no, us, ue                           )
              u(  0, :, po, e) = u(  0, :, po, e) + ue

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
                u(ro + i, j, ro + k, e) = u(ro + i, j, ro + k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 8) % orientation, &
                                    no, us, ue                           )
              u( po, :, po, e) = u( po, :, po, e) + ue

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
                u(i, j, k, e) = u(i, j, k, e) &
                              + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 9) % orientation, &
                                    no, us, ue                           )
              u(  0,  0, :, e) = u(  0,  0, :, e) + ue

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
                u(ro + i, j, k, e) = u(ro + i, j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(10) % orientation, &
                                    no, us, ue                           )
              u( po,  0, :, e) = u( po,  0, :, e) + ue

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
                u(i, ro + j, k, e) = u(i, ro + j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(11) % orientation, &
                                    no, us, ue                           )
              u(  0, po, :, e) = u(  0, po, :, e) + ue

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
                u(ro + i, ro + j, k, e) = u(ro + i, ro + j, k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(12) % orientation, &
                                    no, us, ue                           )
              u( po, po, :, e) = u( po, po, :, e) + ue

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
                u(i, j, k, e) = u(i, j, k, e) &
                              + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u(  0,  0,  0, e) = u(  0,  0,  0, e) + uv

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
                u(ro + i, j, k, e) = u(ro + i, j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u( po,  0,  0, e) = u( po,  0,  0, e) + uv

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
                u(i, ro + j, k, e) = u(i, ro + j, k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u(  0, po,  0, e) = u(  0, po,  0, e) + uv

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
                u(ro + i, ro + j, k, e) = u(ro + i, ro + j, k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u( po, po,  0, e) = u( po, po,  0, e) + uv

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
                u(i, j, ro + k, e) = u(i, j, ro + k, e) &
                                   + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u(  0,  0, po, e) = u(  0,  0, po, e) + uv

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
                u(ro + i, j, ro + k, e) = u(ro + i, j, ro + k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u( po,  0, po, e) = u( po,  0, po, e) + uv

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
                u(i, ro + j, ro + k, e) = u(i, ro + j, ro + k, e) &
                                        + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u(  0, po, po, e) = u(  0, po, po, e) + uv

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
                u(ro + i, ro + j, ro + k, e) = u(ro + i, ro + j, ro + k, e) &
                                             + us(ii, jj, kk, l)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), no, us, uv)
              u( po, po, po, e) = u( po, po, po, e) + uv

            end if
          end if

        end associate
      end do

    end associate

  end subroutine Merge_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(neighbor, element, o, no, us, ue)
    class(MeshElementNeighbor_3D), intent(in)  :: neighbor(:)
    class(MeshElement_3D),         intent(in)  :: element(:)
    integer(IXS),                  intent(in)  :: o
    integer,                       intent(in)  :: no
    real(RWP),                     intent(in)  :: us(:,:,:,:)
    real(RNP),                     intent(out) :: ue(:)

    real(RNP) :: un(size(ue))
    integer :: nn, np
    integer :: i0, ip
    integer :: e, l, n

    nn = size(neighbor)
    np = size(ue)
    i0 = no + 1
    ip = no + np

    ue = 0

    do n = 1, nn

      l = neighbor(n) % id
      e = neighbor(n) % component - 6

      select case(e)
      case( 1)
        un = us( i0:ip, i0, i0, l)
      case( 2)
        un = us( i0:ip, ip, i0, l)
      case( 3)
        un = us( i0:ip, i0, ip, l)
      case( 4)
        un = us( i0:ip, ip, ip, l)
      case( 5)
        un = us( i0, i0:ip, i0, l)
      case( 6)
        un = us( ip, i0:ip, i0, l)
      case( 7)
        un = us( i0, i0:ip, ip, l)
      case( 8)
        un = us( ip, i0:ip, ip, l)
      case( 9)
        un = us( i0, i0, i0:ip, l)
      case(10)
        un = us( ip, i0, i0:ip, l)
      case(11)
        un = us( i0, ip, i0:ip, l)
      case(12)
        un = us( ip, ip, i0:ip, l)
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
