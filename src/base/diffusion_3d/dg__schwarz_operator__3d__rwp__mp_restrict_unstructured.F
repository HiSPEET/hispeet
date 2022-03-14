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

    real(RWP) :: ve(size(v,1)), vv

    integer(IXS) :: orientation
    integer :: c0, c1, e, i, j, k, l, m, n, nn, np
    integer :: s1, s2, s3, t1, t2, t3
    integer :: oi(2,3), ti(2,3)

    associate(po => this%po, no => this%no)

      ! start transfer .........................................................

      call buf_v % Transfer(mesh, v, tag=1000)

      ! assign core regions ....................................................

      np = po + 1    ! number of element points per direction
      c0 = no + 1    ! index of first core point within subdomain
      c1 = c0 + po   ! index of last  core point within subdomain

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

      call buf_v % Merge(v)

      ! fill overlap regions ...................................................

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element)

          ! face 1: ξ = -1 <--> west
          m = element(e) % face(1) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ np - no , po ]  ! original index range
            oi(:,2) = [ 0       , po ]  ! valid for
            oi(:,3) = [ 0       , po ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            ! transformed start indices and strides
            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 0, po
              do j = 0, po
                do i = 1, no
                  vs(i, c0 + j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! face 2:  ξ = 1 <--> east
          m = element(e) % face(2) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ 0 , no - 1 ]  ! original index range
            oi(:,2) = [ 0 , po     ]  ! valid for
            oi(:,3) = [ 0 , po     ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 0, po
              do j = 0, po
                do i = 1, no
                  vs(c1 + i, c0 + j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! face 3: η = -1 <--> south
          m = element(e) % face(3) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ 0       , po ]  ! original index range
            oi(:,2) = [ np - no , po ]  ! valid for
            oi(:,3) = [ 0       , po ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 0, po
              do j = 1, no
                do i = 0, po
                  vs(c0 + i, j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! face 4: η = 1 <--> north
          m = element(e) % face(4) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ 0 , po     ]  ! original index range
            oi(:,2) = [ 0 , no - 1 ]  ! valid for
            oi(:,3) = [ 0 , po     ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 0, po
              do j = 1, no
                do i = 0, po
                  vs(c0 + i, c1 + j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! face 5: ζ = -1 <--> bottom
          m = element(e) % face(5) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ 0       , po ]  ! original index range
            oi(:,2) = [ 0       , po ]  ! valid for
            oi(:,3) = [ np - no , po ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 1, no
              do j = 0, po
                do i = 0, po
                  vs(c0 + i, c0 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! face 6: ζ = 1 <--> top
          m = element(e) % face(6) % i_neighbor
          if (m > 0) then
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            oi(:,1) = [ 0 , po     ]  ! original index range
            oi(:,2) = [ 0 , po     ]  ! valid for
            oi(:,3) = [ 0 , no - 1 ]  ! aligned neighbor

            call TransformIndexRange(orientation, np, 0, oi, 0, ti)

            t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
            t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
            t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

            do k = 1, no
              do j = 0, po
                do i = 0, po
                  vs(c0 + i, c0 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                  t1 = t1 + s1
                end do
                t2 = t2 + s2
              end do
              t3 = t3 + s3
            end do
          end if

          ! edge  1: (x1  , south, bottom)
          nn = element(e) % edge( 1) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 0, po
                    vs(c0 + i, j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 1) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, no, no, e) = real(ve, RWP)

            end if
          end if

          ! edge  2: (x1  , north, bottom)
          nn = element(e) % edge( 2) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 0, po
                    vs(c0 + i, c1 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 2) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, c1 + 1, no, e) = real(ve, RWP)

            end if
          end if

          ! edge  3: (x1  , south, top   )
          nn = element(e) % edge( 3) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 0, po
                    vs(c0 + i, j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 3) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, no, c1 + 1, e) = real(ve, RWP)

            end if
          end if

          ! edge  4: (x1  , north, top   )
          nn = element(e) % edge( 4) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 0, po
                    vs(c0 + i, c1 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 4) % orientation, &
                                    v, ve                                )
              vs(c0 : c0 + po, c1 + 1, c1 + 1, e) = real(ve, RWP)

            end if
          end if

          ! edge  5: (west, x2   , bottom)
          nn = element(e) % edge( 5) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 0, po
                  do i = 1, no
                    vs(i, c0 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 5) % orientation, &
                                    v, ve                                )
              vs(no, c0 : c0 + po, no, e) = real(ve, RWP)

            end if
          end if

          ! edge  6: (east, x2   , bottom)
          nn = element(e) % edge( 6) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 0, po
                  do i = 1, no
                    vs(c1 + i, c0 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 6) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c0 : c0 + po, no, e) = real(ve, RWP)

            end if
          end if

          ! edge  7: (west, x2   , top   )
          nn = element(e) % edge( 7) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 0, po
                  do i = 1, no
                    vs(i, c0 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 7) % orientation, &
                                    v, ve                                )
              vs(no, c0 : c0 + po, c1 + 1, e) = real(ve, RWP)

            end if
          end if

          ! edge  8: (east, x2   , top   )
          nn = element(e) % edge( 8) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 0, po
                  do i = 1, no
                    vs(c1 + i, c0 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 8) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c0 : c0 + po, c1 + 1, e) = real(ve, RWP)

            end if
          end if

          ! edge  9: (west, south, x3    )
          nn = element(e) % edge( 9) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge( 9) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 0, po
                do j = 1, no
                  do i = 1, no
                    vs(i, j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge( 9) % orientation, &
                                    v, ve                                )
              vs(no, no, c0 : c0 + po, e) = real(ve, RWP)

            end if
          end if

          ! edge 10: (east, south, x3    )
          nn = element(e) % edge(10) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(10) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 0, po
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(10) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, no, c0 : c0 + po, e) = real(ve, RWP)

            end if
          end if

          ! edge 11: (west, north, x3    )
          nn = element(e) % edge(11) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(11) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 0, po
                do j = 1, no
                  do i = 1, no
                    vs(i, c1 + j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(11) % orientation, &
                                    v, ve                                )
              vs(no, c1 + 1, c0 : c0 + po, e) = real(ve, RWP)

            end if
          end if

          ! edge 12: (east, north, x3    )
          nn = element(e) % edge(12) % n_neighbor
          if (nn > 0) then
            m = element(e) % edge(12) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 0, po
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, c1 + j, c0 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectEdgeData( element(e) % neighbor(m:n), element, &
                                    element(e) % edge(12) % orientation, &
                                    v, ve                                )
              vs(c1 + 1, c1 + 1, c0 : c0 + po, e) = real(ve, RWP)

            end if
          end if

          ! vertex 1: (west, south, bottom)
          nn = element(e) % vertex(1) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(1) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(i, j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, no, no, e) = real(vv, RWP)

            end if
          end if

          ! vertex 2: (east, south, bottom)
          nn = element(e) % vertex(2) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(2) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, no, no, e) = real(vv, RWP)

            end if
          end if

          ! vertex 3: (west, north, bottom)
          nn = element(e) % vertex(3) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(3) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(i, c1 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, c1 + 1, no, e) = real(vv, RWP)

            end if
          end if

          ! vertex 4: (east, north, bottom)
          nn = element(e) % vertex(4) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(4) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, c1 + j, k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, c1 + 1, no, e) = real(vv, RWP)

            end if
          end if

          ! vertex 5: (west, south, top   )
          nn = element(e) % vertex(5) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(5) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(i, j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, no, c1 + 1, e) = real(vv, RWP)

            end if
          end if

          ! vertex 6: (east, south, top   )
          nn = element(e) % vertex(6) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(6) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, no, c1 + 1, e) = real(vv, RWP)

            end if
          end if

          ! vertex 7: (west, north, top   )
          nn = element(e) % vertex(7) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(7) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(i, c1 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(no, c1 + 1, c1 + 1, e) = real(vv, RWP)

            end if
          end if

          ! vertex 8: (east, north, top   )
          nn = element(e) % vertex(8) % n_neighbor
          if (nn > 0) then
            m = element(e) % vertex(8) % i_neighbor
            orientation = element(e) % neighbor(m) % orientation

            if (nn == 1 .and. orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

              l = element(e) % neighbor(m) % id
              do k = 1, no
                do j = 1, no
                  do i = 1, no
                    vs(c1 + i, c1 + j, c1 + k, e) = real(v(t1, t2, t3, l), RWP)
                    t1 = t1 + s1
                  end do
                  t2 = t2 + s2
                end do
                t3 = t3 + s3
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), v, vv)
              vs(c1 + 1, c1 + 1, c1 + 1, e) = real(vv, RWP)

            end if
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
