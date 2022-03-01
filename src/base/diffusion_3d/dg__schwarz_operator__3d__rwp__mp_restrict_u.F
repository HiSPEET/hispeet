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

  module subroutine Restrict_U(this, mesh, buf_v, v, vs)
    class(DG_SchwarzOperator_3D_RWP), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_v
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: v(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: vs(:,:,:,:)   !< restricted variable

    integer(IXS) :: orientation
    integer :: c0, c1, e, i, j, k, l, m, nn, np
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

!!!      ! offsets of element point indices
!!!      o0 = -1
!!!      o1 = po - no

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
          select case(nn)
          case(1)
            m = element(e) % edge( 1) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
!### identify neighbor edge
!### identify wether to flip or not
!### put info and data into routine
!### extract end possibly flip 1D edge data
!### assign the edge data
            end if
          end select

          ! edge  2: (x1  , north, bottom)
          nn = element(e) % edge( 2) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 2) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  3: (x1  , south, top   )
          nn = element(e) % edge( 3) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 3) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  4: (x1  , north, top   )
          nn = element(e) % edge( 4) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 4) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  5: (west, x2   , bottom)
          nn = element(e) % edge( 5) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 5) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  6: (east, x2   , bottom)
          nn = element(e) % edge( 6) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 6) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  7: (west, x2   , top   )
          nn = element(e) % edge( 7) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 7) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  8: (east, x2   , top   )
          nn = element(e) % edge( 8) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 8) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge  9: (west, south, x3    )
          nn = element(e) % edge( 9) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge( 9) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge 10: (east, south, x3    )
          nn = element(e) % edge(10) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge(10) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge 11: (west, north, x3    )
          nn = element(e) % edge(11) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge(11) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! edge 12: (east, north, x3    )
          nn = element(e) % edge(12) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % edge(12) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 1: (west, south, bottom)
          nn = element(e) % vertex(1) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(1) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 2: (east, south, bottom)
          nn = element(e) % vertex(2) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(2) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 3: (west, north, bottom)
          nn = element(e) % vertex(3) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(3) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 4: (east, north, bottom)
          nn = element(e) % vertex(4) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(4) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ np - no , po     ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 5: (west, south, top   )
          nn = element(e) % vertex(5) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(5) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 6: (east, south, top   )
          nn = element(e) % vertex(6) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(6) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ np - no , po     ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 7: (west, north, top   )
          nn = element(e) % vertex(7) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(7) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ np - no , po     ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

          ! vertex 8: (east, north, top   )
          nn = element(e) % vertex(8) % n_neighbor
          select case(nn)
          case(1)
            m = element(e) % vertex(8) % i_neighbor
            l = element(e) % neighbor(m) % id
            orientation = element(e) % neighbor(m) % orientation

            if (orientation > 0_IXS) then

              oi(:,1) = [ 0       , no - 1 ]  ! original index range
              oi(:,2) = [ 0       , no - 1 ]  ! valid for
              oi(:,3) = [ 0       , no - 1 ]  ! aligned neighbor

              call TransformIndexRange(orientation, np, 0, oi, 0, ti)

              t1 = ti(1,1);  s1 = sign(1, ti(2,1) - ti(1,1))
              t2 = ti(1,2);  s2 = sign(1, ti(2,2) - ti(1,2))
              t3 = ti(1,3);  s3 = sign(1, ti(2,3) - ti(1,3))

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
            end if
          end select

        end associate
      end do

    end associate

  end subroutine Restrict_U
