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
          rs(c0+i, c0+j, c0+k, e) = c_sgn * real(r(i,j,k,e), RWP)
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
              rs(1+i, c0+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
              rs(c1+i, c0+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
              rs(c0+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
              rs(c0+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
              rs(c0+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
              rs(c0+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
            end do
            end do
            end do
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
                rs(c0+i, 1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 1, m, n, r, re)
              rs(c0 : c0+po, no, no, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0, ro, ro ]

                ! additional entries from first neighbor (k = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = lo
                do j = 0, lo-1
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, 1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (j = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = lo
                do k = 0, lo-1
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, 1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c0+i, c1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 2, m, n, r, re)
              rs(c0 : c0+po, c1, no, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0,  0, ro ]

                ! additional entries from first neighbor (j = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = 0
                do k = 0, lo-1
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, c1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (k = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = lo
                do j = 1, lo
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, c1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c0+i, 1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 3, m, n, r, re)
              rs(c0 : c0+po, no, c1, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0, ro,  0 ]

                ! additional entries from first neighbor (j = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = lo
                do k = 1, lo
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, 1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (k = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = 0
                do j = 0, lo-1
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, 1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c0+i, c1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 4, m, n, r, re)
              rs(c0 : c0+po, c1, c1, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0,  0,  0 ]

                ! additional entries from first neighbor (k = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = 0
                do j = 1, lo
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, c1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (j = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = 0
                do k = 1, lo
                do i = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c0+i, c1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 5, m, n, r, re)
              rs(no, c0 : c0+po, no, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [ ro,  0, ro ]

                ! additional entries from first neighbor (i = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = lo
                do k = 0, lo-1
                do j = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (k = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = lo
                do j = 0, po
                do i = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 6, m, n, r, re)
              rs(c1, c0 : c0+po, no, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0,  0, ro ]

                ! additional entries from first neighbor (k = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = lo
                do j = 0, po
                do i = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (i = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = 0
                do k = 0, lo-1
                do j = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c0+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 7, m, n, r, re)
              rs(no, c0 : c0+po, c1, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [ ro,  0,  0 ]

                ! additional entries from first neighbor (k = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = 0
                do j = 0, po
                do i = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (i = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = lo
                do k = 1, lo
                do j = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 8, m, n, r, re)
              rs(c1, c0 : c0+po, c1, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0,  0,  0 ]

                ! additional entries from first neighbor (i = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = 0
                do k = 1, lo
                do j = 0, po
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (k = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                k = 0
                do j = 0, po
                do i = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c0+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 9, m, n, r, re)
              rs(no, no, c0 : c0+po, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [ ro, ro,  0 ]

                ! additional entries from first neighbor (j = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = lo
                do k = 0, po
                do i = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (i = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = lo
                do k = 0, po
                do j = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 10, m, n, r, re)
              rs(c1, no, c0 : c0+po, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0, ro,  0 ]

                ! additional entries from first neighbor (i = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = 0
                do k = 0, po
                do j = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (j = lo)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = lo
                do k = 0, po
                do i = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, 1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 11, m, n, r, re)
              rs(no, c1, c0 : c0+po, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [ ro,  0,  0 ]

                ! additional entries from first neighbor (i = lo)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = lo
                do k = 0, po
                do j = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (j = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = 0
                do k = 0, po
                do i = 0, lo-1
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(c1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectEdgeData(element(e), 12, m, n, r, re)
              rs(c1, c1, c0 : c0+po, e) = c_sgn * real(re, RWP)

              if (no > 1) then

                o = [  0,  0,  0 ]

                ! additional entries from first neighbor (j = 0)
                orientation = element(e) % neighbor(m) % orientation
                l = element(e) % neighbor(m) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                j = 0
                do k = 0, po
                do i = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

                ! additional entries from last neighbor (i = 0)
                orientation = element(e) % neighbor(n) % orientation
                l = element(e) % neighbor(n) % id
                call TransformIndex(orientation, np, 0, o, 0, t, s)
                i = 0
                do k = 0, po
                do j = 1, lo
                  ii = t(1) + s(1,1) * i + s(2,1) * j + s(3,1) * k
                  jj = t(2) + s(1,2) * i + s(2,2) * j + s(3,2) * k
                  kk = t(3) + s(1,3) * i + s(2,3) * j + s(3,3) * k
                  rs(c1+i, c1+j, c0+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
                end do
                end do

              end if

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
                rs(1+i, 1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, no, no, e) = c_sgn * real(rv, RWP)

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
                rs(c1+i, 1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, no, no, e) = c_sgn * real(rv, RWP)

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
                rs(1+i, c1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, c1, no, e) = c_sgn * real(rv, RWP)

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
                rs(c1+i, c1+j, 1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, c1, no, e) = c_sgn * real(rv, RWP)

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
                rs(1+i, 1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, no, c1, e) = c_sgn * real(rv, RWP)

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
                rs(c1+i, 1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, no, c1, e) = c_sgn * real(rv, RWP)

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
                rs(1+i, c1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(no, c1, c1, e) = c_sgn * real(rv, RWP)

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
                rs(c1+i, c1+j, c1+k, e) = c_sgn * real(r(ii, jj, kk, l), RWP)
              end do
              end do
              end do

            else

              n = m + nn - 1
              call CollectVertexData(element(e)%neighbor(m:n), r, rv)
              rs(c1, c1, c1, e) = c_sgn * real(rv, RWP)

            end if
          end if

        end associate
      end do

    end associate

  end subroutine Restrict_Unstructured

  !-----------------------------------------------------------------------------
  !> Collect, align and average neighbor edge data

  pure subroutine CollectEdgeData(element, e, m, n, r, re)
    class(MeshElement_3D), intent(in)  :: element    !< current element
    integer,               intent(in)  :: e          !< element edge ID
    integer,               intent(in)  :: m          !< first neighbor ID
    integer,               intent(in)  :: n          !< last neighbor ID
    real(RNP),             intent(in)  :: r(:,:,:,:) !< element data
    real(RNP),             intent(out) :: re(:)      !< avg neighbor edge data

    real(RNP) :: rk(size(re))
    integer   :: k, l, np

    np = size(re)
    re = 0

    do k = m, n

      l = element % neighbor(k) % id

      select case(ElementEdgeID(element % neighbor(k) % component))
      case( 1)
        rk = r( :, 1, 1,l)
      case( 2)
        rk = r( :,np, 1,l)
      case( 3)
        rk = r( :, 1,np,l)
      case( 4)
        rk = r( :,np,np,l)
      case( 5)
        rk = r( 1, :, 1,l)
      case( 6)
        rk = r(np, :, 1,l)
      case( 7)
        rk = r( 1, :,np,l)
      case( 8)
        rk = r(np, :,np,l)
      case( 9)
        rk = r( 1, 1, :,l)
      case(10)
        rk = r(np, 1, :,l)
      case(11)
        rk = r( 1,np, :,l)
      case(12)
        rk = r(np,np, :,l)
      end select

      if (element % NeigborEdgeIsAligned(e, k)) then
        re = re + rk
      else
        re = re + rk(np:1:-1)
      end if

    end do

    re = re / (n - m + 1)

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
