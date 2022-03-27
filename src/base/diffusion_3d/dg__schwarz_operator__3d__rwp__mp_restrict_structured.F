  !-----------------------------------------------------------------------------
  !> Restrict a structured mesh variable `r` to subdomains
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

  module subroutine Restrict_Structured(this, mesh, buf_r, r, rs, sgn)
    class(DG_SchwarzOperator_3D), intent(in) :: this
    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: buf_r
    class(Mesh_3D), intent(in)    :: mesh          !< mesh partition
    real(RNP),      intent(inout) :: r(0:,0:,0:,:) !< extended mesh variable
    real(RWP),      intent(out)   :: rs(:,:,:,:)   !< restricted variable
    integer, optional, intent(in) :: sgn           !< sign of `r` {+1,-1} [+1]

    real(RWP) :: c_sgn
    integer   :: e, i, j, k, l, m
    integer   :: c0, c1, o0, o1

    if (present(sgn)) then
      c_sgn = sgn
    else
      c_sgn = 1
    end if

    associate(po => this%po, no => this%no)

      ! start transfer .........................................................

      if (no > 0) then
        call buf_r % Transfer(mesh, r, tag=1000)
      end if

      ! assign core regions ....................................................

      ! indices of first and last core points within subdomain
      c0 = no + 1
      c1 = c0 + po

      !$omp do
      do e = 1, mesh % n_elem

        ! ensure that no NANs remain in parts for which no neighbor data exists
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

      ! offsets of element point indices
      o0 = -1
      o1 = po - no

      !$omp do
      do e = 1, mesh % n_elem
        associate(element => mesh % element(e))

          ! face 1: ξ = -1 <--> west
          m = element % face(1) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 0, po
            do i = 1, no
              rs(i, c0 + j, c0 + k, e) = c_sgn * real(r(o1 + i, j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! face 2:  ξ = 1 <--> east
          m = element % face(2) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 0, po
            do i = 1, no
              rs(c1 + i, c0 + j, c0 + k, e) = c_sgn * real(r(o0 + i, j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! face 3: η = -1 <--> south
          m = element % face(3) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 0, po
              rs(c0 + i, j, c0 + k, e) = c_sgn * real(r(i, o1 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! face 4: η = 1 <--> north
          m = element % face(4) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 0, po
              rs(c0 + i, c1 + j, c0 + k, e) = c_sgn * real(r(i, o0 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! face 5: ζ = -1 <--> bottom
          m = element % face(5) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 0, po
              rs(c0 + i, c0 + j, k, e) = c_sgn * real(r(i, j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! face 6: ζ = 1 <--> top
          m = element % face(6) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 0, po
              rs(c0 + i, c0 + j, c1 + k, e) = c_sgn * real(r(i, j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  1: (x1  , south, bottom)
          m = element % edge( 1) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 0, po
              rs(c0 + i, j, k, e) = c_sgn * real(r(i, o1 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  2: (x1  , north, bottom)
          m = element % edge( 2) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 0, po
              rs(c0 + i, c1 + j, k, e) = c_sgn * real(r(i, o0 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  3: (x1  , south, top   )
          m = element % edge( 3) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 0, po
              rs(c0 + i, j, c1 + k, e) = c_sgn * real(r(i, o1 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  4: (x1  , north, top   )
          m = element % edge( 4) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 0, po
              rs(c0 + i, c1 + j, c1 + k, e) = c_sgn * real(r(i, o0 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  5: (west, x2   , bottom)
          m = element % edge( 5) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 1, no
              rs(i, c0 + j, k, e) = c_sgn * real(r(o1 + i, j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  6: (east, x2   , bottom)
          m = element % edge( 6) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 1, no
              rs(c1 + i, c0 + j, k, e) = c_sgn * real(r(o0 + i, j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  7: (west, x2   , top   )
          m = element % edge( 7) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 1, no
              rs(i, c0 + j, c1 + k, e) = c_sgn * real(r(o1 + i, j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  8: (east, x2   , top   )
          m = element % edge( 8) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 0, po
            do i = 1, no
              rs(c1 + i, c0 + j, c1 + k, e) = c_sgn * real(r(o0 + i, j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge  9: (west, south, x3    )
          m = element % edge( 9) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 1, no
              rs(i, j, c0 + k, e) = c_sgn * real(r(o1 + i, o1 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge 10: (east, south, x3    )
          m = element % edge(10) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 1, no
              rs(c1 + i, j, c0 + k, e) = c_sgn * real(r(o0 + i, o1 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge 11: (west, north, x3    )
          m = element % edge(11) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 1, no
              rs(i, c1 + j, c0 + k, e) = c_sgn * real(r(o1 + i, o0 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! edge 12: (east, north, x3    )
          m = element % edge(12) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 0, po
            do j = 1, no
            do i = 1, no
              rs(c1 + i, c1 + j, c0 + k, e) = c_sgn * real(r(o0 + i, o0 + j, k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 1: (west, south, bottom)
          m = element % vertex(1) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(i, j, k, e) = c_sgn * real(r(o1 + i, o1 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 2: (east, south, bottom)
          m = element % vertex(2) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(c1 + i, j, k, e) = c_sgn * real(r(o0 + i, o1 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 3: (west, north, bottom)
          m = element % vertex(3) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(i, c1 + j, k, e) = c_sgn * real(r(o1 + i, o0 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 4: (east, north, bottom)
          m = element % vertex(4) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(c1 + i, c1 + j, k, e) = c_sgn * real(r(o0 + i, o0 + j, o1 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 5: (west, south, top   )
          m = element % vertex(5) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(i, j, c1 + k, e) = c_sgn * real(r(o1 + i, o1 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 6: (east, south, top   )
          m = element % vertex(6) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(c1 + i, j, c1 + k, e) = c_sgn * real(r(o0 + i, o1 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 7: (west, north, top   )
          m = element % vertex(7) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(i, c1 + j, c1 + k, e) = c_sgn * real(r(o1 + i, o0 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

          ! vertex 8: (east, north, top   )
          m = element % vertex(8) % i_neighbor
          if (m > 0) then
            l = element % neighbor(m) % id
            do k = 1, no
            do j = 1, no
            do i = 1, no
              rs(c1 + i, c1 + j, c1 + k, e) = c_sgn * real(r(o0 + i, o0 + j, o0 + k, l), RWP)
            end do
            end do
            end do
          end if

        end associate
      end do

    end associate

  end subroutine Restrict_Structured
