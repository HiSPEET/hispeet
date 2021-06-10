!> summary:  Schwarz method for elliptic equation with equidistant IP/DG
!> author:   Joerg Stiller
!> date:     2018/12/12
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Schwarz method for elliptic equation with equidistant IP/DG
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_SchwarzMethod
  use Execution_Control, only: Error
  use TPO__Schwarz__3D_CA
  use TPO__Schwarz__3D_CI
  use CART__Mesh_Partition
  use CART__Element_Transfer_Buffer
  implicit none

contains

!-------------------------------------------------------------------------------
!> Element-centered overlapping Schwarz method with constant coefficients

module subroutine SchwarzMethod(this, u, f, i_max, r_red, r_max, ni)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations

  ! local variables ............................................................

  !!! save attribute serves to enforce sharing with OpenMP !!!

  real(RNP), allocatable, save :: r (:,:,:,:)  ! residual, including ghosts
  real(RNP), allocatable, save :: u_s(:,:,:,:) ! solution of subsystems
  real(RNP), allocatable, save :: f_s(:,:,:,:) ! RHS of subsystems

  type(ElementTransferBuffer), asynchronous, allocatable, save :: buf_r
  type(ElementTransferBuffer), asynchronous, allocatable, save :: buf_u_s

  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP) :: rr
  integer :: ne, ng, np(3), no(3), ns(3), nc
  integer :: i

  ! initialization .............................................................

  associate(mesh => this%mesh, eop => this%eop, schwarz => this%schwarz)

    ! dimensions
    ne = mesh % ne       ! number of local elements
    ng = mesh % ng       ! number of ghost elements
    np = eop % po + 1    ! number of element points per direction
    no = schwarz % no    ! overlapped point layers
    ns = np + 2*no       ! subdomain extensions
    nc = schwarz % nc    ! number of 1D subdomain configurations

    ! workspace
    !$omp single
    allocate( r   (np(1), np(2), np(3), ne+ng) )
    allocate( u_s (ns(1), ns(2), ns(3), ne+ng) )
    allocate( f_s (ns(1), ns(2), ns(3), ne   ) )
    allocate( buf_r )
    allocate( buf_u_s )
    !$omp end single

    call SetArray(r  (:,:,:,ne+1:ne+ng), ZERO)
    call SetArray(u_s(:,:,:,ne+1:ne+ng), ZERO)

    ! transfer buffers
    call buf_r   % New(mesh, r,  no)
    call buf_u_s % New(mesh, u_s, no)

    ! termination condition
    !$omp single
    if (present(r_max)) then
      rr_term = max(ZERO, r_max)**2
    else
      rr_term = ZERO
    end if
    !$omp end single

    ! Schwarz iterations .......................................................

    do i = 1, i_max

      call this % Residual(u, f, r(:,:,:,:ne))

      ! termination check
      if (present(r_red)) then
        rr = ScalarProduct(r(:,:,:,:ne), r(:,:,:,:ne), mesh%comm)
        !$omp master
        if (mesh%part == 0) then
          if (i == 1) then
            rr_term  = max(rr_term, max(ZERO, sqrt(rr) * r_red)**2)
          end if
          converged = rr <= rr_term
        end if
        call XMPI_Bcast(converged, root=0, comm=mesh%comm)
        !$omp end master
        !$omp barrier
      end if
      if (converged) exit

      call RestrictToSubdomains(mesh, no, buf_r, r, f_s)

      if (schwarz % isotropic) then

         call TPO_Schwarz( schwarz % S1,  schwarz % W1,    &
                           schwarz % cfg, schwarz % D_inv, &
                           f_s, u_s                        )
      else

         call TPO_Schwarz( schwarz % S1,  schwarz % S2,  schwarz % S3, &
                           schwarz % W1,  schwarz % W2,  schwarz % W3, &
                           schwarz % cfg, schwarz % D_inv,             &
                           f_s, u_s                                    )

      end if

      call MergeFromSubdomains(mesh, no, buf_u_s, u_s, u)

    end do

  end associate

  if (present(ni)) ni = min(i, i_max)

  ! clean-up ...................................................................

  !$omp single
  if (allocated( r       )) deallocate( r       )
  if (allocated( u_s     )) deallocate( u_s     )
  if (allocated( f_s     )) deallocate( f_s     )
  if (allocated( buf_r   )) deallocate( buf_r   )
  if (allocated( buf_u_s )) deallocate( buf_u_s )
  !$omp end single

end subroutine SchwarzMethod

!===============================================================================
! Restriction and merging of subdomain data

!-------------------------------------------------------------------------------
!> Restrict a mesh variable to subdomain variable
!>
!> The strategy is as a follows:
!>
!>   1. Send `v` from linked local element to ghosts in remote partitions
!>   2. Assign `v` to subdomain core regions, i.e. 1:1 copy to `vs`
!>   3. Assign data from remote masters to ghosts
!>   4. Assign local and ghost data to overlap regions of `vs`
!>
!> The mesh variable must be dimensioned as `v(0:po,0:po,0:po,mesh%ne+mesh%ng)`.

subroutine RestrictToSubdomains(mesh, no, buf_v, v, vs)

  ! arguments ..................................................................

  class(MeshPartition), intent(in) :: mesh   !< mesh partition
  integer,              intent(in) :: no(3)  !< subdomain overlap
  class(ElementTransferBuffer), asynchronous, intent(inout) :: buf_v

  real(RNP), intent(inout) :: v(0:,0:,0:,:)  !< extended mesh variable
  real(RNP), intent(out)   :: vs(:,:,:,:)    !< restricted variable

  ! local variables ............................................................

  integer :: po
  integer :: e, l, m, n, o
  integer :: i, ie0, ie1, is0, is1
  integer :: j, je0, je1, js0, js1
  integer :: k, ke0, ke1, ks0, ks1

  ! initialization ...........................................................

  po = ubound(v,1)

!opt!  call SetArray(vs, ZERO)

  ! offsets of subdomain point indices
  is0 = no(1) + 1;  is1 = po + is0
  js0 = no(2) + 1;  js1 = po + js0
  ks0 = no(3) + 1;  ks1 = po + ks0

  ! offsets of element point indices
  ie0 = -1;  ie1 = po + 1
  je0 = -1;  je1 = po + 1
  ke0 = -1;  ke1 = po + 1

  ! start transfer .............................................................

  call buf_v % ToGhost_Transfer(mesh, v, tag=1000)

  ! assign core regions ........................................................

  !$omp do
  do e = 1, mesh%ne
    do k = 0, po
    do j = 0, po
    do i = 0, po
      vs(is0 + i, js0 + j, ks0 + k, e) = v(i,j,k,e)
    end do
    end do
    end do
  end do

  ! assign remote data to ghosts ...............................................

  ! for all ghosts: v = v[masters]
  call buf_v % ToGhost_Merge(v, alpha=ZERO, beta=ONE)

  ! fill overlap regions .......................................................

  ! offsets of subdomain point indices
  ! is0, js0, ks0 = 0 : not used
  ! is1, js1, ks1     : unchanged

  ! offsets of element point indices
  ie0 = -1;  ie1 = po - no(1)
  je0 = -1;  je1 = po - no(2)
  ke0 = -1;  ke1 = po - no(3)

  associate(element => mesh % element)

    !$omp do
    do e = 1, mesh % ne

      ! face 1: -x1 <--> west
      o = element(e) % face(1) % neighbor
      if (o > 0) then
        do k = 0, po
        do j = 0, po
        do l = 1, no(1)
          vs(l, js0 + j, ks0 + k, e) = v(ie1 + l, j, k, o)
        end do
        end do
        end do
      end if

      ! face 2:  x1 <--> east
      o = element(e) % face(2) % neighbor
      if (o > 0) then
        do k = 0, po
        do j = 0, po
        do l = 1, no(1)
          vs(is1 + l, js0 + j, ks0 + k, e) = v(ie0 + l, j, k, o)
        end do
        end do
        end do
      end if

      ! face 3: -x2 <--> south
      o = element(e) % face(3) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, m, ks0 + k, e) = v(i, je1 + m, k, o)
        end do
        end do
        end do
      end if

      ! face 4:  x2 <--> north
      o = element(e) % face(4) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, js1 + m, ks0 + k, e) = v(i, je0 + m, k, o)
        end do
        end do
        end do
      end if

      ! face 5: -x3 <--> bottom
      o = element(e) % face(5) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do i = 0, po
          vs(is0 + i, js0 + j, n, e) = v(i, j, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! face 6:  x3 <--> top
      o = element(e) % face(6) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do i = 0, po
          vs(is0 + i, js0 + j, ks1 + n, e) = v(i, j, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! edge  1: (x1  , south, bottom)
      o = element(e) % edge( 1) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, m, n, e) = v(i, je1 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! edge  2: (x1  , north, bottom)
      o = element(e) % edge( 2) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, js1 + m, n, e) = v(i, je0 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! edge  3: (x1  , south, top   )
      o = element(e) % edge( 3) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, m, ks1 + n, e) = v(i, je1 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! edge  4: (x1  , north, top   )
      o = element(e) % edge( 4) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do i = 0, po
          vs(is0 + i, js1 + m, ks1 + n, e) = v(i, je0 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! edge  5: (west, x2   , bottom)
      o = element(e) % edge( 5) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do l = 1, no(1)
          vs(l, js0 + j, n, e) = v(ie1 + l, j, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! edge  6: (east, x2   , bottom)
      o = element(e) % edge( 6) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do l = 1, no(1)
          vs(is1 + l, js0 + j, n, e) = v(ie0 + l, j, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! edge  7: (west, x2   , top   )
      o = element(e) % edge( 7) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do l = 1, no(1)
          vs(l, js0 + j, ks1 + n, e) = v(ie1 + l, j, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! edge  8: (east, x2   , top   )
      o = element(e) % edge( 8) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do j = 0, po
        do l = 1, no(1)
          vs(is1 + l, js0 + j, ks1 + n, e) = v(ie0 + l, j, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! edge  9: (west, south, x3    )
      o = element(e) % edge( 9) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, m, ks0 + k, e) = v(ie1 + l, je1 + m, k, o)
        end do
        end do
        end do
      end if

      ! edge 10: (east, south, x3    )
      o = element(e) % edge(10) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, m, ks0 + k, e) = v(ie0 + l, je1 + m, k, o)
        end do
        end do
        end do
      end if

      ! edge 11: (west, north, x3    )
      o = element(e) % edge(11) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, js1 + m, ks0 + k, e) = v(ie1 + l, je0 + m, k, o)
        end do
        end do
        end do
      end if

      ! edge 12: (east, north, x3    )
      o = element(e) % edge(12) % neighbor
      if (o > 0) then
        do k = 0, po
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, js1 + m, ks0 + k, e) = v(ie0 + l, je0 + m, k, o)
        end do
        end do
        end do
      end if

      ! vertex 1: (west, south, bottom)
      o = element(e) % vertex(1) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, m, n, e) = v(ie1 + l, je1 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 2: (east, south, bottom)
      o = element(e) % vertex(2) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, m, n, e) = v(ie0 + l, je1 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 3: (west, north, bottom)
      o = element(e) % vertex(3) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, js1 + m, n, e) = v(ie1 + l, je0 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 4: (east, north, bottom)
      o = element(e) % vertex(4) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, js1 + m, n, e) = v(ie0 + l, je0 + m, ke1 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 5: (west, south, top   )
      o = element(e) % vertex(5) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, m, ks1 + n, e) = v(ie1 + l, je1 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 6: (east, south, top   )
      o = element(e) % vertex(6) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, m, ks1 + n, e) = v(ie0 + l, je1 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 7: (west, north, top   )
      o = element(e) % vertex(7) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(l, js1 + m, ks1 + n, e) = v(ie1 + l, je0 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

      ! vertex 8: (east, north, top   )
      o = element(e) % vertex(8) % neighbor
      if (o > 0) then
        do n = 1, no(3)
        do m = 1, no(2)
        do l = 1, no(1)
          vs(is1 + l, js1 + m, ks1 + n, e) = v(ie0 + l, je0 + m, ke0 + n, o)
        end do
        end do
        end do
      end if

    end do

  end associate

  ! clean-up ...................................................................

  call buf_v % ToGhost_Finish()

end subroutine RestrictToSubdomains

!-------------------------------------------------------------------------------
!> Merge subdomain contributions into mesh variable
!>
!> The strategy is as a follows:
!>
!>   1. Apply weights to `vs`, taking into account the boundary conditions
!>   2. Send `vs` from masters to ghosts
!>   3. Add `vs` from subdomain core regions to `v`
!>   4. Assign data from remote masters to `vs` ghosts
!>   5. Merge `vs` from local and ghost overlap regions into `v`
!>
!> For transfer from masters to ghosts the subdomain variable must be shaped
!> `vs(ns(1), ns(2), ns(3), mesh%ne + mesh%ng)`.

subroutine MergeFromSubdomains(mesh, no, buf_vs, vs, v)

  ! arguments ..................................................................

  class(MeshPartition), intent(in) :: mesh   !< mesh partition
  integer,              intent(in) :: no(3)  !< subdomain overlap
  class(ElementTransferBuffer), asynchronous, intent(inout) :: buf_vs

  real(RNP), intent(inout) :: vs(:,:,:,:)    !< subdomain variables
  real(RNP), intent(inout) :: v(0:,0:,0:,:)  !< mesh variable

  ! local variables ............................................................

  integer :: po, ns(3)
  integer :: e, l, m, n, o
  integer :: i, ie0, ie1, is0, is1
  integer :: j, je0, je1, js0, js1
  integer :: k, ke0, ke1, ks0, ks1

  ! initialization .............................................................

  ! dimensions
  po = ubound(v,1)
  ns = po + 1 + 2 * no

  ! offsets of subdomain point indices
  is0 = no(1) + 1;  is1 = ns(1) - no(1)
  js0 = no(2) + 1;  js1 = ns(2) - no(2)
  ks0 = no(3) + 1;  ks1 = ns(3) - no(3)

  ! offsets of element point indices
  ie0 = -1;  ie1 = po - no(1)
  je0 = -1;  je1 = po - no(2)
  ke0 = -1;  ke1 = po - no(3)

  ! send vs to ghosts ..........................................................

  call buf_vs % ToGhost_Transfer(mesh, vs, tag=2000)

  ! add vs core regions to v .................................................

  !$omp do
  do e = 1, mesh%ne
    do k = 0, po
    do j = 0, po
    do i = 0, po
      v(i,j,k,e) = v(i,j,k,e) + vs(is0 + i, js0 + j, ks0 + k, e)
    end do
    end do
    end do
  end do

  ! copy received data to ghosts .............................................

  call buf_vs % ToGhost_Merge(vs, alpha=ZERO, beta=ONE)

  ! merge vs from overlap regions into v .....................................

    associate(element => mesh % element)

      !$omp do
      do e = 1, mesh % ne

        ! face 1: -x1 <--> west
        o = element(e) % face(1) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, k, e) = v (ie0 + l,       j,       k, e) &
                                + vs(is1 + l, js0 + j, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 2:  x1 <--> east
        o = element(e) % face(2) % neighbor
        if (o > 0) then
          do k = 0, po
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, k, e) = v (ie1 + l,       j,       k, e) &
                                + vs(      l, js0 + j, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 3: -x2 <--> south
        o = element(e) % face(3) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, k, e) = v (      i, je0 + m,       k, e) &
                                + vs(is0 + i, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 4:  x2 <--> north
        o = element(e) % face(4) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, k, e) = v (      i, je1 + m,       k, e) &
                                + vs(is0 + i,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! face 5: -x3 <--> bottom
        o = element(e) % face(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            v(i, j, ke0 + n, e) = v (      i,       j, ke0 + n, e) &
                                + vs(is0 + i, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! face 6:  x3 <--> top
        o = element(e) % face(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do i = 0, po
            v(i, j, ke1 + n, e) = v (      i,       j, ke1 + n, e) &
                                + vs(is0 + i, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  1: (x1  , south, bottom)
        o = element(e) % edge( 1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, ke0 + n, e) = v (      i, je0 + m, ke0 + n, e) &
                                      + vs(is0 + i, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  2: (x1  , north, bottom)
        o = element(e) % edge( 2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, ke0 + n, e) = v (      i, je1 + m, ke0 + n, e) &
                                      + vs(is0 + i,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  3: (x1  , south, top   )
        o = element(e) % edge( 3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je0 + m, ke1 + n, e) = v (      i, je0 + m, ke1 + n, e) &
                                      + vs(is0 + i, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! edge  4: (x1  , north, top   )
        o = element(e) % edge( 4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do i = 0, po
            v(i, je1 + m, ke1 + n, e) = v (      i, je1 + m, ke1 + n, e) &
                                      + vs(is0 + i,       m,       n, o)
          end do
          end do
          end do
        end if

        ! edge  5: (west, x2   , bottom)
        o = element(e) % edge( 5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, ke0 + n, e) = v (ie0 + l,       j, ke0 + n, e) &
                                      + vs(is1 + l, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  6: (east, x2   , bottom)
        o = element(e) % edge( 6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, ke0 + n, e) = v (ie1 + l,       j, ke0 + n, e) &
                                      + vs(      l, js0 + j, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! edge  7: (west, x2   , top   )
        o = element(e) % edge( 7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie0 + l, j, ke1 + n, e) = v (ie0 + l,       j, ke1 + n, e) &
                                      + vs(is1 + l, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  8: (east, x2   , top   )
        o = element(e) % edge( 8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do j = 0, po
          do l = 1, no(1)
            v(ie1 + l, j, ke1 + n, e) = v (ie1 + l,       j, ke1 + n, e) &
                                      + vs(      l, js0 + j,       n, o)
          end do
          end do
          end do
        end if

        ! edge  9: (west, south, x3    )
        o = element(e) % edge( 9) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, k, e) = v (ie0 + l, je0 + m,       k, e) &
                                      + vs(is1 + l, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 10: (east, south, x3    )
        o = element(e) % edge(10) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, k, e) = v (ie1 + l, je0 + m,       k, e) &
                                      + vs(      l, js1 + m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 11: (west, north, x3    )
        o = element(e) % edge(11) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, k, e) = v (ie0 + l, je1 + m,       k, e) &
                                      + vs(is1 + l,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! edge 12: (east, north, x3    )
        o = element(e) % edge(12) % neighbor
        if (o > 0) then
          do k = 0, po
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, k, e) = v (ie1 + l, je1 + m,       k, e) &
                                      + vs(      l,       m, ks0 + k, o)
          end do
          end do
          end do
        end if

        ! vertex 1: (west, south, bottom)
        o = element(e) % vertex(1) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, ke0 + n, e) = v (ie0 + l, je0 + m, ke0 + n, e) &
                                            + vs(is1 + l, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 2: (east, south, bottom)
        o = element(e) % vertex(2) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, ke0 + n, e) = v (ie1 + l, je0 + m, ke0 + n, e) &
                                            + vs(      l, js1 + m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 3: (west, north, bottom)
        o = element(e) % vertex(3) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, ke0 + n, e) = v (ie0 + l, je1 + m, ke0 + n, e) &
                                            + vs(is1 + l,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 4: (east, north, bottom)
        o = element(e) % vertex(4) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, ke0 + n, e) = v (ie1 + l, je1 + m, ke0 + n, e) &
                                            + vs(      l,       m, ks1 + n, o)
          end do
          end do
          end do
        end if

        ! vertex 5: (west, south, top   )
        o = element(e) % vertex(5) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je0 + m, ke1 + n, e) = v (ie0 + l, je0 + m, ke1 + n, e) &
                                            + vs(is1 + l, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 6: (east, south, top   )
        o = element(e) % vertex(6) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je0 + m, ke1 + n, e) = v (ie1 + l, je0 + m, ke1 + n, e) &
                                            + vs(      l, js1 + m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 7: (west, north, top   )
        o = element(e) % vertex(7) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie0 + l, je1 + m, ke1 + n, e) = v (ie0 + l, je1 + m, ke1 + n, e) &
                                            + vs(is1 + l,       m,       n, o)
          end do
          end do
          end do
        end if

        ! vertex 8: (east, north, top   )
        o = element(e) % vertex(8) % neighbor
        if (o > 0) then
          do n = 1, no(3)
          do m = 1, no(2)
          do l = 1, no(1)
            v(ie1 + l, je1 + m, ke1 + n, e) = v (ie1 + l, je1 + m, ke1 + n, e) &
                                            + vs(      l,       m,       n, o)
          end do
          end do
          end do
        end if

      end do
    end associate

  ! clean-up ...................................................................

  call buf_vs % ToGhost_Finish()

end subroutine MergeFromSubdomains

!===============================================================================

end submodule MP_SchwarzMethod
