
  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real eXplicit

  subroutine Transfer_RX(this, mesh, v, v_ne, tag)

    ! arguments ................................................................

    class(ElementTransferBuffer_3D), asynchronous, intent(inout) :: this
    type(Mesh_3D), intent(in) :: mesh  !< mesh partition
    real(RK),               intent(in) :: v     !< mesh variable
    integer,                intent(in) :: v_ne  !< size of v in element dimension
    integer,                intent(in) :: tag   !< message tag

    dimension :: v(this%np(1), this%np(2), this%np(3), v_ne, this%nc)

    ! internal data ............................................................

    integer :: dest, source
    integer :: i, l, m

    ! sanity check .............................................................

    !$omp master
    if (v_ne < this%ne) then
      call Error('Transfer_RX','size(v,4) < this%ne')
    end if
    !$omp end master

    ! extract and send master data .............................................

    select type (vb => this % master % buf)
    type is (real(RK))

      call CopyToBuffer( nn   = size(this%master%node)  &
                       , nm   = size(v) / this%nc       &
                       , nc   = this%nc                 &
                       , node = this%master%node        &
                       , v    = v                       &
                       , vb   = vb                      )

      !$omp master
      do i = 1, size(mesh%link)
        dest = mesh % link(i) % part
        m = this % master % start(i)
        l = this % master % len(i)
        if (l > 0) then
          call XMPI_Isend( vb(m:m+l-1), dest, tag, mesh%comm &
                         , this%master%request(i)            )
        else
          this%master%request(i) = MPI_REQUEST_NULL
        end if
      end do
      !$omp end master

    end select

    ! receive master data ......................................................

    select type (vb => this % master % buf)
    type is (real(RK))

      !$omp master
      do i = 1, size(mesh%link)
        source = mesh % link(i) % part
        m = this % ghost % start(i)
        l = this % ghost % len(i)
        if (l > 0) then
          call XMPI_Irecv( vb(m:m+l-1), source, tag, mesh%comm &
                         , this%ghost%request(i)              )
        else
          this%ghost%request(i) = MPI_REQUEST_NULL
        end if
      end do
      !$omp end master

    end select

  contains

    subroutine CopyToBuffer(nn, nm, nc, node, v, vb)
      integer,  intent(in)  :: nn
      integer,  intent(in)  :: nm
      integer,  intent(in)  :: nc
      integer,  intent(in)  :: node(nn)
      real(RK), intent(in)  :: v(nm,nc)
      real(RK), intent(out) :: vb(nn,nc)

      integer :: i, j

      !$omp do collapse(2) private(i,j)
      do j = 1, nc
      do i = 1, nn
        vb(i,j) = v(node(i), j)
      end do
      end do

    end subroutine CopyToBuffer

  end subroutine Transfer_RX

  !-----------------------------------------------------------------------------
  !> Complete receive and merge buffer into ghost data -- real eXplicit
  !>
  !> Denoting the buffer with `vb`, the following operation will be executed:
  !>
  !>       v  =  alpha * v  +  beta * vb

  subroutine Merge_RX(this, v, alpha, beta)

    ! arguments ................................................................

    class(ElementTransferBuffer_3D), intent(inout) :: this
    real(RK),           intent(inout) :: v      !< mesh variable
    real(RK), optional, intent(in)    :: alpha  !< coeff of v  [1]
    real(RK), optional, intent(in)    :: beta   !< coeff of vb [1]

    dimension :: v(this%np(1), this%np(2), this%np(3), this%ne+this%ng, this%nc)

    ! internal data ............................................................

    real(RK) :: a, b
    integer  :: ng, nm

    ! wait for receive to complete .............................................

    !$omp master
    ng = size(this % ghost  % request)
    nm = size(this % master % request)
    call MPI_Waitall(ng, this % ghost  % request, MPI_STATUSES_IGNORE)
    call MPI_Waitall(nm, this % master % request, MPI_STATUSES_IGNORE)
    !$omp end master
    !$omp barrier

    ! merge buffer .............................................................

    select type (vb => this % master % buf)
    type is (real(RK))

      if (size(vb) == 0) return

      if (present(alpha)) then
        a = alpha
      else
        a = ZERO
      end if

      if (present(beta)) then
        b = beta
      else
        b = ONE
      end if

      call MergeBuffer( nn   = size(this%ghost%node)  &
                      , nm   = size(v) / this % nc    &
                      , nc   = this % nc              &
                      , node = this % ghost % node    &
                      , a    = a                      &
                      , b    = b                      &
                      , vb   = vb                     &
                      , v    = v                      )
    end select

  contains

    subroutine MergeBuffer(nn, nm, nc, node, a, b, vb, v)
      integer , intent(in)    :: nn
      integer , intent(in)    :: nm
      integer , intent(in)    :: nc
      integer , intent(in)    :: node(nn)
      real(RK), intent(in)    :: a, b
      real(RK), intent(in)    :: vb(nn,nc)
      real(RK), intent(inout) :: v(nm,nc)

      integer :: i, j

      if (a /= ZERO) then

        !$omp do collapse(2) private(i,j)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = a * v(node(i), j)  +  b * vb(i,j)
        end do
        end do

      else

        !$omp do collapse(2) private(i,j)
        do j = 1, nc
        do i = 1, nn
          v(node(i), j) = b * vb(i,j)
        end do
        end do

      end if

    end subroutine MergeBuffer

  end subroutine Merge_RX
