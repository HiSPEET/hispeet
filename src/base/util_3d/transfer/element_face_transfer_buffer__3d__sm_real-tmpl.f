
  !-----------------------------------------------------------------------------
  !> Send master data to ghosts and receive own ghost data -- real eXplicit

  subroutine Transfer_RX(this, mesh, v, v_ne, tag)

    ! arguments ................................................................

    class(ElementFaceTransferBuffer_3D), asynchronous, intent(inout) :: this
    class(Mesh_3D), intent(in) :: mesh  !< mesh partition
    real(RK),       intent(in) :: v     !< element-face variable
    integer,        intent(in) :: v_ne  !< size of v in element dimension
    integer,        intent(in) :: tag   !< message tag

    dimension :: v(this%np, this%np, 6, v_ne, this%nc)

    ! internal data ............................................................

    integer :: c, e, f, f1, f2, i, l

    ! receive remote data ......................................................

    select type (vb => this % recv % buf)
    type is (real(RK))
      !$omp master
      f1 = 1
      do l = 1, size(mesh%link)
        if (mesh % link(l) % n_face > 0) then
          f2 = f1 + mesh % link(l) % n_face - 1
          call XMPI_Irecv( vb(:,:,:,f1:f2), mesh%link(l)%part, tag &
                         , mesh%comm , this%recv%request(l)        )
        else
          this%recv%request(l) = MPI_REQUEST_NULL
        end if
      end do
      !$omp end master
    end select

    ! extract and send local data ..............................................

    select type (vb => this % send % buf)
    type is (real(RK))

      if (this%nc == 1) then

        !$omp do
        do i = 1, this % nf
          vb(:,:,1,i) = v(:,:, this%send%face(i), this%send%element(i), 1)
        end do

      else

        !$omp do
        do i = 1, this % nf
          e = this % send % element(i)
          f = this % send % face(i)
          do c = 1, this % nc
            vb(:,:,c,i) = v(:,:,f,e,c)
          end do
        end do

      end if

      !$omp master
      f1 = 1
      do l = 1, size(mesh%link)
        if (mesh % link(l) % n_face > 0) then
          f2 = f1 + mesh % link(l) % n_face - 1
          call XMPI_Isend( vb(:,:,:,f1:f2), mesh%link(l)%part, tag &
                         , mesh%comm , this%send%request(l)        )
          f1 = f2 + 1
        else
          this%send%request(l) = MPI_REQUEST_NULL
        end if
      end do
      !$omp end master

    end select

  end subroutine Transfer_RX

  !-----------------------------------------------------------------------------
  !> Complete receive and merge remote data -- real eXplicit

  subroutine Merge_RX(this, v)

    ! arguments ................................................................

    class(ElementFaceTransferBuffer_3D), intent(inout) :: this
    real(RK), intent(inout) :: v !< element-face variable, including ghosts

    dimension :: v(this%np, this%np, 6, this%ne+this%ng, this%nc)

    ! internal data ............................................................

    integer :: c, e, f, i, n_req

    ! wait for receive to complete .............................................

    !$omp master
    n_req = size(this % recv % request)
    call MPI_Waitall(n_req, this % recv % request, MPI_STATUSES_IGNORE)
    !$omp end master
    !$omp barrier

    ! merge buffer .............................................................

    select type (vb => this % recv % buf)
    type is (real(RK))

      if (size(vb) == 0) return

      if (this%nc == 1) then
        !$omp do
        do i = 1, this % nf
          v(:,:, this%recv%face(i), this%recv%element(i), 1) = vb(:,:,1,i)
        end do
      else
        !$omp do
        do i = 1, this % nf
          e = this % recv % element(i)
          f = this % recv % face(i)
          do c = 1, this % nc
            v(:,:,f,e,c) = vb(:,:,c,i)
          end do
        end do
      end if

    end select

    !$omp master
    call MPI_Waitall(n_req, this % send % request, MPI_STATUSES_IGNORE)
    !$omp end master

  end subroutine Merge_RX
