module Element_Distribution_Data__3D
  use Kind_Parameters
  use XMPI
  use Mesh_Element__3D
  use Mesh__3D

  implicit none
  private

  public :: ElementDistributionData_3D

  !-----------------------------------------------------------------------------
  !> Type for distributing repartitioned mesh element data

  type ElementDistributionData_3D

    ! MPI
    integer            :: proc = -1   !< MPI rank of source/destination
    type(MPI_Comm)     :: comm        !< MPI communicator
    type(MPI_Request)  :: req_dim(3)  !< MPI request for dimensions
    type(MPI_Request)  :: req_dat(8)  !< MPI request for data arrays

    ! dimensions
    integer :: n_elem     = 0         !< number of elements
    integer :: n_neighbor = 0         !< number of neighbors
    integer :: n_point    = 0         !< number of element points

    ! static element components
    type(MeshElement_3D), allocatable :: element(:)      !< static components

    ! neighbor components, data of element `e` starts at `start_neighbor(e)`
    integer     , allocatable :: start_neighbor(:)       !< neighbor start indices
    integer     , allocatable :: neighbor_id(:)          !< neighbor%id
    integer     , allocatable :: neighbor_part(:)        !< neighbor%part
    integer(IXS), allocatable :: neighbor_component(:)   !< neighbor%component
    integer(IXS), allocatable :: neighbor_orientation(:) !< neighbor%orientation

    ! element points, data of element `e` starts at `start_point(e)`
    integer     , allocatable :: start_point(:)          !< point start indices
    real(RNP)   , allocatable :: x_e(:,:)                !< geometry%x_e

  contains

    procedure :: Send_Dimensions
    procedure :: Send_Data
    procedure :: Send_Finish
    procedure :: Recv_Dimensions
    procedure :: Recv_Data
    procedure :: Recv_Finish
    procedure :: Assign_Data

  end type ElementDistributionData_3D

contains

  !-----------------------------------------------------------------------------
  !> Sending distribution data dimensions

  subroutine Send_Dimensions(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    !$omp master
    associate(dest => this%proc, comm => this%comm, req_dim => this%req_dim)
      if (dest >= 0) then
        call XMPI_Isend(this%n_elem , dest, 101, comm, req_dim(1))
        call XMPI_Isend(this%n_neighbor, dest, 102, comm, req_dim(2))
        call XMPI_Isend(this%n_point   , dest, 103, comm, req_dim(3))
      else
        req_dim = MPI_REQUEST_NULL
      end if
    end associate
    !$omp end master

  end subroutine Send_Dimensions

  !-----------------------------------------------------------------------------
  !> Sending distribution data

  subroutine Send_Data(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    type(MPI_Datatype) :: MPI_MeshElement_3D

    !$omp master
    associate(dest => this%proc, comm => this%comm, req_dat => this%req_dat)

      if (this % n_elem > 0) then

        call Get_MPI_MeshElement_3D(MPI_MeshElement_3D)

        call MPI_Isend( this % element(1) % id  &
                      , this % n_elem        &
                      , MPI_MeshElement_3D      &
                      , dest                    &
                      , 201                     &
                      , comm                    &
                      , req_dat(1)              )

        call XMPI_Isend(this%start_neighbor      , dest, 202, comm, req_dat(2))
        call XMPI_Isend(this%neighbor_id         , dest, 203, comm, req_dat(3))
        call XMPI_Isend(this%neighbor_part       , dest, 204, comm, req_dat(4))
        call XMPI_Isend(this%neighbor_component  , dest, 205, comm, req_dat(5))
        call XMPI_Isend(this%neighbor_orientation, dest, 206, comm, req_dat(6))
        call XMPI_Isend(this%start_point         , dest, 207, comm, req_dat(7))
        call XMPI_Isend(this%x_e                 , dest, 208, comm, req_dat(8))

      else
        req_dat = MPI_REQUEST_NULL
      end if

    end associate
    !$omp end master

  end subroutine Send_Data

  !-----------------------------------------------------------------------------
  !> Finish sending distribution data

  subroutine Send_Finish(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    !$omp master
    associate(req_dim => this%req_dim, req_dat => this%req_dat)
      if (this % proc >= 0) then
        call MPI_Waitall(size(req_dim), req_dim, MPI_STATUSES_IGNORE)
        call MPI_Waitall(size(req_dat), req_dat, MPI_STATUSES_IGNORE)
      end if
    end associate
    !$omp end master

  end subroutine Send_Finish

  !-----------------------------------------------------------------------------
  !> Receiving distribution data dimensions

  subroutine Recv_Dimensions(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    !$omp master
    associate(source => this%proc, comm => this%comm, req_dim => this%req_dim)
      if (source >= 0) then
        call XMPI_Irecv(this%n_elem , source, 101, comm, req_dim(1))
        call XMPI_Irecv(this%n_neighbor, source, 102, comm, req_dim(2))
        call XMPI_Irecv(this%n_point   , source, 103, comm, req_dim(3))
      else
        req_dim = MPI_REQUEST_NULL
      end if
    end associate
    !$omp end master

  end subroutine Recv_Dimensions

  !-----------------------------------------------------------------------------
  !> Receiving distribution data

  subroutine Recv_Data(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    type(MPI_Datatype) :: MPI_MeshElement_3D

    !$omp master
    associate( source => this%proc, req_dim => this%req_dim &
             , comm   => this%comm, req_dat => this%req_dat )

      call Get_MPI_MeshElement_3D(MPI_MeshElement_3D)

      call MPI_Waitall(size(req_dim), req_dim, MPI_STATUSES_IGNORE)

      if (this % n_elem > 0) then

        allocate( this % element              ( this % n_elem  ) )
        allocate( this % start_neighbor       ( this % n_elem  ) )
        allocate( this % neighbor_id          ( this % n_neighbor ) )
        allocate( this % neighbor_part        ( this % n_neighbor ) )
        allocate( this % neighbor_component   ( this % n_neighbor ) )
        allocate( this % neighbor_orientation ( this % n_neighbor ) )
        allocate( this % start_point          ( this % n_elem  ) )
        allocate( this % x_e                  ( this % n_point, 3 ) )

        call MPI_Irecv( this % element(1) % id  &
                      , this % n_elem        &
                      , MPI_MeshElement_3D      &
                      , source                  &
                      , 201                     &
                      , comm                    &
                      , req_dat(1)              )

        call XMPI_Irecv(this%start_neighbor      , source, 202, comm, req_dat(2))
        call XMPI_Irecv(this%neighbor_id         , source, 203, comm, req_dat(3))
        call XMPI_Irecv(this%neighbor_part       , source, 204, comm, req_dat(4))
        call XMPI_Irecv(this%neighbor_component  , source, 205, comm, req_dat(5))
        call XMPI_Irecv(this%neighbor_orientation, source, 206, comm, req_dat(6))
        call XMPI_Irecv(this%start_point         , source, 207, comm, req_dat(7))
        call XMPI_Irecv(this%x_e                 , source, 208, comm, req_dat(8))

      else
        req_dat = MPI_REQUEST_NULL
      end if

    end associate
    !$omp end master

  end subroutine Recv_Data

  !-----------------------------------------------------------------------------
  !> Receiving distribution data

  subroutine Recv_Finish(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    !$omp master
    associate(req_dat => this%req_dat)
      if (this % proc >= 0) then
        call MPI_Waitall(size(req_dat), req_dat, MPI_STATUSES_IGNORE)
      end if
    end associate
    !$omp end master

  end subroutine Recv_Finish

  !-----------------------------------------------------------------------------
  !> Assign element distribution data to new mesh partition

  subroutine Assign_Data(this, mesh)
    class(ElementDistributionData_3D), intent(in) :: this
    class(Mesh_3D), intent(inout) :: mesh

    integer :: cn, cp, nn, po
    integer :: e, i, j, k, l

    !$omp do
    do l = 1, this % n_elem

      e  = this % element(l) % id
      cn = this % start_neighbor(l)
      cp = this % start_point(l)

      ! static element components ..............................................

      mesh % element(e) = this % element(l)

      ! neighbors ..............................................................

      nn = sum( mesh % element(e) % face   % n_neighbor ) &
         + sum( mesh % element(e) % edge   % n_neighbor ) &
         + sum( mesh % element(e) % vertex % n_neighbor )

      allocate(mesh % element(e) % neighbor(nn))

      associate(neighbor => mesh % element(e) % neighbor)
        do i = 1, nn
          neighbor(i) % id          = this % neighbor_id          (cn)
          neighbor(i) % part        = this % neighbor_part        (cn)
          neighbor(i) % component   = this % neighbor_component   (cn)
          neighbor(i) % orientation = this % neighbor_orientation (cn)
          cn = cn + 1
        end do
      end associate

      ! element points .........................................................

      po = mesh % element(e) % geometry % po

      allocate(mesh % element(e) % geometry % x_e(0:po,0:po,0:po,3))

      associate(x_e => mesh % element(e) % geometry % x_e)
        do k = 0, po
        do j = 0, po
        do i = 0, po
          x_e(i,j,k,1) = this % x_e(cp, 1)
          x_e(i,j,k,2) = this % x_e(cp, 2)
          x_e(i,j,k,3) = this % x_e(cp, 3)
          cp = cp + 1
        end do
        end do
        end do
      end associate

    end do

  end subroutine Assign_Data

  !=============================================================================

end module Element_Distribution_Data__3D
