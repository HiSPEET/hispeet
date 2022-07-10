module Element_Distribution_Data__3D
  use Kind_Parameters
  use XMPI
  use Mesh_Element__3D

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
    type(MPI_Request)  :: req_dat(6)  !< MPI request for data arrays

    ! dimensions
    integer :: n_element  = 0         !< number of elements
    integer :: n_neighbor = 0         !< number of neighbors
    integer :: n_point    = 0         !< number of element points

    ! static element components
    type(MeshElement_3D), allocatable :: element(:)      !< static components

    ! concatenated element neighbor components
    integer     , allocatable :: neighbor_id(:)          !< neighbor%id
    integer     , allocatable :: neighbor_part(:)        !< neighbor%part
    integer(IXS), allocatable :: neighbor_component(:)   !< neighbor%component
    integer(IXS), allocatable :: neighbor_orientation(:) !< neighbor%orientation

    ! concatenated element points
    real(RNP), allocatable :: x_e(:,:)                   !< geometry%x_e

  contains

    procedure :: Send_Dimensions
    procedure :: Send_Data
    procedure :: Send_Finish
    procedure :: Recv_Dimensions
    procedure :: Recv_Data
    procedure :: Recv_Finish

  end type ElementDistributionData_3D

contains

  !-----------------------------------------------------------------------------
  !> Sending distribution data dimensions

  subroutine Send_Dimensions(this)
    class(ElementDistributionData_3D), asynchronous, intent(inout) :: this

    !$omp master
    associate(dest => this%proc, comm => this%comm, req_dim => this%req_dim)
      if (dest >= 0) then
        call XMPI_Isend(this%n_element , dest, 101, comm, req_dim(1))
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

      if (this % n_element > 0) then

        call Get_MPI_MeshElement_3D(MPI_MeshElement_3D)

        call MPI_Isend( this % element(1) % id  &
                      , this % n_element        &
                      , MPI_MeshElement_3D      &
                      , dest                    &
                      , 201                     &
                      , comm                    &
                      , req_dat(1)              )

        call XMPI_Isend(this%neighbor_id         , dest, 202, comm, req_dat(2))
        call XMPI_Isend(this%neighbor_part       , dest, 203, comm, req_dat(3))
        call XMPI_Isend(this%neighbor_component  , dest, 204, comm, req_dat(4))
        call XMPI_Isend(this%neighbor_orientation, dest, 205, comm, req_dat(5))
        call XMPI_Isend(this%x_e                 , dest, 206, comm, req_dat(6))

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
        call XMPI_Irecv(this%n_element , source, 101, comm, req_dim(1))
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

      if (this % n_element > 0) then

        allocate( this % element              ( this % n_element  ) )
        allocate( this % neighbor_id          ( this % n_neighbor ) )
        allocate( this % neighbor_part        ( this % n_neighbor ) )
        allocate( this % neighbor_component   ( this % n_neighbor ) )
        allocate( this % neighbor_orientation ( this % n_neighbor ) )
        allocate( this % x_e                  ( this % n_point, 3 ) )

        call MPI_Irecv( this % element(1) % id  &
                      , this % n_element        &
                      , MPI_MeshElement_3D      &
                      , source                  &
                      , 201                     &
                      , comm                    &
                      , req_dat(1)              )

        call XMPI_Irecv(this%neighbor_id         , source, 202, comm, req_dat(2))
        call XMPI_Irecv(this%neighbor_part       , source, 203, comm, req_dat(3))
        call XMPI_Irecv(this%neighbor_component  , source, 204, comm, req_dat(4))
        call XMPI_Irecv(this%neighbor_orientation, source, 205, comm, req_dat(5))
        call XMPI_Irecv(this%x_e                 , source, 206, comm, req_dat(6))

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

  !=============================================================================

end module Element_Distribution_Data__3D
