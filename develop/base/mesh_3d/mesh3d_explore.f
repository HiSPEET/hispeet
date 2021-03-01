program Mesh3d_Explore
  use Kind_Parameters, only: RNP
  use XMPI
  use Mesh_3d__Element
  use Mesh_3d__Partition
  use Mesh_3d__Generate_Regular_Mesh
  use Verify_Mesh_3d
  implicit none

  character(len=*), parameter :: input_file = 'mesh3d_explore.prm'
  real(RNP) :: xo(3) = -1                ! corner closest to -infinity
  real(RNP) :: lx(3) =  2                ! domain extensions
  integer   :: np(3) = [2,1,1]           ! number of partitions per direction
  integer   :: ep(3) = [1,1,1]           ! elements per partition and direction
  integer   :: pg    =  1                ! degree of geometry description
  logical   :: periodic(3) = .false.     ! periodic directions set true
  namelist/input/ xo, lx, np, ep, pg, periodic

  type(MPI_Comm) :: comm                 ! MPI communicator
  integer        :: rank                 ! local MPI rank
  integer        :: n_proc               ! number of MPI processes

  type(Mesh3d_Partition) :: mesh         ! mesh partition
  real(RNP), allocatable :: x(:,:,:,:,:) ! mesh points
  real(RNP)              :: dx(3)        ! element spacing in directions 1:3

  integer :: io, l, part
  logical :: passed, all_passed


  call XMPI_Init()
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  if (rank == 0) then
    open(newunit = io, file = input_file)
    read(io, nml = input)
    close(io)
  end if

  call XMPI_Bcast(xo      , 0, comm)
  call XMPI_Bcast(lx      , 0, comm)
  call XMPI_Bcast(np      , 0, comm)
  call XMPI_Bcast(ep      , 0, comm)
  call XMPI_Bcast(pg      , 0, comm)
  call XMPI_Bcast(periodic, 0, comm)

  dx = lx / (np * ep)

  call GenerateRegularMesh(mesh, np, ep, xo, dx, periodic, comm, pg)
  call VerifyMesh3d(mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)
  if (rank == 0) then
    write(*,'(/,A,G0)') 'VerifyMesh3d: passed = ', all_passed
  end if

  do
    if (rank == 0) then
      write(*,'(/,A,I0,A)',advance='NO') 'partition: ', mesh%n_part,' > part = '
      read(*,*) part
    end if
    call XMPI_Bcast(part, 0, comm)

    if (part < 0) exit

    if (part == rank) then
      call ShowMeshProperties(mesh)
    end if
    call MPI_Barrier(comm)

    do
      if (rank == 0) then
        write(*,'(/,A,I0,A)',advance='NO') 'element: ', mesh%n_elem,' >= l = '
        read(*,*) l
      end if
      call XMPI_Bcast(l, 0, comm)
      if (l < 1) exit
      if (part == rank) then
        if (l <= mesh % n_elem) then
          call ShowMeshElement(mesh % element(l))
        else
          call ShowMeshElement(mesh % ghost(l - mesh%n_elem))
        end if
      end if
      call MPI_Barrier(comm)
    end do

  end do

  call MPI_Finalize()

contains

  subroutine ShowMeshProperties(mesh)
    class(Mesh3d_Partition), intent(in) :: mesh

    integer :: min_vert, max_vert
    integer :: min_edge, max_edge
    integer :: min_face, max_face
    integer :: min_elem, max_elem
    integer :: min_n_nb, max_n_nb
    integer :: i

    write(*,*)
    write(*,'(A,3(1X,G0))') 'n_vert     =', mesh % n_vert
    write(*,'(A,3(1X,G0))') 'n_edge     =', mesh % n_edge
    write(*,'(A,3(1X,G0))') 'n_face     =', mesh % n_face
    write(*,'(A,3(1X,G0))') 'n_elem     =', mesh % n_elem
    write(*,'(A,3(1X,G0))') 'n_ghost    =', mesh % n_ghost
    write(*,'(A,3(1X,G0))') 'p_geom     =', mesh % p_geom
    write(*,'(A,3(1X,G0))') 'n_bound    =', mesh % n_bound
    write(*,'(A,3(1X,G0))') 'n_part     =', mesh % n_part
    write(*,'(A,3(1X,G0))') 'part       =', mesh % part
    write(*,'(A,3(1X,G0))') 'structured =', mesh % structured
    write(*,'(A,3(1X,G0))') 'regular    =', mesh % regular
    write(*,'(A,3(1X,G0))') 'n_elem_1   =', mesh % n_elem_1
    write(*,'(A,3(1X,G0))') 'n_elem_2   =', mesh % n_elem_2
    write(*,'(A,3(1X,G0))') 'n_elem_3   =', mesh % n_elem_3
    write(*,'(A,3(1X,G0))') 'n_face_1   =', mesh % n_face_1
    write(*,'(A,3(1X,G0))') 'n_face_2   =', mesh % n_face_2
    write(*,'(A,3(1X,G0))') 'n_face_3   =', mesh % n_face_3
    write(*,'(A,3(1X,G0))') 'dx         =', mesh % dx

    min_vert = huge(1);  max_vert = -huge(1)
    min_face = huge(1);  max_face = -huge(1)
    min_edge = huge(1);  max_edge = -huge(1)
    min_elem = huge(1);  max_elem = -huge(1)
    min_n_nb = huge(1);  max_n_nb = -huge(1)

    do i = 1, mesh % n_elem
      associate(element => mesh % element(i))
        min_vert = min( min_vert, minval(element % vertex   % id) )
        max_vert = max( max_vert, maxval(element % vertex   % id) )
        min_edge = min( min_edge, minval(element % edge     % id) )
        max_edge = max( max_edge, maxval(element % edge     % id) )
        min_face = min( min_face, minval(element % face     % id) )
        max_face = max( max_face, maxval(element % face     % id) )
        min_elem = min( min_elem, minval(element % neighbor % id) )
        max_elem = max( max_elem, maxval(element % neighbor % id) )
        min_n_nb = min( min_n_nb, size(element % neighbor) )
        max_n_nb = max( max_n_nb, size(element % neighbor) )
      end associate
    end do
    write(*,*)
    write(*,'(A,3(1X,G0))') 'min/max_vert =', min_vert, max_vert
    write(*,'(A,3(1X,G0))') 'min/max_edge =', min_edge, max_edge
    write(*,'(A,3(1X,G0))') 'min/max_face =', min_face, max_face
    write(*,'(A,3(1X,G0))') 'min/max_elem =', min_elem, max_elem
    write(*,'(A,3(1X,G0))') 'min/max_n_nb =', min_n_nb, max_n_nb

  end subroutine ShowMeshProperties

  subroutine ShowMeshElement(element)
    class(Mesh3d_Element), intent(in) :: element

    integer :: i, j, j1, j2, k

    write(*,*)
    write(*,'(A,99(1X,I5))') 'global_id     =', element%global_id
    write(*,'(A,99(1X,I5))') 'local_id      =', element%local_id
    write(*,'(A,99(1X,I5))') 'vertex % id   =', element%vertex%id
    write(*,'(A,99(1X,I5))') 'vertex % rank =', element%vertex%rank
    write(*,'(A,99(1X,I5))') 'vertex % val  =', element%vertex%val
    write(*,'(A,99(1X,I5))') 'edge   % id   =', element%edge%id
    write(*,'(A,99(1X,I5))') 'edge   % rank =', element%edge%rank
    write(*,'(A,99(1X,I5))') 'edge   % val  =', element%edge%val
    write(*,'(A,99(1X,I5))') 'face   % id   =', element%face%id
    write(*,'(A,99(1X,I5))') 'face   % rank =', element%face%rank
    write(*,'(A,99(1X,I5))') 'face   % val  =', element%face%val
    write(*,*)

    if (.not.allocated(element%neighbor)) return

    write(*,'(A,99(1X,I5))') 'n_neighbor  =', size(element%neighbor)
    write(*,*)
    k = 0
    do i = 1, 8
      j1 = element % vertex(i) % i_neighbor
      j2 = element % vertex(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,3I5)') k, '  vertex', i, '  id,part,cc =', &
                                   element % neighbor(j)
      end do
    end do
    do i = 1, 12
      j1 = element % edge(i) % i_neighbor
      j2 = element % edge(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,3I5)') k, '  edge  ', i, '  id,part,cc =', &
                                   element % neighbor(j)
      end do
    end do
    do i = 1, 6
      j1 = element % face(i) % i_neighbor
      j2 = element % face(i) % n_neighbor + j1 - 1
      do j = j1, j2
        k  = k + 1
        write(*,'(I4,A,I3,A,3I5)') k, '  face  ', i, '  id,part,cc =', &
                                   element % neighbor(j)
      end do
    end do

  end subroutine ShowMeshElement

end program Mesh3d_Explore
