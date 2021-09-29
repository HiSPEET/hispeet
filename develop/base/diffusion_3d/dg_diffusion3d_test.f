program DG_Diffusion3D_Test
  use Kind_Parameters
  use Constants
  use XMPI
  use Spectral_Element_Mesh__3D
  use Export_VTK_Volume_Data__3D
  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cylinder
  use Create_Annulus
  implicit none

  ! control parameters .........................................................

  character(len=*), parameter :: input = 'dg_diffusion3d_test.prm'
  ! input file

  integer :: config = 1
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cylindrical domain                                                (u+d)
  !   4  annular domain                                                    (u+d)
  ! configuration > 1 currently available only with one MPI process

  integer :: pg = 1               ! polynomial degree of geometry description
  integer :: po = 7               ! polynomial degree of spectral elements
  logical :: export_vtk = .false. ! generate VTK files

  namelist/control/ config, pg, po, export_vtk

  ! MPI variables ..............................................................

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: mpi_rank
  integer :: mpi_size

  ! problem variables ..........................................................

  type(SpectralElementMesh_3D) :: sem

  ! auxiliary variables ........................................................

  integer :: io

!?!  real(RNP), allocatable, target :: var(:,:,:,:,:)
!?!  real(RNP), allocatable         :: x(:,:,:,:,:)
!?!  real(RNP), allocatable         :: v(:,:,:,:)     ! test variable
!?!  real(RNP), pointer             :: r(:,:,:,:)     ! reference variable
!?!  real(RNP), pointer             :: e(:,:,:,:)     ! error
!?!
!?!  type(ElementTransferBuffer_3D), asynchronous, allocatable :: v_buf
!?!
!?!  real(RNP), allocatable :: area(:)
!?!  real(RNP) :: vol
!?!  real(RNP) :: kappa(3), y(3), err
!?!  logical   :: passed
!?!  integer   :: i, j, k, l
!?!  integer   :: i_err, j_err, k_err, l_err
!?!

  ! initialization .............................................................

  call Init_MPI_Binding()
  call MPI_Comm_rank(comm, mpi_rank)
  call MPI_Comm_size(comm, mpi_size)

  ! read control parameters
  if (mpi_rank == 0) then
    open(newunit = io, file = input)
    read(io, nml = control)
    close(io)

    if (mpi_size > 1) config = 1  ! so far
  end if

  call XMPI_Bcast(config    , 0, comm)
  call XMPI_Bcast(pg        , 0, comm)
  call XMPI_Bcast(po        , 0, comm)
  call XMPI_Bcast(export_vtk, 0, comm)

  ! create spectral element mesh
  select case(config)
  case(2)
    call CreateCuboidDiamonds(comm, input, pg, po, sem)
  case(3)
    call CreateCylinder(comm, input, pg, po, sem)
  case(4)
    call CreateAnnulus(comm, input, pg, po, sem)
  case default
    call CreateCuboidCartesian(comm, input, pg, po, sem)
  end select

  ! finalization ...............................................................

  call MPI_Finalize()


!?!
!?!    se_mesh = SpectralElementMesh_3D(mesh, po)
!?!
!?!    allocate(area(mesh % n_bound))
!?!    call se_mesh % GetVolume(vol)
!?!    call se_mesh % GetSurfaceAreas(area)
!?!
!?!    write(*,'(A)') 'Spectral element mesh'
!?!    write(*,'(2X,A,G0)') 'volume  = ', vol
!?!    do i = 1, mesh % n_bound
!?!      write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
!?!    end do
!?!
!?!    ! set up data ..............................................................
!?!
!?!    if (test_avg .or. export_vtk) then
!?!      call mesh % GetPoints(po, 'L', x)
!?!      allocate(v  (0:po, 0:po, 0:po, mesh%n_elem + mesh%n_ghost))
!?!      allocate(var(0:po, 0:po, 0:po, mesh%n_elem, 2), source = ZERO)
!?!      r(0:,0:,0:,1:) => var(:,:,:,:,1)
!?!      e(0:,0:,0:,1:) => var(:,:,:,:,2)
!?!    end if
!?!
!?!    ! averaging test ...........................................................
!?!
!?!    if (test_avg) then
!?!
!?!      kappa = 2 * PI / h
!?!
!?!      do l = 1, mesh%n_elem
!?!        do k = 0, po
!?!        do j = 0, po
!?!        do i = 0, po
!?!          y(1) = x(i,j,k,l,1)
!?!          y(2) = x(i,j,k,l,2)
!?!          y(3) = x(i,j,k,l,3)
!?!          r(i,j,k,l) = cos(sum(kappa * y))
!?!          v(i,j,k,l) = r(i,j,k,l)
!?!        end do
!?!        end do
!?!        end do
!?!      end do
!?!
!?!      v_buf = ElementTransferBuffer_3D(mesh, v)
!?!      call Assembly_3D(mesh, v, v_buf, avg=.true.)
!?!
!?!      err = 0
!?!      do l = 1, mesh%n_elem
!?!        do k = 0, po
!?!        do j = 0, po
!?!        do i = 0, po
!?!          e(i,j,k,l) = v(i,j,k,l) - r(i,j,k,l)
!?!          if (abs(e(i,j,k,l)) > err) then
!?!            err   = abs(e(i,j,k,l))
!?!            i_err = i
!?!            j_err = j
!?!            k_err = k
!?!            l_err = l
!?!          end if
!?!        end do
!?!        end do
!?!        end do
!?!      end do
!?!      write(*,'(A,ES10.3,A,4I5)') 'Average over element boundaries: err = ', &
!?!                                  err, ' at ', i_err, j_err, k_err, l_err
!?!    end if
!?!
!?!    ! export mesh and data .....................................................
!?!
!?!    if (export_vtk) then
!?!      call ExportVTK_VolumeData( x, var                  &
!?!                               , sname  = ['r','e']      &
!?!                               , file   = 'element_mesh' &
!?!                               , part   = mesh % part    &
!?!                               , n_part = mesh % n_part  )
!?!
!?!      call mesh % GetCuboids(x)
!?!      call ExportVTK_VolumeData( x                       &
!?!                               , file   = 'cuboid_mesh'  &
!?!                               , part   = mesh % part    &
!?!                               , n_part = mesh % n_part  )
!?!
!?!    end if
!?!
!?!  end if

contains

  !-----------------------------------------------------------------------------
  !>
  !=============================================================================


end program DG_Diffusion3D_Test
