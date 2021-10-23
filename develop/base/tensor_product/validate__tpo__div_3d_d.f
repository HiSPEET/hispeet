!> summary:  Validation of the curvilinear tensor-product divergence operator
!> author:   Jerome Michel, Jörg Stiller, Erik Pfister
!> date:     2021/08/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO__Div_3d_D
  use Kind_Parameters
  use Standard_Operators__1D
  use TPO__Div__3D_D
  use TPO__Div__3D_D__Gen
  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Element_Transfer_Buffer__3D
  use Verify_Mesh__3D
  use Assembly__3D
  use Export_VTK_Volume_Data__3D


  implicit none

  !-----------------------------------------------------------------------------
  ! declarations

  ! input parameters ............................................................
  character(len=*), parameter :: &
         input_file     = 'validate__tpo__div_3d_d.prm'
  integer   :: conf     = 1       ! configuration (1 cylinder, 2 annular gap)
  real(RNP) :: r0       = 0.5     ! inner radius  (annular gap only)
  real(RNP) :: r1       = 1       ! outer radius
  real(RNP) :: h        = 2       ! height = axial extension
  integer   :: nr       = 2       ! num elements in radial    direction
  integer   :: np       = 4       ! num elements in azimuthal direction
  integer   :: nz       = 3       ! num elements in axial     direction
  integer   :: po       = 3       ! polynomial order of mesh elements
  integer   :: nt       = 1       ! number of test loops

  logical   :: periodic = .false. ! switch for axial periodicity


  namelist/input/ conf, nt, r0, r1, h, nr, np, nz, po


  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh_3D)         :: generic_mesh
  type(Mesh_3D)                :: mesh
  type(SpectralElementMesh_3D) :: se_mesh


  real(RNP), allocatable :: area(:)
  real(RNP) :: vol
  logical   :: passed
  integer   :: io
  integer   :: e, i, j, k, ne

  ! operators and variables ....................................................

  type(StandardOperators_1D) :: standard_op

  real(RNP), dimension(:,:,:,:,:), allocatable :: u
  real(RNP), dimension(:,:,:,:),   allocatable :: v, w


  real(RNP) :: time
  real(RNP) :: error_gen = 0, mflops_gen = -1, mlups_gen = -1
  real(RNP) :: error_opt = 0, mflops_opt = -1, mlups_opt = -1


  integer :: nflop, npop
  integer :: p, pm1


  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization
  call XMPI_Init()

  call MPI_Comm_rank(comm, rank)

  if (rank == 0) then ! derzeit noch keine Partitionierung => rank = 0

  ! read test parameters .......................................................


  open(newunit = io, file = input_file)
  read(io, nml = input)
  close(io)

  ! create and import generic mesh .............................................

  select case(conf)
  case(2)
    call generic_mesh % CreateAnnulus(r0, r1, h, nr, np, nz, po, periodic)
  case default
    call generic_mesh % CreateCylinder(r1, h, nr, nz, po, periodic)
  end select

  ! operator dimension
  !np = po + 1

  ! problem dimensions
  nflop = (po+1)**3 * (18*(po+1) + 17)
  npop  = (po+1)**3

  ! verification ...............................................................
  call mesh % ImportGenericMesh(generic_mesh, comm = comm)
  call VerifyMesh_3D(mesh, passed)
  write(*,'(/,A,G0,/)') 'VerifyMesh3d: passed = ', passed

  ! spectral element functionality .............................................

  se_mesh = SpectralElementMesh_3D(mesh, po)
  allocate(area(mesh % n_bound))
  call se_mesh % GetVolume(vol)
  call se_mesh % GetSurfaceAreas(area)

  write(*,'(A)') 'Spectral element mesh'
  write(*,'(2X,A,G0)') 'volume  = ', vol
  do i = 1, mesh % n_bound
    write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
  end do

  ! operators ..................................................................

  standard_op = StandardOperators_1D(po)

  ! workspace ..................................................................
  ne =  mesh%n_elem

  allocate( u(0:po,0:po,0:po,ne,3), &
            v(0:po,0:po,0:po,ne),   &
            w(0:po,0:po,0:po,ne)    )


  ! order of test function .....................................................

  p = min(po, 2)

  pm1 = max(p - 1, 0)

  !-----------------------------------------------------------------------------
  ! operand und exact result

  associate( Ms   => standard_op % w              &
           , Ds   => standard_op % D              &
           , x    => se_mesh % metrics % x        &
           , Ji   => se_mesh % metrics % Ji       &
           , Jd   => se_mesh % metrics % Jd)

    do e = 1, ne



      ! operand ................................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po

        u(i,j,k,e,1) =  x(i,j,k,e,1) ** p  *  x(i,j,k,e,2) ** pm1
        u(i,j,k,e,2) =  x(i,j,k,e,2) ** p  *  x(i,j,k,e,3) ** pm1
        u(i,j,k,e,3) =  x(i,j,k,e,3) ** p  *  x(i,j,k,e,1) ** pm1

      end do
      end do
      end do

      ! exact result ...........................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po

        w(i,j,k,e) =  p * ( x(i,j,k,e,1) ** pm1 * x(i,j,k,e,2) ** pm1) &
                   +  p * ( x(i,j,k,e,2) ** pm1 * x(i,j,k,e,3) ** pm1) &
                   +  p * ( x(i,j,k,e,3) ** pm1 * x(i,j,k,e,1) ** pm1)


      end do
      end do
      end do

    end do



  !-----------------------------------------------------------------------------
  ! test generic procedure


    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call TPO_Div_D_Gen(po+1, ne, Ds, Ji, u, v)
    !$acc wait

    call system_clock(count0, rate)

    do i = 1, nt
    call TPO_Div_D_Gen(po+1, ne, Ds, Ji, u, v)
    !$acc wait
    end do

    call system_clock(count)
    !$acc end data
    !$omp end parallel

  !end associate

  time = (count - count0) / real(rate, RNP) / nt

  error_gen  = maxval(abs(v - w))
  mflops_gen = 1E-6 / time * ne * nflop
  mlups_gen  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! test optimized procedure


    !$omp parallel
    !$acc data copyin(u) copyout(v)

    call TPO_Div(Ds, Ji, u, v)
    !$acc wait

    call system_clock(count0, rate)
    do i = 1, nt
      call TPO_Div(Ds, Ji, u, v)
      !$acc wait
    end do

    call system_clock(count)
    !$acc end data
    !$omp end parallel



  end associate


  time = (count - count0) / real(rate, RNP) / nt

  error_opt  = maxval(abs(v - w))
  mflops_opt = 1E-6 / time * ne * nflop
  mlups_opt  = 1E-6 / time * ne * npop


  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                        ',             &
                  '   ------------ generic ------------',  &
                  '   ----------- optimized  ----------'
  write(*,'(3A)') '#  np        ne        nt    ',         &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s'

  write(*,'(I5,2(2X,I8))',  advance='NO') po+1, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

end if

call MPI_Finalize()
!===============================================================================

end program Validate__TPO__Div_3d_D
