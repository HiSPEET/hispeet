!> summary:  Validation of the curvilinear elliptic TPO for constant diffusivity
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/08/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO_Elliptic_DLCI
! validate
  use Kind_Parameters
  use Standard_Element_Operators__1D
  use TPO__Elliptic__3D_DLCI
  use TPO__Elliptic__3D_DLCI__Gen
! mesh
  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Element_Transfer_Buffer__3D
  use Verify_Mesh__3D
  use Assembly__3D
  use Export_VTK_Volume_Data__3D

  !use TPO__Elliptic__3D_DLCI__XSMM

  implicit none

  ! input parameters ...........................................................

  character(len=*), parameter :: &
         input_file     = 'validate__tpo__elliptic_3d_dlci.prm'
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
  logical   :: exact    = .false. ! compare with exact or approximate gradient
  real(RNP) :: lambda   = 1       ! Helmholtz parameter
  real(RNP) :: nu       = 1       ! diffusivity

  namelist/input/ conf, nt, r0, r1, h, nr, np, nz, po, exact, lambda, nu

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh_3D)         :: generic_mesh
  type(Mesh_3D)                :: mesh
  type(SpectralElementMesh_3D) :: se_mesh

  real(RNP), allocatable :: area(:)
  real(RNP) :: vol
  logical   :: passed
  integer   :: io
  integer   :: e, i, j, k, l, ne

  ! operators and variables ....................................................

  type(StandardElementOperators_1D) :: standard_op

  real(RNP), dimension(:,:,:,:), allocatable :: u, v, w
  real(RNP), dimension(:,:,:,:), allocatable :: Jd_Ji_grad_u

  real(RNP) :: dx_u, dy_u, dz_u, d1_u, d2_u, d3_u

  real(RNP) :: time  = 0
  real(RNP) :: error_gen = 0, mflops_gen = -1, mlups_gen = -1
  real(RNP) :: error_opt = 0, mflops_opt = -1, mlups_opt = -1

  integer :: npop, nflop
  integer :: p, pm1, pm2

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

  ! test parameters ............................................................

  ! problem dimensions
  npop  = (po+1)**3                 ! number of operands per element
  nflop = npop * (12*(po+1) + 29)   ! number of FLOPs per element

  ! verification ...............................................................

  call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  call VerifyMesh_3D(mesh, passed)
  write(*,'(/,A,G0,/)') 'VerifyMesh3d: passed = ', passed

  ! spectral element functionality .............................................

  se_mesh = SpectralElementMesh_3D(mesh, po)
  allocate(area(mesh % n_bound))
  call se_mesh % Get_Volume(vol)
  call se_mesh % Get_SurfaceAreas(area)

  write(*,'(A)') 'Spectral element mesh'
  write(*,'(2X,A,G0)') 'volume  = ', vol
  do i = 1, mesh % n_bound
    write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
  end do

  ! operators ..................................................................

  standard_op = StandardElementOperators_1D(po)

  ! workspace ..................................................................

  ne =  mesh%n_elem ! number of elements

  allocate( u(0:po,0:po,0:po,ne), &
            v(0:po,0:po,0:po,ne), &
            w(0:po,0:po,0:po,ne)  )

  allocate( Jd_Ji_grad_u(0:po,0:po,0:po,3) )

  ! order of test function .....................................................

  p = min(po, 2)

  pm1 = max(p - 1, 0)
  pm2 = max(p - 2, 0)

  !-----------------------------------------------------------------------------
  ! operand und exact result

  associate( Ms => standard_op % w        &
           , Ds => standard_op % D        &
           , x  => se_mesh % metrics % x  &
           , Ji => se_mesh % metrics % Ji &
           , Jd => se_mesh % metrics % Jd &
           , G  => se_mesh % metrics % G  )

    do e = 1, ne

      do k = 0, po
      do j = 0, po
      do i = 0, po

        ! operand ...............................................................

        u(i,j,k,e) =  x(i,j,k,e,1) ** p  *  x(i,j,k,e,2) ** pm1  &
                   +  x(i,j,k,e,2) ** p  *  x(i,j,k,e,3) ** pm1  &
                   +  x(i,j,k,e,3) ** p  *  x(i,j,k,e,1) ** pm1

      end do
      end do
      end do

      do k = 0, po
      do j = 0, po
      do i = 0, po

        ! gradient .............................................................

        if (exact) then  ! exact

          dx_u  =  p    *  x(i,j,k,e,1) ** pm1  *  x(i,j,k,e,2) ** pm1  &
                +  pm1  *  x(i,j,k,e,3) ** p    *  x(i,j,k,e,1) ** pm2

          dy_u  =  p    *  x(i,j,k,e,2) ** pm1  *  x(i,j,k,e,3) ** pm1  &
                +  pm1  *  x(i,j,k,e,1) ** p    *  x(i,j,k,e,2) ** pm2

          dz_u  =  p    *  x(i,j,k,e,3) ** pm1  *  x(i,j,k,e,1) ** pm1  &
                +  pm1  *  x(i,j,k,e,2) ** p    *  x(i,j,k,e,3) ** pm2

        else  ! approximate

          d1_u = sum(Ds(i,:) * u(:,j,k,e))
          d2_u = sum(Ds(j,:) * u(i,:,k,e))
          d3_u = sum(Ds(k,:) * u(i,j,:,e))

          dx_u = d1_u * Ji(i,j,k,e,1,1) &
               + d2_u * Ji(i,j,k,e,2,1) &
               + d3_u * Ji(i,j,k,e,3,1)

          dy_u = d1_u * Ji(i,j,k,e,1,2) &
               + d2_u * Ji(i,j,k,e,2,2) &
               + d3_u * Ji(i,j,k,e,3,2)

          dz_u = d1_u * Ji(i,j,k,e,1,3) &
               + d2_u * Ji(i,j,k,e,2,3) &
               + d3_u * Ji(i,j,k,e,3,3)

        end if

        Jd_Ji_grad_u(i,j,k,1)  =  Jd(i,j,k,e) * (  Ji(i,j,k,e,1,1) * dx_u  &
                                                +  Ji(i,j,k,e,1,2) * dy_u  &
                                                +  Ji(i,j,k,e,1,3) * dz_u  )

        Jd_Ji_grad_u(i,j,k,2)  =  Jd(i,j,k,e) * (  Ji(i,j,k,e,2,1) * dx_u  &
                                                +  Ji(i,j,k,e,2,2) * dy_u  &
                                                +  Ji(i,j,k,e,2,3) * dz_u  )

        Jd_Ji_grad_u(i,j,k,3)  =  Jd(i,j,k,e) * (  Ji(i,j,k,e,3,1) * dx_u  &
                                                +  Ji(i,j,k,e,3,2) * dy_u  &
                                                +  Ji(i,j,k,e,3,3) * dz_u  )

      end do
      end do
      end do

      ! reference ..............................................................

      do k = 0, po
      do j = 0, po
      do i = 0, po

        w(i,j,k,e) = lambda * Ms(i) * Ms(j) * Ms(k) * Jd(i,j,k,e) * u(i,j,k,e)

        do l = 0, po
          w(i,j,k,e) = w(i,j,k,e)                                           &
            + nu * Ms(l) * Ms(j) * Ms(k) * Ds(l,i) * Jd_Ji_grad_u(l,j,k,1)  &
            + nu * Ms(i) * Ms(l) * Ms(k) * Ds(l,j) * Jd_Ji_grad_u(i,l,k,2)  &
            + nu * Ms(i) * Ms(j) * Ms(l) * Ds(l,k) * Jd_Ji_grad_u(i,j,l,3)
        end do

      end do
      end do
      end do

    end do

    !---------------------------------------------------------------------------
    ! test generic procedure

    !$omp parallel

    call TPO_Elliptic_DLCI_Gen(Ms, Ds, Jd, G, lambda, nu, u, v)

    call system_clock(count0, rate)

    do i = 1, nt
      call TPO_Elliptic_DLCI_Gen(Ms, Ds, Jd, G, lambda, nu, u, v)
    end do

    call system_clock(count)
    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / max(nt, 1)

    error_gen  = maxval(abs(v - w))
    mflops_gen = 1E-6 / time * ne * nflop
    mlups_gen  = 1E-6 / time * ne * npop

    !---------------------------------------------------------------------------
    ! test optimized procedure

    !$omp parallel

    call TPO_Elliptic_DLCI(Ms, Ds, Jd, G, lambda, nu, u, v)

    call system_clock(count0, rate)

    do i = 1, nt
      call TPO_Elliptic_DLCI(Ms, Ds, Jd, G, lambda, nu, u, v)
    end do

    call system_clock(count)

    !$omp end parallel

    time = (count - count0) / real(rate, RNP) / nt

    error_opt  = maxval(abs(v - w))
    mflops_opt = 1E-6 / time * ne * nflop
    mlups_opt  = 1E-6 / time * ne * npop

  end associate

  !-----------------------------------------------------------------------------
  ! print results

  write(*,*)
  write(*,'(3A)') '#                        ',             &
                  '   ------ generic operator --------',   &
                  '   ------ optimized operator ------'
  write(*,'(3A)') '#  np        ne        nt    ',         &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s'
  write(*,'(I5,2(2X,I8))',  advance='NO') po+1, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_gen, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

  end if

  call MPI_Finalize()

  !=============================================================================

end program Validate__TPO_Elliptic_DLCI
