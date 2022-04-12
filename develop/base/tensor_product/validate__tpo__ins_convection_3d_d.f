!> summary:  Validation of the curvilinear TPO for 3D convection
!> author:   Jerome Michel, Joerg Stiller
!> date:     2021/01/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program Validate__TPO_INS_Convection_D
  use Kind_Parameters
  use TPO__INS_Convection__3D_D
  use TPO__INS_Convection__3D_D__Gen
  !use TPO__INS_Convection__3D_D__XSMM_RDP

  use Execution_Control
  use Array_Assignments
  use Array_Reductions

  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Mesh__3D
  use Spectral_Element_Mesh__3D
  use Convert_Traces__3D
  use Spectral_Element_Variable__3D

  use Verify_Mesh__3D
  use Export_VTK_Volume_Data__3D

  use INS__Operator__3D

  implicit none

  ! input parameters ...........................................................

  character(len=*), parameter :: &
         input_file     = 'validate__tpo__ins_convection_3d_d.prm'
  integer   :: conf     = 1         ! configuration (1 cylinder, 2 annular gap)
  real(RNP) :: r0       = 0.5       ! inner radius  (annular gap only)
  real(RNP) :: r1       = 1         ! outer radius
  real(RNP) :: h        = 2         ! height = axial extension
  integer   :: nr       = 2         ! num elements in radial    direction
  integer   :: np       = 4         ! num elements in azimuthal direction
  integer   :: nz       = 3         ! num elements in axial     direction
  integer   :: po_g     = 3         ! polynomial order of geometry
  integer   :: po_v     = 8         ! polynomial order of velocity
  integer   :: po_q     = 12        ! polynomial order for convective operator
  integer   :: nt       = 1         ! number of test loops
  logical   :: periodic = .false.   ! switch for axial periodicity

  namelist/input/ conf, nt, r0, r1, h, nr, np, nz, po_g, po_v, po_q

  type(MPI_Comm) :: comm = MPI_COMM_WORLD
  integer :: rank

  type(GenericMesh_3D)         :: generic_mesh
  type(Mesh_3D)                :: mesh

  real(RNP), allocatable :: area(:)
  real(RNP) :: vol
  logical   :: passed
  integer   :: io
  integer   :: e, i, j, k, l, c, nv, nq, ne, s1, s2, s3, s4

  ! operators and variables ....................................................

  real(RNP), dimension(:,:,:,:,:), allocatable :: F_c_gen, F_c_libxsmm
  real(RNP), dimension(:,:,:,:,:), allocatable :: v

  real(RNP), dimension(:,:,:,:,:), allocatable :: vm
  real(RNP), dimension(:,:,:,:,:), allocatable :: vp

  type(SpectralElementVariable_3D) :: sev_v

  type(INS_Operator_3D) :: ins_op

  type(INS_Options_3D) :: ins_options

  real(RNP) :: time  = 0
  real(RNP) :: mflops_gen = -1, mlups_gen = -1
  real(RNP) :: error_opt = 0, mflops_opt = -1, mlups_opt = -1

  integer :: npop, nflop, nflop_f
  integer :: nvvv, nvvq, nvqq, nqqq
  integer :: nvv, nvq, nqq
  integer :: p, pm1, pm2

  logical :: interpolate
  logical :: lobatto

  integer(IXL) :: count, count0, rate

  !-----------------------------------------------------------------------------
  ! initialization
  call XMPI_Init()

  call MPI_Comm_rank(comm, rank)

  SO_FAR_ONLY_ONE_PROCESS: if (rank == 0) then

  ! read test parameters .......................................................

  open(newunit = io, file = input_file)
  read(io, nml = input)
  close(io)

  ! create and import generic mesh .............................................

  select case(conf)
  case(2)
    call generic_mesh % CreateAnnulus(r0, r1, h, nr, np, nz, po_g, periodic)
  case default
    call generic_mesh % CreateCylinder(r1, h, nr, nz, po_g, periodic)
  end select

  ! test parameters ............................................................

  ! verification ...............................................................

  ! import u-mesh
  call mesh % ImportGenericMesh(generic_mesh, comm = comm)

  ! verify u-mesh
  call VerifyMesh_3D(mesh, passed)

  write(*,'(/,A,G0,/)') 'VerifyMesh3d: passed = ', passed

  ! ins operators...............................................................

  ins_options % opt_v % po = po_v
  ins_options % opt_p % po = po_v - 1
  ins_options % opt_q % po = po_q

  call ins_op % Init(ins_options, mesh)

  ! workspace ..................................................................

  nv = po_v + 1
  nq = po_q + 1
  ne = mesh%n_elem ! number of elements

  ! order of test function .....................................................

  p = min(po_q, 2)

  pm1 = max(p - 1, 0)
  pm2 = max(p - 2, 0)

  ! spectral element functionality .............................................

  allocate( F_c_gen(0:po_v,0:po_v,0:po_v,ne,3),     &
            F_c_libxsmm(0:po_v,0:po_v,0:po_v,ne,3), &
            v(0:po_v,0:po_v,0:po_v,ne,3),           &
            vm(0:po_v,0:po_v,6,ne,3),               &
            vp(0:po_v,0:po_v,6,ne,3)                )

  associate( se_mesh_v => ins_op % sem_v                 &
           , x         => ins_op % sem_v % metrics % x   &
           , ins_mesh  => ins_op % mesh                  )

    do e = 1, ne
      do k = 0, po_v
      do j = 0, po_v
      do i = 0, po_v

        v(i,j,k,e,1) =  x(i,j,k,e,1) ** p  *  x(i,j,k,e,2) ** pm1
        v(i,j,k,e,2) =  x(i,j,k,e,2) ** p  *  x(i,j,k,e,3) ** pm1
        v(i,j,k,e,3) =  x(i,j,k,e,3) ** p  *  x(i,j,k,e,1) ** pm1

      end do
      end do
      end do
    end do

    allocate(area(ins_mesh % n_bound))
    call se_mesh_v % Get_Volume(vol)
    call se_mesh_v % Get_SurfaceAreas(area)

    write(*,'(A)') 'Spectral element mesh of order u'
    write(*,'(2X,A,G0)') 'volume  = ', vol
    do i = 1, ins_mesh % n_bound
      write(*,'(2X,A,I0,A,G0)') 'area(',i,') = ', area(i)
    end do

    sev_v = SpectralElementVariable_3D(se_mesh_v,v)

    ! traces
    call sev_v % GetTraces(vm)
    call ConvertTraces(ins_mesh, vm, vp)

  end associate

  !-----------------------------------------------------------------------------
  ! Operand und exact result

  associate ( I_vq => ins_op % iop_vq % A )
    lobatto = I_vq(1,1) == 1 .and. I_vq(nq,nv) == 1
    interpolate = ((nq /= nv)  .or. (.not. lobatto))
  end associate

  ! problem dimensions..........................................................

  npop = nv**3         ! number of operands per element

  nvvv = nv**3
  nvvq = nv**2 * nq
  nvqq = nv * nq**2
  nqqq = nq**3

  nvv = nv**2
  nvq = nv * nq
  nqq = nq**2

  ! calculate number of flops...................................................

  ! v-interpolation
  if (interpolate) then
    nflop = 3 * (nvvq * 2 * nv + nvqq * 2 * nv + nqqq * 2 * nv)
  else
    nflop = 0
  end if

  ! mass matrix
  nflop = nflop + 3 * nqqq

  ! MJ⁻¹⋅vv
  nflop = nflop + nqqq * (6 + 9 * 6)

  ! volume integral
  if (interpolate) then
    nflop = nflop + 3 * (nvqq * 2 * nq + nvvq * 2 * nq + nvvv * 2 * nq) &
                  + 6 * (nvqq * 2 * nq + nvvq * 2 * nq + nvvv * (2 * nq + 1))
  else
    nflop = nflop + 3 * (nvvv * 2 * nv) + 6 * (nvvv * (2 * nv + 1))
  end if

  ! face integrals
  nflop_f = 8             ! iv/iq
  if (interpolate) then
    if (lobatto) then
      nflop_f = nflop_f + 0
    else
      nflop_f = nflop_f + 3 * (nvq * 2 * nv + nqq * 2 * nv)
    end if
    nflop_f = nflop_f + 3 * (nvq * 2 * nv + nqq * 2 * nv)
  else
    nflop_f = nflop_f + 0
  end if
  nflop_f = nflop_f + nqq * (5 + 5 + 2 + 3 * 8)
  if (interpolate) then
    nflop_f = nflop_f + 3 * (nvq * 2 * nq + nvv * (2 * nq + 1))
  else
    nflop_f = nflop_f + 3 * nvv
  end if

  ! number of FLOPs per element
  nflop = nflop + 6 * nflop_f

  !-----------------------------------------------------------------------------
  ! Test generic procedure

  if (interpolate) then
    associate(  w_q  => ins_op % sop_q  % w               &
              , I_vq => ins_op % iop_vq % A               &
              , D_v  => ins_op % eop_v  % D               &
              , Ji_q => ins_op % sem_q  % metrics % Ji    &
              , Jd_q => ins_op % sem_q  % metrics % Jd    &
              , a_q  => ins_op % sem_q  % metrics % a     &
              , n_q  => ins_op % sem_q  % metrics % n     )

    !$omp parallel

    call TPO_INS_Convection_D_Gen(nv, nq, ne, &
             D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_gen)

    call system_clock(count0, rate)

    do i = 1, nt
      call TPO_INS_Convection_D_Gen(nv, nq, ne, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_gen)
    end do

    call system_clock(count)
    !$omp end parallel

    end associate
  else
    associate(  w_q  => ins_op % sop_q  % w               &
              , I_vq => ins_op % iop_vq % A               &
              , D_v  => ins_op % eop_v  % D               &
              , Ji_q => ins_op % sem_v  % metrics % Ji    &
              , Jd_q => ins_op % sem_v  % metrics % Jd    &
              , a_q  => ins_op % sem_v  % metrics % a     &
              , n_q  => ins_op % sem_v  % metrics % n     )

      !$omp parallel

      call TPO_INS_Convection_D_Gen(nv, nq, ne, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_gen)

      call system_clock(count0, rate)

      do i = 1, nt
        call TPO_INS_Convection_D_Gen(nv, nq, ne, &
                 D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_gen)
      end do

      call system_clock(count)
      !$omp end parallel
    end associate
  end if

  time = (count - count0) / real(rate, RNP) / max(nt, 1)

  mflops_gen = 1E-6 / time * ne * nflop
  mlups_gen  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! Test optimized procedure

  if (interpolate) then
    associate(  w_q  => ins_op % sop_q  % w               &
              , I_vq => ins_op % iop_vq % A               &
              , D_v  => ins_op % eop_v  % D               &
              , Ji_q => ins_op % sem_q  % metrics % Ji    &
              , Jd_q => ins_op % sem_q  % metrics % Jd    &
              , a_q  => ins_op % sem_q  % metrics % a     &
              , n_q  => ins_op % sem_q  % metrics % n     )

      !$omp parallel

      call TPO_INS_Convection(nv, nq, ne, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_libxsmm)

      call system_clock(count0, rate)

      do i = 1, nt
        call TPO_INS_Convection(nv, nq, ne, &
                 D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_libxsmm)
      end do

      call system_clock(count)

      !$omp end parallel
    end associate
  else
    associate(  w_q  => ins_op % sop_q  % w               &
              , I_vq => ins_op % iop_vq % A               &
              , D_v  => ins_op % eop_v  % D               &
              , Ji_q => ins_op % sem_v  % metrics % Ji    &
              , Jd_q => ins_op % sem_v  % metrics % Jd    &
              , a_q  => ins_op % sem_v  % metrics % a     &
              , n_q  => ins_op % sem_v  % metrics % n     )
      !$omp parallel

      call TPO_INS_Convection(nv, nq, ne, &
               D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_libxsmm)

      call system_clock(count0, rate)

      do i = 1, nt
        call TPO_INS_Convection(nv, nq, ne, &
                 D_v, I_vq, w_q, Jd_q, Ji_q, a_q, n_q, v, vp, F_c_libxsmm)
      end do

      call system_clock(count)

      !$omp end parallel
    end associate
  end if

  time = (count - count0) / real(rate, RNP) / nt

  error_opt  = maxval(abs(F_c_gen - F_c_libxsmm))
  mflops_opt = 1E-6 / time * ne * nflop
  mlups_opt  = 1E-6 / time * ne * npop

  !-----------------------------------------------------------------------------
  ! Print results

  write(*,*)
  write(*,'(3A)') '                                ',      &
                  '   ------ generic operator ---------',  &
                  '   ------ optimized operator -------'
  write(*,'(3A)') '   np_v   np_q      ne        nt    ',  &
                  '   error     MFLOP/s      MLUP/s    ',  &
                  '   error     MFLOP/s      MLUP/s'
  write(*,'(2(2X,I5),I8,2X,I8)',  advance='NO') nv, nq, ne, nt
  write(*,'(3(2X,ES10.3))', advance='NO') error_opt, mflops_gen, mlups_gen
  write(*,'(3(2X,ES10.3))') error_opt, mflops_opt, mlups_opt
  write(*,*)

  end if SO_FAR_ONLY_ONE_PROCESS

  call MPI_Finalize()

  !=============================================================================

end program Validate__TPO_INS_Convection_D
