!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Testing elliptic multilevel solvers with static refinement
!> author:   Joerg Stiller
!> date:     2024/09/16
!===============================================================================

program ML_Elliptic_Test_Static
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Execution_Control
  use Logging_Levels
  use Array_Assignments
  use Array_Reductions

  use Elliptic_Problem__3D
  use Elliptic_Problem__Knotty__3D
  use Elliptic_Problem__Simple__3D
  use Elliptic_Problem__Sphere__3D
  use Elliptic_Problem__TGV_Pressure__3D

  use Create_Cuboid_Cartesian
  use Create_Cuboid_Diamonds
  use Create_Cuboid_OneRotated
  use Create_Cylinder
  use Create_Annulus
  use Create_Spherical_Shell_Segment

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use Volume_Integrals__3D
  use VTK__Export_Mesh_SFC__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__Array_Reductions__3D
  use ML__DG__Elliptic_Solver__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP .............................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control .....................................................................

  ! NOTE
  ! The case name is defined by the first command line argument.
  ! If no argument is given, the default case is assumed.

  character(len=*), parameter :: default_case = 'ml_elliptic_test_static'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(case_name).prm

  integer :: solution_method = 21
  ! 11:  CS-MG
  ! 12:  CS-MGCG
  ! 21:  FAS-MG

  logical :: hdf5_check = .false. ! write and re-read ML mesh before solving
  logical :: vtk_export = .false. ! switch for VTK export
  integer :: vtk_mode   =  3      ! 1/2/3: all/active/leaf elements

  namelist/control_prm/ solution_method, hdf5_check, vtk_export, vtk_mode

  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  ! domain .....................................................................

  integer :: test_domain = 1 ! computational domain
  ! configuration (u/s = un/structured, r = regular, d = deformed)
  !   1  cuboidal domain with Cartesian mesh                               (s+r)
  !   2  cuboidal domain with unstructured "diamond" mesh                  (u+d)
  !   3  cuboidal domain with  3x3x3 elements and rotated center           (u+r)
  !   4  cylindrical domain                                                (u+d)
  !   5  annular domain                                                    (u+d)
  !   6  spherical shell segment                                           (u+d)

  character(len=80) :: gmsh_file = '../../gmsh/cylinder_2d'

  namelist/domain_prm/ test_domain, gmsh_file

  ! mesh .......................................................................

  type(GenericMesh_3D)     , save :: generic_mesh
  type(Mesh_3D)            , save :: base_mesh
  type(ML_Mesh_Options_3D) , save :: ml_mesh_opt
  type(ML_Mesh_3D)         , save :: ml_mesh

  ! problem ....................................................................

  integer :: test_problem = 3 ! 1...6: simple_{1/2/3}d/knotty/TGV/sphere
  integer :: start_values = 0 ! 0/1/2: zero, exact, random

  namelist/problem_prm/ test_problem, start_values

  ! common
  real(RNP) :: lambda =  0    ! Helmholtz parameter
  real(RNP) :: nu_0   =  1    ! constant or mean diffusivity ν₀
  real(RNP) :: nu_1   =  0    ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: d_nu   =  0    ! diffusivity fluctuation phase shift
  integer   :: k_nu   =  1    ! diffusivity fluctuation wave number
  integer   :: k_u    =  1    ! solution wave number

  namelist/problem_prm/ lambda, nu_0, nu_1, d_nu, k_nu, k_u

  ! sphere
  real(RNP) :: alpha  =  200  ! radial scaling factor
  real(RNP) :: x_c(3) = -0.05 ! sphere center
  real(RNP) :: r_f(1) =  0.7  ! sphere radius

  namelist/problem_prm/ alpha, x_c, r_f

  ! boundary conditions {'D','N','P'} ['D']
  character, allocatable :: bc(:)

  namelist/problem_prm/  bc

  class(EllipticProblem_3D), allocatable, save :: problem
  character(len=80) :: problem_name = ''

  ! solvers .....................................................................

  integer, allocatable, save :: po(:) ! sequence of polynomial orders
  integer :: fc_smooth = 0 ! fine-to-coarse discontinuity handling {0,1,2}
  integer :: fc_filter = 0 ! fine-to-coarse filter order, 0: no filtering

  type(ML_MeshOperators_3D),      save :: ml_op
  type(ML_DG_EllipticSolver_3D),  save :: ml_elliptic
  type(ML_DG_EllipticOptions_3D), save :: ml_elliptic_opt

  namelist/solver_prm/ po, fc_smooth, fc_filter, ml_elliptic_opt

  ! variables ...................................................................

  character(len=:), allocatable, save :: name_var(:) ! names of variables

  type(ML_MeshVariable_3D), save :: var ! container of variables
  type(ML_MeshVariable_3D), save :: mm  ! diagonal mass matrix
  type(ML_MeshVariable_3D), save :: nu  ! diffusivity
  type(ML_MeshVariable_3D), save :: f   ! sources
  type(ML_MeshVariable_3D), save :: s   ! exact solution
  type(ML_MeshVariable_3D), save :: u   ! numerical solution
  type(ML_MeshVariable_3D), save :: e   ! error
  type(ML_MeshVariable_3D), save :: r   ! residual or just workspace

  type(ML_BoundaryVariable_3D), save :: bv ! boundary values

  ! auxiliaries ................................................................

  character(:), allocatable, save :: domain_name
  real(RNP), dimension(:), allocatable, save :: ea0_0, ra0_2
  real(RNP), dimension(:), allocatable, save :: ean_0, ran_2
  real(RNP), save :: e0_0, en_0, e0_1, en_1, r0_2, rn_2
  real(RNP), save :: t_start, t_solve
  integer,   save :: n_i

  logical :: exists, passed, all_passed, singular
  integer :: io, stat, dim
  integer :: l_top

  !-----------------------------------------------------------------------------
  ! Initialization

  ! MPI and OpenMP .............................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  !$omp parallel
  n_thread = OMP_Num_Threads()
  !$omp end parallel

  ! control and domain parameters ..............................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Validation of elliptic multilevel solvers'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) then
      case_name = default_case
    end if
    case_file = trim(case_name) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      read(io, nml = domain_prm)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      close(io)
    else
      call Error('ML_Elliptic_Test_Static', &
                 'file "'// trim(case_file) //'" not found')
    end if

  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize remaining parameters
  call XMPI_Bcast(case_name       , 0, comm)
  call XMPI_Bcast(hdf5_check      , 0, comm)
  call XMPI_Bcast(vtk_export      , 0, comm)
  call XMPI_Bcast(vtk_mode        , 0, comm)
  call XMPI_Bcast(solution_method , 0, comm)
  call XMPI_Bcast(test_domain     , 0, comm)
  call XMPI_Bcast(gmsh_file       , 0, comm)

  ! base mesh ..................................................................

  select case(test_domain)
  case(0)
    if (rank == 0) then
      call ImportGMSH_3D(gmsh_file, generic_mesh)
    end if
    call base_mesh % ImportGenericMesh(generic_mesh, comm)
    domain_name = gmsh_file
  case(1)
    call CreateCuboidCartesian(comm, case_file, base_mesh)
    domain_name = 'Cuboidal domain with Cartesian mesh'
  case(2)
    call CreateCuboidDiamonds(comm, case_file, base_mesh)
    domain_name = 'Cuboidal domain with unstructured "diamond" mesh'
  case(3)
    call CreateCuboidOneRotated(comm, case_file, base_mesh)
    domain_name = 'Cuboidal domain with 3x3x3 elements and rotated center'
  case(4)
    call CreateCylinder(comm, case_file, base_mesh)
    domain_name = 'Cylindrical domain with unstructured mesh'
  case(5)
    call CreateAnnulus(comm, case_file, base_mesh)
    domain_name = 'Annular domain with unstructured mesh'
  case(6)
    call CreateSphericalShellSegment(comm, case_file, base_mesh)
    domain_name = 'Spherical shell segment with (un)structured mesh'
  end select

  if (rank == 0) then
    write(*,'(/,2A)') 'verifying base mesh for ', trim(domain_name)
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'passed = ', all_passed
  end if

  call MPI_Barrier(comm)

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
  l_top   = size(ml_mesh%mesh)

  if (hdf5_check) then
    call ml_mesh % WriteHDF5(case_name)
    deallocate(ml_mesh % mesh)
    call ml_mesh % ReadHDF5(case_name, comm)
  end if

  associate(mesh => ml_mesh%mesh)
    block
      integer :: l

      ! verification
      do l = 1, l_top
        call VerifyMesh_3D(mesh(l), passed)
        call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
        if (.not. all_passed) exit
      end do
      call MPI_Barrier(comm)
      if (rank == 0) then
        if (all_passed) then
          write(*,'(2X,A,/)') 'verification: all levels passed'
        else
          write(*,'(2X,A,I0,A,/)') 'verification of level ',l,' failed'
        end if
      end if

    end block
  end associate

  ! problem ....................................................................

  allocate(bc(base_mesh % n_bound), source = 'D')

  if (rank == 0) then
    write(*,'(/,A)') 'initializing elliptic problem'
    open(newunit = io, file = case_file)
    read(io, nml = problem_prm)
    close(io)
  end if

  ! globalize problem parameters
  call XMPI_Bcast( test_problem   , 0, comm )
  call XMPI_Bcast( start_values   , 0, comm )
  call XMPI_Bcast( lambda         , 0, comm )
  call XMPI_Bcast( nu_0           , 0, comm )
  call XMPI_Bcast( nu_1           , 0, comm )
  call XMPI_Bcast( d_nu           , 0, comm )
  call XMPI_Bcast( k_nu           , 0, comm )
  call XMPI_Bcast( k_u            , 0, comm )
  call XMPI_Bcast( alpha          , 0, comm )
  call XMPI_Bcast( x_c            , 0, comm )
  call XMPI_Bcast( r_f            , 0, comm )
  call XMPI_Bcast( bc             , 0, comm )

  select case(test_problem)
  case(1:3)
    dim = test_problem
    write(problem_name,'(A,I0,A)') 'Simple_', dim, 'D'
    problem = EllipticProblem_Simple_3D(lambda, nu_0, nu_1, d_nu, k_nu, k_u, dim)
  case(4)
    problem_name = 'Knotty'
    problem = EllipticProblem_Knotty_3D(lambda, nu_0, nu_1, d_nu, k_nu, k_u)
  case(5)
    problem_name = 'TGV_Pressure'
    problem = EllipticProblem_TGV_Pressure_3D(lambda, nu_0, nu_1, d_nu, k_nu)
  case(6)
    problem_name = 'Sphere'
    problem = EllipticProblem_Sphere_3D(lambda, alpha, x_c, r_f)
  end select

  ! enforce periodicity at coupled boundaries
  where(base_mesh % boundary % coupled > 0) bc = 'P'

  singular = lambda == ZERO .and. all(bc == 'P' .or. bc == 'N')

  if (rank == 0) then
    write(*,'(T3,A,T26,9(G0,X))') 'problem name:'        , trim(problem_name)
    write(*,'(T3,A,T26,9(G0,X))') 'variable diffusivity:', problem % nu_1 > 0
    write(*,'(T3,A,T26,9(G0,X))') 'boundary conditions:' , bc
  end if

  ! solver ......................................................................

  if (rank == 0) then
    allocate(po(1000), source = -1)
    open(newunit = io, file = case_file)
    read(io, nml = solver_prm)
    close(io)
    write(*,'(/,A)') 'initializing multilevel operators'
    po = po(1:l_top)
  else
    allocate(po(l_top), source = -1)
  end if

  call XMPI_Bcast(po       , 0, comm)
  call XMPI_Bcast(fc_smooth, 0, comm)
  call XMPI_Bcast(fc_filter, 0, comm)
  ml_op = ML_MeshOperators_3D(ml_mesh, po, smooth=fc_smooth, filter=fc_filter)

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel solver'
  end if

  call ml_elliptic_opt % Bcast(0, comm)
  ml_elliptic = ML_DG_EllipticSolver_3D(ml_op, ml_elliptic_opt)

  ! print info
  call ml_op % Print_MeshCharacteristics()

  ! variables ..................................................................

  if (rank == 0) then
    write(*,'(/,A)') 'initializing multilevel variables'
  end if

  name_var = [ 'nu', 'f ', 's ', 'u ', 'e ', 'r ' ]

  call var % Init(ml_op, size(name_var), name_var)

  ! handles for accessing individual variables
  call var % GetSlice(nu, first = 1, last = 1)
  call var % GetSlice(f , first = 2, last = 2)
  call var % GetSlice(s , first = 3, last = 3)
  call var % GetSlice(u , first = 4, last = 4)
  call var % GetSlice(e , first = 5, last = 5)
  call var % GetSlice(r , first = 6, last = 6)

  ! structure for keeping diagonal mass matrices
  call mm % Init(ml_op, 1)

  ! structure for keeping boundary values
  call bv % Init(ml_op, nc = 1)

  allocate(ea0_0(l_top), source = ZERO)
  allocate(ean_0, ra0_2, ran_2, source = ea0_0)

  block
    real(RNP), allocatable, save :: q(:,:,:,:,:)
    integer :: i, l

    do l = 1, l_top
      associate( sem_l => ml_op % sem(l)                 &
               , x_l   => ml_op % sem(l) % metrics % x   &
               , mm_l  => mm % level(l) % val(:,:,:,:,1) &
               , nu_l  => nu % level(l) % val(:,:,:,:,1) &
               , f_l   => f  % level(l) % val(:,:,:,:,1) &
               , r_l   => r  % level(l) % val(:,:,:,:,1) &
               , s_l   => s  % level(l) % val(:,:,:,:,1) &
               , u_l   => u  % level(l) % val(:,:,:,:,1) &
               , bu_l  => bv % level(l) % var            )

        call sem_l % Get_DG_DiagonalMassMatrix(mm_l)

        call problem % GetExactSolution (x_l, s_l)
        call problem % GetDiffusivity   (x_l, nu_l)
        call problem % GetSource        (x_l, r_l)     ! r = λ u - ∇·(ν ∇u)
        f_l = mm_l * r_l                               ! f = M r

        ! start values
        select case(start_values)
        case(0)
          call SetArray(u_l, ZERO)
        case(1)
          call SetArray(u_l, s_l)
        case default
          call random_number(u_l)
          u_l = 2*u_l - 1
        end select

        if (any(bc == 'N')) then
          allocate(q, mold = x_l)
          call problem % GetExactGradient(x_l, q)
          q(:,:,:,:,1) = nu_l * q(:,:,:,:,1)
          q(:,:,:,:,2) = nu_l * q(:,:,:,:,2)
          q(:,:,:,:,3) = nu_l * q(:,:,:,:,3)
        end if

        ! extract boundary conditions
        do i = 1, size(bc)
          select case(bc(i))
          case('D')
            call bu_l(i) % Extract(s_l)
          case('N')
            call bu_l(i) % ExtractNormalComponent(sem_l, q)
          end select
        end do

        if (allocated(q)) deallocate(q)

      end associate
    end do

  end block

  !-----------------------------------------------------------------------------
  ! Initial error and residual norms

  call Evaluation(e0_0, e0_1, r0_2, ea0_0, ra0_2)

  !-----------------------------------------------------------------------------
  ! Solution

  if (rank == 0) then
    write(*,'(/,A)') 'executing multilevel solver'
    t_start = MPI_Wtime()
  end if

  if (problem % nu_1 > 0) then
    select case(solution_method)
    case(11)
      call ml_elliptic % CS_MG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    case(12)
      call ml_elliptic % CS_MGCG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    case(21)
      call ml_elliptic % FAS_MG_Solver(bc, lambda, nu, bv, f, u, ni=n_i)
    end select
  else
    select case(solution_method)
    case(11)
      call ml_elliptic % CS_MG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    case(12)
      call ml_elliptic % CS_MGCG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    case(21)
      call ml_elliptic % FAS_MG_Solver(bc, lambda, nu_0, bv, f, u, ni=n_i)
    end select
  end if

  if (rank == 0) then
    t_solve = MPI_Wtime() - t_start
  end if

  !-----------------------------------------------------------------------------
  ! Evaluation

  block
    real(RNP) :: a, ck, cs, q, rate, t10, tau10, v10
    integer   :: l, l_min, n10, np_leaf
    integer, allocatable :: n_leaf(:,:)

    if (rank == 0) then
      write(*,'(/,A)') 'evaluation'
    end if

    ! get number of leaf elements
    call ml_op % Get_MeshCharacteristics(n_leaf)
    np_leaf = sum(n_leaf(:,4) * (po+1)**3)

    ! final error and residual norms
    call Evaluation(en_0, en_1, rn_2, ean_0, ran_2)

    if (rank == 0) then
      write(*,*)
      write(*,'(A,I0)')     'num iterations:  n_i = ', n_i
      write(*,'(A,ES10.3)') 'solution time:   t_s =' , t_solve

      write(*,'(/,A,/)') 'composite error and residual'

      write(*,'(2X,A,ES10.3)') ' e0_0   = ', e0_0
      write(*,'(2X,A,ES10.3)') ' e0_1   = ', e0_1
      write(*,'(2X,A,ES10.3)') ' r0_2   = ', r0_2
      write(*,*)
      write(*,'(2X,A,ES10.3)') ' en_0   = ', en_0
      write(*,'(2X,A,ES10.3)') ' en_1   = ', en_1
      write(*,'(2X,A,ES10.3)') ' rn_2   = ', rn_2
      write(*,*)

      rate  = log10(r0_2 / rn_2) / n_i
      n10   = ceiling(10 / rate)
      t10   = t_solve / n_i * n10
      tau10 = t10 * 1E6 * n_proc / np_leaf

      ! equivalent number of V-cycles
      associate( ns_1 => ml_elliptic % ns_1 &
               , ns_2 => ml_elliptic % ns_2 &
               , ns_c => ml_elliptic % ns_c )

        ! Krylov acceleration cost
        if (solution_method == 12) then
          ck = 1
        else
          ck = 0
        end if

        ! smoother cost
        if (ml_elliptic%smooth_method >= 3) then
          cs = 2
        else
          cs = 1
        end if

        ! refinement rate
        q = 8

        ! cost ratio
        if (po(l_top) > po(1)) then
          ! assume p-refinement
          a = 1 / (2*q)
        else
          ! assume h-refinement
          a = 1 / q
        end if
        a = (1 - a) / (1 - a**l_top)

        v10 = ( ((ns_1 + ns_2) * cs + 1) * n10          &
              - (ns_1 + ns_2 - ns_c) * cs * a * (n10-1) &
              + ck * a * n10                            &
              ) / 3

      end associate

      write(*,'(X,A7,X,2(7X,A3),6X,A5,7X,A8)') &
          '-lg ρ', 'n10', 'v10', 't10/s', 'τ10/μs'
      write(*,'(G11.3,I7,F10.1,3X,ES10.3,3X,ES9.2)') rate, n10, v10, t10, tau10

      write(*,'(/,A,/)') 'errors and residuals over active elements per level'

      select case(solution_method)
      case(11,12)
        l_min = l_top
      case default
        l_min = 1
      end select

      write(*,'(4X,A,3X,7(A7,5X))') &
        'l', '  e0_0 ', '  en_0 ', '  r0_2 ', '  rn_2 ', '-lg ρ'
      do l = l_min, l_top
        write(*,'(I5,7(2X,ES10.3))') &
           l, ea0_0(l), ean_0(l), ra0_2(l), ran_2(l), &
           log10(ra0_2(l) / ran_2(l)) / n_i
      end do

    end if
  end block

  !-----------------------------------------------------------------------------
  ! VTK export

  if (vtk_export) then
    vtk_mode = max(1, min(3, vtk_mode))
    call var % ExportVTK(ml_op, trim(case_name), vtk_mode)
    if (vtk_mode < 3) then
      associate(mesh => ml_mesh%mesh)
        block
          character(len=9) :: tag
          integer :: l
          do l = 1, size(mesh)
            if (mesh(l) % has_sfc) then
              write(tag,'(A,I0,A)') '_sfc_l', l
              call VTK_ExportMeshSFC_3D(mesh(l), file = trim(case_name)//trim(tag))
            end if
          end do
        end block
      end associate
    end if
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Evaluation

  subroutine Evaluation(e_0, e_1, r_2, ea_0, ra_2)
    real(RNP),             intent(out) :: e_0    !< H0 composite leaf error
    real(RNP),             intent(out) :: e_1    !< H1 composite leaf error
    real(RNP),             intent(out) :: r_2    !< E2 composite leaf residual
    real(RNP), contiguous, intent(out) :: ea_0(:) !< H0 active error    @ level
    real(RNP), contiguous, intent(out) :: ra_2(:) !< E2 active residual @ level

    real(RNP) :: avg, vol
    integer   :: l

    ! error ε ..................................................................

    do l = 1, l_top
      associate( e_l => e % level(l) % val(:,:,:,:,1) &
               , r_l => r % level(l) % val(:,:,:,:,1) &
               , s_l => s % level(l) % val(:,:,:,:,1) &
               , u_l => u % level(l) % val(:,:,:,:,1) )

        e_l = u_l - s_l
        r_l = 1

      end associate
    end do

    ! calibration
    if (singular) then
      vol = ML_WeightedScalarProduct_3D(mm, r, r, leaf = .true.)
      avg = ML_WeightedScalarProduct_3D(mm, e, r, leaf = .true.) / vol
      do l = 1, l_top
        associate(e_l => e % level(l) % val(:,:,:,:,1))
          e_l = e_l - avg
        end associate
      end do
    end if

    ! H0 error ε₀ ..............................................................

    ! composite H0 norm
    e_0 = sqrt(ML_WeightedScalarProduct_3D(mm, e, e, leaf = .true.))

    ! H0 norm over active elements of present level
    do l = 1, l_top
      associate( e_l => e % level(l) % val(:,:,:,:,1) &
               , r_l => r % level(l) % val(:,:,:,:,1) )

        r_l = e_l ** 2
        call GetVolumeIntegral(ml_op%sem(l), r_l, ea_0(l))
        ea_0(l) = sqrt(ea_0(l))

      end associate
    end do

    ! H1 error ε₁ ..............................................................

    do l = 1, l_top
      associate( ell_l => ml_elliptic % elliptic_op(l)   &
               , nu_l  => nu % level(l) % val(:,:,:,:,1) &
               , e_l   => e  % level(l) % val(:,:,:,:,1) &
               , r_l   => r  % level(l) % val(:,:,:,:,1) )

        ! r = Aε
        if (problem % nu_1 > 0) then
          call ell_l % Apply(bc, lambda, nu_l, u = e_l, r = r_l)
        else
          call ell_l % Apply(bc, lambda, nu_0, u = e_l, r = r_l)
        end if

      end associate
    end do

    ! ε₁ = √ εAε
    e_1 = sqrt(ML_ScalarProduct_3D(e, r, leaf = .true.))

    ! residual .................................................................

    if (problem % nu_1 > 0) then
      call ml_elliptic % FAS_MG_Residual(bc, lambda, nu, bv, f, u, r)
    else
      call ml_elliptic % FAS_MG_Residual(bc, lambda, nu_0, bv, f, u, r)
    end if

    ! E2 composite leaf residual
    r_2 = sqrt(ML_ScalarProduct_3D(r, r, leaf = .true.))

    ! E2 active residual @ level
    do l = 1, l_top
      associate( mesh_l => ml_op % sem(l) % mesh         &
               , r_l    => r % level(l) % val(:,:,:,:,1) )

        if (mesh_l % part >= 0) then
          ra_2(l) = sqrt(ScalarProduct(r_l, r_l, mesh_l%comm_parts))
        end if

      end associate
    end do

  end subroutine Evaluation

  !=============================================================================

end program ML_Elliptic_Test_Static
