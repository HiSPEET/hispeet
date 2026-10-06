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

!> summary:  Test of multilevel INS projection step
!> author:   Joerg Stiller
!> date:     2026/07/01
!>
!> If present, the first argument of the invoking command will be interpreted
!> as the base name of the control file. If omitted, the program looks for
!> `ml_ins_solver_3d.prm`.
!===============================================================================

program Test_ML_INS_Projection_3D
  use Kind_Parameters
  use Constants
  use OpenMP_Binding
  use XMPI
  use Array_Assignments
  use Execution_Control
  use Logging_Levels

  use Create_Cuboid_Cartesian

  use Mesh__3D
  use Generic_Mesh__3D
  use Verify_Mesh__3D

  use INS__Problem__3D
  use INS__Problem__Test_Suite__3D

  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  use ML__Boundary_Variable__3D
  use ML__INS__Operator__3D
  use ML__INS__Projection__3D
  use ML__INS__Flow_Characteristics__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! MPI and OpenMP variables ...................................................

  type(MPI_Comm) :: comm      ! MPI communicator
  integer        :: rank      ! local MPI rank
  integer        :: n_proc    ! number of MPI processes
  integer        :: n_thread  ! number of OpenMP threads

  ! control parameters .........................................................

  ! control of logging levels
  namelist/control_prm/ log_level
  namelist/control_prm/ log_level_inner_iteration
  namelist/control_prm/ log_level_outer_iteration
  namelist/control_prm/ log_level_multigrid_cycle

  character(len=*), parameter :: default_case = 'ml_ins_projection_3d'
  character(len=80) :: flow_case ! flow case name
  character(len=80) :: case_file ! flow case input file: trim(flow_case).prm

  character(len=80) :: flow_problem = 'VariableViscosity'  ! problem name
  character(len=80) :: problem_file = 'variable_viscosity' ! parameters file

  integer :: vtk_mode = 0 ! VTK export mode, 0/1/2/3: none/all/active/leafs

  namelist/control_prm/ flow_problem, problem_file, vtk_mode

  ! mesh .......................................................................

  type(Mesh_3D), allocatable, save :: base_mesh
  type(ML_Mesh_3D),           save :: ml_mesh

  type(ML_Mesh_Options_3D), save :: ml_mesh_opt
    ! multilevel mesh options, defining top level, partitioning and refinement,
    ! stacic components defined in namelist `ml_mesh_options_3d__static` and
    ! dynamic components in `ml_mesh_options_3d__dynamic`

  ! problem ....................................................................

  class(INS_Problem_3D), allocatable, save :: problem

  ! spatial operators ..........................................................

  integer, allocatable, save :: po(:)
  type(ML_INS_OperatorOptions_3D)  , save :: ml_ins_opt
  type(ML_INS_ProjectionOptions_3D), save :: ml_proj_opt

  namelist/spatial_prm/ po, ml_ins_opt, ml_proj_opt

  type(ML_MeshOperators_3D) , save :: ml_op
  type(ML_INS_Operator_3D)  , save :: ml_ins
  type(ML_INS_Projection_3D), save :: ml_proj

  ! time stepping ..............................................................

  real(RNP) :: t  = 0 ! problem time
  real(RNP) :: dt = 1 ! step size

  namelist/temporal_prm/ dt

  ! pressure disturbance .......................................................

  ! p' = a_w * sin(φ₁) sin(φ₂) sin(φ₃)
  ! φᵢ = 2π(xᵢ - x_wᵢ)/l_wᵢ

  integer   :: k_w    = 1 ! wave number
  real(RNP) :: x_w(3) = 0 ! wave starting position xw
  real(RNP) :: l_w(3) = 1 ! wave lengths
  real(RNP) :: a_w    = 1 ! wave amplitude

  namelist/disturbance_prm/ k_w, x_w, l_w, a_w

  ! variables ..................................................................

  type(ML_MeshVariable_3D),     save :: u   ! solution variables
  type(ML_MeshVariable_3D),     save :: q   ! disturbed pressure
  type(ML_MeshVariable_3D),     save :: mu  ! bulk viscosity  (not really used)
  type(ML_MeshVariable_3D),     save :: nu  ! shear viscosity (not really used)
  type(ML_MeshVariable_3D),     save :: vtk ! variables exported to VTK
  type(ML_BoundaryVariable_3D), save :: bv  ! boundary values

  ! auxiliaries ................................................................

  type(ML_INS_FlowCharacteristics_3D) :: ml_flow_char

  character(len=20), allocatable :: var_name(:)
  character(len=20), allocatable :: vtk_name(:)

  real(RNP) :: volume
  logical   :: exists
  integer   :: io, stat
  integer   :: l_top, l_max, n_bound
  integer   :: n_comp  ! number of solution components
  integer   :: n_var   ! number of variables
  integer   :: n_vtk   ! number of VTK quantities
  integer   :: l

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

  ! parameters .................................................................

  ! read control parameters
  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Multilevel projection step for incompressible flow'
    write(*,*)
    write(*,'(T3,A,T30,9(G0,X))') 'number of processes:', n_proc
    write(*,'(T3,A,T30,9(G0,X))') 'number of threads:'  , n_thread
    write(*,*)

    call get_command_argument(1, flow_case, status=stat)
    if (stat /= 0 .or. len_trim(flow_case) == 0) then
      flow_case = default_case
    end if
    case_file = trim(flow_case) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      read(io, nml = disturbance_prm)
      close(io)
    else
       call Error( 'ML_INS_Solver_3D', &
                   'input file "' // trim(case_file) // '" not found' )
    end if

  end if

  ! globalize logging levels
  call XMPI_Bcast_LoggingLevels(0, comm)

  ! globalize control parameters
  call XMPI_Bcast( flow_case   , 0, comm)
  call XMPI_Bcast( case_file   , 0, comm)
  call XMPI_Bcast( flow_problem, 0, comm)
  call XMPI_Bcast( problem_file, 0, comm)
  call XMPI_Bcast( vtk_mode    , 0, comm)

  ! globalize disturbance parameters
  call XMPI_Bcast( k_w, 0, comm)
  call XMPI_Bcast( x_w, 0, comm)
  call XMPI_Bcast( l_w, 0, comm)
  call XMPI_Bcast( a_w, 0, comm)

  ! create mesh ................................................................

  allocate(base_mesh)

  call CreateCuboidCartesian(comm, case_file, base_mesh)

  if (rank == 0) then
    write(*,'(2X,A)') 'creating multilevel mesh'
    ! read multilevel mesh options
    open(newunit = io, file = case_file)
    ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
    close(io)
  end if

  call ml_mesh_opt % Bcast(0, comm)
  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
  deallocate(base_mesh)

  l_top   = size(ml_mesh % mesh)
  l_max   = max(l_top, ml_mesh_opt % l_max)
  n_bound = ml_mesh % mesh(1) % n_bound

  ! problem ....................................................................

  call Set_INS_TestProblem_3D(problem, flow_problem, problem_file, n_bound, comm)

  ! check & fix boundary conditions
  block
    integer :: b
    do b = 1, n_bound
      if (ml_mesh % mesh(1) % boundary(b) % coupled > 0) then
        if (problem % bc_v(b) /= 'P') then
          if (rank == 0) then
            write(*,'(2X,A,X,I0)') '*** enforcing periodic BC on boundary', b
          end if
          problem % bc_v(b) = 'P'
          problem % bc_p(b) = 'P'
        end if
      end if
    end do
  end block

  ! spatial ....................................................................

  allocate(po(l_max), source = -1)

  ! read
  if (rank == 0) then
    write(*,'(/,A)') 'creating spatial operators'
    open(newunit = io, file = case_file)
    read(io, nml = spatial_prm)
    close(io)
  end if

  ! globalize
  call XMPI_Bcast(po, 0, comm)
  call ml_ins_opt  % Bcast(0, comm)
  call ml_proj_opt % Bcast(0, comm)

  ml_op = ML_MeshOperators_3D(ml_mesh, po)
  ml_ins = ML_INS_Operator_3D(ml_ins_opt, problem, ml_mesh, ml_op)
  ml_proj = ML_INS_Projection_3D(ml_proj_opt, ml_ins)

  call ml_ins % ml_op_u % Get_Volume(volume)

  ! variables ..................................................................

  n_comp = problem % nc
  n_var  = n_comp

  var_name = [ 'v_x', 'v_y', 'v_z', 'p  ']

  call u  % Init(ml_ins % ml_op_u, nc = n_var, name = var_name)
  call q  % Init(ml_ins % ml_op_u, nc = 1)
  call mu % Init(ml_ins % ml_op_u, nc = 1)
  call nu % Init(ml_ins % ml_op_u, nc = 1)
  call bv % Init(ml_ins % ml_op_u, nc = ml_ins % problem % nc)

  !-----------------------------------------------------------------------------
  ! Initial conditions and boundary values

  do l = 1, l_top
    associate( x_l    => ml_ins % ml_op_u % sem(l) % metrics % x  &
             , u_l    => u  % level(l) % val                      &
             , q_l    => q  % level(l) % val                      &
             , bv_l   => bv % level(l) % var                      )
      block
        real(RNP) :: alpha(3), kappa(3), phi(3)
        integer   :: b, i, j, k, e

        ! exact initial conditions .............................................

        call problem % GetInitialValues(x_l, u_l)

        ! add pressure disturbance .............................................

        kappa = 2 * k_w * PI / l_w
        alpha = kappa * a_w * dt

        do e = 1, ubound(u_l, 4)
          do k = 0, po(l)
          do j = 0, po(l)
          do i = 0, po(l)

            phi(1) = kappa(1) * (x_l(i,j,k,e,1) - x_w(1))
            phi(2) = kappa(2) * (x_l(i,j,k,e,2) - x_w(2))
            phi(3) = kappa(3) * (x_l(i,j,k,e,3) - x_w(3))

            ! disturbed velocity
            u_l(i,j,k,e,1) = u_l(i,j,k,e,1)  &
                           + alpha(1) * cos(phi(1)) * sin(phi(2)) * sin(phi(3))
            u_l(i,j,k,e,2) = u_l(i,j,k,e,2)  &
                           + alpha(2) * sin(phi(1)) * cos(phi(2)) * sin(phi(3))
            u_l(i,j,k,e,3) = u_l(i,j,k,e,3)  &
                           + alpha(3) * sin(phi(1)) * sin(phi(2)) * cos(phi(3))

            ! disturbed pressure
            q_l(i,j,k,e,1) = u_l(i,j,k,e,4)  &
                           - a_w * sin(phi(1)) * sin(phi(2)) * sin(phi(3))

            ! reset pressure
            u_l(i,j,k,e,4) = 0

          end do
          end do
          end do
        end do

        ! boundary values ......................................................

        do b = 1, size(bv_l)
          call bv_l(b) % Extract(u_l)
        end do

      end block
    end associate
  end do

  !-----------------------------------------------------------------------------
  ! Print info

  ! mesh
  if (rank == 0) then
    write(*,'(/,A)') 'multilevel mesh characteristics'
  end if
  call ml_ins % ml_op_u % Print_MeshCharacteristics()

  !-----------------------------------------------------------------------------
  ! Projection

  if (rank == 0) then
    write(*,*)
  end if

  ! initial flow characteristics
  call ml_flow_char % Evaluate(ml_ins, t, mu, nu, u, dt, volume, leaf = .true.)
  call ml_flow_char % PrintHeader()
  call ml_flow_char % PrintValues()

  ! projection step
  call ml_proj % ProjectionStep(dt, bv, u)

  ! add disturbed to computed pressure
  do l = 1, l_top
    associate( p_l => u % level(l) % val(:,:,:,:,4) &
             , q_l => q % level(l) % val(:,:,:,:,1) )

      p_l = p_l + q_l

    end associate
  end do

  ! remove mean pressure
  !call ml_ins % CalibratePressure(u)

  call ml_flow_char % Evaluate(ml_ins, t, mu, nu, u, dt, volume, leaf=.true.)
  call ml_flow_char % PrintValues()

  if (rank == 0) then
    write(*,*)
  end if

  !-----------------------------------------------------------------------------
  ! Write plot files

  EXPORT_VTK: if (vtk_mode > 0) then

    if (rank == 0) then
      write(*,'(A)') 'writing VTK file'
    end if

    if (problem % HasExactSolution()) then

      ! export data including exact solution and error .........................

      block
        integer :: i, k

        n_vtk = 3 * n_comp

        allocate(vtk_name(n_vtk))
        vtk_name(1:n_var) = var_name
        k = n_var
        do i = 1, n_comp
          write(vtk_name(k + i         ),'(2A)') trim(var_name(i)), '__exact'
          write(vtk_name(k + i + n_comp),'(2A)') trim(var_name(i)), '__error'
        end do

        call vtk % Init(ml_ins%ml_op_u, n_vtk, vtk_name)

        do l = 1, l_top
          call SetArray(vtk%level(l)%val(:,:,:,:,1:n_var), u%level(l)%val)
          associate( x_l => ml_ins % ml_op_u % sem(l) % metrics % x           &
                   , u_l => vtk % level(l) % val(:,:,:,:,1+0*n_comp:1*n_comp) &
                   , s_l => vtk % level(l) % val(:,:,:,:,1+1*n_comp:2*n_comp) &
                   , e_l => vtk % level(l) % val(:,:,:,:,1+2*n_comp:3*n_comp) )

            call problem % GetExactSolution(x_l, t, s_l)
            e_l = u_l - s_l
          end associate
        end do

        call vtk % ExportVTK(ml_ins%ml_op_u, file = flow_case, mode = vtk_mode)

      end block

    else

      ! export solution and averaged variables .................................

      call u % ExportVTK(ml_ins%ml_op_u, file = flow_case, mode = vtk_mode)

    end if

  end if EXPORT_VTK

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

  !=============================================================================

end program Test_ML_INS_Projection_3D
