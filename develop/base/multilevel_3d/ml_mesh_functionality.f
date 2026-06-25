!> summary:  Program for testing basic multilevel mesh functionality
!> author:   Joerg Stiller
!> date:     2024/07/02
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Mesh_Functionality
  use Kind_Parameters
  use Constants
  use Logging_Levels
  use XMPI
  use Create_Cuboid_Cartesian
  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D
  use Export_VTK_Mesh_SFC__3D
  use Child_To_Parent_Projection__3D
  use Child_To_Parent_Restriction__3D
  use Parent_To_Child_Interpolation__3D
  use ML__Mesh__3D
  use ML__Mesh_Operators__3D
  use ML__Mesh_Variable__3D
  implicit none

  ! variables ..................................................................

  integer :: domain = 0 ! 0: import GMSH, 1: create cuboidal

  character(len=*), parameter :: default_case = 'ml_mesh_functionality'
  character(len=80) :: case_name ! case name
  character(len=80) :: case_file ! case input file: trim(case_name).prm
  character(len=80) :: gmsh_file = '../../gmsh/cylinder_2d'

  logical :: restart    = .false. ! check restart functionality
  logical :: vtk_export = .false. ! switch for VTK export
  integer :: vtk_mode   =  3      ! 1/2/3: all/active/leaf elements

  integer, allocatable :: po(:) ! sequence of polynomial orders

  namelist /control/ log_level, domain, gmsh_file, restart, vtk_export, vtk_mode
  namelist /operators/ po

  type(GenericMesh_3D)     , allocatable, save :: generic_mesh
  type(Mesh_3D)            , allocatable, save :: base_mesh
  type(ML_Mesh_Options_3D) , allocatable, save :: ml_mesh_opt
  type(ML_Mesh_3D)         , allocatable, save :: ml_mesh
  type(ML_MeshOperators_3D), allocatable, save :: ml_op
  type(ML_MeshVariable_3D) , allocatable, save :: ml_var

  type(MPI_Comm) :: comm = MPI_COMM_WORLD

  character(len=9), save :: var_name(14), tag

  real(RNP), allocatable, save :: delta(:)
  real(RNP), allocatable, save :: delta_loc(:)
  real(RNP), save :: kappa = real(2 * PI, RNP)

  character(:), allocatable :: mesh_file
  character(:), allocatable :: data_file
  character(:), allocatable :: plot_file

  real(RNP) :: v
  logical :: passed, all_passed
  integer :: n_level, n_proc, rank, prm
  integer :: ne_max, ne_min, ne_tot
  integer :: stat
  integer :: e, i, j, k, l, nc

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  allocate(generic_mesh)
  allocate(base_mesh)
  allocate(ml_mesh_opt)
  allocate(ml_mesh)
  allocate(ml_op)
  allocate(ml_var)

  ! read parameters
  if (rank == 0) then

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) then
      case_name = default_case
    end if
    case_file = trim(case_name) // '.prm'

    write(*,'(/,A,/)') 'Testing basic multilevel mesh functionality'
    write(*,'(2X,2A)') 'reading input parameters from ', trim(case_file)

    open(newunit = prm, file = case_file)
    read(prm, nml = control)
    ml_mesh_opt = ML_Mesh_Options_3D(prm, n_proc)
    allocate(po(ml_mesh_opt%l_top))
    read(prm, nml = operators)
    close(prm)
  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  if (rank > 0) then
    allocate(po(ml_mesh_opt%l_top), source = -1)
  end if

  ! globalize remaining parameters
  call XMPI_Bcast(domain     , 0, comm)
  call XMPI_Bcast(gmsh_file  , 0, comm)
  call XMPI_Bcast(restart    , 0, comm)
  call XMPI_Bcast(vtk_export , 0, comm)
  call XMPI_Bcast(vtk_mode   , 0, comm)
  call XMPI_Bcast(po         , 0, comm)

  call XMPI_Bcast_LoggingLevels(0, comm)


  ! mesh import ................................................................

  select case(domain)
  case(0)
    if (rank == 0) then
      call ImportGMSH_3D(gmsh_file, generic_mesh)
    end if
    call base_mesh % ImportGenericMesh(generic_mesh, comm)
  case(1)
    call CreateCuboidCartesian(comm, case_file, base_mesh)
  end select

  if (rank == 0) then
    write(*,'(/,A)') 'initial mesh'
  end if

  call VerifyMesh_3D(base_mesh, passed)
  call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

  if (rank == 0) then
    write(*,'(2X,A,G0)') 'passed = ', all_passed
  end if

  call MPI_Barrier(comm)

  ! multilevel mesh ............................................................

  ml_mesh = ML_Mesh_3D(base_mesh, ml_mesh_opt)
  n_level = size(ml_mesh%mesh)

  associate(mesh => ml_mesh%mesh)

    ! verification
    do l = 1, n_level
      call VerifyMesh_3D(mesh(l), passed)
      call XMPI_Allreduce(passed, all_passed, MPI_LAND, comm)
      if (.not. all_passed) exit
    end do
    call MPI_Barrier(comm)
    if (rank == 0) then
      if (all_passed) then
        write(*,'(2X,9G0)') 'verification: all levels passed'
      else
        write(*,'(2X,9G0)') 'verification of level ',l,' failed'
      end if
    end if

    ! print info
    do l = 1, n_level
      if (mesh(l)%part >= 0) then
        call XMPI_Reduce(mesh(l)%n_elem, ne_min, MPI_MIN, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_max, MPI_MAX, 0, mesh(l)%comm_parts)
        call XMPI_Reduce(mesh(l)%n_elem, ne_tot, MPI_SUM, 0, mesh(l)%comm_parts)
      end if
      if (mesh(l)%part == 0) then
        write(*,'(2X,A,I4,A,I5,A,3(A,I6))') &
          'level ',l,': n_parts =',mesh(l)%n_parts,',  ', &
          'min/max/sum(n_elem) = ',ne_min,' / ',ne_max,' / ',ne_tot
      end if
    end do

  end associate

  ! multilevel operators .......................................................

  ml_op = ML_MeshOperators_3D(ml_mesh, po)

  ! multilevel variable ........................................................

  nc = size(var_name)
  var_name( 1) = 'mesh_part'
  var_name( 2) = 'elem_id'
  var_name( 3) = 'elem_type'
  var_name( 4) = 'elem_q_Js'
  var_name( 5) = 'var_order'
  var_name( 6) = 'v'
  var_name( 7) = 'Mv'
  var_name( 8) = 'Iv_p'
  var_name( 9) = 'Iv_p - v'
  var_name(10) = 'Iv_c'
  var_name(11) = 'Iv_c - v'
  var_name(12) = 'Pv_c'
  var_name(13) = 'Pv_c - v'
  var_name(14) = 'R(Mv_c)'

  call ml_var % Init(ml_op, nc, var_name)

  do l = 1, n_level
    associate( sem => ml_op  % sem(l)                &
             , eop => ml_op  % sem(l) % std_op       &
             , x   => ml_op  % sem(l) % metrics % x  &
             , Jd  => ml_op  % sem(l) % metrics % Jd &
             , var => ml_var % level(l) % val        )

      do e = 1, sem%mesh%n_elem
        var(:,:,:,e,1) = sem%mesh%part
        var(:,:,:,e,2) = e
        if (sem%mesh%element(e)%frozen) then
          var(:,:,:,e,3) = 0
        else if (sem%mesh%element(e)%adaptation%refinement < 1000) then
          var(:,:,:,e,3) = 1
        else
          var(:,:,:,e,3) = 2
        end if
        var(:,:,:,e,4) = minval(Jd(:,:,:,e)) / maxval(Jd(:,:,:,e))
        var(:,:,:,e,5) = po(l)

        ! test function values and mass-weighted values
        do k = 0, po(l)
        do j = 0, po(l)
        do i = 0, po(l)
          v = sin(kappa * (x(i,j,k,e,1) + x(i,j,k,e,2) + x(i,j,k,e,3)))
          var(i,j,k,e,6) = v
          var(i,j,k,e,7) = v * eop%w(i) * eop%w(j) * eop%w(k) * Jd(i,j,k,e)
        end do
        end do
        end do

        var(:,:,:,e, 8) = 0
        var(:,:,:,e, 9) = 0
        var(:,:,:,e,10) = 0
        var(:,:,:,e,11) = 0
        var(:,:,:,e,12) = 0
        var(:,:,:,e,13) = 0
        var(:,:,:,e,14) = 0

      end do

    end associate
  end do

  ! coarse-to-fine interpolation ...............................................

  allocate(delta(n_level-1), source = ZERO)
  allocate(delta_loc, source = delta)

  do l = 1, n_level - 1
    associate( parent => ml_op  % sem(l)     % mesh &
             , child  => ml_op  % sem(l+1)   % mesh &
             , var_p  => ml_var % level(l  ) % val  &
             , var_c  => ml_var % level(l+1) % val  &
             , iop    => ml_op  % iop_cf_x(l)       )

      call ParentToChildInterpolation_3D( parent, child, iop     &
                                        , v_p = var_p(:,:,:,:,6) &
                                        , v_c = var_c(:,:,:,:,8) )

      do e = 1, child%n_elem
        var_c(:,:,:,e,9) = var_c(:,:,:,e,8) - var_c(:,:,:,e,6)
        delta_loc(l) = max(delta_loc(l), maxval(abs(var_c(:,:,:,e,9))))
      end do

    end associate
  end do

  call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'coarse-to-fine interpolation error'
    do l = 1, n_level-1
      write(*,'(2X,5G0,ES10.3)') '|I v_',l,' - v_',l+1,'| =', delta(l)
    end do
  end if

  ! fine-to-coarse interpolation ...............................................

  do l = 1, n_level - 1
    associate( parent => ml_op  % sem(l)     % mesh &
             , child  => ml_op  % sem(l+1)   % mesh &
             , var_p  => ml_var % level(l  ) % val  &
             , var_c  => ml_var % level(l+1) % val  &
             , pop    => ml_op  % iop_fc_x(l+1)     )

      call ChildToParentProjection_3D( child, parent, pop &
                                     , var_c(:,:,:,:,6)   &
                                     , var_p(:,:,:,:,10)  )

      delta_loc(l) = 0
      do e = 1, parent%n_elem
        if (parent%element(e)%adaptation%refinement < 1000) cycle
        var_p(:,:,:,e,11) = var_p(:,:,:,e,10) - var_p(:,:,:,e,6)
        delta_loc(l) = max(delta_loc(l), maxval(abs(var_p(:,:,:,e,11))))
      end do

    end associate
  end do

  call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'fine-to-coarse interpolation error'
    do l = 1, n_level-1
      write(*,'(2X,5G0,ES10.3)') '|I v_',l+1,' - v_',l,'| =', delta(l)
    end do
  end if

  ! fine-to-coarse L2-projection ...............................................

  do l = 1, n_level - 1
    associate( parent => ml_op  % sem(l)     % mesh &
             , child  => ml_op  % sem(l+1)   % mesh &
             , var_p  => ml_var % level(l  ) % val  &
             , var_c  => ml_var % level(l+1) % val  &
             , pop    => ml_op  % pop_fc_x(l+1)     )

      call ChildToParentProjection_3D( child, parent, pop &
                                     , var_c(:,:,:,:,6)   &
                                     , var_p(:,:,:,:,12)  )

      delta_loc(l) = 0
      do e = 1, parent%n_elem
        if (parent%element(e)%adaptation%refinement < 1000) cycle
        var_p(:,:,:,e,13) = var_p(:,:,:,e,12) - var_p(:,:,:,e,6)
        delta_loc(l) = max(delta_loc(l), maxval(abs(var_p(:,:,:,e,13))))
      end do

    end associate
  end do

  call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'fine-to-coarse L²-projection error'
    do l = 1, n_level-1
      write(*,'(2X,5G0,ES10.3)') '|P v_',l+1,' - v_',l,'| =', delta(l)
    end do
  end if

  ! fine-to-coarse restriction .................................................

  do l = 1, n_level - 1
    associate( parent => ml_op  % sem(l)     % mesh &
             , child  => ml_op  % sem(l+1)   % mesh &
             , var_p  => ml_var % level(l  ) % val  &
             , var_c  => ml_var % level(l+1) % val  &
             , iop    => ml_op  % iop_cf_x(l)       )

      call ChildToParentRestriction_3D( child, parent, iop &
                                      , var_c(:,:,:,:,7)   &
                                      , var_p(:,:,:,:,14)  )

      delta_loc(l) = 0
      do e = 1, parent%n_elem
        if (parent%element(e)%adaptation%refinement < 1000) cycle
        v = abs(sum(var_p(:,:,:,e,14)) - sum(var_p(:,:,:,e,7)))
        delta_loc(l) = max(delta_loc(l), v)
      end do

    end associate
  end do

  call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'fine-to-coarse restriction element integral error'
    do l = 1, n_level-1
      write(*,'(2X,5G0,ES10.3)') '|Σ R(Mv_',l+1,') - Σ Mv_',l,'| =', delta(l)
    end do
  end if

  ! VTK export .................................................................

  if (vtk_export) then

    ! mesh and variables
    call ml_var % ExportVTK(ml_op, trim(plot_file), mode=vtk_mode)

    ! space filling curve
    do l = 1, n_level
      associate(mesh => ml_op % sem(l) % mesh)
        if (mesh % has_sfc) then
          write(tag,'(A,I0,A)') '_sfc_l', l
          call ExportVTK_MeshSFC(mesh, file = trim(plot_file)//trim(tag))
        end if
      end associate
    end do
  end if

  ! HDF5 write/read ............................................................

  if (restart) then

    if (rank == 0) then
      write(*,'(/,A)') 'writing and re-reading multilevel mesh and variables'
    end if

    mesh_file = trim(case_name) // '_mesh'
    data_file = trim(case_name) // '_data'

    ! write
    call ml_mesh % WriteHDF5(mesh_file)
    call ml_var  % WriteHDF5(data_file)

    ! destroy, allocate and read mesh
    deallocate(ml_mesh)
    allocate(ml_mesh)
    call ml_mesh % ReadHDF5(mesh_file, comm)

    ! rebuild operators
    ml_op = ML_MeshOperators_3D(ml_mesh, po)

    ! destroy, allocate and read variables
    deallocate(ml_var )
    allocate(ml_var )
    call ml_var % Init(ml_op, nc, var_name)
    call ml_var % ReadHDF5(data_file)

    ! verification
    passed = .true.
    do l = 1, n_level
      associate( sem => ml_op  % sem(l)                &
               , eop => ml_op  % sem(l) % std_op       &
               , x   => ml_op  % sem(l) % metrics % x  &
               , Jd  => ml_op  % sem(l) % metrics % Jd &
               , var => ml_var % level(l) % val        )

        do e = 1, sem%mesh%n_elem

          passed = passed .and. all(var(:,:,:,e,1) == sem%mesh%part)
          passed = passed .and. all(var(:,:,:,e,2) == e)

          if (sem%mesh%element(e)%frozen) then
            passed = passed .and. all(var(:,:,:,e,3) == 0)
          else if (sem%mesh%element(e)%adaptation%refinement < 1000) then
            passed = passed .and. all(var(:,:,:,e,3) == 1)
          else
            passed = passed .and. all(var(:,:,:,e,3) == 2)
          end if

          passed = passed .and. all( abs( var(:,:,:,e,4)      &
                                        - minval(Jd(:,:,:,e)) &
                                        / maxval(Jd(:,:,:,e)) &
                                        ) < epsilon(ONE) )

          passed = passed .and. all(var(:,:,:,e,5) == po(l))

          do k = 0, po(l)
          do j = 0, po(l)
          do i = 0, po(l)
            v = sin(kappa * (x(i,j,k,e,1) + x(i,j,k,e,2) + x(i,j,k,e,3)))
            passed = passed .and. abs( var(i,j,k,e,6) - v ) < epsilon(ONE)
            passed = passed .and. abs( var(i,j,k,e,7)      &
                                     - v * eop%w(i)        &
                                         * eop%w(j)        &
                                         * eop%w(k)        &
                                         * Jd(i,j,k,e)     &
                                     ) < epsilon(ONE)
          end do
          end do
          end do

        end do

      end associate
    end do

    call XMPI_Reduce(passed, all_passed, MPI_LAND, 0, comm)

    if (rank == 0) then
      write(*,'(2X,A,G0)') 'passed = ', all_passed
    end if

  end if

  ! finalization ...............................................................

  if ( allocated(delta_loc) ) deallocate(delta_loc)
  if ( allocated(delta)     ) deallocate(delta)

  call MPI_Finalize()

  !=============================================================================

end program ML_Mesh_Functionality
