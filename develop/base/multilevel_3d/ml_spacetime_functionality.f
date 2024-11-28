!> summary:  Program for testing basic multilevel spacetime functionality
!> author:   Joerg Stiller
!> date:     2024/11/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

program ML_Spacetime_Functionality
  use Kind_Parameters
  use Constants
  use Logging_Levels
  use Array_Assignments
  use XMPI

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D

  use ML__Mesh__3D
  use ML__Spacetime_Operators__3D
  use ML__Spacetime_Variable__3D

! use CF__Spacetime_Interpolation__3D
! use FC__Spacetime_Projection__3D
! use FC__Spacetime_Restriction__3D

  implicit none

  ! variables ..................................................................

  character(len=100) :: gmsh_file = '../gmsh_3d/cylinder_2d'
  character(len=100) :: plot_file = ''

  integer, allocatable :: po_x(:)    ! polynomial orders in x
  integer, allocatable :: po_t(:)    ! polynomial orders in t
  integer, allocatable :: ne_t(:)    ! elements per slice in t

  character(len=2) :: nodes_x = 'L'  ! type of space nodes
  character(len=2) :: nodes_t = 'RR' ! type of time nodes {E,L,RR}
  integer          :: smooth  =  0   ! discontinuity smoothing {0,1,2}

  namelist /control/   log_level, gmsh_file, plot_file
  namelist /operators/ po_x, po_t, ne_t, nodes_t, smooth

  type(GenericMesh_3D)          , save :: generic_mesh
  type(Mesh_3D)                 , save :: base_mesh
  type(ML_Mesh_Options_3D)      , save :: ml_mesh_opt
  type(ML_Mesh_3D)              , save :: ml_mesh
  type(ML_SpacetimeOperators_3D), save :: ml_op
  type(ML_SpacetimeVariable_3D) , save :: ml_var

  type(MPI_Comm) :: comm = MPI_COMM_WORLD

  character(len=9), save :: var_name(9)

!!   real(RNP), allocatable, save :: delta(:)
!!   real(RNP), allocatable, save :: delta_loc(:)
  real(RNP), save :: kappa = real(2 * PI, RNP)

  real(RNP), allocatable :: t(:)
  real(RNP) :: dt
  logical :: passed, all_passed
  integer :: n_level, n_proc, rank, prm
  integer :: ne_max, ne_min, ne_tot

  integer ::  e, i, j, k, l, m, n, nc

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  ! read parameters
  if (rank == 0) then
    write(*,'(/,A,/)') 'Testing basic multilevel spacetime functionality'
    write(*,'(2X,A)') 'reading input parameters'
    open(newunit = prm, file = 'ml_spacetime_functionality.prm')
    read(prm, nml = control)
    ml_mesh_opt = ML_Mesh_Options_3D(prm, n_proc)
    allocate(po_x(ml_mesh_opt%l_top), source = -1)
    allocate(po_t, ne_t, source = po_x)
    read(prm, nml = operators)
    close(prm)
  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  if (rank > 0) then
    allocate(po_x(ml_mesh_opt%l_top), source = -1)
    allocate(po_t, ne_t, source = po_x)
  end if

  ! globalize remaining parameters
  call XMPI_Bcast(gmsh_file , 0, comm)
  call XMPI_Bcast(plot_file , 0, comm)
  call XMPI_Bcast(po_x      , 0, comm)
  call XMPI_Bcast(po_t      , 0, comm)
  call XMPI_Bcast(ne_t      , 0, comm)
  call XMPI_Bcast(nodes_x   , 0, comm)
  call XMPI_Bcast(nodes_t   , 0, comm)
  call XMPI_Bcast(smooth    , 0, comm)

  ! mesh import ................................................................

  if (rank == 0) then
    call ImportGMSH_3D(gmsh_file, generic_mesh)
  end if

  call base_mesh % ImportGenericMesh(generic_mesh, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'verifying imported mesh'
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

  ! multilevel spacetime operators .............................................

  ml_op = ML_SpacetimeOperators_3D( ml_mesh, po_x, po_t, ne_t &
                                  , nodes_x, nodes_t, smooth  )

  ! multilevel spacetime variable ..............................................

  nc = size(var_name)
  var_name(1) = 'v'
  var_name(2) = 'Mv'
  var_name(3) = 'Iv_p'
  var_name(4) = 'Iv_p - v'
  var_name(5) = 'Iv_c'
  var_name(6) = 'Iv_c - v'
  var_name(7) = 'Pv_c'
  var_name(8) = 'Pv_c - v'
  var_name(9) = 'R(Mv_c)'

  ml_var = ML_SpacetimeVariable_3D(ml_op, nc, var_name)

  do l = 1, n_level
    associate( sem   => ml_op  % sem(l)                &
             , eop   => ml_op  % sem(l) % std_op       &
             , x     => ml_op  % sem(l) % metrics % x  &
             , Jd    => ml_op  % sem(l) % metrics % Jd &
             , sdc   => ml_op  % sdc(l)                &
             , level => ml_var % level(l)              )

      allocate(t(0:level%po_t))
      dt    = ONE / level % ne_t
      t(0:) = sdc % CollocationPoints(ZERO, dt)

      do n = 1, level % ne_t
        do m = 0, level % po_t
          associate(v => level % var(m,n) % val)

            do e = 1, sem%mesh%n_elem
              do k = 0, eop%po
              do j = 0, eop%po
              do i = 0, eop%po

                v(i,j,k,e,1) = sin(kappa * ( x(i,j,k,e,1)    &
                                           + x(i,j,k,e,2)    &
                                           + x(i,j,k,e,3)    &
                                           + t(m)         ))

                v(i,j,k,e,2) = v(i,j,k,e,1) &
                             * eop % w(i)   &
                             * eop % w(j)   &
                             * eop % w(k)   &
                             * Jd(i,j,k,e)
              end do
              end do
              end do
            end do

            call SetArray(v(:,:,:,:,3:), ZERO, multi=.true.)

          end associate
        end do

        t = t + dt
      end do

      deallocate(t)

    end associate
  end do

!!   ! coarse-to-fine interpolation ...............................................
!!
!!   allocate(delta(n_level-1), source = ZERO)
!!   allocate(delta_loc, source = delta)
!!
!!   do l = 1, n_level - 1
!!     associate( parent => ml_op  % sem(l)     % mesh &
!!              , child  => ml_op  % sem(l+1)   % mesh &
!!              , var_p  => ml_var % level(l  ) % val  &
!!              , var_c  => ml_var % level(l+1) % val  &
!!              , iop    => ml_op  % iop_cf_x(l)       )
!!
!!       call ParentToChildInterpolation_3D( parent, child, iop     &
!!                                         , v_p = var_p(:,:,:,:,6) &
!!                                         , v_c = var_c(:,:,:,:,8) )
!!
!!       do e = 1, child%n_elem
!!         var_c(:,:,:,e,9) = var_c(:,:,:,e,8) - var_c(:,:,:,e,6)
!!         delta_loc(l) = max(delta_loc(l), maxval(abs(var_c(:,:,:,e,9))))
!!       end do
!!
!!     end associate
!!   end do
!!
!!   call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)
!!
!!   if (rank == 0) then
!!     write(*,'(/,A)') 'coarse-to-fine interpolation error'
!!     do l = 1, n_level-1
!!       write(*,'(2X,5G0,ES10.3)') '|I v_',l,' - v_',l+1,'| =', delta(l)
!!     end do
!!   end if
!!
!!   ! fine-to-coarse interpolation ...............................................
!!
!!   do l = 1, n_level - 1
!!     associate( parent => ml_op  % sem(l)     % mesh &
!!              , child  => ml_op  % sem(l+1)   % mesh &
!!              , var_p  => ml_var % level(l  ) % val  &
!!              , var_c  => ml_var % level(l+1) % val  &
!!              , pop    => ml_op  % iop_fc_x(l+1)     )
!!
!!       call ChildToParentProjection_3D( child, parent, pop &
!!                                      , var_c(:,:,:,:,6)   &
!!                                      , var_p(:,:,:,:,10)  )
!!
!!       delta_loc(l) = 0
!!       do e = 1, parent%n_elem
!!         if (parent%element(e)%adaptation%refinement < 1000) cycle
!!         var_p(:,:,:,e,11) = var_p(:,:,:,e,10) - var_p(:,:,:,e,6)
!!         delta_loc(l) = max(delta_loc(l), maxval(abs(var_p(:,:,:,e,11))))
!!       end do
!!
!!     end associate
!!   end do
!!
!!   call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)
!!
!!   if (rank == 0) then
!!     write(*,'(/,A)') 'fine-to-coarse interpolation error'
!!     do l = 1, n_level-1
!!       write(*,'(2X,5G0,ES10.3)') '|I v_',l+1,' - v_',l,'| =', delta(l)
!!     end do
!!   end if
!!
!!   ! fine-to-coarse L2-projection ...............................................
!!
!!   do l = 1, n_level - 1
!!     associate( parent => ml_op  % sem(l)     % mesh &
!!              , child  => ml_op  % sem(l+1)   % mesh &
!!              , var_p  => ml_var % level(l  ) % val  &
!!              , var_c  => ml_var % level(l+1) % val  &
!!              , pop    => ml_op  % pop_fc_x(l+1)     )
!!
!!       call ChildToParentProjection_3D( child, parent, pop &
!!                                      , var_c(:,:,:,:,6)   &
!!                                      , var_p(:,:,:,:,12)  )
!!
!!       delta_loc(l) = 0
!!       do e = 1, parent%n_elem
!!         if (parent%element(e)%adaptation%refinement < 1000) cycle
!!         var_p(:,:,:,e,13) = var_p(:,:,:,e,12) - var_p(:,:,:,e,6)
!!         delta_loc(l) = max(delta_loc(l), maxval(abs(var_p(:,:,:,e,13))))
!!       end do
!!
!!     end associate
!!   end do
!!
!!   call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)
!!
!!   if (rank == 0) then
!!     write(*,'(/,A)') 'fine-to-coarse L²-projection error'
!!     do l = 1, n_level-1
!!       write(*,'(2X,5G0,ES10.3)') '|P v_',l+1,' - v_',l,'| =', delta(l)
!!     end do
!!   end if
!!
!!   ! fine-to-coarse restriction .................................................
!!
!!   do l = 1, n_level - 1
!!     associate( parent => ml_op  % sem(l)     % mesh &
!!              , child  => ml_op  % sem(l+1)   % mesh &
!!              , var_p  => ml_var % level(l  ) % val  &
!!              , var_c  => ml_var % level(l+1) % val  &
!!              , iop    => ml_op  % iop_cf_x(l)       )
!!
!!       call ChildToParentRestriction_3D( child, parent, iop &
!!                                       , var_c(:,:,:,:,7)   &
!!                                       , var_p(:,:,:,:,14)  )
!!
!!       delta_loc(l) = 0
!!       do e = 1, parent%n_elem
!!         if (parent%element(e)%adaptation%refinement < 1000) cycle
!!         v = abs(sum(var_p(:,:,:,e,14)) - sum(var_p(:,:,:,e,7)))
!!         delta_loc(l) = max(delta_loc(l), v)
!!       end do
!!
!!     end associate
!!   end do
!!
!!   call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)
!!
!!   if (rank == 0) then
!!     write(*,'(/,A)') 'fine-to-coarse restriction element integral error'
!!     do l = 1, n_level-1
!!       write(*,'(2X,5G0,ES10.3)') '|Σ R(Mv_',l+1,') - Σ Mv_',l,'| =', delta(l)
!!     end do
!!   end if
!!
!!   ! VTK export .................................................................
!!
!!   if (len_trim(plot_file) > 0) then
!!     call ml_var % ExportVTK(ml_op, trim(plot_file)//'_full', mode=1)
!!     call ml_var % ExportVTK(ml_op, trim(plot_file)//'_leaf', mode=3)
!!   end if

  ! finalization ...............................................................

!!   if ( allocated(delta_loc) ) deallocate(delta_loc)
!!   if ( allocated(delta)     ) deallocate(delta)

  call MPI_Finalize()

  !=============================================================================

end program ML_Spacetime_Functionality
