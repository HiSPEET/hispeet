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

!> summary:  Program for testing basic multilevel spacetime functionality
!> author:   Joerg Stiller
!> date:     2024/11/28
!===============================================================================

program ML_Spacetime_Functionality
  use Kind_Parameters
  use Constants
  use Logging_Levels
  use Execution_Control
  use Array_Assignments
  use XMPI

  use Import_GMSH__3D
  use Generic_Mesh__3D
  use Mesh__3D
  use Verify_Mesh__3D

  use Spacetime_Variable__3D

  use ML__Mesh__3D
  use ML__Spacetime_Operators__3D
  use ML__Spacetime_Variable__3D

  use CF__Spacetime_Interpolation__3D
  use FC__Spacetime_Projection__3D
  use FC__Spacetime_Restriction__3D

  implicit none

  ! variables ..................................................................

  character(len=*), parameter :: default_case = 'ml_spacetime_functionality'
  character(len= 80) :: test_case ! test case name
  character(len=100) :: case_file ! test case input file: trim(test_case).prm
  character(len=100) :: gmsh_file = '../../gmsh/cylinder_2d'

  integer, allocatable :: po_x(:)    ! polynomial orders in x
  integer, allocatable :: po_t(:)    ! polynomial orders in t
  integer, allocatable :: ne_t(:)    ! elements per slice in t

  character(len=2) :: nodes_x = 'L'  ! type of space nodes
  character(len=2) :: nodes_t = 'RR' ! type of time nodes {E,L,RR}
  integer          :: smooth  =  0   ! discontinuity smoothing {0,1,2}
  real(RNP)        :: dt_slab =  1

  namelist /control/ log_level, gmsh_file
  namelist /discretization/ po_x, po_t, ne_t, nodes_t, smooth, dt_slab

  type(GenericMesh_3D)          , save :: generic_mesh
  type(Mesh_3D)                 , save :: base_mesh
  type(ML_Mesh_Options_3D)      , save :: ml_mesh_opt
  type(ML_Mesh_3D)              , save :: ml_mesh
  type(ML_SpacetimeOperators_3D), save :: ml_op
  type(ML_SpacetimeVariable_3D) , save :: ml_var

  type(SpacetimeVariable_3D), save :: st_v
  type(SpacetimeVariable_3D), save :: st_r

  type(MPI_Comm) :: comm = MPI_COMM_WORLD

  character(len=11), save :: var_name(9)

  real(RNP), allocatable, save :: delta(:)
  real(RNP), allocatable, save :: delta_loc(:)
  real(RNP), save :: kappa = real(2 * PI, RNP)

  real(RNP), allocatable :: t(:), wt(:)
  real(RNP) :: dt
  logical :: exists, passed, all_passed
  integer :: n_level, n_proc, rank
  integer :: io, stat
  integer :: ne_max, ne_min, ne_tot

  integer :: e, i, j, k, l, m, n, nc

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  ! read parameters
  if (rank == 0) then
    write(*,'(/,A,/)') 'Testing basic multilevel spacetime functionality'

    call get_command_argument(1, test_case, status=stat)
    if (stat /= 0 .or. len_trim(test_case) == 0) then
      test_case = default_case
    end if
    case_file = trim(test_case) // '.prm'

    inquire(file=case_file, exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control)
      ml_mesh_opt = ML_Mesh_Options_3D(io, n_proc)
      allocate(po_x(ml_mesh_opt%l_top), source = -1)
      allocate(po_t, ne_t, source = po_x)
      read(io, nml = discretization)
      close(io)
    else
       call Error( 'ML_Spacetime_Functionality', &
                   'input file "' // trim(case_file) // '" not found' )
    end if
  end if

  ! globalize multilevel mesh options
  call ml_mesh_opt % Bcast(0, comm)

  if (rank > 0) then
    allocate(po_x(ml_mesh_opt%l_top), source = -1)
    allocate(po_t, ne_t, source = po_x)
  end if

  ! globalize remaining parameters
  call XMPI_Bcast(test_case , 0, comm)
  call XMPI_Bcast(gmsh_file , 0, comm)
  call XMPI_Bcast(po_x      , 0, comm)
  call XMPI_Bcast(po_t      , 0, comm)
  call XMPI_Bcast(ne_t      , 0, comm)
  call XMPI_Bcast(nodes_x   , 0, comm)
  call XMPI_Bcast(nodes_t   , 0, comm)
  call XMPI_Bcast(smooth    , 0, comm)
  call XMPI_Bcast(dt_slab   , 0, comm)

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
  var_name(3) = 'I_cf(v)'
  var_name(4) = 'I_cf(v) - v'
  var_name(5) = 'I_fc(v)'
  var_name(6) = 'I_fc(v) - v'
  var_name(7) = 'P_fc(v)'
  var_name(8) = 'P_fc(v) - v'
  var_name(9) = 'R_fc(Mv)'

  call ml_var % Init(ml_op, nc, var_name)

  do l = 1, n_level
    associate( sem   => ml_op  % sem(l)                &
             , eop   => ml_op  % sem(l) % std_op       &
             , x     => ml_op  % sem(l) % metrics % x  &
             , Jd    => ml_op  % sem(l) % metrics % Jd &
             , sdc   => ml_op  % sdc(l)                &
             , level => ml_var % level(l)              )

      allocate(t(0:level%po_t))
      allocate(wt, mold = t)
      dt     = dt_slab / level % ne_t
      t(0:)  = sdc % CollocationPoints(ZERO, dt)
      wt(0:) = sdc % CollocationWeights()

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
                             * Jd(i,j,k,e)  &
                             * dt * wt(m)
              end do
              end do
              end do
            end do

            call SetArray(v(:,:,:,:,3:), ZERO, multi=.true.)

          end associate
        end do

        t = t + dt
      end do

      deallocate(t, wt)

    end associate
  end do

  ! coarse-to-fine interpolation ...............................................

  allocate(delta(n_level-1), source = ZERO)
  allocate(delta_loc, source = delta)

  do l = 1, n_level - 1
    call ml_var % level(l  ) % GetSlice(st_v, first=1, last=1)
    call ml_var % level(l+1) % GetSlice(st_r, first=3, last=3)

    call CF_SpacetimeInterpolation_3D(ml_op, l, st_v, st_r)

    delta_loc(l) = 0
    do n = 1, ml_var % level(l+1) % ne_t
    do m = 0, ml_var % level(l+1) % po_t
      associate(v => ml_var % level(l+1) % var(m,n) % val)
        do e = 1, size(v,4)
          v(:,:,:,e,4) = v(:,:,:,e,3) - v(:,:,:,e,1) ! =  I_cf(v_c) - v_f
          delta_loc(l) = max(delta_loc(l), maxval(abs(v(:,:,:,e,4))))
        end do
      end associate
    end do
    end do
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
    call ml_var % level(l+1) % GetSlice(st_v, first=1, last=1)
    call ml_var % level(l  ) % GetSlice(st_r, first=5, last=5)

    call FC_SpacetimeProjection_3D(ml_op, 'I', l+1, st_v, st_r)

    delta_loc(l) = 0
    do n = 1, ml_var % level(l) % ne_t
    do m = 0, ml_var % level(l) % po_t
      associate( mesh => ml_var % level(l) % var(m,n) % mesh &
               , v    => ml_var % level(l) % var(m,n) % val  )
        do e = 1, mesh % n_elem_active
          if (mesh%element(e)%adaptation%refinement < 1000) cycle
          v(:,:,:,e,6) = v(:,:,:,e,5) - v(:,:,:,e,1) ! =  I_fc(v_f) - v_c
          delta_loc(l) = max(delta_loc(l), maxval(abs(v(:,:,:,e,6))))
        end do
      end associate
    end do
    end do
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
    call ml_var % level(l+1) % GetSlice(st_v, first=1, last=1)
    call ml_var % level(l  ) % GetSlice(st_r, first=7, last=7)

    call FC_SpacetimeProjection_3D(ml_op, 'P', l+1, st_v, st_r)

    delta_loc(l) = 0
    do n = 1, ml_var % level(l) % ne_t
    do m = 0, ml_var % level(l) % po_t
      associate( mesh => ml_var % level(l) % var(m,n) % mesh &
               , v    => ml_var % level(l) % var(m,n) % val  )
        do e = 1, mesh % n_elem_active
          if (mesh%element(e)%adaptation%refinement < 1000) cycle
          v(:,:,:,e,8) = v(:,:,:,e,7) - v(:,:,:,e,1) ! =  P_fc(v_f) - v_c
          delta_loc(l) = max(delta_loc(l), maxval(abs(v(:,:,:,e,8))))
        end do
      end associate
    end do
    end do
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
    call ml_var % level(l+1) % GetSlice(st_v, first=2, last=2)
    call ml_var % level(l  ) % GetSlice(st_r, first=9, last=9)

    call FC_SpacetimeRestriction_3D(ml_op, l+1, st_v, st_r)

    delta_loc(l) = 0
    do n = 1, ml_var % level(l) % ne_t
    do m = 0, ml_var % level(l) % po_t
      associate( mesh  => ml_var % level(l) % var(m,n) % mesh           &
               , Mv_c  => ml_var % level(l) % var(m,n) % val(:,:,:,:,2) &
               , RMv_f => ml_var % level(l) % var(m,n) % val(:,:,:,:,9) )

        do e = 1, mesh % n_elem_active
          if (mesh%element(e)%adaptation%refinement < 1000) cycle
          delta_loc(l) = max( delta_loc(l)                             &
                            , abs(sum(RMv_f(:,:,:,e) - Mv_c(:,:,:,e))) )
        end do

      end associate
    end do
    end do
  end do

  call XMPI_Reduce(delta_loc, delta, MPI_MAX, 0, comm)

  if (rank == 0) then
    write(*,'(/,A)') 'fine-to-coarse restriction element integral error'
    do l = 1, n_level-1
      write(*,'(2X,5G0,ES10.3)') '|Σ R(Mv_',l+1,') - Σ Mv_',l,'| =', delta(l)
    end do
  end if

  ! finalization ...............................................................

  !$omp master
  deallocate(delta_loc, delta)
  !$omp end master

  call MPI_Finalize()

  !=============================================================================

end program ML_Spacetime_Functionality
