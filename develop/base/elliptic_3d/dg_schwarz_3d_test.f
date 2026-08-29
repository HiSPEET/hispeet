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

!> summary:  Investigation of 3D DG Schwarz methods
!> author:   Joerg Stiller
!> date:     2025/10/01
!===============================================================================

program DG_Schwarz_3D_Test
  use Kind_Parameters
  use Constants
  use Execution_Control
  use XMPI

  use Mesh__3D
  use Generate_Regular_Mesh__3D
  use Spectral_Element_Mesh__3D
  use Boundary_Variable__3D
  use Element_Transfer_Buffer__3D
  use VTK__Export_Mesh_Data__3D

  use DG__Elliptic_Operator__3D
  use DG__Schwarz_Operator__3D
  use TPO__Schwarz__3D

  use Elliptic_Problem__3D
  use Elliptic_Problem__Knotty__3D

  implicit none

  !-----------------------------------------------------------------------------
  ! Declarations

  ! control parameters .........................................................

  ! input file (*.prm)
  character(len=*), parameter :: default_case = 'dg_schwarz_3d_test'
  character(len=80) :: case_name = ''
  character(len=80) :: case_file = ''

  logical :: export_vtk = .false.  ! switch for VTK export
  logical :: subdiv_vtk = .false.

  namelist/control_prm/ export_vtk, subdiv_vtk

  ! problem parameters .........................................................

  real(RNP) :: lambda       =  0     ! Helmholtz parameter
  real(RNP) :: nu_0         =  1     ! diffusivity mean value ν₀
  real(RNP) :: nu_1         =  0     ! diffusivity fluctuation amplitude ν₁
  real(RNP) :: d_nu         =  0     ! diffusivity fluctuation phase shift
  integer   :: k_nu         =  1     ! diffusivity fluctuation wave number
  integer   :: k_u          =  2     ! solution wave number
  integer   :: start_values =  1     ! 0/1/2: zero, exact, random

  namelist/problem_prm/ lambda, start_values

  character, allocatable :: bc(:)    ! boundary conditions

  ! discretization parameters ..................................................

  integer   :: po         = 7        ! polynomial degree of spectral elements
  real(RNP) :: penalty    = 5        ! penalty parameter > 1
  logical   :: regular    = .true.
  logical   :: structured = .true.

  namelist/discretization_prm/ po, penalty, regular, structured

  ! solution ...................................................................

  type(DG_SchwarzOptions_3D) :: schwarz_opt
  logical :: skip(27) = .false.     ! elements for which r = 0 is set

  namelist /solver_prm/ schwarz_opt, skip

  ! MPI ........................................................................

  type(MPI_Comm) :: comm             ! MPI communicator
  integer        :: rank             ! local MPI rank
  integer        :: n_proc           ! number of MPI processes

  ! mesh and variables .........................................................

  type(Mesh_3D)                :: mesh
  type(SpectralElementMesh_3D) :: sem

  ! discrete operators
  type(DG_EllipticOperator_3D) :: elliptic_op

  ! problem
  class(EllipticProblem_3D), allocatable :: problem

  ! space and names for variables
  real(RNP), allocatable, target :: var(:,:,:,:,:)
  character(len=80), allocatable :: var_names(:)

  ! variables
  real(RNP), pointer, contiguous :: s   (:,:,:,:) ! exact solution
  real(RNP), pointer, contiguous :: u   (:,:,:,:) ! approximate solution
  real(RNP), pointer, contiguous :: f   (:,:,:,:) ! RHS
  real(RNP), pointer, contiguous :: r   (:,:,:,:) ! residual
  real(RNP), pointer, contiguous :: e   (:,:,:,:) ! error
  real(RNP), pointer, contiguous :: rs_c(:,:,:,:) ! residual @ central subdomain
  real(RNP), pointer, contiguous :: us_c(:,:,:,:) ! result @ central subdomain
  real(RNP), pointer, contiguous :: es_c(:,:,:,:) ! error @ central subdomain

  ! diagonal mass matrix
  real(RNP), allocatable :: mm(:,:,:,:)

  ! Schwarz
  real(RDP), allocatable, save :: vs(:)        ! subdomain viscosities
  real(RDP), allocatable, save :: rs(:,:,:,:)  ! RHS of subsystems
  real(RDP), allocatable, save :: us(:,:,:,:)  ! solution of subsystems
  integer,   allocatable, save :: cfg(:,:)     ! subdomain configurations
  type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_r
  type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_us

  type(BoundaryVariable_3D), allocatable :: bv_u(:)

  ! auxiliary variables ........................................................

  logical :: exists
  integer :: io, stat
  integer :: n_bound, n_elem, n_var
  integer :: na, np, ne, ng, nl(3), no, ns, wp
  integer :: i, ic = 14

  !-----------------------------------------------------------------------------
  ! Initialization

  ! MPI ........................................................................

  call XMPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank)
  call MPI_Comm_size(comm, n_proc)

  if (n_proc > 1) then
    call Error('DG_Schwarz_Test_3D', 'must be run with only 1 process')
  end if

  ! control parameters .........................................................

  if (rank == 0) then

    write(*,'(/,A)') repeat('=',80)
    write(*,'(A)') 'Investigation of 3D DG Schwarz methods'
    write(*,*)

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) then
      case_name = default_case
    end if
    case_file = trim(case_name) // '.prm'

    inquire(file=trim(case_file), exist=exists)
    if (exists) then
      write(*,'(2X,A)') 'reading control parameters from ' // trim(case_file)
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      close(io)
    else
       call Warning( 'DG_Schwarz_3D_Test', 'using default parameters' )
    end if

  end if

  ! mesh generation ............................................................

  block

    real(RNP) :: xo(3) = -1
    real(RNP) :: dx(3) =  2 * THIRD
    logical   :: periodic(3) = .false.

    call GenerateRegularMesh(mesh, [1,1,1], [3,3,3], xo, dx, periodic, comm)

    mesh % regular    = regular
    mesh % structured = structured

    n_elem  = mesh % n_elem
    n_bound = mesh % n_bound

  end block

  ! problem ....................................................................

  allocate(bc(n_bound), source = 'D')

  ! read problem parameters
  if (rank == 0) then
    write(*,'(/,A)') 'initializing elliptic problem'
    open(newunit = io, file = case_file)
    read(io, nml = problem_prm)
    close(io)
  end if

  problem = EllipticProblem_Knotty_3D(lambda, nu_0, nu_1, d_nu, k_nu, k_u)

  ! spectral element mesh ......................................................

  ! read problem parameters
  if (rank == 0) then
    write(*,'(/,A)') 'initializing spectral element operators and solver'
    open(newunit = io, file = case_file)
    read(io, nml = discretization_prm)
    read(io, nml = solver_prm)
    schwarz_opt % wp = RDP ! enforce double precision
    close(io)
  end if

  sem = SpectralElementMesh_3D(mesh, po)

  ! variables ..................................................................

  n_var = 8

  allocate(var(0:po,0:po,0:po,1:n_elem,1:n_var), source = ZERO)

  s   (0:,0:,0:,1:) => var(:,:,:,:,1)
  u   (0:,0:,0:,1:) => var(:,:,:,:,2)
  f   (0:,0:,0:,1:) => var(:,:,:,:,3)
  r   (0:,0:,0:,1:) => var(:,:,:,:,4)
  e   (0:,0:,0:,1:) => var(:,:,:,:,5)
  rs_c(0:,0:,0:,1:) => var(:,:,:,:,6)
  us_c(0:,0:,0:,1:) => var(:,:,:,:,7)
  es_c(0:,0:,0:,1:) => var(:,:,:,:,8)

  var_names = [ 's   ', 'u   ', 'f   ', 'r   ', 'e   ', 'rs_c', 'us_c', 'es_c' ]

  allocate(mm(0:po,0:po,0:po,1:n_elem))
  call sem % Get_DG_DiagonalMassMatrix(mm)

  allocate(bv_u(n_bound))
  do i = 1, n_bound
    call bv_u(i) % Init(sem%mesh%boundary(i), po, nc = 1)
  end do

  ! solution and RHS ...........................................................

  associate(x => sem % metrics % x)

    ! exact solution, diffusivity and flux vector
    call problem % GetExactSolution(x, s)

    ! r = λ u - ∇·(ν ∇u)
    call problem % GetSource(x, f)

    ! project source
    f = mm * f

    ! extract boundary conditions
    do i = 1, n_bound
      call bv_u(i) % Extract(s)
    end do

    ! start values
    select case(start_values)
    case(0)
      u = 0
    case(1)
      u = s
    case default
      call random_number(u)
      u = 2*u - 1
    end select

  end associate

  ! operators ..................................................................

  elliptic_op = DG_EllipticOperator_3D(sem, schwarz_opt, penalty)

  ! Schwarz method .............................................................

  np = po + 1
  na = mesh % n_elem_active
  ne = mesh % n_elem
  ng = mesh % n_ghost
  no = elliptic_op % schwarz % no
  wp = elliptic_op % schwarz % wp
  ns = np + 2*no

  nl = no

  allocate(cfg(3,na), vs(na), rs(ns,ns,ns,ne), us(ns,ns,ns,ne))

  vs = real(nu_0, RDP)
  rs = 0
  us = 0

  call elliptic_op % schwarz % GetSubdomainConfigurations(mesh, bc, cfg)

  buf_r  = ElementTransferBuffer_3D(mesh, r , nl)
  buf_us = ElementTransferBuffer_3D(mesh, us, nl)

  !-----------------------------------------------------------------------------
  ! Schwarz method for the central element

  ! residual computation .......................................................

  call elliptic_op % Residual(bc, lambda, nu_0, f, bv_u, u, r)

  ! where requested set r = 0 ..................................................

  do i = 1, size(skip)
    if (skip(i)) r(:,:,:,i) = 0
  end do

  ! residual restriction .......................................................

  call elliptic_op % schwarz % RestrictResidual(mesh, buf_r, r, rs)

  ! subdomain corrections ......................................................

  call TPO_Schwarz( elliptic_op % schwarz % ops_dp % S  &
                  , elliptic_op % schwarz % ops_dp % V  &
                  , elliptic_op % schwarz % ops_dp % W  &
                  , elliptic_op % schwarz % ops_dp % g  &
                  , cfg(:,:na)                          &
                  , real(lambda, RDP)                   &
                  , vs(:na)                             &
                  , rs(:,:,:,:na)                       &
                  , us(:,:,:,:na)                       )

  call InjectCentralSubdomain(rs_c, rs)

  ! solution assembly ..........................................................

  do i = 1, size(skip)
    if (skip(i)) us(:,:,:,i) = 0
  end do

  call elliptic_op % schwarz % MergeCorrections(mesh, buf_us, us, u)

  us_c = s
  call InjectCentralSubdomain(us_c, us)

  block
    integer :: k
    print *
    i = 0
    do k = 0, po
      print '(2(A,I0),A,99(X,ES10.3))', 'u(',i,',:,',k,',ic) =', u(i,:,k,ic)
    end do
    print *
    print '(A,3I3)','cfg(ic) =',cfg(:,ic)
    print *
    print '(A,99(X,ES10.3))','W1 =',elliptic_op % schwarz % ops_dp % W(:,1)
    print *
  end block

  !-----------------------------------------------------------------------------
  ! Evaluation

  e    = u    - s
  es_c = us_c - s

  write(*,*)
  write(*,'(A,ES10.3)') 'r_l2   =', sqrt(sum(mm * r**2))
  write(*,'(A,ES10.3)') 'r_max  =', maxval(abs(r))
  write(*,*)
  write(*,'(A,ES10.3)') 'e_l2   =', sqrt(sum(mm * e**2))
  write(*,'(A,ES10.3)') 'e_max  =', maxval(abs(e))
  write(*,*)
  write(*,'(A,ES10.3)') 'ec_l2  =', sqrt(sum(mm(:,:,:,ic) * e(:,:,:,ic)**2))
  write(*,'(A,ES10.3)') 'ec_max =', maxval(abs(e(:,:,:,ic)))
  write(*,*)
  write(*,'(A,ES10.3)') 'es_l2  =', sqrt(sum(mm(:,:,:,ic) * es_c(:,:,:,ic)**2))
  write(*,'(A,ES10.3)') 'es_max =', maxval(abs(es_c(:,:,:,ic)))

  !-----------------------------------------------------------------------------
  ! Export

  if (export_vtk .and. mesh%part >= 0) then

    call VTK_ExportMeshData_3D( sem % metrics % x         &
                              , s       = var             &
                              , sname   = var_names       &
                              , file    = trim(case_name) &
                              , part    = mesh % part     &
                              , n_parts = mesh % n_parts  &
                              , subdiv  = subdiv_vtk      )
  end if

  !-----------------------------------------------------------------------------
  ! Finalization

  call MPI_Finalize()

contains

  subroutine InjectCentralSubdomain(a, as)
    real(RNP), intent(inout) :: a(0:,0:,0:,:)
    real(RDP), intent(in)    :: as(:,:,:,:)

    integer :: m0, m1, m2

    m0 = no + 1
    m1 = ns - no
    m2 = m1 + 1

    a(np-no:po  , np-no:po  , np-no:po  ,  1) = as( 1:no,  1:no,  1:no, ic)
    a(    0:po  , np-no:po  , np-no:po  ,  2) = as(m0:m1,  1:no,  1:no, ic)
    a(    0:no-1, np-no:po  , np-no:po  ,  3) = as(m2:ns,  1:no,  1:no, ic)

    a(np-no:po  ,     0:po  , np-no:po  ,  4) = as( 1:no, m0:m1,  1:no, ic)
    a(    0:po  ,     0:po  , np-no:po  ,  5) = as(m0:m1, m0:m1,  1:no, ic)
    a(    0:no-1,     0:po  , np-no:po  ,  6) = as(m2:ns, m0:m1,  1:no, ic)

    a(np-no:po  ,     0:no-1, np-no:po  ,  7) = as( 1:no, m2:ns,  1:no, ic)
    a(    0:po  ,     0:no-1, np-no:po  ,  8) = as(m0:m1, m2:ns,  1:no, ic)
    a(    0:no-1,     0:no-1, np-no:po  ,  9) = as(m2:ns, m2:ns,  1:no, ic)

    a(np-no:po  , np-no:po  ,     0:po  , 10) = as( 1:no,  1:no, m0:m1, ic)
    a(    0:po  , np-no:po  ,     0:po  , 11) = as(m0:m1,  1:no, m0:m1, ic)
    a(    0:no-1, np-no:po  ,     0:po  , 12) = as(m2:ns,  1:no, m0:m1, ic)

    a(np-no:po  ,     0:po  ,     0:po  , 13) = as( 1:no, m0:m1, m0:m1, ic)
    a(    0:po  ,     0:po  ,     0:po  , 14) = as(m0:m1, m0:m1, m0:m1, ic)
    a(    0:no-1,     0:po  ,     0:po  , 15) = as(m2:ns, m0:m1, m0:m1, ic)

    a(np-no:po  ,     0:no-1,     0:po  , 16) = as( 1:no, m2:ns, m0:m1, ic)
    a(    0:po  ,     0:no-1,     0:po  , 17) = as(m0:m1, m2:ns, m0:m1, ic)
    a(    0:no-1,     0:no-1,     0:po  , 18) = as(m2:ns, m2:ns, m0:m1, ic)

    a(np-no:po  , np-no:po  ,     0:no-1, 19) = as( 1:no,  1:no, m2:ns, ic)
    a(    0:po  , np-no:po  ,     0:no-1, 20) = as(m0:m1,  1:no, m2:ns, ic)
    a(    0:no-1, np-no:po  ,     0:no-1, 21) = as(m2:ns,  1:no, m2:ns, ic)

    a(np-no:po  ,     0:po  ,     0:no-1, 22) = as( 1:no, m0:m1, m2:ns, ic)
    a(    0:po  ,     0:po  ,     0:no-1, 23) = as(m0:m1, m0:m1, m2:ns, ic)
    a(    0:no-1,     0:po  ,     0:no-1, 24) = as(m2:ns, m0:m1, m2:ns, ic)

    a(np-no:po  ,     0:no-1,     0:no-1, 25) = as( 1:no, m2:ns, m2:ns, ic)
    a(    0:po  ,     0:no-1,     0:no-1, 26) = as(m0:m1, m2:ns, m2:ns, ic)
    a(    0:no-1,     0:no-1,     0:no-1, 27) = as(m2:ns, m2:ns, m2:ns, ic)

  end subroutine InjectCentralSubdomain

  !=============================================================================

end program DG_Schwarz_3D_Test
