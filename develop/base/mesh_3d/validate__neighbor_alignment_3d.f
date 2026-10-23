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

!> summary:  Geometric validation of neighbor orientation and alignment
!> author:   Joerg Stiller
!> date:     2026/10/18
!>
!> Usage: validate__neighbor_alignment_3d [case]
!>
!> Reads the namelist `control_prm` (po, lx, xo, max_print, n_config,
!> periodic) from `case.prm` (default `validate__neighbor_alignment_3d.prm`),
!> if present.
!>
!> For each periodicity configuration `periodic(:,c)`, c = 1..n_config, and
!> all 4³ rotations of the center element of a 3x3x3 Cartesian mesh
!> (`CreateOneRotated`), the mesh is imported with `ImportGenericMesh`, such
!> that the orientations are determined by the production code.
!> All checks are based on the physical coordinates of the element points
!> and do not use `TransformIndex` or the orientation tables:
!>
!>   * face neighbors
!>       - orientation: neighbor directions d1, d2 match ξ, η of the element
!>       - AlignFromNeighborFace, AlignWithNeighborFace (real and integer)
!>       - face and vertex IDs of the coupled faces
!>   * edge neighbors
!>       - orientation as above
!>       - AlignFromNeighborEdge, AlignWithNeighborEdge (real and integer)
!>       - edge and vertex IDs of the coupled edges
!>   * vertex neighbors
!>       - orientation as above (only regular vertices, others are skipped)
!>       - coordinates and IDs of the coupled vertices
!>   * mesh faces (element face normal and rotation)
!>       - AlignWithMesh: all element faces referring to the same mesh face
!>         yield the same mesh face data (real and integer)
!>       - AlignFromMesh reverses AlignWithMesh
!>   * mesh edges (element edge orientation)
!>       - as for mesh faces, with all element edges of a mesh edge
!>
!> Real data are the point coordinates, integer data encode the point index.
!> Coordinates are compared modulo the domain extensions in periodic
!> directions.
!> Output per configuration: the number of wrong cases per check (must be 0),
!> the distinct pairs (element component, orientation) that failed and the
!> coverage of component/orientation pairs. For mesh faces the orientation
!> is given as `nr` = n+r or n-r (normal ±1, rotation r), for mesh edges as
!> +1 or -1. The final line `RES total wrong` sums over all configurations.
!> With `mutate = T`, face 1 and edge 1 of the center element get a wrong
!> rotation and orientation as negative control. Runs on rank 0 only.
!===============================================================================

program Validate__Neighbor_Alignment_3D

  use Kind_Parameters
  use Constants
  use XMPI
  use Generic_Mesh__3D
  use Mesh__3D
  use Mesh_Element__3D
  use Mesh_Element_Indexing__3D

  implicit none

  ! control parameters
  character(len=*), parameter :: DEFAULT_CASE = &
      'validate__neighbor_alignment_3d'
  character(len=80) :: case_name = '', case_file = ''

  integer   :: po        = 3                 ! polynomial order
  real(RNP) :: lx(3)     = [ 1.0, 1.3, 1.7 ] ! domain extensions
  real(RNP) :: xo(3)     = [ 0.1, 0.2, 0.3 ] ! domain origin
  integer   :: max_print = 20                ! max printed failures
  integer   :: n_config  = 5                 ! number of configurations
  logical   :: periodic(3,8) = reshape( &    ! periodic directions
      [ .false., .false., .false.       &
      , .true. , .false., .false.       &
      , .false., .true. , .false.       &
      , .false., .false., .true.        &
      , .true. , .true. , .true.        &
      , .false., .false., .false.       &
      , .false., .false., .false.       &
      , .false., .false., .false.       ], [3, 8])

  logical   :: mutate    = .false.                ! negative control

  namelist/control_prm/ po, lx, xo, max_print, n_config, periodic, mutate

  ! neighbor orientations
  integer(IXS), parameter :: ORIENTATION(24) =                            &
      [ 12_IXS, 13_IXS, 15_IXS, 16_IXS, 21_IXS, 23_IXS, 24_IXS, 26_IXS    &
      , 31_IXS, 32_IXS, 34_IXS, 35_IXS, 42_IXS, 43_IXS, 45_IXS, 46_IXS    &
      , 51_IXS, 53_IXS, 54_IXS, 56_IXS, 61_IXS, 62_IXS, 64_IXS, 65_IXS ]

  ! checks
  integer, parameter :: N_CHECK = 21
  character(len=*), parameter :: CHECK_NAME(N_CHECK) =  &
      [ 'face   orientation          '                  &
      , 'face   AlignFrom RNP        '                  &
      , 'face   AlignWith RNP        '                  &
      , 'face   AlignFrom IDK        '                  &
      , 'face   AlignWith IDK        '                  &
      , 'face   component IDs        '                  &
      , 'edge   orientation          '                  &
      , 'edge   AlignFrom RNP        '                  &
      , 'edge   AlignWith RNP        '                  &
      , 'edge   AlignFrom IDK        '                  &
      , 'edge   AlignWith IDK        '                  &
      , 'edge   component IDs        '                  &
      , 'vertex orientation          '                  &
      , 'vertex coordinates          '                  &
      , 'vertex IDs                  '                  &
      , 'mesh face AlignWith RNP     '                  &
      , 'mesh face AlignFrom RNP     '                  &
      , 'mesh face AlignWith IDK     '                  &
      , 'mesh edge AlignWith RNP     '                  &
      , 'mesh edge AlignFrom RNP     '                  &
      , 'mesh edge AlignWith IDK     ' ]

  integer, parameter :: C_FACE_ORIENT = 1, C_FACE_FROM = 2, C_FACE_WITH = 3
  integer, parameter :: C_FACE_FROM_I = 4, C_FACE_WITH_I = 5
  integer, parameter :: C_FACE_VID = 6
  integer, parameter :: C_EDGE_ORIENT = 7, C_EDGE_FROM = 8, C_EDGE_WITH = 9
  integer, parameter :: C_EDGE_FROM_I = 10, C_EDGE_WITH_I = 11
  integer, parameter :: C_EDGE_VID = 12
  integer, parameter :: C_VERT_ORIENT = 13, C_VERT_X = 14, C_VERT_VID = 15
  integer, parameter :: C_MF_WITH = 16, C_MF_FROM = 17, C_MF_WITH_I = 18
  integer, parameter :: C_ME_WITH = 19, C_ME_FROM = 20, C_ME_WITH_I = 21

  ! results
  ! results, reset per configuration
  integer :: n_case(N_CHECK)         ! checked cases
  integer :: n_bad(N_CHECK)          ! failed cases
  logical :: bad_pair(N_CHECK,12,0:24) ! failed (component, o)
  logical :: seen_face(6,24)
  logical :: seen_edge(12,24)
  logical :: seen_vert(8,24)
  logical :: seen_mface(6,8)         ! (face, normal/rotation)
  logical :: seen_medge(12,2)        ! (edge, orientation)
  integer :: n_invalid               ! orientation codes not in list
  integer :: n_skip_vert             ! vertices without orientation
  integer :: n_print
  integer :: n_wrong_total = 0       ! failed cases of all configurations

  ! current periodicity
  logical :: per(3)

  ! first contribution to mesh faces and edges
  real(RNP), allocatable :: x_mface(:,:,:,:), x_medge(:,:,:)
  logical,   allocatable :: set_mface(:), set_medge(:)

  ! mesh
  type(GenericMesh_3D)         :: generic_mesh
  type(Mesh_3D), allocatable   :: mesh
  real(RNP),     allocatable   :: x(:,:,:,:,:)

  integer   :: rank, io, stat, r1, r2, r3, c
  logical   :: exists
  real(RNP) :: tol

  ! initialization .............................................................

  call XMPI_Init()
  call MPI_Comm_rank(MPI_COMM_WORLD, rank)

  if (rank == 0) then

    call get_command_argument(1, case_name, status=stat)
    if (stat /= 0 .or. len_trim(case_name) == 0) case_name = DEFAULT_CASE
    case_file = trim(case_name) // '.prm'
    inquire(file=trim(case_file), exist=exists)
    if (exists) then
      open(newunit = io, file = case_file)
      read(io, nml = control_prm)
      close(io)
    end if

    tol = 1e-12_RNP * maxval(abs(xo) + lx)

    ! all configurations and rotations of the center element .................

    do c = 1, n_config
      per = periodic(:,c)
      call ResetResults()
      do r3 = 0, 3
      do r2 = 0, 3
      do r1 = 0, 3
        call generic_mesh % CreateOneRotated( xo, lx, po, [r1, r2, r3]    &
                                            , periodic = per )
        allocate(mesh)
        call mesh % ImportGenericMesh(generic_mesh, comm = MPI_COMM_SELF)
        call mesh % GetPoints(po, 'L', x)
        if (mutate) then ! corrupt face 1 and edge 1 of the center element
          associate( face => mesh % element(14) % face(1)                &
                   , edge => mesh % element(14) % edge(1) )
            face % rotation = mod(face % rotation + 1, 4)
            edge % orientation = -edge % orientation
          end associate
        end if
        call CheckMesh([r1, r2, r3])
        deallocate(mesh, x)
      end do
      end do
      end do
      call PrintResults()
    end do

    write(*,'(/,A,I0)') 'RES total wrong ', n_wrong_total

  end if

  call MPI_Finalize()

contains

  !-----------------------------------------------------------------------------
  !> Checks all neighbor relations of the current mesh

  subroutine CheckMesh(rot)
    integer, intent(in) :: rot(3) !< rotation of the center element

    integer   :: e, en, f, fc, ed, ec, v, vc, n, n1, n2
    integer(IXS) :: o
    real(RNP) :: ae(3,3), an(3,3)

    ! mesh faces and edges
    allocate(x_mface(0:po,0:po,3,mesh % n_face), set_mface(mesh % n_face))
    allocate(x_medge(0:po,3,mesh % n_edge), set_medge(mesh % n_edge))
    set_mface = .false.
    set_medge = .false.
    do e = 1, mesh % n_elem
      do f = 1, 6
        call CheckMeshFace(rot, e, f)
      end do
      do ed = 1, 12
        call CheckMeshEdge(rot, e, ed)
      end do
    end do
    deallocate(x_mface, set_mface, x_medge, set_medge)

    ! neighbors
    do e = 1, mesh % n_elem
      associate(element => mesh % element(e))

        ae = Axes(e)

        ! face neighbors .......................................................

        do f = 1, 6
          if (element % face(f) % n_neighbor /= 1) cycle
          n  = element % face(f) % i_neighbor
          en = element % neighbor(n) % id
          fc = ElementFaceID(element % neighbor(n) % component)
          o  = element % neighbor(n) % orientation
          an = Axes(en)
          call CheckFace(rot, e, f, n, en, fc, o, ae, an)
        end do

        ! edge neighbors .......................................................

        do ed = 1, 12
          if (element % edge(ed) % n_neighbor < 1) cycle
          n1 = element % edge(ed) % i_neighbor
          n2 = element % edge(ed) % n_neighbor + n1 - 1
          do n = n1, n2
            en = element % neighbor(n) % id
            ec = ElementEdgeID(element % neighbor(n) % component)
            o  = element % neighbor(n) % orientation
            an = Axes(en)
            call CheckEdge(rot, e, ed, n, en, ec, o, ae, an)
          end do
        end do

        ! vertex neighbors .....................................................

        do v = 1, 8
          if (element % vertex(v) % n_neighbor < 1) cycle
          n1 = element % vertex(v) % i_neighbor
          n2 = element % vertex(v) % n_neighbor + n1 - 1
          do n = n1, n2
            en = element % neighbor(n) % id
            vc = ElementVertexID(element % neighbor(n) % component)
            o  = element % neighbor(n) % orientation
            an = Axes(en)
            call CheckVertex(rot, e, v, en, vc, o, ae, an)
          end do
        end do

      end associate
    end do

  end subroutine CheckMesh

  !-----------------------------------------------------------------------------
  !> Checks face f of element e coupled to face fc of neighbor en

  subroutine CheckFace(rot, e, f, n, en, fc, o, ae, an)
    integer,      intent(in) :: rot(3), e, f, n, en, fc
    integer(IXS), intent(in) :: o
    real(RNP),    intent(in) :: ae(3,3), an(3,3)

    real(RNP) :: xe(0:po,0:po,3), xn(0:po,0:po,3), xa(0:po,0:po,3)
    integer   :: ie(0:po,0:po), in_(0:po,0:po), ia(0:po,0:po)
    integer   :: d, i, j, p, q, k
    logical   :: ok

    associate(element => mesh % element(e), neighbor => mesh % element(en))

      k = OrientationIndex(o)
      if (k > 0) seen_face(f,k) = .true.

      xe = FacePoints(e, f)
      xn = FacePoints(en, fc)

      ! orientation
      call Record(C_FACE_ORIENT, f, k, OrientationMatches(o, ae, an)      &
                 , 'face', rot, e, f, en, fc, o)

      ! AlignFromNeighborFace, real: neighbor coordinates to element
      do d = 1, 3
        call element % AlignFromNeighborFace(f, n, xn(:,:,d), xa(:,:,d))
      end do
      call Record(C_FACE_FROM, f, k, Gap(size(xa)/3, xa, xe) <= tol       &
                 , 'face', rot, e, f, en, fc, o)

      ! AlignWithNeighborFace, real: element coordinates to neighbor
      do d = 1, 3
        call element % AlignWithNeighborFace(f, n, xe(:,:,d), xa(:,:,d))
      end do
      call Record(C_FACE_WITH, f, k, Gap(size(xa)/3, xa, xn) <= tol       &
                 , 'face', rot, e, f, en, fc, o)

      ! AlignFromNeighborFace, integer: neighbor point index to element
      do q = 0, po
      do p = 0, po
        in_(p,q) = 1000 * p + q
      end do
      end do
      call element % AlignFromNeighborFace(f, n, in_, ia)
      ok = .true.
      do j = 0, po
      do i = 0, po
        p = ia(i,j) / 1000
        q = ia(i,j) - 1000 * p
        ok = ok .and. Inside(p, q) .and. SamePoint(xn(p,q,:), xe(i,j,:))
      end do
      end do
      call Record(C_FACE_FROM_I, f, k, ok, 'face', rot, e, f, en, fc, o)

      ! AlignWithNeighborFace, integer: element point index to neighbor
      do j = 0, po
      do i = 0, po
        ie(i,j) = 1000 * i + j
      end do
      end do
      call element % AlignWithNeighborFace(f, n, ie, ia)
      ok = .true.
      do q = 0, po
      do p = 0, po
        i = ia(p,q) / 1000
        j = ia(p,q) - 1000 * i
        ok = ok .and. Inside(i, j) .and. SamePoint(xe(i,j,:), xn(p,q,:))
      end do
      end do
      call Record(C_FACE_WITH_I, f, k, ok, 'face', rot, e, f, en, fc, o)

      ! face and vertex IDs of the coupled faces
      ok = element % face(f) % id == neighbor % face(fc) % id .and.       &
           SameSet( element  % vertex(V_FACE(:,f )) % id                  &
                  , neighbor % vertex(V_FACE(:,fc)) % id )
      call Record(C_FACE_VID, f, k, ok, 'face', rot, e, f, en, fc, o)

    end associate

  end subroutine CheckFace

  !-----------------------------------------------------------------------------
  !> Checks edge ed of element e coupled to edge ec of neighbor en

  subroutine CheckEdge(rot, e, ed, n, en, ec, o, ae, an)
    integer,      intent(in) :: rot(3), e, ed, n, en, ec
    integer(IXS), intent(in) :: o
    real(RNP),    intent(in) :: ae(3,3), an(3,3)

    real(RNP) :: xe(0:po,3), xn(0:po,3), xa(0:po,3)
    integer   :: ie(0:po), in_(0:po), ia(0:po)
    integer   :: d, i, p, k
    logical   :: ok

    associate(element => mesh % element(e), neighbor => mesh % element(en))

      k = OrientationIndex(o)
      if (k > 0) seen_edge(ed,k) = .true.

      xe = EdgePoints(e, ed)
      xn = EdgePoints(en, ec)

      call Record(C_EDGE_ORIENT, ed, k, OrientationMatches(o, ae, an)     &
                 , 'edge', rot, e, ed, en, ec, o)

      do d = 1, 3
        call element % AlignFromNeighborEdge(ed, n, xn(:,d), xa(:,d))
      end do
      call Record(C_EDGE_FROM, ed, k, Gap(size(xa)/3, xa, xe) <= tol       &
                 , 'edge', rot, e, ed, en, ec, o)

      do d = 1, 3
        call element % AlignWithNeighborEdge(ed, n, xe(:,d), xa(:,d))
      end do
      call Record(C_EDGE_WITH, ed, k, Gap(size(xa)/3, xa, xn) <= tol       &
                 , 'edge', rot, e, ed, en, ec, o)

      in_ = [(p, p = 0, po)]
      call element % AlignFromNeighborEdge(ed, n, in_, ia)
      ok = .true.
      do i = 0, po
        ok = ok .and. Inside(ia(i), 0) .and. SamePoint(xn(ia(i),:), xe(i,:))
      end do
      call Record(C_EDGE_FROM_I, ed, k, ok, 'edge', rot, e, ed, en, ec, o)

      ie = [(i, i = 0, po)]
      call element % AlignWithNeighborEdge(ed, n, ie, ia)
      ok = .true.
      do p = 0, po
        ok = ok .and. Inside(ia(p), 0) .and. SamePoint(xe(ia(p),:), xn(p,:))
      end do
      call Record(C_EDGE_WITH_I, ed, k, ok, 'edge', rot, e, ed, en, ec, o)

      ok = element % edge(ed) % id == neighbor % edge(ec) % id .and.  &
           SameSet( element  % vertex(V_EDGE(:,ed)) % id              &
                  , neighbor % vertex(V_EDGE(:,ec)) % id )
      call Record(C_EDGE_VID, ed, k, ok, 'edge', rot, e, ed, en, ec, o)

    end associate

  end subroutine CheckEdge

  !-----------------------------------------------------------------------------
  !> Checks vertex v of element e coupled to vertex vc of neighbor en

  subroutine CheckVertex(rot, e, v, en, vc, o, ae, an)
    integer,      intent(in) :: rot(3), e, v, en, vc
    integer(IXS), intent(in) :: o
    real(RNP),    intent(in) :: ae(3,3), an(3,3)

    integer :: k
    logical :: ok

    associate(element => mesh % element(e), neighbor => mesh % element(en))

      if (o > 0) then
        k = OrientationIndex(o)
        if (k > 0) seen_vert(v,k) = .true.
        call Record(C_VERT_ORIENT, v, k, OrientationMatches(o, ae, an)  &
                   , 'vertex', rot, e, v, en, vc, o)
      else
        k = 0
        n_skip_vert = n_skip_vert + 1
      end if

      ok = SamePoint(VertexPoint(e, v), VertexPoint(en, vc))
      call Record(C_VERT_X, v, max(k,0), ok, 'vertex', rot, e, v, en, vc, o)

      ok = element % vertex(v) % id == neighbor % vertex(vc) % id
      call Record(C_VERT_VID, v, max(k,0), ok, 'vertex', rot, e, v, en, vc, o)

    end associate

  end subroutine CheckVertex

  !-----------------------------------------------------------------------------
  !> Checks the alignment of face f of element e with its mesh face

  subroutine CheckMeshFace(rot, e, f)
    integer, intent(in) :: rot(3), e, f

    real(RNP) :: xe(0:po,0:po,3), xm(0:po,0:po,3), xa(0:po,0:po,3)
    integer   :: ie(0:po,0:po), im(0:po,0:po)
    integer   :: d, i, j, p, q, k, m
    integer(IXS) :: o
    logical   :: ok

    associate(face => mesh % element(e) % face(f))

      m = face % id
      k = 1 + face % rotation
      if (face % normal /= 1) k = k + 4
      o = 10 * merge(1, 2, face % normal == 1) + face % rotation
      seen_mface(f,k) = .true.

      xe = FacePoints(e, f)

      ! real data
      do d = 1, 3
        call face % AlignWithMesh(xe(:,:,d), xm(:,:,d))
        call face % AlignFromMesh(xm(:,:,d), xa(:,:,d))
      end do
      call Record(C_MF_FROM, f, k, Gap(size(xa)/3, xa, xe) <= tol  &
                 , 'face', rot, e, f, 0, m, o)

      ! integer data: element point index
      do j = 0, po
      do i = 0, po
        ie(i,j) = 1000 * i + j
      end do
      end do
      call face % AlignWithMesh(ie, im)

      if (.not. set_mface(m)) then
        set_mface(m) = .true.
        x_mface(:,:,:,m) = xm
      else
        call Record( C_MF_WITH, f, k                               &
                   , Gap(size(xm)/3, xm, x_mface(:,:,:,m)) <= tol  &
                   , 'face', rot, e, f, 0, m, o )
        ok = .true.
        do q = 0, po
        do p = 0, po
          i = im(p,q) / 1000
          j = im(p,q) - 1000 * i
          ok = ok .and. Inside(i, j) .and.          &
               SamePoint(xe(i,j,:), x_mface(p,q,:,m))
        end do
        end do
        call Record(C_MF_WITH_I, f, k, ok, 'face', rot, e, f, 0, m, o)
      end if

    end associate

  end subroutine CheckMeshFace

  !-----------------------------------------------------------------------------
  !> Checks the alignment of edge ed of element e with its mesh edge

  subroutine CheckMeshEdge(rot, e, ed)
    integer, intent(in) :: rot(3), e, ed

    real(RNP) :: xe(0:po,3), xm(0:po,3), xa(0:po,3)
    integer   :: ie(0:po), im(0:po)
    integer   :: d, i, p, k, m
    integer(IXS) :: o
    logical   :: ok

    associate(edge => mesh % element(e) % edge(ed))

      m = edge % id
      o = edge % orientation
      k = merge(1, 2, o == 1)
      seen_medge(ed,k) = .true.

      xe = EdgePoints(e, ed)

      do d = 1, 3
        call edge % AlignWithMesh(xe(:,d), xm(:,d))
        call edge % AlignFromMesh(xm(:,d), xa(:,d))
      end do
      call Record(C_ME_FROM, ed, k, Gap(size(xa)/3, xa, xe) <= tol  &
                 , 'edge', rot, e, ed, 0, m, o)

      ie = [(i, i = 0, po)]
      call edge % AlignWithMesh(ie, im)

      if (.not. set_medge(m)) then
        set_medge(m) = .true.
        x_medge(:,:,m) = xm
      else
        call Record( C_ME_WITH, ed, k                            &
                   , Gap(size(xm)/3, xm, x_medge(:,:,m)) <= tol  &
                   , 'edge', rot, e, ed, 0, m, o )
        ok = .true.
        do p = 0, po
          ok = ok .and. Inside(im(p), 0) .and.      &
               SamePoint(xe(im(p),:), x_medge(p,:,m))
        end do
        call Record(C_ME_WITH_I, ed, k, ok, 'edge', rot, e, ed, 0, m, o)
      end if

    end associate

  end subroutine CheckMeshEdge

  !-----------------------------------------------------------------------------
  !> Resets the results

  subroutine ResetResults()

    n_case      = 0
    n_bad       = 0
    bad_pair    = .false.
    seen_face   = .false.
    seen_edge   = .false.
    seen_vert   = .false.
    seen_mface  = .false.
    seen_medge  = .false.
    n_invalid   = 0
    n_skip_vert = 0
    n_print     = 0

  end subroutine ResetResults

  !-----------------------------------------------------------------------------
  !> Records the result of a check

  subroutine Record(c, comp, k, ok, kind, rot, e, ce, en, cn, o)
    integer,          intent(in) :: c        !< check
    integer,          intent(in) :: comp     !< element component (index)
    integer,          intent(in) :: k        !< orientation index, 0 if none
    logical,          intent(in) :: ok       !< check passed
    character(len=*), intent(in) :: kind     !< component kind
    integer,          intent(in) :: rot(3)   !< rotation of center element
    integer,          intent(in) :: e, ce    !< element and its component
    integer,          intent(in) :: en, cn   !< neighbor and its component,
                                             !< 0 and mesh component ID
    integer(IXS),     intent(in) :: o        !< orientation

    n_case(c) = n_case(c) + 1
    if (ok) return

    n_bad(c) = n_bad(c) + 1
    bad_pair(c,comp,k) = .true.

    if (n_print < max_print) then
      n_print = n_print + 1
      if (en > 0) then
        write(*,'(A,A,A,3I1,A,I2,1X,A,I2,A,I2,1X,A,I2,A,I0)')           &
            'FAIL ', trim(CHECK_NAME(c)), '  rot=', rot                 &
            , '  element ', e, kind, ce, ' -> neighbor ', en, kind, cn  &
            , '  o=', o
      else
        write(*,'(A,A,A,3I1,A,I2,1X,A,I2,A,A,I0,A,I0)')                 &
            'FAIL ', trim(CHECK_NAME(c)), '  rot=', rot                 &
            , '  element ', e, kind, ce, ' -> mesh ', kind, cn          &
            , '  o=', o
      end if
    end if

  end subroutine Record

  !-----------------------------------------------------------------------------
  !> Prints the summary

  subroutine PrintResults()

    integer :: c, comp, k, n_comp, n_miss
    character(len=400) :: line
    character(len=4)   :: label

    n_wrong_total = n_wrong_total + sum(n_bad)

    write(*,'(/,A,3(1X,L1))') '### periodic =', per
    write(*,'(A,I0)') 'po = ', po
    write(*,'(A,I0)')   'invalid orientation codes:      ', n_invalid
    write(*,'(A,I0)')   'vertices without orientation:   ', n_skip_vert

    n_miss = count(.not. seen_face)
    write(*,'(A,I0,A)') 'face   (face, o) pairs seen:    '  &
                       , 144 - n_miss, ' / 144'
    n_miss = count(.not. seen_edge)
    write(*,'(A,I0,A)') 'edge   (edge, o) pairs seen:    '  &
                       , 288 - n_miss, ' / 288'
    n_miss = count(.not. seen_vert)
    write(*,'(A,I0,A)') 'vertex (vertex, o) pairs seen:  '  &
                       , 192 - n_miss, ' / 192'
    n_miss = count(.not. seen_mface)
    write(*,'(A,I0,A)') 'mesh face (face, nr) seen:      '  &
                       , 48 - n_miss, ' / 48'
    n_miss = count(.not. seen_medge)
    write(*,'(A,I0,A)') 'mesh edge (edge, ±1) seen:      '  &
                       , 24 - n_miss, ' / 24'

    write(*,'(/,A)') 'RES check                        cases      wrong'
    do c = 1, N_CHECK
      write(*,'(A,A,2I11)') 'RES ', CHECK_NAME(c), n_case(c), n_bad(c)
    end do

    do c = 1, N_CHECK
      if (n_bad(c) == 0) cycle
      select case(c)
      case(C_FACE_ORIENT:C_FACE_VID)
        n_comp = 6
      case(C_EDGE_ORIENT:C_EDGE_VID, C_ME_WITH:C_ME_WITH_I)
        n_comp = 12
      case(C_MF_WITH:C_MF_WITH_I)
        n_comp = 6
      case default
        n_comp = 8
      end select
      write(*,'(/,2A)') 'wrong (component: orientations) for '  &
                       , trim(CHECK_NAME(c))
      do comp = 1, n_comp
        if (.not. any(bad_pair(c,comp,:))) cycle
        write(line,'(I2,A)') comp, ':'
        do k = 0, 24
          if (.not. bad_pair(c,comp,k)) cycle
          select case(c)
          case(C_MF_WITH:C_MF_WITH_I)
            write(label,'(A,I1)') merge('n+r', 'n-r', k <= 4), mod(k-1, 4)
          case(C_ME_WITH:C_ME_WITH_I)
            label = merge('+1  ', '-1  ', k == 1)
          case default
            if (k == 0) then
              label = 'none'
            else
              write(label,'(I2)') ORIENTATION(k)
            end if
          end select
          line = trim(line) // ' ' // trim(label)
        end do
        write(*,'(2X,A)') trim(line)
      end do
    end do

  end subroutine PrintResults

  !-----------------------------------------------------------------------------
  !> Index of orientation o in ORIENTATION, 0 if invalid

  integer function OrientationIndex(o) result(k)
    integer(IXS), intent(in) :: o

    do k = 1, 24
      if (ORIENTATION(k) == o) return
    end do
    k = 0
    n_invalid = n_invalid + 1

  end function OrientationIndex

  !-----------------------------------------------------------------------------
  !> Physical unit vectors of the element directions ξ, η, ζ

  function Axes(e) result(a)
    integer, intent(in) :: e
    real(RNP) :: a(3,3)

    integer :: d

    a(:,1) = x(po,0,0,e,:) - x(0,0,0,e,:)
    a(:,2) = x(0,po,0,e,:) - x(0,0,0,e,:)
    a(:,3) = x(0,0,po,e,:) - x(0,0,0,e,:)
    do d = 1, 3
      a(:,d) = a(:,d) / norm2(a(:,d))
    end do

  end function Axes

  !-----------------------------------------------------------------------------
  !> True if neighbor directions d1, d2 of o = 10 d1 + d2 match ξ, η

  logical function OrientationMatches(o, ae, an) result(ok)
    integer(IXS), intent(in) :: o
    real(RNP),    intent(in) :: ae(3,3) !< element directions
    real(RNP),    intent(in) :: an(3,3) !< neighbor directions

    integer :: d1, d2

    d1 = o / 10
    d2 = o - 10 * d1
    ok = d1 >= 1 .and. d1 <= 6 .and. d2 >= 1 .and. d2 <= 6
    if (.not. ok) return
    ok = norm2(Direction(d1, an) - ae(:,1)) < 1e-8_RNP .and.  &
         norm2(Direction(d2, an) - ae(:,2)) < 1e-8_RNP

  end function OrientationMatches

  !-----------------------------------------------------------------------------
  !> Direction d (1..3: +ξ,+η,+ζ, 4..6: -ξ,-η,-ζ) of given axes

  function Direction(d, a) result(v)
    integer,   intent(in) :: d
    real(RNP), intent(in) :: a(3,3)
    real(RNP) :: v(3)

    if (d <= 3) then
      v =  a(:,d)
    else
      v = -a(:,d-3)
    end if

  end function Direction

  !-----------------------------------------------------------------------------
  !> Coordinates of the points on face f of element e

  function FacePoints(e, f) result(xf)
    integer, intent(in) :: e, f
    real(RNP) :: xf(0:po,0:po,3)

    select case(f)
    case(1); xf = x( 0,  :,  :, e, :)
    case(2); xf = x(po,  :,  :, e, :)
    case(3); xf = x( :,  0,  :, e, :)
    case(4); xf = x( :, po,  :, e, :)
    case(5); xf = x( :,  :,  0, e, :)
    case(6); xf = x( :,  :, po, e, :)
    end select

  end function FacePoints

  !-----------------------------------------------------------------------------
  !> Coordinates of the points on edge ed of element e

  function EdgePoints(e, ed) result(xe)
    integer, intent(in) :: e, ed
    real(RNP) :: xe(0:po,3)

    select case(ed)
    case( 1); xe = x( :,  0,  0, e, :)
    case( 2); xe = x( :, po,  0, e, :)
    case( 3); xe = x( :,  0, po, e, :)
    case( 4); xe = x( :, po, po, e, :)
    case( 5); xe = x( 0,  :,  0, e, :)
    case( 6); xe = x(po,  :,  0, e, :)
    case( 7); xe = x( 0,  :, po, e, :)
    case( 8); xe = x(po,  :, po, e, :)
    case( 9); xe = x( 0,  0,  :, e, :)
    case(10); xe = x(po,  0,  :, e, :)
    case(11); xe = x( 0, po,  :, e, :)
    case(12); xe = x(po, po,  :, e, :)
    end select

  end function EdgePoints

  !-----------------------------------------------------------------------------
  !> Coordinates of vertex v of element e

  function VertexPoint(e, v) result(xv)
    integer, intent(in) :: e, v
    real(RNP) :: xv(3)

    integer :: i, j, k

    i = po * mod(v - 1, 2)
    j = po * mod((v - 1) / 2, 2)
    k = po * ((v - 1) / 4)
    xv = x(i, j, k, e, :)

  end function VertexPoint

  !-----------------------------------------------------------------------------
  !> True if two points coincide

  logical function SamePoint(a, b)
    real(RNP), intent(in) :: a(3), b(3)

    SamePoint = Gap(1, a, b) <= tol

  end function SamePoint

  !-----------------------------------------------------------------------------
  !> Maximum distance of n points, modulo lx in periodic directions

  real(RNP) function Gap(n, a, b)
    integer,   intent(in) :: n
    real(RNP), intent(in) :: a(n,3), b(n,3)

    real(RNP) :: d
    integer   :: i, j

    Gap = 0
    do j = 1, 3
      do i = 1, n
        d = a(i,j) - b(i,j)
        if (per(j)) d = d - lx(j) * anint(d / lx(j))
        Gap = max(Gap, abs(d))
      end do
    end do

  end function Gap

  !-----------------------------------------------------------------------------
  !> True if 0 <= i, j <= po

  logical function Inside(i, j)
    integer, intent(in) :: i, j

    Inside = i >= 0 .and. i <= po .and. j >= 0 .and. j <= po

  end function Inside

  !-----------------------------------------------------------------------------
  !> True if a and b contain the same values

  logical function SameSet(a, b)
    integer, intent(in) :: a(:), b(:)

    integer :: i

    SameSet = size(a) == size(b)
    do i = 1, size(a)
      SameSet = SameSet .and. any(b == a(i)) .and. any(a == b(i))
    end do

  end function SameSet

end program Validate__Neighbor_Alignment_3D
