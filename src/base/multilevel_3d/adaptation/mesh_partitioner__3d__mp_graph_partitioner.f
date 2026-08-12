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

!> summary:  Graph-based mesh partitioning using Metis or ParMetis
!> author:   Joerg Stiller
!> date:     2026/01/27
!===============================================================================

submodule(Mesh_Partitioner__3D) MP_Graph_Partitioner
  use ParMETIS_Binding
  use Element_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Graph-based mesh partitioning using Metis or ParMetis

  module subroutine Graph_Partitioner(opt, mesh, tp_elem, n_parts)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh !< current mesh partition
    integer, intent(out) :: tp_elem(:) !< new/child element target partitions
    integer, intent(out) :: n_parts    !< actual number of partitions

    ! internal variables .......................................................

    integer(METIS_IDX_T), allocatable :: vtxdist(:), vwgt(:,:)
    integer(METIS_IDX_T), allocatable :: adjncy(:), adjwgt(:), xadj(:)

    integer, allocatable :: vtx_elem(:), vtx_part(:), vtx_part_loc(:)
    integer :: n_parts_loc, n_proc, proc
    integer :: n_vtx
    integer :: i, e, tp

    type(MPI_Comm) :: comm

    ! prerequisites ............................................................

    call IdentifyGraphVertices(opt, mesh, vtx_elem, n_vtx)
    call GetGraphVertexWeights(opt, mesh, n_vtx, vtx_elem, vwgt)
    call GlobalizeGraphVertices(opt, mesh, n_vtx, vtx_elem, comm, vtxdist)
    call GetAdjacency(opt, mesh, n_vtx, vtx_elem, xadj, adjncy, adjwgt)

    ! partitioning .............................................................

    if (n_vtx > 0) then

      if (log_level > 0) then
        call MPI_Comm_rank(comm, proc)
        call MPI_Comm_size(comm, n_proc)
        if (proc == 0) then
          write(*,'(/,A)') 'partitioner input'
          write(*,'(2X,A,X,I0)')   'opt%n_parts   =', opt%n_parts
          write(*,'(2X,A,X,I0)')   'n_proc        =', n_proc
          write(*,'(2X,99(G0,X))') 'vtxdist       =', vtxdist
        end if
        if (proc == 0 .or. log_level > 1) then
          write(*,'(2X,999(G0,X))') 'n_elem  [',proc,'] =', mesh%n_elem
          write(*,'(2X,999(G0,X))') 'n_ghost [',proc,'] =', mesh%n_ghost
          if (log_level > 2) then
            write(*,'(2X,999(G0,X))') 'vtx_elem[',proc,'] =', vtx_elem
            do i = lbound(vwgt,1),ubound(vwgt,1)
              write(*,'(2X,999(G0,X))') 'vwgt(',i,')[',proc,'] =', vwgt(i,:)
            end do
          end if
        end if
      end if

      if (opt%n_parts == 1) then
        where(vtx_elem >= 0)
          tp_elem = 0
        elsewhere
          tp_elem = -1
        end where
        n_parts_loc = 1
      else if (opt%child .and. opt%split .or. mesh%n_parts == 1) then
        call METIS_Partitioner( opt, mesh, vtx_elem, xadj, adjncy        &
                              , vwgt, adjwgt, comm, tp_elem, n_parts_loc )
      else
        call ParMETIS_Partitioner( opt, mesh, vtx_elem, vtxdist, xadj, adjncy &
                                 , vwgt, adjwgt, comm, tp_elem, n_parts_loc   )
      end if

    else
      tp_elem     = -1
      n_parts_loc = -1
    end if

    ! globalize n_parts
    call XMPI_Allreduce(n_parts_loc, n_parts, MPI_MAX, mesh%comm_parts)

    if (log_level > 0 .and. n_vtx > 0) then
      if (proc == 0) then
        write(*,'(/,A)') 'partitioning results'
        write(*,'(2X,A,I0)') 'n_parts    = ', n_parts
      end if
      if (log_level > 2) then
        write(*,'(2X,999(G0,X))') 'tp_elem[',proc,'] =', tp_elem(1:mesh%n_elem)
      end if
      allocate(vtx_part_loc(0:n_parts-1), source = 0)
      allocate(vtx_part, mold = vtx_part_loc)
      do e = 1, mesh%n_elem
        tp = tp_elem(e)
        if (tp >= 0) then
          vtx_part_loc(tp) = vtx_part_loc(tp) + 1
        end if
      end do
      call XMPI_Reduce(vtx_part_loc, vtx_part, MPI_SUM, 0, comm)
      if (proc == 0) then
        write(*,'(2X,A,I0)') 'sum(n_vtx) = ', sum(vtx_part)
        write(*,'(2X,A,I0)') 'min(n_vtx) = ', minval(vtx_part)
        write(*,'(2X,A,I0)') 'max(n_vtx) = ', maxval(vtx_part)
      end if
    end if

    ! clean-up .................................................................

    call MPI_Comm_free(comm)

  end subroutine Graph_Partitioner

  !-----------------------------------------------------------------------------
  !> METIS partitioner

  subroutine METIS_Partitioner( opt, mesh, vtx_elem, xadj, adjncy    &
                              , vwgt, adjwgt, comm, tp_elem, n_parts )

    ! arguments ...............................................................

    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(in) :: vtx_elem(:)
    integer(METIS_IDX_T), intent(in) :: xadj(0:), adjncy(0:)
    integer(METIS_IDX_T), intent(in) :: vwgt(0:,0:), adjwgt(0:)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: tp_elem(:) !< child element target partitions
    integer, intent(out) :: n_parts    !< actual number of partitions

    ! internal variables .......................................................

    real(METIS_REAL_T),   allocatable :: tpwgts(:,:)
    real(METIS_REAL_T),   allocatable :: ubvec(:)
    integer(METIS_IDX_T), allocatable :: part(:)

    integer(METIS_IDX_T) :: ncon, nvtx, nparts, edgecut

    integer, allocatable :: n_parts_proc(:)
    integer :: proc, n_proc, o_parts
    integer :: e, i

    ! prerequisites ............................................................

    call MPI_Comm_rank(comm, proc)
    call MPI_Comm_size(comm, n_proc)

    ncon = size(vwgt,1)
    nvtx = size(vwgt,2)

    ! number of partitions contributed by local parent
    nparts = max(0, min(opt%n_parts/n_proc, int(nvtx)))

    ! METIS array arguments, using C-style numbering
    allocate( tpwgts ( 0:ncon-1, 0:nparts-1 ) )
    allocate( ubvec  ( 0:ncon-1             ) )
    allocate( part   ( 0:nvtx-1             ) )

    ! partitioning .............................................................

    if (nparts > 1) then

      ! tolerance for multi-constraint weighting
      if (ncon == 1) then
        ubvec = 1.001
      else
        ubvec = 1.01
      end if

      ! fractions of vertex weight per partition
      tpwgts = 1.00 / nparts

      call METIS_PartGraphRecursive( nvtx, ncon, xadj, adjncy, vwgt, adjwgt &
                                   , nparts, tpwgts, ubvec, edgecut, part   )

    else
      part = 0
    end if

    ! result ...................................................................

    ! distribution of new partitions over processes
    allocate(n_parts_proc(0:n_proc-1))
    call MPI_Allgather( int(nparts) , 1, MPI_INTEGER       &
                      , n_parts_proc, 1, MPI_INTEGER, comm )

    ! total number of new partitions
    n_parts = sum(n_parts_proc)

    ! offset for locally generated partitions
    o_parts = sum(n_parts_proc(0:proc-1))

    i = 0
    do e = 1, mesh % n_elem
      if (vtx_elem(e) >= 0) then
        tp_elem(e) = part(i) + o_parts
        i = i + 1
      else
        tp_elem(e) = -1
      end if
    end do

  end subroutine METIS_Partitioner

  !-----------------------------------------------------------------------------
  !> ParMETIS partitioner

  subroutine ParMETIS_Partitioner( opt, mesh, vtx_elem, vtxdist, xadj, adjncy &
                                 , vwgt, adjwgt, comm, tp_elem, n_parts       )

    ! arguments ...............................................................

    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(in) :: vtx_elem(:)
    integer(METIS_IDX_T), intent(in) :: vtxdist(0:), xadj(0:), adjncy(0:)
    integer(METIS_IDX_T), intent(in) :: vwgt(0:,0:), adjwgt(0:)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(out) :: tp_elem(:) !< child element target partitions
    integer, intent(out) :: n_parts    !< actual number of partitions

    ! internal variables .......................................................

    real(METIS_REAL_T),   allocatable :: tpwgts(:,:)
    real(METIS_REAL_T),   allocatable :: ubvec(:)
    integer(METIS_IDX_T), allocatable :: part(:)

    integer(METIS_IDX_T) :: wgtflag, numflag, options(3)
    integer(METIS_IDX_T) :: ncon, nvtx, nparts, edgecut

    integer :: proc, n_proc
    integer :: e, i

    ! prerequisites ............................................................

    call MPI_Comm_rank(comm, proc)
    call MPI_Comm_size(comm, n_proc)

    numflag = 0
    nparts  = min(opt%n_parts, int(vtxdist(n_proc)))
    ncon    = size(vwgt,1)
    nvtx    = size(vwgt,2)

    options = 0

    if (any(opt % w_adj > 0)) then
      wgtflag = 3   ! graph vertex and edge constraints
    else
      wgtflag = 2   ! graph vertex constraints only
    end if

    ! ParMETIS array arguments, using C-style numbering
    allocate( tpwgts ( 0:ncon-1, 0:nparts-1 ) )
    allocate( ubvec  ( 0:ncon-1             ) )
    allocate( part   ( 0:nvtx-1             ) )

    ! tolerance for multi-constraint weighting
    ubvec  = 1.05

    ! fractions of vertex weight per partition
    tpwgts = 1.00 / nparts

    ! partitioning .............................................................

    call ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt, adjwgt, wgtflag  &
                             , numflag, ncon, nparts, tpwgts, ubvec, options &
                             , edgecut, part, comm % MPI_VAL                 )

    if (log_level > 0) then
      print '(999(G0,X))', '#PMKway# proc[',proc,'] vtxdist =',vtxdist
      print '(999(G0,X))', '#PMKway# proc[',proc,'] xadj    =',xadj
      print '(999(G0,X))', '#PMKway# proc[',proc,'] adjncy  =',adjncy
      print '(999(G0,X))', '#PMKway# proc[',proc,'] vwgt    =',vwgt
      print '(999(G0,X))', '#PMKway# proc[',proc,'] adjwgt  =',adjwgt
      print '(999(G0,X))', '#PMKway# proc[',proc,'] wgtflag =',wgtflag
      print '(999(G0,X))', '#PMKway# proc[',proc,'] numflag =',numflag
      print '(999(G0,X))', '#PMKway# proc[',proc,'] ncon    =',ncon
      print '(999(G0,X))', '#PMKway# proc[',proc,'] nparts  =',nparts
      print '(999(G0,X))', '#PMKway# proc[',proc,'] tpwgts  =',tpwgts
      print '(999(G0,X))', '#PMKway# proc[',proc,'] ubvec   =',ubvec
      print '(999(G0,X))', '#PMKway# proc[',proc,'] options =',options
      print '(999(G0,X))', '#PMKway# proc[',proc,'] edgecut =',edgecut
      print '(999(G0,X))', '#PMKway# proc[',proc,'] part    =',part
    end if

    ! result ...................................................................

    n_parts = nparts

    i = 0
    do e = 1, mesh % n_elem
      if (vtx_elem(e) >= 0) then
        tp_elem(e) = part(i)
        i = i + 1
      else
        tp_elem(e) = -1
      end if
    end do

  end subroutine ParMETIS_Partitioner

  !-----------------------------------------------------------------------------
  !> Identification of local graph vertices

  subroutine IdentifyGraphVertices(opt, mesh, vtx_elem, n_vtx)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
      !< local mesh partition
    integer, allocatable, intent(out) :: vtx_elem(:)
      !< map from elements to graph vertices
    integer, intent(out)   :: n_vtx
      !< number of local graph vertices

    integer :: e, n

    allocate(vtx_elem(mesh%n_elem + mesh%n_ghost), source = -1)

    n = 0
    do e = 1, mesh%n_elem
      if (opt%child .and. mesh%element(e)%adaptation%mark < 100) cycle
      vtx_elem(e) = n
      n = n + 1
    end do
    n_vtx = n

  end subroutine IdentifyGraphVertices

  !-----------------------------------------------------------------------------
  !> Determination of graph vertex weights

  subroutine GetGraphVertexWeights(opt, mesh, n_vtx, vtx_elem, vwgt)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in)  :: mesh
    integer, intent(in) :: n_vtx
    integer, contiguous, intent(in) :: vtx_elem(:)
    integer(METIS_IDX_T), allocatable, intent(out) :: vwgt(:,:)
      !< graph vertex weights

    integer, allocatable, save :: vwgt_raw(:,:)
    integer, allocatable, save :: max_vwgt(:), max_vwgt_loc(:)
    integer :: n_con, n_con_raw
    integer :: c_active, c_frozen
    integer :: c, g, q, w
    integer :: e, i, l

    ! initialization ...........................................................

    if (opt%child) then
      n_con_raw = max(1, min(opt%n_con_child,3))
    else
      n_con_raw = max(1, min(opt%n_con_root,4))
    end if

    allocate(vwgt_raw(0:n_con_raw-1, 0:n_vtx-1))
    allocate(max_vwgt(0:n_con_raw-1), source = 0)
    allocate(max_vwgt_loc, source = max_vwgt)

    c_active = opt % c_active
    c_frozen = opt % c_frozen

    ! raw weights ..............................................................

    do e = 1, mesh%n_elem
      associate(element => mesh % element(e))

        i = vtx_elem(e)
        if (i < 0) cycle

        ! children
        select case(element % adaptation % mark)
        case(100:108)
          c = 1 * c_frozen  ! 1 frozen child @ vertex or cloned
        case(201:212)
          c = 2 * c_frozen  ! 2 frozen children @ edge
        case(401:406)
          c = 4 * c_frozen  ! 4 frozen children @ face
        case(800)
          c = 8 * c_frozen  ! 8 frozen children @ element
        case(1000,8000)
          c = 8 * c_active  ! 8 active children @ element or cloned
        end select

        ! grandchildren
        if (element%adaptation%sublevels > 1) then
          g = 1
        else
          g = 0
        end if

        ! great-grandchildren and further descendants
        q = 0
        w = 1
        do l = 3, min(element%adaptation%sublevels, opt%n_con_sub)
          q = q + w
          w = 8 * w
        end do

        if (opt%child) then
          select case(n_con_raw)
          case(1)
            vwgt_raw(0,i) = c + 64 * (g + 8*q)
          case(2)
            vwgt_raw(0,i) = c
            vwgt_raw(1,i) = g + 8 * q
          case(3)
            vwgt_raw(0,i) = c
            vwgt_raw(1,i) = g
            vwgt_raw(2,i) = q
          end select
        else
          vwgt_raw(0,i) = 1
          select case(n_con_raw)
          case(2)
            vwgt_raw(1,i) = c + 64 * (g + 8*q)
          case(3)
            vwgt_raw(1,i) = c
            vwgt_raw(2,i) = g + 8 * q
          case(4)
            vwgt_raw(1,i) = c
            vwgt_raw(2,i) = g
            vwgt_raw(3,i) = q
          end select
        end if

      end associate
    end do

    ! remove constraints with zero sum weight ..................................

    do i = 0, n_con_raw-1
      max_vwgt_loc(i) = maxval(vwgt_raw(i,:))
    end do
    call XMPI_Allreduce(max_vwgt_loc, max_vwgt, MPI_MAX, mesh%comm_parts)

    n_con = 0
    do i = 0, n_con_raw-1
      if (max_vwgt(i) == 0) exit
      n_con = n_con + 1
    end do

    if (n_con == 0) then
      call Error( 'GetGraphVertexWeights'                     &
                , 'invalid vertex weights'                    &
                , 'Mesh_Partitioner__3D:MP_Graph_Partitioner' )
    else if (n_con == n_con_raw) then
      call move_alloc(vwgt_raw, vwgt)
    else
      allocate(vwgt(0:n_con-1, 0:n_vtx-1), source = vwgt_raw(0:n_con-1,:))
      deallocate(vwgt_raw)
    end if

    ! finalization .............................................................

    deallocate(max_vwgt, max_vwgt_loc)

  end subroutine GetGraphVertexWeights

  !-----------------------------------------------------------------------------
  !> Globalization of graph vertices

  subroutine GlobalizeGraphVertices(opt, mesh, n_vtx, vtx_elem, comm, vtxdist)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
      !< local mesh partition
    integer, intent(in)   :: n_vtx
      !< number of local graph vertices
    integer, contiguous, target, intent(inout) :: vtx_elem(:)
      !< map from elements to graph vertices, will be converted to global IDs
    type(MPI_Comm) :: comm
      !< communicator between mesh partitions contributing to the graph
    integer(METIS_IDX_T), allocatable, intent(out) :: vtxdist(:)
      !< distribution graph vertices over processes in comm

    ! internal variables .......................................................

    type(ElementTransferBuffer_3D), asynchronous :: buf_vtx_elem
    integer, contiguous, pointer :: var_vtx_elem(:,:,:,:)
    integer, allocatable :: nvtx_proc(:)
    integer :: i, proc, n_proc

    ! MPI communicator for contributors ........................................

    if (n_vtx > 0) then
      i = 1
    else
      i = 0
    end if
    call MPI_Comm_split(mesh%comm_parts, i, mesh%part, comm)

    ! graph vertex distribution ................................................

    if (n_vtx > 0) then

      call MPI_Comm_rank(comm, proc)
      call MPI_Comm_size(comm, n_proc)

      allocate(vtxdist(0:n_proc), nvtx_proc(0:n_proc-1))

      call MPI_Allgather(n_vtx, 1, MPI_INTEGER, nvtx_proc, 1, MPI_INTEGER, comm)

      vtxdist(0) = 0
      do i = 1, n_proc
        vtxdist(i) = vtxdist(i-1) + nvtx_proc(i-1)
      end do

    else
      allocate(vtxdist(0), nvtx_proc(0), source = 0)
    end if

    ! globalization of graph vertex IDs ........................................

    if (opt%child .and. opt%split) return

    ! global ID of local vertices
    if (n_vtx > 0) then
      where(vtx_elem >= 0)
        vtx_elem = vtx_elem + vtxdist(proc)
      end where
    end if

    ! transfer vertex IDs to ghosts
    var_vtx_elem(1:1, 1:1, 1:1, 1:size(vtx_elem)) => vtx_elem
    buf_vtx_elem = ElementTransferBuffer_3D(mesh, var_vtx_elem)
    call buf_vtx_elem % Transfer(mesh, var_vtx_elem, tag=1000)
    call buf_vtx_elem % Merge(var_vtx_elem)

  end subroutine GlobalizeGraphVertices

  !-----------------------------------------------------------------------------
  !> Determine graph adjacency and corresponding weights

  subroutine GetAdjacency(opt, mesh, n_vtx, vtx_elem, xadj, adjncy, adjwgt)
    class(MeshPartitionerOptions_3D), intent(in) :: opt
    class(Mesh_3D), intent(in) :: mesh
    integer, intent(in)   :: n_vtx
    integer, contiguous, intent(in) :: vtx_elem(:)
    integer(METIS_IDX_T), allocatable, intent(out) :: xadj(:)
    integer(METIS_IDX_T), allocatable, intent(out) :: adjncy(:)
    integer(METIS_IDX_T), allocatable, intent(out) :: adjwgt(:)

    integer :: e, i, j, k, l, m, n

    if (n_vtx == 0) return

    associate(w_adj => opt % w_adj)

      ! adjacency offsets ......................................................

      allocate(xadj(0:n_vtx))

      xadj(0) = 0

      m = 0
      do e = 1, mesh%n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

          m = m + 1

          xadj(m) = xadj(m-1)

          ! count graph edges contributed by element faces
          do k = 1, 6
            n = element % face(k) % n_neighbor - 1
            if (n < 0) cycle
            i = element % face(k) % i_neighbor
            do j = i, i+n
              if (vtx_elem(element % neighbor(j) % id) < 0) cycle
              xadj(m) = xadj(m) + 1
            end do
          end do

          if (w_adj(2) > 0) then
            ! count graph edges contributed by element edges
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                if (vtx_elem(element % neighbor(j) % id) < 0) cycle
                xadj(m) = xadj(m) + 1
              end do
            end do
          end if

          if (w_adj(3) > 0) then
            ! count graph edges contributed by element vertices
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                if (vtx_elem(element % neighbor(j) % id) < 0) cycle
                xadj(m) = xadj(m) + 1
              end do
            end do
          end if

        end associate
      end do

      ! adjacency and weights ..................................................

      m = xadj(n_vtx) - 1 ! number of graph edges
      allocate(adjncy(0:m))
      if (any(w_adj > 0)) then
        allocate(adjwgt(0:m), source = 0)
      else
        allocate(adjwgt(0:0))
      end if

      m = 0

      do e = 1, mesh%n_elem
        if (vtx_elem(e) < 0) cycle
        associate(element => mesh % element(e))

          ! element faces
          do k = 1, 6
            n = element % face(k) % n_neighbor - 1
            if (n < 0) cycle
            i = element % face(k) % i_neighbor
            do j = i, i+n
              l = vtx_elem(element % neighbor(j) % id)
              if (l < 0) cycle
              adjncy(m) = l
              if (w_adj(1) > 0) then
                adjwgt(m) = w_adj(1)
              end if
              m = m + 1
            end do
          end do

          ! element edges
          if (w_adj(2) > 0) then
            do k = 1, 12
              n = element % edge(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % edge(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_adj(2)
                m = m + 1
              end do
            end do
          end if

          ! element vertices
          if (w_adj(3) > 0) then
            do k = 1, 8
              n = element % vertex(k) % n_neighbor - 1
              if (n < 0) cycle
              i = element % vertex(k) % i_neighbor
              do j = i, i+n
                l = vtx_elem(element % neighbor(j) % id)
                if (l < 0) cycle
                adjncy(m) = l
                adjwgt(m) = w_adj(3)
                m = m + 1
              end do
            end do
          end if

        end associate
      end do

    end associate

  end subroutine GetAdjacency

  !=============================================================================

end submodule MP_Graph_Partitioner
