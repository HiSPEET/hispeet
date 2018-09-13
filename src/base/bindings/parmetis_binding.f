!> summary:  Interface to ParMETIS
!> author:   Joerg Stiller
!> date:     2014/10/09
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Interface to ParMETIS
!>
!> @todo
!>
!>   *  Add interfaces for further ParMETIS routines
!===============================================================================

module ParMETIS_Binding
  use, intrinsic :: ISO_C_Binding, only: C_INT, C_FLOAT
  implicit none

  private :: C_INT, C_FLOAT

  !------------------------------------------------------------------------------
  ! Kind parameters

  !> Integer kind parameter matching Metis' idx_t
  integer, parameter, public :: METIS_IDX_T = C_INT

  !> Real kind parameter matching Metis' real_t
  integer, parameter, public :: METIS_REAL_T = C_FLOAT

  !------------------------------------------------------------------------------
  !> Interface to ParMETIS_V3_PartKway.
  !>
  !> This interface provides a one-to-one binding to ParMETIS_V3_PartKway.
  !> Though possible, it is not advisable to access ParMETIS directly through
  !> the interface. We recommend to encapsulate the call in an intermediate
  !> layer which translates the application data into ParMETIS form.

  interface

    subroutine ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt, &
        wgtflag, numflag, ncon, nparts, tpwgts, ubvec, options, edgecut, &
        part, comm) bind(C, name='ParMETIS_V3_PartKway_F2C')

      import :: METIS_IDX_T, METIS_REAL_T

      !> Distribution of the graph vertices among processors. The size of
      !> vtxdist equals the number of the processors calling the routine,
      !> i.e., the size of the MPI communicator comm).
      integer(METIS_IDX_T), intent(in) :: vtxdist(*)

      !> Indices for accessing the adjacency list adjncy. The size of xadj
      !> equals n+1, where n is the number of locally-stored vertices.
      integer(METIS_IDX_T), intent(in) :: xadj(*)

      !> Adjacency list: adjncy(xadj(i):xadj(i+1)-1) identifies the vertices
      !> that are connected with vertex i. The size of adjncy equals 2m, where
      !> m is the number of the edges of the graph.
      integer(METIS_IDX_T), intent(in) :: adjncy(*)

      !> Number of weights (constraints) per vertex.
      integer(METIS_IDX_T), intent(in) :: ncon

      !> Specifies the vertex weights. The dimensions of vwgt are (ncon, n),
      !> where n is the number of locally-stored vertices. However, if not
      !> required (wgtflag = 0 or 1), a zero-sized array can be passed instead.
      integer(METIS_IDX_T), intent(in) :: vwgt(ncon,*)

      !> Edge weights. The size of adjwgt equals that of adjncy. However, if not
      !> required (wgtflag = 0 or 2), a zero-sized array can be passed instead.
      integer(METIS_IDX_T), intent(in) :: adjwgt(*)

      !> Weighting method:
      !>   * 0  None (vwgt and adjwgt not beeing accessed).
      !>   * 1  Edges only (vwgt not beeing accessed).
      !>   * 2  Vertices only (adjwgt not beeing accessed).
      !>   * 3  Vertices and edges.
      integer(METIS_IDX_T), intent(in) :: wgtflag

      !> Numbering scheme:
      !>   * 0  C-style numbering that starts from 0.
      !>   * 1  Fortran-style numbering that starts from 1.
      integer(METIS_IDX_T), intent(in) :: numflag

      !> Specifies the number of partitions that are desired. Note that the
      !> number of partitions is independent of the number of processors that
      !> call this routine.
      integer(METIS_IDX_T), intent(in) :: nparts

      !> Specifies the fraction of vertex weight that should be distributed to
      !> each partition for each balance constraint. If all of the partitions
      !> are to be of the same size for every vertex weight, then tpwgts should
      !> be set to a value of 1/nparts. The sum of all of the tpwgts for a given
      !> vertex weight should be one.
      real(METIS_REAL_T), intent(in) :: tpwgts(ncon, nparts)

      !> Array of size ncon that specifying the imbalance tolerance for each
      !> vertex weight, with 1 being perfect balance and nparts being perfect
      !> imbalance. A value of 1.05 for each of the ncon weights is recommended.
      real(METIS_REAL_T), intent(in) :: ubvec(ncon)

      !> Additional parameters for the routine. Pass options = 0 for the default
      !> behavior. For other choices confer to the ParMETIS manual.
      integer(METIS_IDX_T), intent(in) :: options(0:2)

      !> Number of edges that are cut by the partitioning.
      integer(METIS_IDX_T), intent(out) :: edgecut

      !> Partition vector of the locally-stored vertices. The size of part
      !> equals the number of locally-stored vertices.
      integer(METIS_IDX_T), intent(out) :: part(*)

      !> MPI communicator value. With the MPI-3 F08 binding use the MPI_VAL
      !> component of the type MPI_Comm
      integer, intent(in) ::  comm

    end subroutine ParMETIS_V3_PartKway

  end interface

!===============================================================================

end module ParMETIS_Binding
