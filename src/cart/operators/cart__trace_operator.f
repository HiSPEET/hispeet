!> summary:  Extraction of face traces from scalar or array variables
!> author:   Joerg Stiller
!> date:     2018/03/13
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Extraction of face traces from scalar or array variables
!>
!> @todo
!> OpenACC
!>
!>   *  assumption: mesh (element and trace) variables on GPU
!>   *  test/decide what can be done on GPU
!>   *  implications when using with MPI via TraceTransferBuffer
!>   *  synchronisation between CPU/GPU
!>   *  OpenACC binding --> function to determine the vector length
!>
!> @endtodo
!===============================================================================

module CART__Trace_Operator

  use Kind_Parameters,         only: RNP
  use Constants,               only: ONE, ZERO
  use Array_Assignments,       only: SetArray
  use CART__Mesh_Partition,    only: MeshPartition
  use CART__Boundary_Variable, only: BoundaryVariable
  use CART__Trace_Transfer_Buffer

  implicit none
  private

  public :: TraceOperator

  !-----------------------------------------------------------------------------
  !> Operator for generating, keeping and handling trace data

  type TraceOperator

    type(TraceTransferBuffer) :: buf  !< trace transfer buffer

  contains

    generic,   public  :: GetTrace_Start =>    &
                            GetTrace_Start_SB, &
                            GetTrace_Start_SO, &
                            GetTrace_Start_AB, &
                            GetTrace_Start_AO
    procedure, private :: GetTrace_Start_SB
    procedure, private :: GetTrace_Start_SO
    procedure, private :: GetTrace_Start_AB
    procedure, private :: GetTrace_Start_AO

    generic,   public  :: GetTrace_Finish =>   &
                            GetTrace_Finish_S, &
                            GetTrace_Finish_A
    procedure, private :: GetTrace_Finish_S
    procedure, private :: GetTrace_Finish_A

  end type TraceOperator

  !-----------------------------------------------------------------------------
  ! private variables

  !> mapping of boundary face orientation to inner trace side
  integer, parameter :: inner_side(-3:3) = [ 2, 2, 2, 0, 1, 1, 1 ]

contains

!===============================================================================
! GetTrace_Start variants and auxiliary routines

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- scalar with BC

subroutine GetTrace_Start_SB(this, mesh, u, bv_u, tr_u, tag)
  class(TraceOperator), asynchronous, intent(inout) :: this !< trace operator
  type(MeshPartition),    intent(in)  :: mesh            !< mesh partition
  real(RNP),              intent(in)  :: u(0:,0:,0:,:)   !< mesh variable
  type(BoundaryVariable), intent(in)  :: bv_u(:)         !< boundary cond.
  real(RNP),              intent(out) :: tr_u(0:,0:,:,:) !< trace of u
  integer,                intent(in)  :: tag             !< message tag

  integer :: po, ne, nf, nc

  po = ubound(u,1)
  ne = size(u,4)
  nf = size(tr_u,4)
  nc = 1

  call GetTrace_Start_X(this, mesh, po, ne, nf, nc, u, bv_u, tr_u, tag)

end subroutine GetTrace_Start_SB

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- scalar, open (no BC)

subroutine GetTrace_Start_SO(this, mesh, u, tr_u, tag)
  class(TraceOperator), asynchronous, intent(inout) :: this !< trace operator
  type(MeshPartition), intent(in)  :: mesh            !< mesh partition
  real(RNP),           intent(in)  :: u(0:,0:,0:,:)   !< mesh variable
  real(RNP),           intent(out) :: tr_u(0:,0:,:,:) !< trace of u
  integer,             intent(in)  :: tag             !< message tag

  integer :: po, ne, nf, nc

  po = ubound(u,1)
  ne = size(u,4)
  nf = size(tr_u,4)
  nc = 1

  call GetTrace_Start_X(this, mesh, po, ne, nf, nc, u, tr_u=tr_u, tag=tag)

end subroutine GetTrace_Start_SO

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- array with BC

subroutine GetTrace_Start_AB(this, mesh, u, bv_u, tr_u, tag)
  class(TraceOperator), asynchronous, intent(inout) :: this !< trace operator
  type(MeshPartition),    intent(in)  :: mesh              !< mesh partition
  real(RNP),              intent(in)  :: u(0:,0:,0:,:,:)   !< mesh variable
  type(BoundaryVariable), intent(in)  :: bv_u(:)           !< boundary cond.
  real(RNP),              intent(out) :: tr_u(0:,0:,:,:,:) !< trace of u
  integer,                intent(in)  :: tag               !< message tag

  integer :: po, ne, nf, nc

  po = ubound(u,1)
  ne = size(u,4)
  nf = size(tr_u,4)
  nc = size(tr_u,5)

  call GetTrace_Start_X(this, mesh, po, ne, nf, nc, u, bv_u, tr_u, tag)

end subroutine GetTrace_Start_AB

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- array, open (no BC)

subroutine GetTrace_Start_AO(this, mesh, u, tr_u, tag)
  class(TraceOperator), asynchronous, intent(inout) :: this !< trace operator
  type(MeshPartition), intent(in)  :: mesh              !< mesh partition
  real(RNP),           intent(in)  :: u(0:,0:,0:,:,:)   !< mesh variable
  real(RNP),           intent(out) :: tr_u(0:,0:,:,:,:) !< trace of u
  integer,             intent(in)  :: tag               !< message tag

  integer :: po, ne, nf, nc

  po = ubound(u,1)
  ne = size(u,4)
  nf = size(tr_u,4)
  nc = size(tr_u,5)

  call GetTrace_Start_X(this, mesh, po, ne, nf, nc, u, tr_u=tr_u, tag=tag)

end subroutine GetTrace_Start_AO

!-------------------------------------------------------------------------------
!> Extraction of face traces from scalar or array variables -- eXplicit

subroutine GetTrace_Start_X(this, mesh, po, ne, nf, nc, u, bv_u, tr_u, tag)

  ! arguments ..................................................................

  class(TraceOperator), asynchronous, intent(inout) :: this !< trace operator
  type(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                intent(in)    :: po    !< polynomial order
  integer,                intent(in)    :: ne    !< number of elements
  integer,                intent(in)    :: nf    !< number of faces
  integer,                intent(in)    :: nc    !< number of components
  real(RNP),              intent(in)    :: u     !< mesh variable
  type(BoundaryVariable), intent(in)    :: bv_u  !< boundary conditions & values
  real(RNP),              intent(out)   :: tr_u  !< trace of u on faces
  integer,                intent(in)    :: tag   !< message tag

  dimension :: u    (0:po, 0:po, 0:po, ne, nc)
  dimension :: tr_u (0:po, 0:po, 1:2 , nf, nc)
  dimension :: bv_u (mesh%n_boundary)
  optional  :: bv_u

  ! initialization .............................................................

  ! prevent delayed threads from interfering
  !$omp master
  this % buf = TraceTransferBuffer(mesh, tr_u)
  !$omp end master
  !$omp barrier

  ! tr_u = 0
  call SetArray(tr_u, ZERO, multi=.true.)

  ! get local trace ............................................................

  call GetLocalTrace(mesh, po, ne, nf, nc, u, tr_u)

  ! start transfer of shared trace .............................................

  call this % buf % Transfer(mesh, tr_u, tag)

  ! apply boundary conditions ..................................................

  if (present(bv_u)) then
    call ApplyBoundaryConditions(mesh, po, nf, nc, bv_u, tr_u)
  else
    call ExtendToBoundary(mesh, po, nf, nc, tr_u)
  end if

end subroutine GetTrace_Start_X

!-------------------------------------------------------------------------------
!> Extract local trace

subroutine GetLocalTrace(mesh, po, ne, nf, nc, u, tr_u)

  ! arguments ..................................................................

  class(MeshPartition),  intent(in)    :: mesh    !< mesh partition
  integer,               intent(in)    :: po      !< polynomial order
  integer,               intent(in)    :: ne      !< number of elements
  integer,               intent(in)    :: nf      !< number of faces
  integer,               intent(in)    :: nc      !< number of components
  real(RNP),             intent(in)    :: u       !< mesh variable
  real(RNP),             intent(inout) :: tr_u    !< trace of u

  dimension :: u    (0:po, 0:po, 0:po, ne, nc)
  dimension :: tr_u (0:po, 0:po, 1:2 , nf, nc)

  ! internal variables .........................................................

  integer :: c, e, f(6), i, j, k

  ! extract local trace ........................................................

  !$acc data present(u, tr_u)
  !$acc parallel
  !$acc loop collapse(2) gang worker private(f)

  !$omp do collapse(2) private(c,e,f)
  do c = 1, nc
  do e = 1, mesh%ne

    f = mesh % element(e) % face % id

    ! west
    !$acc loop collapse(2) vector async(1)
    do k = 0, po
    do j = 0, po
      tr_u(j,k,2,f(1),c) = u( 0,j,k,e,c)
    end do
    end do

    ! east
    !$acc loop collapse(2) vector async(2)
    do k = 0, po
    do j = 0, po
      tr_u(j,k,1,f(2),c) = u(po,j,k,e,c)
    end do
    end do

    ! south
    !$acc loop collapse(2) vector async(3)
    do k = 0, po
    do i = 0, po
      tr_u(i,k,2,f(3),c) = u(i, 0,k,e,c)
    end do
    end do

    ! north
    !$acc loop collapse(2) vector async(4)
    do k = 0, po
    do i = 0, po
      tr_u(i,k,1,f(4),c) = u(i,po,k,e,c)
    end do
    end do

    ! bottom
    !$acc loop collapse(2) vector async(5)
    do j = 0, po
    do i = 0, po
      tr_u(i,j,2,f(5),c) = u(i,j, 0,e,c)
    end do
    end do

    ! top
    !$acc loop collapse(2) vector async(6)
    do j = 0, po
    do i = 0, po
      tr_u(i,j,1,f(6),c) = u(i,j,po,e,c)
    end do
    end do

    !$acc wait
  end do
  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine GetLocalTrace

!-------------------------------------------------------------------------------
!> Application of boundary conditions to trace data

subroutine ApplyBoundaryConditions(mesh, po, nf, nc, bv_u, tr_u)

  ! arguments ..................................................................

  type(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                intent(in)    :: po    !< polynomial order
  integer,                intent(in)    :: nf    !< number of faces
  integer,                intent(in)    :: nc    !< number of components
  type(BoundaryVariable), intent(in)    :: bv_u  !< boundary conditions & values
  real(RNP),              intent(inout) :: tr_u  !< trace of u

  dimension :: bv_u (mesh%n_boundary)
  dimension :: tr_u (0:po, 0:po, 1:2, nf, nc)

  ! internal variables .........................................................

  real(RNP), contiguous, pointer, save :: ub(:,:,:) => null()
  integer :: b, c, f, i, k, o, p, q

  ! set boundary conditions ....................................................

  do b = 1, mesh%n_boundary
    associate(boundary => mesh%boundary(b))

      do c = 1, nc
        select case(bv_u(b) % BoundaryCondition(c))

        case('D')
          !$omp single
          ub(0:,0:,1:) => bv_u(b) % Component(c)
          !$omp end single

          !$omp do private(k)
          do k = 1, boundary%nf
            f = boundary % face(k) % mesh_face % id
            i = inner_side(boundary % face(k) % orientation)
            o = 3 - i
            do q = 0, po
            do p = 0, po
              tr_u(p,q,o,f,c) = 2 * ub(p,q,k) - tr_u(p,q,i,f,c)
            end do
            end do
          end do

        case('N')
          !$omp do private(k)
          do k = 1, boundary%nf
            f = boundary % face(k) % mesh_face % id
            i = inner_side(boundary % face(k) % orientation)
            o = 3 - i
            do q = 0, po
            do p = 0, po
              tr_u(p,q,o,f,c) = tr_u(p,q,i,f,c)
            end do
            end do
          end do

        end select
      end do

    end associate
  end do

end subroutine ApplyBoundaryConditions

!-------------------------------------------------------------------------------
!>

subroutine ExtendToBoundary(mesh, po, nf, nc, tr_u)

  ! arguments ..................................................................

  type(MeshPartition), intent(in)    :: mesh  !< mesh partition
  integer,             intent(in)    :: po    !< polynomial order
  integer,             intent(in)    :: nf    !< number of faces
  integer,             intent(in)    :: nc    !< number of components
  real(RNP),           intent(inout) :: tr_u  !< trace of u

  dimension :: tr_u (0:po, 0:po, 1:2, nf, nc)

  ! internal variables .........................................................

  integer :: b, c, f, i, k, o, p, q

  ! extend to boundary .........................................................

  do b = 1, mesh%n_boundary
    associate(boundary => mesh%boundary(b))

      !$omp do collapse(2) private(c,k)
      do c = 1, nc
        do k = 1, boundary%nf

          f = boundary % face(k) % mesh_face % id
          i = inner_side(boundary % face(k) % orientation)
          o = 3 - i

          do q = 0, po
          do p = 0, po
            tr_u(p,q,o,f,c) = tr_u(p,q,i,f,c)
          end do
          end do

        end do
      end do

    end associate
  end do

end subroutine ExtendToBoundary

!===============================================================================
! GetTrace_Finish variants

!-------------------------------------------------------------------------------
!> Merge remote traces and finish transfer -- single variable

subroutine GetTrace_Finish_S(this, mesh, tr_u)
  class(TraceOperator), intent(inout) :: this            !< trace operator
  type(MeshPartition),  intent(in)    :: mesh            !< mesh partition
  real(RNP),            intent(inout) :: tr_u(:,:,:,:)   !< trace of u

  call this % buf % Merge(mesh, tr_u, alpha=ZERO, beta=ONE)

end subroutine GetTrace_Finish_S

!-------------------------------------------------------------------------------
!> Merge remote traces and finish transfer -- array of variables

subroutine GetTrace_Finish_A(this, mesh, tr_u)
  class(TraceOperator), intent(inout) :: this            !< trace operator
  type(MeshPartition),  intent(in)    :: mesh            !< mesh partition
  real(RNP),            intent(inout) :: tr_u(:,:,:,:,:) !< trace of u

  call this % buf % Merge(mesh, tr_u, alpha=ZERO, beta=ONE)

end subroutine GetTrace_Finish_A

!===============================================================================

end module CART__Trace_Operator
