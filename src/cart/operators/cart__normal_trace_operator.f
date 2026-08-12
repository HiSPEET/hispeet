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

!> summary:  Extraction of face normal traces from vector-valued variables
!> author:   Joerg Stiller
!> date:     2018/09/17
!>
!>### Extraction of face normal traces from vector-valued variables
!>
!> Boundary conditions
!>
!>   *  Dirichlet conditions are enforced if boundary values are supplied.
!>
!>   *  At Neumann boundaries the normal component is extrapolated.
!>      Inhomogeneous Neumann conditions may require additional treatment,
!>      which is deliberately skipped here.
!>
!===============================================================================

module CART__Normal_Trace_Operator

  use Kind_Parameters,         only: RNP
  use Constants,               only: ONE, ZERO
  use Array_Assignments,       only: SetArray
  use CART__Mesh_Partition,    only: MeshPartition
  use CART__Boundary_Variable, only: BoundaryVariable
  use CART__Trace_Transfer_Buffer

  implicit none
  private

  public :: NormalTraceOperator

  !-----------------------------------------------------------------------------
  !> Operator for generating, keeping and handling trace data

  type NormalTraceOperator

    type(TraceTransferBuffer) :: buf  !< trace transfer buffer

  contains

    generic,   public  :: GetTrace_Start => GetTrace_Start_B, GetTrace_Start_O
    procedure, private :: GetTrace_Start_B
    procedure, private :: GetTrace_Start_O

    procedure, public  :: GetTrace_Finish

  end type NormalTraceOperator

  !-----------------------------------------------------------------------------
  ! private variables

  !> mapping of boundary face orientation to inner trace side
  integer, parameter :: inner_side(-3:3) = [ 2, 2, 2, 0, 1, 1, 1 ]

contains

!===============================================================================
! GetTrace_Start variants and auxiliary routines

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- with BC

subroutine GetTrace_Start_B(this, mesh, u, bv_u, tr_un, tag)
  class(NormalTraceOperator), asynchronous, intent(inout) :: this
  type(MeshPartition),    intent(in)  :: mesh             !< mesh partition
  real(RNP),              intent(in)  :: u(0:,0:,0:,:,:)  !< vector variable
  type(BoundaryVariable), intent(in)  :: bv_u(:)          !< boundary conditions
  real(RNP),              intent(out) :: tr_un(0:,0:,:,:) !< normal trace of u
  integer,                intent(in)  :: tag              !< message tag

  integer :: po, ne, nf

  po = ubound(u,1)
  ne = mesh % ne
  nf = mesh % nf

  call GetTrace_Start_X(this, mesh, po, ne, nf, u, bv_u, tr_un, tag)

end subroutine GetTrace_Start_B

!-------------------------------------------------------------------------------
!> Extract and start transferring traces -- open (no BC)

subroutine GetTrace_Start_O(this, mesh, u, tr_un, tag)
  class(NormalTraceOperator), asynchronous, intent(inout) :: this
  type(MeshPartition), intent(in)    :: mesh             !< mesh partition
  real(RNP),           intent(in)    :: u(0:,0:,0:,:,:)  !< vector variable
  real(RNP),           intent(out)   :: tr_un(0:,0:,:,:) !< normal trace of u
  integer,             intent(in)    :: tag              !< message tag

  integer :: po, ne, nf

  po = ubound(u,1)
  ne = mesh % ne
  nf = mesh % nf

  call GetTrace_Start_X(this, mesh, po, ne, nf, u, tr_un=tr_un, tag=tag)

end subroutine GetTrace_Start_O

!-------------------------------------------------------------------------------
!> Extraction of face traces from scalar or array variables -- eXplicit

subroutine GetTrace_Start_X(this, mesh, po, ne, nf, u, bv_u, tr_un, tag)

  ! arguments ..................................................................

  class(NormalTraceOperator), asynchronous, intent(inout) :: this
  type(MeshPartition),        intent(in)    :: mesh  !< mesh partition
  integer,                    intent(in)    :: po    !< polynomial order
  integer,                    intent(in)    :: ne    !< number of elements
  integer,                    intent(in)    :: nf    !< number of faces
  real(RNP),                  intent(in)    :: u     !< vector variable
  type(BoundaryVariable),     intent(in)    :: bv_u  !< boundary conditions
  real(RNP),                  intent(out)   :: tr_un !< trace of u on faces
  integer,                    intent(in)    :: tag   !< message tag

  dimension :: u     (0:po, 0:po, 0:po, ne, 3)
  dimension :: tr_un (0:po, 0:po, 1:2 , nf)
  dimension :: bv_u  (mesh%n_boundary)
  optional  :: bv_u

  ! initialization .............................................................

  ! prevent delayed threads from interfering
  !$omp master
  this % buf = TraceTransferBuffer(mesh, tr_un)
  !$omp end master
  !$omp barrier

  ! tr_un = 0
  call SetArray(tr_un, ZERO)

  ! get local trace ............................................................

  call GetLocalTrace(mesh, po, ne, nf, u, tr_un)

  ! start transfer of shared trace .............................................

  call this % buf % Transfer(mesh, tr_un, tag)

  ! apply boundary conditions ..................................................

  if (present(bv_u)) then
    call ApplyBoundaryConditions(mesh, po, nf, bv_u, tr_un)
  else
    call ExtendToBoundary(mesh, po, nf, tr_un)
  end if

end subroutine GetTrace_Start_X

!-------------------------------------------------------------------------------
!> Extract local trace

subroutine GetLocalTrace(mesh, po, ne, nf, u, tr_un)

  ! arguments ..................................................................

  class(MeshPartition),  intent(in)    :: mesh    !< mesh partition
  integer,               intent(in)    :: po      !< polynomial order
  integer,               intent(in)    :: ne      !< number of elements
  integer,               intent(in)    :: nf      !< number of faces
  real(RNP),             intent(in)    :: u       !< mesh variable
  real(RNP),             intent(inout) :: tr_un   !< normal trace of u

  dimension :: u     (0:po, 0:po, 0:po, ne, 3)
  dimension :: tr_un (0:po, 0:po, 1:2 , nf)

  ! internal variables .........................................................

  integer :: e, f(6), i, j, k

  ! extract local trace ........................................................

  !$acc data present(u, tr_un)
  !$acc parallel
  !$acc loop gang worker private(f)

  !$omp do private(e,f)
  do e = 1, mesh%ne

    f = mesh % element(e) % face % id

    ! west
    !$acc loop collapse(2) vector async(1)
    do k = 0, po
    do j = 0, po
      tr_un(j,k,2,f(1)) = -u(0,j,k,e,1)
    end do
    end do

    ! east
    !$acc loop collapse(2) vector async(2)
    do k = 0, po
    do j = 0, po
      tr_un(j,k,1,f(2)) = u(po,j,k,e,1)
    end do
    end do

    ! south
    !$acc loop collapse(2) vector async(3)
    do k = 0, po
    do i = 0, po
      tr_un(i,k,2,f(3)) = -u(i,0,k,e,2)
    end do
    end do

    ! north
    !$acc loop collapse(2) vector async(4)
    do k = 0, po
    do i = 0, po
      tr_un(i,k,1,f(4)) = u(i,po,k,e,2)
    end do
    end do

    ! bottom
    !$acc loop collapse(2) vector async(5)
    do j = 0, po
    do i = 0, po
      tr_un(i,j,2,f(5)) = -u(i,j,0,e,3)
    end do
    end do

    ! top
    !$acc loop collapse(2) vector async(6)
    do j = 0, po
    do i = 0, po
      tr_un(i,j,1,f(6)) = u(i,j,po,e,3)
    end do
    end do

    !$acc wait
  end do
  !$omp end do

  !$acc end parallel
  !$acc end data

end subroutine GetLocalTrace

!-------------------------------------------------------------------------------
!> Application of boundary conditions to trace data

subroutine ApplyBoundaryConditions(mesh, po, nf, bv_u, tr_un)

  ! arguments ..................................................................

  type(MeshPartition),    intent(in)    :: mesh  !< mesh partition
  integer,                intent(in)    :: po    !< polynomial order
  integer,                intent(in)    :: nf    !< number of faces
  type(BoundaryVariable), intent(in)    :: bv_u  !< boundary conditions & values
  real(RNP),              intent(inout) :: tr_un !< normal trace of u

  dimension :: bv_u (mesh%n_boundary)
  dimension :: tr_un (0:po, 0:po, 1:2, nf)

  ! internal variables .........................................................

  real(RNP), contiguous, pointer, save :: ub(:,:,:,:) => null()

  character :: bc(3)
  integer   :: b, d, f, i, k, n, o, p, q, s

  ! set boundary conditions ....................................................

  do b = 1, mesh%n_boundary
    associate(boundary => mesh%boundary(b))

      bc = bv_u(b) % bc(1:3)

      !$omp single
      ub(0:,0:,1:,1:) => bv_u(b) % Components(1,3)
      !$omp end single

      !$omp do private(k)
      do k = 1, boundary%nf

        f = boundary % face(k) % mesh_face % id  ! mesh face
        n = boundary % face(k) % orientation     ! normal orientation
        d = abs(n)                               ! normal direction
        s = sign(1,n)                            ! normal sense
        i = inner_side(n)                        ! inner side of face
        o = 3 - i                                ! outer side of face

        select case(bc(d))

        case('D') ! Dirichlet: enforce boundary condition
          do q = 0, po
          do p = 0, po
            tr_un(p,q,o,f) = tr_un(p,q,i,f) - s * 2 * ub(p,q,k,d)
          end do
          end do

        case('N') ! Neumann: extrapolate normal derivative
          do q = 0, po
          do p = 0, po
            tr_un(p,q,o,f) = -tr_un(p,q,i,f)
          end do
          end do

        end select
        end do

    end associate
  end do

end subroutine ApplyBoundaryConditions

!-------------------------------------------------------------------------------
!> Extrapolation to outer side

subroutine ExtendToBoundary(mesh, po, nf, tr_un)

  ! arguments ..................................................................

  type(MeshPartition), intent(in)    :: mesh  !< mesh partition
  integer,             intent(in)    :: po    !< polynomial order
  integer,             intent(in)    :: nf    !< number of faces
  real(RNP),           intent(inout) :: tr_un !< trace of u

  dimension :: tr_un (0:po, 0:po, 1:2, nf)

  ! internal variables .........................................................

  integer :: b, f, i, k, o, p, q

  ! extend to boundary .........................................................

  do b = 1, mesh%n_boundary
    associate(boundary => mesh%boundary(b))

      !$omp do private(f)
      do k = 1, boundary%nf

        f = boundary % face(k) % mesh_face % id
        i = inner_side(boundary % face(k) % orientation)
        o = 3 - i

        do q = 0, po
        do p = 0, po
          tr_un(p,q,o,f) = -tr_un(p,q,i,f)
        end do
        end do

      end do

    end associate
  end do

end subroutine ExtendToBoundary

!===============================================================================
! GetTrace_Finish

!-------------------------------------------------------------------------------
!> Merge remote traces and finish transfer

subroutine GetTrace_Finish(this, mesh, tr_un)
  class(NormalTraceOperator), intent(inout) :: this
  type(MeshPartition), intent(in)    :: mesh            !< mesh partition
  real(RNP),           intent(inout) :: tr_un(:,:,:,:)  !< normal trace of u

  call this % buf % Merge(mesh, tr_un, alpha=ZERO, beta=ONE)

end subroutine GetTrace_Finish

!===============================================================================

end module CART__Normal_Trace_Operator
