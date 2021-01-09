!> summary:  Application of the IP/DG elliptic operator
!> author:   Joerg Stiller
!> date:     2018/11/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CART__Elliptic_Operator_IP) MP_Apply
  use CART__Mesh_Partition
  implicit none

  !> mapping of boundary face orientation to inner trace side
  integer, parameter :: inner_side(-3:3) = [ 2, 2, 2, 0, 1, 1, 1 ]

  interface

    !> Application of the operator with constant isotropic diffusivity
    module subroutine Apply_CI(this, u, v)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply_CI

    !> Application of the operator with variable isotropic diffusivity
    module subroutine Apply_VI(this, u, v)
      class(EllipticOperator3D_IP), intent(in) :: this
      real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
      real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result
    end subroutine Apply_VI

  end interface

contains

!-------------------------------------------------------------------------------
!> Application of the IP/DG elliptic operator

module subroutine Apply(this, u, v)
  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)  :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)  !< result

  if (allocated(this % nu_ci)) then
    call Apply_CI(this, u, v)
  else
    call Apply_VI(this, u, v)
  end if

end subroutine Apply

!===============================================================================
! Procedures shared between Apply_CI and Apply_VI

!-------------------------------------------------------------------------------
!> Modify boundary traces to yield correct contribution to the operator

subroutine ApplyBoundaryConditions(mesh, bc, tr_u, tr_qn)
  class(MeshPartition), intent(in)    :: mesh             !< mesh partition
  character,            intent(in)    :: bc(:)            !< boundary conditions
  real(RNP),            intent(inout) :: tr_u (0:,0:,:,:) !< trace of u
  real(RNP),            intent(inout) :: tr_qn(0:,0:,:,:) !< trace of ν du/dn

  integer :: b, f, i, l, o

  do b = 1, size(bc)
    associate(face => mesh % boundary(b) % face)

      select case(bc(b))

      case('D')
        ! Dirichlet
        !   – interior solution contributes twice to [u]
        !   - du/dn is mirrored from interior
        !$omp do
        do l = 1, size(face)
          f = face(l) % mesh_face % id          ! mesh face
          i = inner_side(face(l) % orientation) ! inner side
          o = 3 - i                             ! outer side
          tr_u (:,:,o,f) = -tr_u (:,:,i,f)
          tr_qn(:,:,o,f) = -tr_qn(:,:,i,f)
        end do

      case('N')
        ! Neumann: du/dn does not contribute, [u] = 0 due to extrapolation
        !$omp do
        do l = 1, size(face)
          f = face(l) % mesh_face % id          ! mesh face
          i = inner_side(face(l) % orientation) ! inner side
          o = 3 - i                             ! outer side
          tr_u (:,:,o,f) = tr_u (:,:,i,f)
          tr_qn(:,:,:,f) = 0
        end do

      end select
    end associate
  end do

end subroutine ApplyBoundaryConditions

!-------------------------------------------------------------------------------
!> Jump across faces
!>
!> Pass `normal=.true.` if the trace `tr_u` contains the normal components of a
!> vector variable

subroutine ComputeJumps(tr_u, J_u, normal)
  real(RNP), intent(in)  :: tr_u(0:,0:,:,:) !< trace of u
  real(RNP), intent(out) :: J_u (0:,0:,:)   !< [u]ᵢ
  logical, optional, intent(in) :: normal   !< switch for normal traces [F]

  integer   :: f

  if (present(normal)) then
    !$omp do
    do f = 1, size(J_u, 3)
      J_u(:,:,f) = tr_u(:,:,1,f) + tr_u(:,:,2,f)
    end do
  else
    !$omp do
    do f = 1, size(J_u, 3)
      J_u(:,:,f) = tr_u(:,:,1,f) - tr_u(:,:,2,f)
    end do
  end if

end subroutine ComputeJumps

!-------------------------------------------------------------------------------
!> Average over faces
!>
!> Pass `normal=.true.` if the trace `tr_u` contains the normal components of a
!> vector variable

subroutine ComputeAverages(tr_u, A_u, normal)
  real(RNP), intent(in)  :: tr_u(0:,0:,:,:) !< trace of u
  real(RNP), intent(out) :: A_u (0:,0:,:)   !< {u}
  logical, optional, intent(in) :: normal   !< switch for normal traces [F]

  integer   :: f

  if (present(normal)) then
    !$omp do
    do f = 1, size(A_u, 3)
      A_u(:,:,f) = HALF * (tr_u(:,:,1,f) - tr_u(:,:,2,f))
    end do
  else
    !$omp do
    do f = 1, size(A_u, 3)
      A_u(:,:,f) = HALF * (tr_u(:,:,1,f) + tr_u(:,:,2,f))
    end do
  end if

end subroutine ComputeAverages

!===============================================================================

end submodule MP_Apply
