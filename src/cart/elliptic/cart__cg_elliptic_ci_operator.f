!> summary:  Cartesian elliptic CG-SE operator, constant isotropic, unstructured
!> author:   Joerg Stiller
!> date:     2018/10/22
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Cartesian elliptic CG-SE operator, constant isotropic, unstructured
!===============================================================================

module CART__CG_Elliptic_CI_Operator

  use Kind_Parameters,   only: RNP
  use Constants,         only: ZERO, ONE, HALF
  use Execution_Control, only: Error
  use Array_Assignments

  use TPO__Elliptic_3d_RLCI
  use CART__Mesh_Partition
  use CART__Assembly_Operator
  use CART__CG_Element_Operators
  use CART__CG_Elliptic_CI_BC

  implicit none
  private

  public :: EllipticOperator

contains

!-------------------------------------------------------------------------------
!> Elliptic operator, v = lambda M u + nu L u
!>
!> Unless the assembly is switched off, the result array must include the ghost
!> entries, i.e.
!>
!>     size(v,4) == mesh%ne + mesh%ng

subroutine EllipticOperator(mesh, eop, lambda, nu, bc, u, v, assemble)

  ! arguments ..................................................................

  class(MeshPartition),         intent(in)  :: mesh !< mesh partition
  class(CG_ElementOperators3D), intent(in)  :: eop  !< CG element operators

  real(RNP), intent(in)  :: lambda          !< Helmholtz parameter
  real(RNP), intent(in)  :: nu              !< diffusivity
  character, intent(in)  :: bc(:)           !< boundary conditions {P,D,N}
  real(RNP), intent(in)  :: u(0:,0:,0:,:)   !< approximate solution
  real(RNP), intent(out) :: v(0:,0:,0:,:)   !< result

  logical, optional, intent(in) :: assemble !< switch on/off assembly [T]

  ! local variables ............................................................

  type(AssemblyOperator), allocatable, save :: assembly_op

  logical :: assemble_
  integer :: ne, ng, np
  integer, parameter :: nl(3) = 1

  ! initialization .............................................................

  ne = mesh % ne
  ng = mesh % ng
  np = eop  % po + 1

  if (present(assemble)) then
    assemble_ = assemble
  else
    assemble_ = .true.
  end if

  ! apply element operators ....................................................

  call TPO_Elliptic_RLCI(eop%w, eop%L, lambda, nu, eop%dx, u, v)

  if (size(v,4) > ne) then
    call SetArray(v(:,:,:,ne+1:), ZERO)
  end if

  ! assembly ...................................................................

  if (assemble_) then

    if (size(v,4) /= ne + ng) then
      call Error( 'EllipticOperator' &
                , 'size(v,4) /= ne + ng' &
                , 'CART__CG_Elliptic_CI_Operator' )
    end if

    !$omp single
    allocate(assembly_op)
    !$omp end single

    call assembly_op % New(mesh, v, nl)
    call assembly_op % StartAssembly(mesh, v, 1000)
    call assembly_op % FinishAssembly(mesh, v)

    !$omp wait
    !$omp master
    deallocate(assembly_op)
    !$omp end master

  end if

  ! zero Dirichlet entries .....................................................

  call ZeroDirichletEntries(mesh, bc, v)

end subroutine EllipticOperator

!===============================================================================

end module CART__CG_Elliptic_CI_Operator
