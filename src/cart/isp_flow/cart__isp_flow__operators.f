!> summary:  DGM with implicit-explicit Euler pressure correction scheme
!> author:   Joerg Stiller
!> date:     2018/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### DGM with implicit-explicit Euler pressure correction scheme
!===============================================================================

module CART__ISP_Flow__Operators

  use Kind_Parameters,           only: RNP
  use Constants,                 only: ZERO, ONE
  use Embedded_Interpolation_3D, only: InterpolationOperator3D
  use Projection_Operator_3D,    only: ProjectionOperator3D
  use IP_Element_Operators_1D

  use ISP_Flow_Problem

  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Elliptic_Operator_IP
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Boundary_Values

  implicit none
  private

  public :: FlowOperators

  !-----------------------------------------------------------------------------
  !>

  type FlowOperators

    integer :: po_u = -1 !< polynomial order of variables except for pressure
    integer :: po_p = -1 !< polynomial order of pressure
    integer :: po_q = -1 !< polynomial order for quadrature of nonlinear terms

    type(MeshPartition), pointer :: mesh => null() !< local mesh partition

    real(RNP), allocatable :: x(:,:,:,:,:)         !< mesh points
    type(BoundaryVariable), allocatable :: bv_x(:) !< boundary points
    type(BoundaryVariable), allocatable :: bv_u(:) !< boundary values

    type(EllipticOperator3D_IP) :: eop_u !< element operators for u \ p
    type(EllipticOperator3D_IP) :: eop_p !< element operators for p
    type(IP_ElementOperators1D) :: eop_q !< element operators for nonlinear terms

    type(ProjectionOperator3D)    :: pop_up !< projection    from po_u to po_p
    type(InterpolationOperator3D) :: iop_uq !< interpolation from po_u to po_q
    type(InterpolationOperator3D) :: iop_pu !< interpolation from po_p to po_u

    type(PMG_Method3D) :: pmg_u !< p-MG/element operators for u \ p
    type(PMG_Method3D) :: pmg_p !< p-MG/element operators for p

    integer :: monitor_level = 0 !< no, essential or full monitoring {0,1,2}

  contains

    procedure :: New => New_FlowOperators

  end type FlowOperators

contains

!-------------------------------------------------------------------------------
!>

subroutine New_FlowOperators( this                       &
                            , problem                    &
                            , mesh                       &
                            , po_u, po_p, po_q, penalty  &
                            , pmg_u_opt                  &
                            , pmg_p_opt                  &
                            , monitor_level              &
                            )

  ! arguments ..................................................................

  class(FlowOperators),         intent(inout) :: this
  class(FlowProblem),           intent(in)    :: problem !< flow problem
  class(MeshPartition), target, intent(in)    :: mesh    !< mesh partition

  ! discretization parameters
  integer,   intent(in) :: po_u    !< order of variables except for pressure
  integer,   intent(in) :: po_p    !< order of pressure
  integer,   intent(in) :: po_q    !< order for quadrature of nonlinear terms
  real(RNP), intent(in) :: penalty !< penalty parameter op SIP method > 1

  ! multigrid parameters
  class(PMG_Options3D), intent(in) :: pmg_u_opt !< options for u
  class(PMG_Options3D), intent(in) :: pmg_p_opt !< options for p

  ! control
  integer, optional, intent(in) :: monitor_level !< monitor level [0] {0,1,2}

  ! local variables ............................................................

  integer :: b, l, l_top_u, l_top_p
  character(len=20) :: line = ''

  ! initialization of components ...............................................

  this % po_u = po_u
  this % po_p = po_p
  this % po_q = po_q

  this % mesh => mesh
  call mesh % GetPoints(po_u, 'GLL', this % x)
  allocate(this % bv_x(mesh%n_boundary))
  allocate(this % bv_u(mesh%n_boundary))
  call GetBoundaryPoints(mesh, this % x, this % bv_x)
  do b = 1, mesh%n_boundary
    call this % bv_u(b) % New(mesh, po_u, b, problem%bc(b,:))
  end do

!################################ HIER WEITER !################################!
!### lambda = ONE -- dt noch unbekannt ?
!### nu -- was tun im variablen Fall?
!### bc = problem % bc ?
!### schwarz_opt  --  Schwarz Options in PMG Options integrieren ?
!###
! call this % pmg_u % New(pmg_u_opt, mesh, penalty)
  this % pmg_u = PMG_Method3D(mesh, lambda, nu, bc, ip_opt, schwarz_opt, pmg_u_opt )

  associate(level => this % pmg_u % level)
    this % eop_u = level(ubound(level,1)) % eop
  end associate

  call this % pmg_p % New(pmg_p_opt, mesh, penalty)

  associate(level => this % pmg_p % level)
    this % eop_p = level(ubound(level,1)) % eop
  end associate

  if (po_p /= po_u) then
    this % pop_up = ProjectionOperator3D(    this%eop_u,   &
                                             this%eop_p%x, &
                                             this%eop_p%w, &
                                             mesh%dx       )
    this % iop_pu = InterpolationOperator3D( this%eop_p,   &
                                             this%eop_u%x  )
  end if

  if (po_q /= po_u) then
    this % eop_q  = IP_ElementOperators1D( IP_ElementOptions1D(eop_u, po_q) )
    this % iop_uq = InterpolationOperator3D( this%eop_u, this%eop_q%x )
  end if

  if (present(monitor_level)) this % monitor_level = monitor_level

  ! control output .............................................................

  !$omp single
  if (this%monitor_level > 0 .and. mesh%part == 0) then
    associate( pl_u => this % pmg_u % level % po &
             , pl_p => this % pmg_p % level % po )

      write(*,'(/,2X,A)') 'Polynomial levels:'
      write(*,'(A5,A5,A5)') 'l', 'u', 'p'
      l_top_u = ubound(pl_u,1)
      l_top_p = ubound(pl_p,1)
      do l = 1, max(l_top_u, l_top_p)
        write(line,'(I5)') l - 1
        if (l <= l_top_u) write(line( 6:),'(I5)') pl_u(l)
        if (l <= l_top_p) write(line(11:),'(I5)') pl_p(l)
        write(*,'(A)') line
        line = ''
      end do
      write(*,*)

    end associate
  end if
  !$omp end single

 end subroutine New_FlowOperators

!-------------------------------------------------------------------------------
!>

subroutine SetPolynomialLevels(po, rc, pl)
  integer,              intent(in)  :: po    !< polynomial order of top level
  real(RNP),            intent(in)  :: rc    !< coarsening ratio > 1
  integer, allocatable, intent(out) :: pl(:) !< polynomial order of top level

  integer, allocatable :: p(:)
  integer :: i, l

  allocate(p(0:po-1), source = po)

  l = 0
  do
    l    = l + 1
    p(l) = max( min( nint( p(l-1)/rc), p(l-1)-1 ), 1 )
    if (p(l) == 1) exit
  end do

  allocate(pl(0:l))
  do i = 0, l
    pl(i) = p(l-i)
  end do

end subroutine SetPolynomialLevels

!===============================================================================

end module CART__ISP_Flow__Operators


