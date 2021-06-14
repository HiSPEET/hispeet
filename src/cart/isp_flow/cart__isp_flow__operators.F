!> summary:  Operators for semi-implicit incompressible flow solvers
!> author:   Joerg Stiller
!> date:     2018/04/10
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CART__ISP_Flow__Operators

  use Kind_Parameters, only: RNP
  use Constants,       only: ZERO, ONE
  use XMPI

  use Standard_Operators__1D
  use IP_Element_Operators__1D
  use Embedded_Interpolation__3D
  use Projection_Operator__3D

  use ISP_Flow_Problem

  use CART__Mesh_Partition
  use CART__Boundary_Variable
  use CART__Elliptic_Operator_IP
  use CART__Elliptic_PMG
  use CART__ISP_Flow__Boundary_Values

  implicit none
  private

  public :: FlowOperators
  public :: FlowOpControl

  !-----------------------------------------------------------------------------
  !>

  type FlowOpControl

    real(RNP) :: chi           = -2.0    !< factor of ∇[ν∇·v] such that with
                                         !! constant viscosity
                                         !!   χ = -2.   ⟺ ∇·τ = -ν∇×∇×v
                                         !!   χ = -1.   ⟺ ∇·τ = ν∇²v
                                         !!   χ = -2./3 ⟺ ∇·τ = ν(∇²v + ∇∇·v/3)

    integer   :: div_penalty   =  0      !< divergence penalty type
                                         !!   0: no div penalty
                                         !!   1: penalty scales with viscosity
                                         !!   2: penalty scales with (Δx/P)²/Δt
    real(RNP) :: div_coeff     =  1      !< divergence penalty coefficient
    real(RNP) :: div_r_red     =  1E-10  !< min divergence residual reduction
    integer   :: div_i_max     =  10     !< max num div correction iterations
    logical   :: div_final     = .true.  !< perform final projection step

    integer   :: monitor       =  0      !< no, basic or full monitoring {0,1,2}

    logical   :: svv           = .false. !< switch for the SVV model

  contains

    procedure ::  Bcast => Bcast_FlowOpControl

  end type FlowOpControl

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

    type(StandardOperators_1D) :: eop_u !< element operators for u \ p
    type(StandardOperators_1D) :: eop_p !< element operators for p
    type(StandardOperators_1D) :: eop_q !< element operators for nonlinear terms

    type(ProjectionOperator_3D)    :: pop_up !< projection    from po_u to po_p
    type(EmbeddedInterpolation_3D) :: iop_up !< interpolation from po_u to po_p
    type(EmbeddedInterpolation_3D) :: iop_uq !< interpolation from po_u to po_q
    type(EmbeddedInterpolation_3D) :: iop_pu !< interpolation from po_p to po_u

    type(PMG_Method3D) :: pmg_u !< p-MG/element operators for u \ p
    type(PMG_Method3D) :: pmg_p !< p-MG/element operators for p

    type(FlowOpControl) :: control !< control parameters

  contains

    procedure :: Init_FlowOperators

  end type FlowOperators

  ! constructor interface
  interface FlowOperators
    module procedure New_FlowOperators
  end interface

contains

  !-----------------------------------------------------------------------------
  !>

  subroutine Bcast_FlowOpControl(this, root, comm)
    class(FlowOpControl), intent(inout) :: this
    integer,              intent(in)    :: root !< rank of broadcast root
    type(MPI_Comm),       intent(in)    :: comm !< MPI communicator

    type(MPI_Request) :: request(8)
    type(MPI_Status)  :: stat(size(request))
    integer :: n

    n = 1
    call XMPI_Ibcast( this % chi        , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % div_penalty, root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % div_coeff  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % div_r_red  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % div_i_max  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % div_final  , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % monitor    , root, comm, request(n) );  n = n + 1
    call XMPI_Ibcast( this % svv        , root, comm, request(n) )

    call MPI_Waitall(n, request, stat)

  end subroutine Bcast_FlowOpControl

  !-----------------------------------------------------------------------------
  !> Constructor

  function New_FlowOperators( problem                    &
                            , mesh                       &
                            , po_u, po_p, po_q, penalty  &
                            , pmg_u_opt                  &
                            , pmg_p_opt                  &
                            , control                    &
                            ) result(this)

    type(FlowOperators) :: this

    class(FlowProblem),           intent(in) :: problem !< flow problem
    class(MeshPartition), target, intent(in) :: mesh    !< mesh partition

    integer,   intent(in) :: po_u    !< order of variables except for pressure
    integer,   intent(in) :: po_p    !< order of pressure
    integer,   intent(in) :: po_q    !< order for quadrature of nonlinear terms
    real(RNP), intent(in) :: penalty !< penalty parameter op SIP method > 1

    class(PMG_Options3D), intent(in) :: pmg_u_opt !< PMG options for u
    class(PMG_Options3D), intent(in) :: pmg_p_opt !< PMG options for p
    class(FlowOpControl), intent(in) :: control   !< control parameters

    call Init_FlowOperators( this, problem, mesh, po_u, po_p, po_q, penalty, &
                             pmg_u_opt, pmg_p_opt, control                   )

  end function New_FlowOperators

  !-----------------------------------------------------------------------------
  !> Initialization of flow operators

  subroutine Init_FlowOperators( this                       &
                               , problem                    &
                               , mesh                       &
                               , po_u, po_p, po_q, penalty  &
                               , pmg_u_opt                  &
                               , pmg_p_opt                  &
                               , control                    &
                               )

    ! arguments ................................................................

    class(FlowOperators),         intent(inout) :: this
    class(FlowProblem),           intent(in)    :: problem !< flow problem
    class(MeshPartition), target, intent(in)    :: mesh    !< mesh partition

    ! discretization parameters
    integer,   intent(in) :: po_u    !< order of variables except for pressure
    integer,   intent(in) :: po_p    !< order of pressure ≤ po_u
    integer,   intent(in) :: po_q    !< order for quadrature of nonlinear terms
    real(RNP), intent(in) :: penalty !< penalty parameter op SIP method > 1

    ! multigrid parameters
    class(PMG_Options3D), intent(in) :: pmg_u_opt !< options for u
    class(PMG_Options3D), intent(in) :: pmg_p_opt !< options for p

    ! control
    class(FlowOpControl), intent(in) :: control   !< control parameters

    ! local variables ..........................................................

    integer :: b, l, l_top_u, l_top_p
    character(len=20) :: line = ''

    ! initialization of components .............................................

    this % po_u = po_u
    this % po_p = min(po_p, po_u)
    this % po_q = po_q

    this % eop_u = StandardOperators_1D(this % po_u, svv = control%svv)
    this % eop_p = StandardOperators_1D(this % po_p)
    if (po_q /= po_u) then
      this % eop_q = StandardOperators_1D(po_q, no_vdm = .true.)
    end if

    this % mesh => mesh
    call mesh % GetPoints(po_u, 'L', this % x)
    allocate(this % bv_x(mesh%n_boundary))
    allocate(this % bv_u(mesh%n_boundary))
    call GetBoundaryPoints(mesh, this % x, this % bv_x)
    do b = 1, mesh%n_boundary
      this % bv_u(b) = BoundaryVariable(mesh, po_u, b, problem%bc(b,:))
    end do

    this%pmg_u = PMG_Method3D( mesh,                                       &
                               IP_ElementOptions_1D( po      = this%po_u,   &
                                                    penalty = penalty,     &
                                                    svv     = control%svv  &
                                                  ),                       &
                               pmg_u_opt                                   )

    this%pmg_p = PMG_Method3D( mesh,                                       &
                               IP_ElementOptions_1D( po      = this%po_p,   &
                                                    penalty = penalty      &
                                                  ),                       &
                               pmg_p_opt                                   )

    call this % pmg_p % SetProblem(lambda=ZERO, nu=ONE, bc=problem%bc(:,4))

    if (this % po_p /= this % po_u) then
      this % pop_up = ProjectionOperator_3D( this%eop_p,   &
                                            this%eop_u%x, &
                                            this%eop_u%w, &
                                            mesh%dx       )

      this % iop_up = EmbeddedInterpolation_3D( this%eop_u, this%eop_p%x )
      this % iop_pu = EmbeddedInterpolation_3D( this%eop_p, this%eop_u%x )
    end if

    if (po_q /= po_u) then
      this % iop_uq = EmbeddedInterpolation_3D( this%eop_u, this%eop_q%x )
    end if

    this % control = control

    ! control output ...........................................................

    !$omp single
    if (control%monitor > 0 .and. mesh%part == 0) then
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

   end subroutine Init_FlowOperators

  !=============================================================================

end module CART__ISP_Flow__Operators
