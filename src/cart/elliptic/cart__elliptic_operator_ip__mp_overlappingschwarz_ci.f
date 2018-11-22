
submodule(CART__Elliptic_Operator_IP) MP_OverlappingSchwarz_CI
  use Execution_Control, only: Error

contains

!-------------------------------------------------------------------------------
!> Element-centered overlapping Schwarz method with constant coefficients

module subroutine OverlappingSchwarz_CI( this, lambda, nu, bc, u, f, &
                                         i_max, r_red, r_max, ni     )

  class(EllipticOperator3D_IP), intent(in) :: this
  real(RNP), intent(in)    :: lambda         !< Helmholtz parameter
  real(RNP), intent(in)    :: nu             !< diffusivity
  character, intent(in)    :: bc(:)          !< BC types {P,D,N}
  real(RNP), intent(inout) :: u(0:,0:,0:,:)  !< approximate solution
  real(RNP), intent(in)    :: f(0:,0:,0:,:)  !< right hand side
  integer,   intent(in)    :: i_max          !< max num iterations
  real(RNP), optional, intent(in)  :: r_red  !< min residual reduction
  real(RNP), optional, intent(in)  :: r_max  !< max admissible residual
  integer,   optional, intent(out) :: ni     !< exec num iterations

  ! local variables ............................................................

  !!! save attribute serves to enforce sharing with OpenMP !!!

  real(RNP), allocatable, save :: r (:,:,:,:)  ! residual, including ghosts
  real(RNP), allocatable, save :: u_s(:,:,:,:) ! solution of subsystems
  real(RNP), allocatable, save :: f_s(:,:,:,:) ! RHS of subsystems

  type(ElementTransferBuffer), allocatable, save :: buf_r
  type(ElementTransferBuffer), allocatable, save :: buf_u_s

! procedure(TPO_Schwarz_Proc),     pointer :: SchwarzGenOP
  procedure(TPO_Schwarz_Iso_Proc), pointer :: SchwarzIsoOP

  real(RNP), save :: rr_term
  logical  , save :: converged

  real(RNP) :: rr
  integer :: ne, ng, np(3), no(3), ns(3), nc
  integer :: i

  ! initialization .............................................................

  ! dimensions
  ne = mesh % ne              ! number of local elements
  ng = mesh % ng              ! number of ghost elements
  np = this % eop % po + 1    ! number of element points per direction
  no = this % no              ! overlapped point layers
  ns = np + 2*no              ! subdomain extensions
  nc = size(boundary_configuration, 2)

  ! workspace
  !$omp single
  allocate( r   (np(1), np(2), np(3), ne+ng) )
  allocate( u_s (ns(1), ns(2), ns(3), ne+ng) )
  allocate( f_s (ns(1), ns(2), ns(3), ne   ) )
  allocate( buf_r )
  allocate( buf_u_s )
  !$omp end single

  call AssignScalar(r  , ZERO)
  call AssignScalar(u_s, ZERO)

  ! transfer buffers
  call buf_r   % New(mesh, r,  no)
  call buf_u_s % New(mesh, u_s, no)

  ! subdomain operator
  if (this % isotropic) then
    call TPO_Schwarz_Iso_Assign(ns(1), SchwarzIsoOP)
  else
!   call TPO_Schwarz_Assign(ns(1), ns(2), ns(3), SchwarzGenOP)
    call Error('Iteration_C', &
               'Anisotropic Schwarz operator not implemented yet', &
               'CART__DG_Elliptic_CI_Schwarz')
  end if

  ! termination condition
  !$omp single
  if (present(r_max)) then
    rr_term = max(ZERO, r_max)**2
  else
    rr_term = ZERO
  end if
  !$omp end single

  ! Schwarz iterations .........................................................

  do i = 1, i_max

    call EllipticResidual(mesh, this%eop, lambda, nu, bc, u, f, r(:,:,:,:ne))

    ! termination check
    if (present(r_red)) then
      rr = ScalarProduct(r(:,:,:,:ne), r(:,:,:,:ne), mesh%comm)
      !$omp master
      if (mesh%part == 0) then
        if (i == 1) then
          rr_term  = max(rr_term, max(ZERO, sqrt(rr) * r_red)**2)
        end if
        converged = rr <= rr_term
      end if
      call XMPI_Bcast(converged, root=0, comm=mesh%comm)
      !$omp end master
      !$omp barrier
    end if
    if (converged) exit

    call RestrictToSubdomains(mesh, this, buf_r, r, f_s)

    if (this % isotropic) then
      call SchwarzIsoOP( ns(1), nc, ne,                   &
                         this % S1,                       &
                         this % V1, this % V2, this % V3, &
                         this % W1,                       &
                         ec, mesh%dx, lambda,             &
                         nu_s, f_s, u_s                   )

    else
!     call SchwarzGenOP( ns(1), ns(2), ns(3), nc, ne,     &
!                        this % S1, this % S2, this % S3, &
!                        this % D1, this % D2, this % D3, &
!                        this % W1, this % W2, this % W3, &
!                        ec, f_s, u_s                       )
    end if

    call MergeFromSubdomains(mesh, this, buf_u_s, u_s, u)

  end do

  ! clean-up ...................................................................

  !$omp single
  deallocate(r, u_s, f_s, nu_s, ec)
  deallocate(buf_r, buf_u_s)
  !$omp end single


end subroutine OverlappingSchwarz_CI

!===============================================================================

end submodule MP_OverlappingSchwarz_CI
