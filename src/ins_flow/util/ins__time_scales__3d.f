!> summary:  Incompressible Navier-Stokes time scales
!> author:   Joerg Stiller
!> date:     2023/02/19
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module INS__Time_Scales__3D
  use Kind_Parameters, only: RDP, RNP
  use Constants,       only: ONE, HALF, THIRD
  use Eigenproblems
  use XMPI
  use INS__Problem__3D
  use INS__Operator__3D

  implicit none
  private

  public :: INS_TimeScales_3D

  !-----------------------------------------------------------------------------
  !> Time scales based on 1D model problems and mean point spacing, respectively

  type INS_TimeScales_3D
    real(RNP) :: tau_conv_ve = -1 !< convective, based on eigenvalues  and v
    real(RNP) :: tau_conv_vm = -1 !< convective, based on mean spacing and v
    real(RNP) :: tau_conv_re = -1 !< convective, based on eigenvalues  and v_ref
    real(RNP) :: tau_conv_rm = -1 !< convective, based on mean spacing and v_ref
    real(RNP) :: tau_diff_re = -1 !< diffusive,  based on eigenvalues  and ν_ref
    real(RNP) :: tau_diff_rm = -1 !< diffusive,  based on mean spacing and ν_ref
  contains
    procedure :: Evaluate
  end type INS_TimeScales_3D

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of time scales

  subroutine Evaluate(this, problem, ins_op, u)
    class(INS_TimeScales_3D), intent(out) :: this
    class(INS_Problem_3D),    intent(in)  :: problem      !< INS problem
    class(INS_Operator_3D),   intent(in)  :: ins_op       !< INS operator
    real(RNP), contiguous,    intent(in)  :: u(:,:,:,:,:) !< flow variables

    real(RNP), parameter :: eps = epsilon(ONE)

    real(RNP), save :: tau_conv_v_loc = huge(ONE)
    real(RNP), save :: tau_conv_r_loc = huge(ONE)
    real(RNP), save :: tau_diff_r_loc = huge(ONE)

    complex(RDP), allocatable :: lmb_conv(:)
    real(RDP),    allocatable :: lmb_diff(:), A_conv(:,:), A_diff(:,:)

    real(RNP) :: lmb_conv_max, tau_conv_r, tau_conv_v
    real(RNP) :: lmb_diff_max, tau_diff_r
    real(RNP) :: h1, h2, h3, he, hh, vv
    integer   :: e, i, j, k, np, po

    associate( mesh => ins_op % mesh  &
             , eop  => ins_op % eop_u )

      if (mesh % part < 0) return

      po = eop % po
      np = po + 1

      ! scales neglecting polynomial order .....................................

      !$omp do reduction(min: tau_conv_v_loc, tau_conv_r_loc, tau_diff_r_loc)
      do e = 1, mesh % n_elem
        associate(x_c => mesh % element(e) % geometry % x_c)
          h1 = 2 * sqrt(x_c(1,1)**2 + x_c(1,2)**2 + x_c(1,3)**2)  ! 2|∂x/∂ξ|
          h2 = 2 * sqrt(x_c(2,1)**2 + x_c(2,2)**2 + x_c(2,3)**2)  ! 2|∂x/∂η|
          h3 = 2 * sqrt(x_c(3,1)**2 + x_c(3,2)**2 + x_c(3,3)**2)  ! 2|∂x/∂ζ|
          he = (h1 * h2 * h3)**THIRD
          hh = h1**2 + h2**2 + h3**2
          vv = 0
          do k = 1, np
          do j = 1, np
          do i = 1, np
            vv = max(vv, u(i,j,k,e,1)**2 + u(i,j,k,e,2)**2 + u(i,j,k,e,3)**2)
          end do
          end do
          end do
          tau_conv_v_loc = min(tau_conv_v_loc, he / max(eps, sqrt(vv))         )
          tau_conv_r_loc = min(tau_conv_r_loc, he / max(eps, problem % v_ref)  )
          tau_diff_r_loc = min(tau_diff_r_loc, hh / max(eps, problem % nu_ref) )
        end associate
      end do

      !$omp master
      call XMPI_Reduce(tau_conv_v_loc, tau_conv_v, MPI_MIN, 0, mesh%comm_parts)
      call XMPI_Reduce(tau_conv_r_loc, tau_conv_r, MPI_MIN, 0, mesh%comm_parts)
      call XMPI_Reduce(tau_diff_r_loc, tau_diff_r, MPI_MIN, 0, mesh%comm_parts)
      !$omp end master

      ! scales adjusted to polynomial order .......................................

      !$omp master

      if (mesh % part == 0) then

        ! eigenvalues of 1D convection and diffusion problems
        select case(po)
        case(0:1)
          lmb_conv_max = HALF
          lmb_diff_max = ONE
        case default
          ! convection problem with unit velocity in standard element
          allocate(A_conv(po,po), lmb_conv(po))
          A_conv = eop % D(1:,1:)
          call SolveNonsymmetricEigenproblem(A_conv, lmb_conv)
          lmb_conv_max = real(maxval(abs(lmb_conv)), RNP)
          ! diffusion problem with Dirichlet conditions in standard element
          allocate(A_diff(po-1,po-1), lmb_diff(po-1))
          A_diff = eop % L(1:po-1,1:po-1)
          call SolveSymmetricEigenproblem(A_diff, lmb_diff)
          lmb_diff_max = real(maxval(abs(lmb_diff)), RNP)
        end select

        ! resulting time scales
        this % tau_conv_ve = tau_conv_v / (2 * lmb_conv_max)
        this % tau_conv_vm = tau_conv_v / max(po, 1)
        this % tau_conv_re = tau_conv_r / (2 * lmb_conv_max)
        this % tau_conv_rm = tau_conv_r / max(po, 1)
        this % tau_diff_re = tau_diff_r / (4 * lmb_diff_max)
        this % tau_diff_rm = tau_diff_r / (4 * max(po, 1)**2)

      end if

      ! globalize
      call XMPI_Bcast(this % tau_conv_ve, 0, mesh % comm_parts)
      call XMPI_Bcast(this % tau_conv_vm, 0, mesh % comm_parts)
      call XMPI_Bcast(this % tau_conv_re, 0, mesh % comm_parts)
      call XMPI_Bcast(this % tau_conv_rm, 0, mesh % comm_parts)
      call XMPI_Bcast(this % tau_diff_re, 0, mesh % comm_parts)
      call XMPI_Bcast(this % tau_diff_rm, 0, mesh % comm_parts)

      !$omp end master

      ! finalization ...........................................................

      !$omp master
      tau_conv_v_loc = huge(ONE)
      tau_conv_r_loc = huge(ONE)
      tau_diff_r_loc = huge(ONE)
      !$omp end master

    end associate

  end subroutine Evaluate

  !=============================================================================

end module INS__Time_Scales__3D
