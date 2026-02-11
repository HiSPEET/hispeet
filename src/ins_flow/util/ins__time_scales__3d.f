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
  !> Time scales based on 1D model problems

  type INS_TimeScales_3D
    real(RNP) :: tau_conv_v !< convective, based on local velocity
    real(RNP) :: tau_conv_r !< convective, based on reference velocity
    real(RNP) :: tau_diff_r !< diffusive,  based on reference viscosity
  contains
    procedure :: Evaluate
    procedure :: Globalize
  end type INS_TimeScales_3D

  !-----------------------------------------------------------------------------
  !> Max abs eigenvalues of 1D single element convection problem for P ≤ 32

  real(RNP), parameter ::                                        &
     LMB_CONV(2:32) = [ 1.00000000000000D+0, 1.66481452033766D+0 &
                      , 2.47159803597853D+0, 3.40883748057081D+0 &
                      , 4.47073758437598D+0, 5.65599485263299D+0 &
                      , 6.96732548733230D+0, 8.41089833350806D+0 &
                      , 9.99507344978968D+0, 1.17283555836581D+1 &
                      , 1.36174001000091D+1, 1.56661755899168D+1 &
                      , 1.78764016027800D+1, 2.02484698509996D+1 &
                      , 2.27821848530042D+1, 2.54771710732096D+1 &
                      , 2.83330470912398D+1, 3.13494832098709D+1 &
                      , 3.45262112697018D+1, 3.78630182222503D+1 &
                      , 4.13597359477187D+1, 4.50162316204434D+1 &
                      , 4.88323997653336D+1, 5.28081560475890D+1 &
                      , 5.69434325202662D+1, 6.12381740124714D+1 &
                      , 6.56923353834880D+1, 7.03058794266417D+1 &
                      , 7.50787752594481D+1, 8.00109970785380D+1 &
                      , 8.51025231895497D+1                      ]

  !-----------------------------------------------------------------------------
  !> Maximum eigenvalues of 1D single element diffusion problem for P ≤ 32

  real(RNP), parameter ::                                        &
     LMB_DIFF(2:32) = [ 2.00000000000000D+0, 7.50000000000000D+0 &
                      , 1.70344011421667D+1, 3.22249721603218D+1 &
                      , 5.60616739241431D+1, 9.27060546644337D+1 &
                      , 1.46811569105794D+2, 2.23347883992649D+2 &
                      , 3.27773246116899D+2, 4.66136982753496D+2 &
                      , 6.45103438870549D+2, 8.71952866452095D+2 &
                      , 1.15457888274874D+3, 1.50148640181074D+3 &
                      , 1.92179033173631D+3, 2.42521479671065D+3 &
                      , 3.02209267218295D+3, 3.72336530190282D+3 &
                      , 4.54058232137468D+3, 5.48590154484064D+3 &
                      , 6.57208889111166D+3, 7.81251833375077D+3 &
                      , 9.22117186688794D+3, 1.08126394813015D+4 &
                      , 1.26021191473786D+4, 1.46054168027795D+4 &
                      , 1.68389463433729D+4, 1.93197296164838D+4 &
                      , 2.20653964158007D+4, 2.50941844774896D+4 &
                      , 2.84249394771984D+4                      ]

contains

  !-----------------------------------------------------------------------------
  !> Evaluation of time scales
  !>
  !> If the communicator `comm` is specified, it will be used to globalize the
  !> results.

  subroutine Evaluate(this, ins_op, u, comm)
    class(INS_TimeScales_3D), intent(out) :: this
    class(INS_Operator_3D),   intent(in)  :: ins_op       !< INS operator
    real(RNP), contiguous,    intent(in)  :: u(:,:,:,:,:) !< flow variables
    type(MPI_Comm), optional, intent(in)  :: comm

    real(RNP), parameter :: eps = epsilon(ONE)

    real(RNP), save :: tau_conv_v = huge(ONE)
    real(RNP), save :: tau_conv_r = huge(ONE)
    real(RNP), save :: tau_diff_r = huge(ONE)

    real(RNP) :: lmb_conv_max, lmb_diff_max
    real(RNP) :: dx_c(3), he, hh, vv
    integer   :: e, i, j, k, np, po

    associate( problem => ins_op % problem &
             , mesh    => ins_op % mesh    &
             , eop     => ins_op % eop_u   )


      if (mesh % part < 0) then

        this % tau_conv_v = huge(ONE)
        this % tau_conv_r = huge(ONE)
        this % tau_diff_r = huge(ONE)

      else

        po = eop % po
        np = po + 1

        ! scales for linear unit element .......................................

        !$omp do reduction(min: tau_conv_v, tau_conv_r, tau_diff_r)
        do e = 1, mesh % n_elem_active
          call mesh % element(e) % GetCuboidDimensions(dx_c)
          he = product(dx_c)**THIRD
          hh = sum(dx_c**2)
          vv = 0
          do k = 1, np
          do j = 1, np
          do i = 1, np
            vv = max(vv, u(i,j,k,e,1)**2 + u(i,j,k,e,2)**2 + u(i,j,k,e,3)**2)
          end do
          end do
          end do
          tau_conv_v = min(tau_conv_v, he / max(eps, sqrt(vv))         )
          tau_conv_r = min(tau_conv_r, he / max(eps, problem % v_ref)  )
          tau_diff_r = min(tau_diff_r, hh / max(eps, problem % nu_ref) )
        end do

        ! scales adjusted to polynomial order ..................................

        !$omp master
        select case(po)

        case(0)
          lmb_conv_max = HALF
          lmb_diff_max = ONE

        case(1:32)
          lmb_conv_max = LMB_CONV(po)
          lmb_diff_max = LMB_DIFF(po)

        case default
          block

            complex(RDP), allocatable :: lmb_conv(:)
            real(RDP),    allocatable :: lmb_diff(:), A_conv(:,:), A_diff(:,:)

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

          end block
        end select
        !$omp end master
        !$omp barrier

        this % tau_conv_v = tau_conv_v / (2 * lmb_conv_max)
        this % tau_conv_r = tau_conv_r / (2 * lmb_conv_max)
        this % tau_diff_r = tau_diff_r / (4 * lmb_diff_max)

      end if

      ! globalization ..........................................................

      if (present(comm)) then
        call this % Globalize(comm)
      else
        !$omp barrier
      end if

      ! finalization ...........................................................

      !$omp master
      tau_conv_v = huge(ONE)
      tau_conv_r = huge(ONE)
      tau_diff_r = huge(ONE)
      !$omp end master

    end associate

  end subroutine Evaluate

  !-----------------------------------------------------------------------------
  !> Globalization

  subroutine Globalize(this, comm)
    class(INS_TimeScales_3D), intent(inout) :: this
    type(MPI_Comm),           intent(in)    :: comm

    real(RNP), save :: tau_conv_v_loc, tau_conv_v
    real(RNP), save :: tau_conv_r_loc, tau_conv_r
    real(RNP), save :: tau_diff_r_loc, tau_diff_r

    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    !$omp master

    tau_conv_v_loc = this % tau_conv_v
    tau_conv_r_loc = this % tau_conv_r
    tau_diff_r_loc = this % tau_diff_r

    call XMPI_Allreduce(tau_conv_v_loc, tau_conv_v, MPI_MIN, comm)
    call XMPI_Allreduce(tau_conv_r_loc, tau_conv_r, MPI_MIN, comm)
    call XMPI_Allreduce(tau_diff_r_loc, tau_diff_r, MPI_MIN, comm)

    !$omp end master
    !$omp barrier
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    this % tau_conv_v = tau_conv_v
    this % tau_conv_r = tau_conv_r
    this % tau_diff_r = tau_diff_r

  end subroutine Globalize

  !=============================================================================

end module INS__Time_Scales__3D
