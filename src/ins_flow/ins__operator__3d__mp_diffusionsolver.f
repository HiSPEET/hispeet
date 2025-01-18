!> summary:  Incompressible Navier-Stokes DG-SEM diffusion solver
!> author:   Joerg Stiller
!> date:     2022/09/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
!>   - add standby option?
!>   - ensure consistent handling of empty partitions
!===============================================================================

submodule(INS__Operator__3D) MP_DiffusionSolver
  use Logging_Levels , only: log_level_inner_iteration
  use Array_Reductions
  use TPO__Average__3D
  use TPO__Schwarz__3D
  use Element_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !>  IPCG Diffusion solver with Schwarz preconditioner

  module subroutine DiffusionSolver &
      (this, tau, mu, nu, f, bv_u, v, i_max, r_red, r_max, ni)

    class(INS_Operator_3D), intent(in) :: this
    !< incompressible Navier-Stokes operator

    real(RNP), intent(in) :: tau
    !< τ, effective time step width

    real(RNP), contiguous, optional, intent(in) :: mu(:,:,:,:)
    !< kinematic bulk viscosity, μ(np,np,np,ne,3)

    real(RNP), contiguous, optional, intent(in) :: nu(:,:,:,:)
    !< kinematic shear viscosity, ν(np,np,np,ne,3)

    real(RNP), contiguous, intent(in) :: f(:,:,:,:,:)
    !< sources, f(np,np,np,ne,3)

    class(BoundaryVariable_3D), intent(in) :: bv_u(:)
    !< boundary values at final time t
    !!   - Γᴰ :  vᵇ    in components 1:3
    !!   - Γᴼ :  τ_nn  in component    4

    real(RNP), contiguous, intent(inout) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    integer,               intent(in)    :: i_max  !< max num iterations
    real(RNP),   optional, intent(in)    :: r_red  !< min residual reduction
    real(RNP),   optional, intent(in)    :: r_max  !< max admissible residual
    integer,     optional, intent(out)   :: ni     !< executed num iterations

    ! internal variables .......................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: p, q, r, s, z
    real(RNP), dimension(:),         allocatable, save :: nu_avg
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, pq, rr
    logical   :: check_convergence
    integer   :: i, i_max_

    if (i_max < 1) return

    ! skip empty partition
    if (this % mesh % part < 0) return

    associate(mesh => this % mesh)

      ! initialization .........................................................

      check_convergence = log_level_inner_iteration > 0
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      !$omp master
      allocate(r, mold = f)
      allocate(p, mold = f)
      allocate(q, mold = f)
      allocate(s, mold = f)
      allocate(z, mold = f)
      allocate(nu_avg(mesh%n_elem))
      !$omp end master
      !$omp barrier

      ! initial residual
      call this % GetDiffusionResidual(tau, mu, nu, f, bv_u, v, r)

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r, r, mesh%comm_parts)
        !$omp master
        if (present(r_red)) then
          rr_term  = max(ZERO, sqrt(rr) * r_red)**2
        else
          rr_term = 0
        end if
        if (present(r_max)) then
          rr_term = max(rr_term, max(ZERO, r_max)**2)
        end if
        converged = rr <= rr_term
        call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
        if (log_level_inner_iteration > 1 .and. mesh%part == 0) then
          print '(A,T25,A,I5,A,ES12.5)', &
                '#INS:DiffusionSolver','>>>  i  =',0,',  |r| =', sqrt(rr)
        end if
        !$omp end master
        !$omp barrier
      else
        !$omp master
        converged = .false.
        !$omp end master
        !$omp barrier
      end if

      if (converged) then
        i_max_ = 0
        i      = 0
      else
        i_max_ = i_max
      end if

      ! element-averaged viscosity .............................................

      if (present(nu)) then
        call TPO_Average(this%eop_u%w, nu, nu_avg)
      else
        call SetArray(nu_avg, this % nu_0)
      end if

      ! iteration ..............................................................

      do i = 1, i_max_

        call Schwarz_Preconditioner(this, tau, nu_avg, r, z, standby = i < i_max_)

        ! set/update search vector
        if (i == 1) then
          call SetArray(p, z, multi = .true.)                 ! p = z
        else
          call SetArray(q, r, multi = .true.)                 ! q = r
          call MergeArrays(ONE, q, -ONE, s, multi = .true.)   ! q = r - s
          beta = ScalarProduct(q, z, mesh%comm_parts) / delta
          call MergeArrays(beta, p, ONE, z, multi = .true.)   ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r, multi = .true.)

        ! correction
        call this % ApplyDiffusionOperator(tau, mu, nu, p, q)  ! q = Ap

        delta = ScalarProduct(r, z, mesh%comm_parts)
        pq    = ScalarProduct(p, q, mesh%comm_parts)
        alpha = delta / pq
        call MergeArrays(ONE, v, alpha, p, multi = .true.)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          call this % GetDiffusionResidual(tau, mu, nu, f, bv_u, v, r)
        else
          call MergeArrays(ONE, r, -alpha, q, multi = .true.)
        end if

        if (check_convergence) then
          rr = ScalarProduct(r, r, mesh%comm_parts)
          !$omp master
          converged = rr <= rr_term
          call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
          !$omp end master
          !$omp barrier
        end if

        if (converged .or. i == i_max_) exit

        !$omp master
        if (log_level_inner_iteration > 1 .and. mesh%part == 0) then
          print '(A,T25,A,I5,A,ES12.5)', &
                '#INS:DiffusionSolver','>>>  i  =',i,',  |r| =', sqrt(rr)
        end if
        !$omp end master

      end do

      !$omp master
      if (log_level_inner_iteration > 0 .and. mesh%part == 0) then
        print '(A,T25,A,I5,A,ES12.5)', &
              '#INS:DiffusionSolver','>>>  ni =',i,',  |r| =', sqrt(rr)
      end if
      !$omp end master

      ! finalization ...........................................................

      if (present(ni)) ni = i

      !$omp master
      deallocate(p, q, r, s, z)
      deallocate(nu_avg)
      !$omp end master

    end associate

  end subroutine DiffusionSolver

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz preconditioner
  !>
  !> Performs a single Schwarz sweep for each velocity component. In order to
  !> accomodate different boundary conditions, an individual Schwarz operator
  !> is used for each component. However, for efficiency reasons, the overlap
  !> is assumed to be the same in all cases.

  subroutine Schwarz_Preconditioner(ins_op, tau, nu_avg, r, z, standby)

    class(INS_Operator_3D), intent(in)  :: ins_op       !< INS DG operator
    real(RNP),              intent(in)  :: tau          !< τ
    real(RNP), contiguous,  intent(in)  :: nu_avg(:)    !< ν element mean values
    real(RNP), contiguous,  intent(in)  :: r(:,:,:,:,:) !< residual
    real(RNP), contiguous,  intent(out) :: z(:,:,:,:,:) !< correction
    logical, intent(in) :: standby !< switch for reusing workspace

    ! local variables ..........................................................

    real(RNP), allocatable, save :: rg(:,:,:,:)    ! residual with ghost entries

    ! auxiliaries for the DP Schwarz preconditioner
    real(RDP) :: lambda_dp
    real(RDP), allocatable, save :: nu_dp(:)       ! subdomain viscosities
    real(RDP), allocatable, save :: zs_dp(:,:,:,:) ! solution of subsystems
    real(RDP), allocatable, save :: rs_dp(:,:,:,:) ! RHS of subsystems

    ! auxiliaries for the SP Schwarz preconditioner
    real(RSP) :: lambda_sp
    real(RSP), allocatable, save :: nu_sp(:)       ! subdomain viscosities
    real(RSP), allocatable, save :: zs_sp(:,:,:,:) ! solution of subsystems
    real(RSP), allocatable, save :: rs_sp(:,:,:,:) ! RHS of subsystems

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_rg
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_zs

    ! parameters saved for reuse
    integer :: np = -1  ! number of element points per direction
    integer :: ne = -1  ! number of elements
    integer :: ng = -1  ! number of ghosts
    integer :: no = -1  ! overlap
    integer :: wp = -1  ! working precision

    integer :: d, nl(3), ns
    logical :: reuse

    associate(mesh => ins_op % mesh, schwarz_v => ins_op % schwarz_v)

      ! initialization .........................................................

      !$omp master

      reuse = np == size(z,1)           .and. &
              ne == mesh % n_elem       .and. &
              ng == mesh % n_ghost      .and. &
              no == schwarz_v(1) % no   .and. &
              wp == schwarz_v(1) % wp

      if (.not. reuse) then

        if (allocated( rg     )) deallocate( rg     )
        if (allocated( nu_dp  )) deallocate( nu_dp  )
        if (allocated( zs_dp  )) deallocate( zs_dp  )
        if (allocated( rs_dp  )) deallocate( rs_dp  )
        if (allocated( nu_sp  )) deallocate( nu_sp  )
        if (allocated( zs_sp  )) deallocate( zs_sp  )
        if (allocated( rs_sp  )) deallocate( rs_sp  )
        if (allocated( buf_rg )) deallocate( buf_rg )
        if (allocated( buf_zs )) deallocate( buf_zs )

        np = size(z,1)
        ne = mesh % n_elem
        ng = mesh % n_ghost
        no = schwarz_v(1) % no
        wp = schwarz_v(1) % wp

        nl = no
        ns = np + 2*no

        allocate(rg(np, np, np, ne+ng), source = ZERO)
        buf_rg = ElementTransferBuffer_3D(mesh, rg, nl)

        if (wp == RSP) then
          allocate(nu_sp(ne))
          allocate(zs_sp(ns, ns, ns, ne+ng), source = 0E0)
          allocate(rs_sp(ns, ns, ns, ne)   , source = 0E0)
          buf_zs = ElementTransferBuffer_3D(mesh, zs_sp, nl)
        else
          allocate(nu_dp(ne))
          allocate(zs_dp(ns, ns, ns, ne+ng), source = 0D0)
          allocate(rs_dp(ns, ns, ns, ne),    source = 0D0)
          buf_zs = ElementTransferBuffer_3D(mesh, zs_dp, nl)
        end if

      end if

      !$omp end master

      ! coefficients
      select case(wp)
      case(RSP)
        lambda_sp = real(1/tau, RSP)
        !$omp workshare
        nu_sp = real(nu_avg, RSP)
        !$omp workshare nowait
      case default
        lambda_dp = real(1 / tau, RDP)
        !$omp workshare
        nu_dp = real(nu_avg, RDP)
        !$omp workshare nowait
      end select

      call SetArray(z, ZERO, multi = .true.)

      ! Schwarz sweeps .........................................................

      do d = 1, 3

        call SetArray(rg(:,:,:,:ne), r(:,:,:,:,d))

        select case(wp)

        case(RSP)
          call schwarz_v(d) % RestrictResidual(mesh, buf_rg, rg, rs_sp)
          call TPO_Schwarz( schwarz_v(d) % ops_sp % S  &
                          , schwarz_v(d) % ops_sp % V  &
                          , schwarz_v(d) % ops_sp % W  &
                          , schwarz_v(d) % ops_sp % g  &
                          , schwarz_v(d) % cfg         &
                          , lambda_sp                  &
                          , nu_sp                      &
                          , rs_sp                      &
                          , zs_sp                      )
          call schwarz_v(d) % MergeCorrections(mesh, buf_zs, zs_sp, z(:,:,:,:,d))

        case default
          call schwarz_v(d) % RestrictResidual(mesh, buf_rg, rg, rs_dp)
          call TPO_Schwarz( schwarz_v(d) % ops_dp % S  &
                          , schwarz_v(d) % ops_dp % V  &
                          , schwarz_v(d) % ops_dp % W  &
                          , schwarz_v(d) % ops_dp % g  &
                          , schwarz_v(d) % cfg         &
                          , lambda_dp                  &
                          , nu_dp                      &
                          , rs_dp                      &
                          , zs_dp                      )
          call schwarz_v(d) % MergeCorrections(mesh, buf_zs, zs_dp, z(:,:,:,:,d))

        end select

      end do

      ! cleanup ................................................................

      ! keep workspace in case of standby
      if (standby) return

      !$omp master

      if (allocated( rg     )) deallocate( rg     )
      if (allocated( nu_dp  )) deallocate( nu_dp  )
      if (allocated( zs_dp  )) deallocate( zs_dp  )
      if (allocated( rs_dp  )) deallocate( rs_dp  )
      if (allocated( nu_sp  )) deallocate( nu_sp  )
      if (allocated( zs_sp  )) deallocate( zs_sp  )
      if (allocated( rs_sp  )) deallocate( rs_sp  )
      if (allocated( buf_rg )) deallocate( buf_rg )
      if (allocated( buf_zs )) deallocate( buf_zs )

      np = -1
      ne = -1
      ng = -1
      no = -1
      wp = -1

      !$omp end master

    end associate

  end subroutine Schwarz_Preconditioner

  !=============================================================================

end submodule MP_DiffusionSolver
