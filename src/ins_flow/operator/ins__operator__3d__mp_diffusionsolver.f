!> summary:  Incompressible Navier-Stokes DG-SEM diffusion solver
!> author:   Joerg Stiller
!> date:     2022/09/28
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!> @todo
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

  module subroutine DiffusionSolver(this, tau, mu, nu, f, bv, v, precon, ni)

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

    class(BoundaryVariable_3D), intent(in) :: bv(:)
    !< boundary values at final time t
    !!   - Γᴰ :  vᵇ    in components 1:3
    !!   - Γᴼ :  τ_nn  in component    4

    real(RNP), contiguous, intent(inout) :: v(:,:,:,:,:)
    !< velocity, v(np,np,np,ne,3)

    logical, optional, intent(in)  :: precon !< switch preconditioner mode
    integer, optional, intent(out) :: ni     !< executed num iterations

    ! internal variables .......................................................

    real(RNP), dimension(:,:,:,:,:), allocatable, save :: p, q, r, s, z
    real(RNP), dimension(:),         allocatable, save :: nu_avg
    real(RNP), save :: rr_term
    logical  , save :: converged

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, pq, rr
    real(RNP) :: r_max, r_red
    integer   :: i, i_max
    logical   :: check_convergence

    associate(mesh => this % mesh)

      ! initialization .........................................................

      i_max = this % i_max_v
      r_red = this % r_red
      r_max = this % r_max
      if (present(precon)) then
        if (precon) i_max = this % k_pre_v
      end if

      if (i_max < 1 .or. mesh % part < 0) then
        if (present(ni)) ni = 0
        return
      end if

      check_convergence = r_red > 0 .or. &
                          r_max > 0 .or. &
                          log_level_inner_iteration > 0

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
      call this % GetDiffusionResidual(tau, mu, nu, f, bv, v, r)

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r, r, mesh%comm_parts)
        !$omp master
        rr_term  = max(ZERO, sqrt(rr) * r_red, r_max)**2
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
        i_max = 0
        i     = 0
      else
        i_max = i_max
      end if

      ! element-averaged viscosity .............................................

      if (this % HasVariableViscosity()) then
        call TPO_Average(this%eop_u%w, nu, nu_avg)
      else
        call SetArray(nu_avg, this % nu_0)
      end if

      ! iteration ..............................................................

      do i = 1, i_max

        select case(this % diffusion_solver)
        case('DPCG')
          call Diagonal_Preconditioner(this, tau, r, z, standby = i < i_max)
        case('SPCG')
          call Schwarz_Preconditioner(this, tau, nu_avg, r, z, standby = i < i_max)
        case default
          call Error( 'DiffusionSolver'              &
                    , 'Preconditioner "'             &
                      // trim(this%diffusion_solver) &
                      // '" not supported'           &
                    , 'INS__Operator__3D'            )
        end select

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
          call this % GetDiffusionResidual(tau, mu, nu, f, bv, v, r)
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

        if (converged .or. i == i_max) exit

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
  !> Diagonal preconditioner based on the mass matrix: z = τ M⁻¹ r

  subroutine Diagonal_Preconditioner(ins_op, tau, r, z, standby)

    class(INS_Operator_3D), intent(in)  :: ins_op       !< INS DG operator
    real(RNP),              intent(in)  :: tau          !< τ
    real(RNP), contiguous,  intent(in)  :: r(:,:,:,:,:) !< residual
    real(RNP), contiguous,  intent(out) :: z(:,:,:,:,:) !< correction
    logical, intent(in) :: standby !< switch for reusing workspace

    ! local variables ..........................................................

    real(RNP), allocatable, save :: mm_inv(:,:,:,:) ! inverse mass matrix

    ! parameters saved for reuse
    integer :: np = -1  ! number of element points per direction
    integer :: na = -1  ! number of active elements
    integer :: ne = -1  ! number of elements

    logical :: reuse
    integer :: e

    ! initialization .........................................................

    !$omp master

    reuse = np == size(z,1)                     .and. &
            na == ins_op % mesh % n_elem_active .and. &
            ne == ins_op % mesh % n_elem

    if (.not. reuse) then
      if (allocated(mm_inv)) deallocate(mm_inv)
      np = size(z,1)
      na = ins_op % mesh % n_elem_active
      ne = ins_op % mesh % n_elem
      allocate(mm_inv(np,np,np,ne))
    end if

    !$omp end master

    if (.not. reuse) then
      call ins_op % sem_u % Get_DG_DiagonalMassMatrix(mm_inv)

      !$omp do
      do e = 1, na
        mm_inv(:,:,:,e) = ONE / mm_inv(:,:,:,e)
      end do

    end if

    ! preconditioning ..........................................................

    !$omp do
    do e = 1, na
      z(:,:,:,e,1) = tau * mm_inv(:,:,:,e) * r(:,:,:,e,1)
      z(:,:,:,e,2) = tau * mm_inv(:,:,:,e) * r(:,:,:,e,2)
      z(:,:,:,e,3) = tau * mm_inv(:,:,:,e) * r(:,:,:,e,3)
    end do

    !$omp do
    do e = na+1, ne
      z(:,:,:,e,1) = ZERO
      z(:,:,:,e,2) = ZERO
      z(:,:,:,e,3) = ZERO
    end do

    ! cleanup ................................................................

    ! keep workspace in case of standby
    if (standby) return

    !$omp master
    if (allocated(mm_inv)) deallocate(mm_inv)
    np = -1
    na = -1
    ne = -1
    !$omp end master

  end subroutine Diagonal_Preconditioner

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

    character, allocatable, save :: bc_schwarz(:)  ! BC adapted for Schwarz
    integer,   allocatable, save :: cfg(:,:,:)     ! subdomain configurations

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_rg
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_zs

    ! parameters saved for reuse
    integer :: np = -1  ! number of element points per direction
    integer :: na = -1  ! number of active elements
    integer :: ne = -1  ! number of elements
    integer :: ng = -1  ! number of ghosts
    integer :: no = -1  ! overlap
    integer :: wp = -1  ! working precision

    integer :: d, nl(3), ns
    logical :: reuse

    associate(mesh => ins_op % mesh, schwarz => ins_op % schwarz_u)

      ! initialization .........................................................

      !$omp master

      reuse = np == size(z,1)            .and. &
              na == mesh % n_elem_active .and. &
              ne == mesh % n_elem        .and. &
              ng == mesh % n_ghost       .and. &
              no == schwarz % no         .and. &
              wp == schwarz % wp

      if (.not. reuse) then

        if (allocated( rg         )) deallocate( rg         )
        if (allocated( nu_dp      )) deallocate( nu_dp      )
        if (allocated( zs_dp      )) deallocate( zs_dp      )
        if (allocated( rs_dp      )) deallocate( rs_dp      )
        if (allocated( nu_sp      )) deallocate( nu_sp      )
        if (allocated( zs_sp      )) deallocate( zs_sp      )
        if (allocated( rs_sp      )) deallocate( rs_sp      )
        if (allocated( buf_rg     )) deallocate( buf_rg     )
        if (allocated( buf_zs     )) deallocate( buf_zs     )
        if (allocated( bc_schwarz )) deallocate( bc_schwarz )
        if (allocated( cfg        )) deallocate( cfg        )

        np = size(z,1)
        na = mesh % n_elem_active
        ne = mesh % n_elem
        ng = mesh % n_ghost
        no = schwarz % no
        wp = schwarz % wp

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

        allocate(bc_schwarz, source = ins_op % problem % bc_v)
        where(bc_schwarz == 'O')
          bc_schwarz = 'N'
        end where

        allocate(cfg(3,na,3))

      end if

      !$omp end master

      if (.not. reuse) then
        ! subdomain configurations -- so far identical for all directions
        call schwarz % GetSubdomainConfigurations(mesh, bc_schwarz, cfg(:,:,1))
        cfg(:,:,2) = cfg(:,:,1)
        cfg(:,:,3) = cfg(:,:,1)
      end if

      ! coefficients
      select case(wp)
      case(RSP)
        lambda_sp = real(1/tau, RSP)
        !$omp workshare
        nu_sp = real(nu_avg, RSP)
        !$omp end workshare nowait
      case default
        lambda_dp = real(1 / tau, RDP)
        !$omp workshare
        nu_dp = real(nu_avg, RDP)
        !$omp end workshare nowait
      end select

      call SetArray(z, ZERO, multi = .true.)

      ! Schwarz sweeps .........................................................

      do d = 1, 3

        call SetArray(rg(:,:,:,:ne), r(:,:,:,:,d))

        select case(wp)

        case(RSP)
          call schwarz % RestrictResidual(mesh, buf_rg, rg, rs_sp)
          call TPO_Schwarz( schwarz % ops_sp % S  &
                          , schwarz % ops_sp % V  &
                          , schwarz % ops_sp % W  &
                          , schwarz % ops_sp % g  &
                          , cfg(:,:,d)            &
                          , lambda_sp             &
                          , nu_sp(:na)            &
                          , rs_sp(:,:,:,:na)      &
                          , zs_sp(:,:,:,:na)      )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_sp, z(:,:,:,:,d))

        case default
          call schwarz % RestrictResidual(mesh, buf_rg, rg, rs_dp)
          call TPO_Schwarz( schwarz % ops_dp % S  &
                          , schwarz % ops_dp % V  &
                          , schwarz % ops_dp % W  &
                          , schwarz % ops_dp % g  &
                          , cfg(:,:,d)            &
                          , lambda_dp             &
                          , nu_dp(:na)            &
                          , rs_dp(:,:,:,:na)      &
                          , zs_dp(:,:,:,:na)      )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_dp, z(:,:,:,:,d))

        end select

      end do

      ! cleanup ................................................................

      ! keep workspace in case of standby
      if (standby) return

      !$omp master

      if (allocated( rg         )) deallocate( rg         )
      if (allocated( nu_dp      )) deallocate( nu_dp      )
      if (allocated( zs_dp      )) deallocate( zs_dp      )
      if (allocated( rs_dp      )) deallocate( rs_dp      )
      if (allocated( nu_sp      )) deallocate( nu_sp      )
      if (allocated( zs_sp      )) deallocate( zs_sp      )
      if (allocated( rs_sp      )) deallocate( rs_sp      )
      if (allocated( buf_rg     )) deallocate( buf_rg     )
      if (allocated( buf_zs     )) deallocate( buf_zs     )
      if (allocated( bc_schwarz )) deallocate( bc_schwarz )
      if (allocated( cfg        )) deallocate( cfg        )

      np = -1
      na = -1
      ne = -1
      ng = -1
      no = -1
      wp = -1

      !$omp end master

    end associate

  end subroutine Schwarz_Preconditioner

  !=============================================================================

end submodule MP_DiffusionSolver
