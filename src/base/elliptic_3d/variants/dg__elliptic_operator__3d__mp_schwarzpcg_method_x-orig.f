!> summary:  3D DG elliptic operator: IPCG with Schwarz preconditioner
!> author:   Joerg Stiller
!> date:     2022/10/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_SchwarzPCG_Method_X
  use Array_Assignments
  use Array_Reductions
  use TPO__Average__3D
  use TPO__Schwarz__3D
  use Element_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Inexact Schwarz-preconditioned CG method with either constant or variable ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent
  !>
  !> This routine implements the IPCG method proposed in: G. Golub & Q. Ye,
  !> SIAM J. Sci. Comput. 21(4):1305-1320, 1999

  module subroutine SchwarzPCG_Method_X( this, bc, lambda, nu_c, nu_v, u, f &
                                       , bv, i_max, r_red, r_max, ni        )

    class(DG_EllipticOperator_3D),        intent(in)    :: this
    character,                            intent(in)    :: bc(:)
    real(RNP),                            intent(in)    :: lambda        !< λ
    real(RNP),                  optional, intent(in)    :: nu_c          !< νᵖ+νˢ
    real(RNP), contiguous,      optional, intent(in)    :: nu_v(:,:,:,:) !< νᵖ
    real(RNP), contiguous,                intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous,                intent(in)    :: f(:,:,:,:)
    class(BoundaryVariable_3D), optional, intent(in)    :: bv(:)
    integer,                              intent(in)    :: i_max
    real(RNP),                  optional, intent(in)    :: r_red
    real(RNP),                  optional, intent(in)    :: r_max
    integer,                    optional, intent(out)   :: ni

    ! internal variables .......................................................

    ! residual with ghost entries
    real(RNP), allocatable, target, save :: rg(:,:,:,:)

    ! residual with no ghost entries
    real(RNP), contiguous, pointer, save :: r(:,:,:,:)

    ! work arrays
    real(RNP), dimension(:,:,:,:), allocatable, save :: p, q, s, z
    real(RNP), dimension(:),       allocatable, save :: nu_avg

    ! control
    real(RNP), save :: rr_term
    logical  , save :: converged

    ! auxiliaries for the DP Schwarz preconditioner
    real(RDP) :: lambda_dp                         ! Helmholtz parameter
    real(RDP), allocatable, save :: nu_dp(:)       ! subdomain diffusivities
    real(RDP), allocatable, save :: fs_dp(:,:,:,:) ! RHS of subsystems
    real(RDP), allocatable, save :: zs_dp(:,:,:,:) ! solution of subsystems

    ! auxiliaries for the SP Schwarz preconditioner
    real(RSP) :: lambda_sp                         ! Helmholtz parameter
    real(RSP), allocatable, save :: nu_sp(:)       ! subdomain diffusivities
    real(RSP), allocatable, save :: fs_sp(:,:,:,:) ! RHS of subsystems
    real(RSP), allocatable, save :: zs_sp(:,:,:,:) ! solution of subsystems

    integer, allocatable, save :: cfg(:,:) ! subdomain configurations

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_rg
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_zs

    logical, save :: singular, singular_loc

    real(RNP) :: alpha, beta, delta, rr
    logical   :: check_convergence
    integer   :: na, ne, ng, nl(3), no, np, ns, wp
    integer   :: i, i_max_

    ! skip empty partition
    if (this % sem % mesh % part < 0) then
      if (present(ni)) ni = -1
      return
    end if

    associate( mesh    => this % sem % mesh &
             , eop     => this % eop        &
             , schwarz => this % schwarz    )

      ! initialization .........................................................

      na = mesh % n_elem_active
      ne = mesh % n_elem
      ng = mesh % n_ghost
      wp = schwarz % wp
      no = schwarz % no
      np = this % eop % po + 1
      nl = no
      ns = np + 2*no

      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      !$omp master

      allocate(p(np, np, np, ne), source = ZERO)
      allocate(q, s, z, source = p)

      allocate(nu_avg(ne), source = ZERO)

      allocate(rg(np, np, np, ne+ng), source = ZERO)
      buf_rg = ElementTransferBuffer_3D(mesh, rg, nl)
      r => rg(:,:,:,1:ne)

      select case(wp)
      case(RSP)
        allocate(nu_sp(ne))
        allocate(fs_sp(ns, ns, ns, ne)   , source = 0E0)
        allocate(zs_sp(ns, ns, ns, ne+ng), source = 0E0)
        buf_zs = ElementTransferBuffer_3D(mesh, zs_sp, nl)
      case default
        allocate(nu_dp(ne))
        allocate(fs_dp(ns, ns, ns, ne)   , source = 0D0)
        allocate(zs_dp(ns, ns, ns, ne+ng), source = 0D0)
        buf_zs = ElementTransferBuffer_3D(mesh, zs_dp, nl)
      end select

      allocate(cfg(3,na))

      singular = abs(lambda) < epsilon(ONE) .and. all(bc /= 'D')
      if (singular .and. .not. mesh%is_root)  then
        singular_loc = mesh%n_elem_frozen == 0
        call XMPI_Allreduce(singular_loc, singular, MPI_LAND, mesh%comm_parts)
      end if

      !$omp end master
      !$omp barrier
      !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      call schwarz % GetSubdomainConfigurations(mesh, bc, cfg)

      check_convergence = log_level_inner_iteration > 0
      if (present(r_red)) check_convergence = r_red > 0 .or. check_convergence
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      ! coefficients ...........................................................

      ! element-averaged diffusivity
      if (present(nu_c)) then
        call SetArray(nu_avg, nu_c)
      else
        call TPO_Average(eop%w, nu_v, nu_avg)
      end if

      select case(wp)
      case(RSP)
        lambda_sp = real(lambda, RSP)
        !$omp workshare
        nu_sp = real(nu_avg, RSP)
        !$omp end workshare nowait
      case default
        lambda_dp = real(lambda, RDP)
        !$omp workshare
        nu_dp = real(nu_avg, RDP)
        !$omp end workshare nowait
      end select
      ! omp barrier not needed because Residual is blocking

      ! initial residual .......................................................

      ! r = f - Au
      if (present(nu_c)) then
        call this % Residual(bc, lambda, nu_c, f, bv, u, r)
      else
        call this % Residual(bc, lambda, nu_v, f, bv, u, r)
      end if
      if (singular) then
        call CalibrateArray(r(:,:,:,:na), mesh%comm_parts)
      end if

      ! termination conditions
      if (check_convergence) then
        rr = ScalarProduct(r(:,:,:,:na), r(:,:,:,:na), mesh%comm_parts)
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
                '#Elliptic:SchwarzPCG','>>>  i  =',0,',  |r| =', sqrt(rr)
        end if
        !$omp end master
      else
        !$omp master
        converged = .false.
        !$omp end master
      end if
      !$omp barrier

      if (converged) then
        i_max_ = 0
        i      = 0
      else
        i_max_ = i_max
      end if

      ! iteration ...............................................................

      do i = 1, i_max_

        ! Schwarz preconditioner, z = (Aˢ)⁻¹ r
        call SetArray(z, ZERO)
        select case(wp)
        case(RSP)
          call schwarz % RestrictResidual(mesh, buf_rg, rg, fs_sp)
          call TPO_Schwarz( schwarz % ops_sp % S  &
                          , schwarz % ops_sp % V  &
                          , schwarz % ops_sp % W  &
                          , schwarz % ops_sp % g  &
                          , cfg                   &
                          , lambda_sp             &
                          , nu_sp(:na)            &
                          , fs_sp(:,:,:,:na)      &
                          , zs_sp(:,:,:,:na)      )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_sp, z)
        case default
          call schwarz % RestrictResidual(mesh, buf_rg, rg, fs_dp)
          call TPO_Schwarz( schwarz % ops_dp % S  &
                          , schwarz % ops_dp % V  &
                          , schwarz % ops_dp % W  &
                          , schwarz % ops_dp % g  &
                          , cfg                   &
                          , lambda_dp             &
                          , nu_dp(:na)            &
                          , fs_dp(:,:,:,:na)      &
                          , zs_dp(:,:,:,:na)      )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_dp, z)
        end select

        ! set/update search vector
        if (i == 1) then
          if (singular) then
            call CalibrateArray(z(:,:,:,:na), mesh%comm_parts)
          end if
          call SetArray(p(:,:,:,:na), z(:,:,:,:na))                ! p = z
        else
          call SetArray(q(:,:,:,:na), r(:,:,:,:na))                ! q = r
          call MergeArrays(ONE, q(:,:,:,:na), -ONE, s(:,:,:,:na))  ! q = r - s
          beta = ScalarProduct( q(:,:,:,:na) &
                              , z(:,:,:,:na) &
                              , mesh%comm_parts ) / delta
          call MergeArrays(beta, p(:,:,:,:na), ONE, z(:,:,:,:na))  ! p = β p + z
        end if

        ! save old residual
        call SetArray(s, r)

        ! operator application with no source and homogeneous BC
        if (present(nu_c)) then
          call this % Apply(bc, lambda, nu_c, u=p, r=q)
        else
          call this % Apply(bc, lambda, nu_v, u=p, r=q)
        end if

        ! correction
        delta = ScalarProduct(r, z, mesh%comm_parts)
        alpha = delta / ScalarProduct(p, q, mesh%comm_parts)
        call MergeArrays(ONE, u(:,:,:,:na), alpha, p(:,:,:,:na))

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          if (present(nu_c)) then
            call this % Residual(bc, lambda, nu_c, f, bv, u, r)
          else
            call this % Residual(bc, lambda, nu_v, f, bv, u, r)
          end if
          if (singular) then
            call CalibrateArray(r(:,:,:,:na), mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r(:,:,:,:na), -alpha, q(:,:,:,:na))
        end if

        if (check_convergence) then
          rr = ScalarProduct(r(:,:,:,:na), r(:,:,:,:na), mesh%comm_parts)
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
                '#Elliptic:SchwarzPCG','>>>  i  =',i,',  |r| =', sqrt(rr)
        end if
        !$omp end master

      end do

      !$omp master
      if (log_level_inner_iteration > 0 .and. mesh%part == 0) then
        print '(A,T25,A,I5,A,ES12.5)', &
              '#Elliptic:SchwarzPCG','>>>  ni =',i,',  |r| =', sqrt(rr)
      end if
      !$omp end master

      ! finalization ...........................................................

      if (present(ni)) ni = i

      !$omp master
      deallocate(p, q, s, z)
      deallocate(nu_avg, rg)
      deallocate(buf_rg, buf_zs)
      select case(wp)
      case(RSP)
        deallocate(nu_sp, fs_sp, zs_sp)
      case default
        deallocate(nu_dp, fs_dp, zs_dp)
      end select
      deallocate(cfg)
      !$omp end master

    end associate

  end subroutine SchwarzPCG_Method_X

  !=============================================================================

end submodule MP_SchwarzPCG_Method_X
