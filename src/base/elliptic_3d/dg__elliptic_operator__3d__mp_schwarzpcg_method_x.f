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
  !> This routine implements the IPCG method proposed in: G. Golub & Q. Ye,
  !> SIAM J. Sci. Comput. 21(4):1305-1320, 1999

  module subroutine SchwarzPCG_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                       , i_max, r_red, r_max, ni            )

    class(DG_EllipticOperator_3D),             intent(in)    :: this
    real(RNP),                                 intent(in)    :: lambda
    real(RNP),                       optional, intent(in)    :: nu_c
    real(RNP), contiguous,           optional, intent(in)    :: nu_v(:,:,:,:)
    real(RNP), contiguous,                     intent(inout) :: u(:,:,:,:)
    real(RNP), contiguous,                     intent(in)    :: f(:,:,:,:)
    class(SpectralElementBoundaryVariable_3D), intent(in)    :: bv(:)
    integer,                                   intent(in)    :: i_max
    real(RNP),                       optional, intent(in)    :: r_red
    real(RNP),                       optional, intent(in)    :: r_max
    integer,                         optional, intent(out)   :: ni

    ! internal variables .......................................................

    ! residual with ghost entries
    real(RNP), allocatable, target, save :: rg(:,:,:,:)

    ! residual with no ghost entries
    real(RNP), pointer, save :: r(:,:,:,:)

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

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_rg
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_zs

    real(RNP), parameter :: eps = epsilon(ONE) * 1e-3
    real(RNP) :: alpha, beta, delta, rr
    logical   :: check_convergence, singular
    integer   :: ne, ng, nl(3), no, np, ns, wp
    integer   :: i, i_max_

    ! skip empty partition
    if (this % sem % mesh % part < 0) return


    associate( mesh    => this % sem % mesh &
             , eop     => this % eop        &
             , schwarz => this % schwarz    )

      ! initialization .........................................................

      ne = mesh % n_elem
      ng = mesh % n_ghost
      wp = schwarz % wp
      no = schwarz % no
      np = this % eop % po + 1
      nl = no
      ns = np + 2*no

      !$omp master

      allocate(p, mold = u)
      allocate(q, mold = u)
      allocate(s, mold = u)
      allocate(z, mold = u)

      allocate(nu_avg(ne))

      allocate(rg(np, np, np, ne+ng), source = ZERO)
      buf_rg = ElementTransferBuffer_3D(mesh, r, nl)
      r => rg(:,:,:,1:ne)

      select case(wp)
      case(RSP)
        allocate(fs_sp(ns, ns, ns, ne)   , source = 0E0)
        allocate(zs_sp(ns, ns, ns, ne+ng), source = 0E0)
        buf_zs = ElementTransferBuffer_3D(mesh, zs_sp, nl)
      case default
        allocate(fs_dp(ns, ns, ns, ne)   , source = 0D0)
        allocate(zs_dp(ns, ns, ns, ne+ng), source = 0D0)
        buf_zs = ElementTransferBuffer_3D(mesh, zs_dp, nl)
      end select

      !$omp end master
      !$omp barrier

      check_convergence = present(r_red) .or. present(r_max)
      singular = abs(lambda) < epsilon(ONE) .and. all(this%bc /= 'D')

      ! coefficients ...........................................................

      ! element-averaged diffusivity
      if (present(nu_c)) then
        call SetArray(nu_avg, this % TotalDiffusivity(nu_c))
      else
        call TPO_Average(ONE/8, eop%w, nu_v, nu_avg)
      end if

      !$omp master
      select case(wp)
      case(RSP)
        lambda_sp = real(lambda, RSP)
        nu_sp     = real(nu_avg, RSP)
      case default
        lambda_dp = real(lambda, RDP)
        nu_dp     = real(nu_avg, RDP)
      end select
      !$omp end master
      ! omp barrier not needed because Apply is blocking

      ! initial residual .......................................................

      ! r = f - Au
      if (present(nu_c)) then
        call this % Residual(lambda, nu_c, f, bv, u, r)
      else
        call this % Residual(lambda, nu_v, f, bv, u, r)
      end if
      if (singular) then
        call CalibrateArray(r, mesh%comm_parts)
      end if

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

      ! iteration ...............................................................

      do i = 1, i_max_

        ! Schwarz preconditioner, z = (Aˢ)⁻¹ r
        call SetArray(z, ZERO)
        select case(wp)
        case(RSP)
          call schwarz % RestrictResidual(mesh, buf_rg, rg, fs_sp)
          call TPO_Schwarz( schwarz % ops_sp % S      &
                          , schwarz % ops_sp % V      &
                          , schwarz % ops_sp % W      &
                          , schwarz % ops_sp % g      &
                          , schwarz % cfg             &
                          , lambda_sp                 &
                          , nu_sp                     &
                          , fs_sp                     &
                          , zs_sp                     )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_sp, z)
        case default
          call schwarz % RestrictResidual(mesh, buf_rg, rg, fs_dp)
          call TPO_Schwarz( schwarz % ops_dp % S      &
                          , schwarz % ops_dp % V      &
                          , schwarz % ops_dp % W      &
                          , schwarz % ops_dp % g      &
                          , schwarz % cfg             &
                          , lambda_dp                 &
                          , nu_dp                     &
                          , fs_dp                     &
                          , zs_dp                     )
          call schwarz % MergeCorrections(mesh, buf_zs, zs_dp, z)
        end select

        ! set/update search vector
        if (i == 1) then
          if (singular) then
            call CalibrateArray(z, mesh%comm_parts)
          end if
          call SetArray(p, z)                                 ! p = z
        else
          call SetArray(q, r)                                 ! q = r
          call MergeArrays(ONE, q, -ONE, s)                   ! q = r - s
          beta = ScalarProduct(q, z, mesh%comm_parts) / delta
          call MergeArrays(beta, p, ONE, z)                   ! p = beta p + z
        end if

        ! save old residual
        call SetArray(s, r)

        ! operator application with no source and homogeneous BC
        if (present(nu_c)) then
          call this % Apply(lambda, nu_c, p, q)
        else
          call this % Apply(lambda, nu_v, p, q)
        end if

        ! correction
        delta = ScalarProduct(r, z, mesh%comm_parts)
        alpha = delta / ScalarProduct(p, q, mesh%comm_parts)
        call MergeArrays(ONE, u, alpha, p)

        if (mod(i,50) == 0) then
          ! compute true residual to get rid of round-off errors
          if (present(nu_c)) then
            call this % Residual(lambda, nu_c, f, bv, u, r)
          else
            call this % Residual(lambda, nu_v, f, bv, u, r)
          end if
          if (singular) then
            call CalibrateArray(r, mesh%comm_parts)
          end if
        else
          call MergeArrays(ONE, r, -alpha, q)
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

      end do

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
      !$omp end master

    end associate

  end subroutine SchwarzPCG_Method_X

  !=============================================================================

end submodule MP_SchwarzPCG_Method_X
