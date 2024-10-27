!> summary:  3D DG elliptic operator: Overlapping additive Schwarz method
!> author:   Joerg Stiller
!> date:     2022/10/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(DG__Elliptic_Operator__3D) MP_Schwarz_Method_X
  use Array_Assignments
  use Array_Reductions
  use TPO__Average__3D
  use TPO__Schwarz__3D
  use Element_Transfer_Buffer__3D
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Element-centered overlapping Schwarz method with constant or variable ν
  !>
  !> Homogeneous conditions are used if boundary values `bv` are absent
  !>
  !> @remark
  !> One and only one of the parameters `nu_c` and `nu_v` is to be passed

  module subroutine Schwarz_Method_X( this, lambda, nu_c, nu_v, u, f, bv &
                                    , i_max, r_red, r_max, ni            )

    class(DG_EllipticOperator_3D),        intent(in)    :: this
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

    ! local variables ..........................................................

    real(RNP), allocatable, save :: r(:,:,:,:)     ! residual, including ghosts
    real(RNP), allocatable, save :: nu_avg(:)      ! element averaged ν

    real(RDP) :: lambda_dp                         ! Helmholtz parameter
    real(RDP), allocatable, save :: nu_dp(:)       ! subdomain diffusivities
    real(RDP), allocatable, save :: fs_dp(:,:,:,:) ! RHS of subsystems
    real(RDP), allocatable, save :: us_dp(:,:,:,:) ! solution of subsystems

    real(RSP) :: lambda_sp                         ! Helmholtz parameter
    real(RSP), allocatable, save :: nu_sp(:)       ! subdomain diffusivities
    real(RSP), allocatable, save :: fs_sp(:,:,:,:) ! RHS of subsystems
    real(RSP), allocatable, save :: us_sp(:,:,:,:) ! solution of subsystems

    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_r
    type(ElementTransferBuffer_3D), asynchronous, allocatable, save :: buf_us

    real(RNP), save :: rr_term, rr_red
    logical  , save :: converged

    real(RNP) :: rr
    integer   :: i, ne, ng, nl(3), no, np, ns, wp
    logical   :: check_convergence

    ! skip empty partition
    if (this % sem % mesh % part < 0) then
      if (present(ni)) ni = -1
      return
    end if

    associate( mesh    => this % sem % mesh &
             , eop     => this % eop        &
             , schwarz => this % schwarz    )

      ! initialization .........................................................

      check_convergence = .false.
      if (present(r_red)) check_convergence = r_red > 0
      if (present(r_max)) check_convergence = r_max > 0 .or. check_convergence

      ne = mesh % n_elem
      ng = mesh % n_ghost
      wp = schwarz % wp
      no = schwarz % no
      np = this % eop % po + 1
      nl = no
      ns = np + 2*no

      !$omp master

      allocate(r(np, np, np, ne+ng), source = ZERO)
      allocate(nu_avg(ne))

      select case(wp)
      case(RSP)
        allocate(fs_sp(ns, ns, ns, ne)   , source = 0E0)
        allocate(us_sp(ns, ns, ns, ne+ng), source = 0E0)
        buf_us = ElementTransferBuffer_3D(mesh, us_sp, nl)
      case default
        allocate(fs_dp(ns, ns, ns, ne)   , source = 0D0)
        allocate(us_dp(ns, ns, ns, ne+ng), source = 0D0)
        buf_us = ElementTransferBuffer_3D(mesh, us_dp, nl)
      end select
      buf_r = ElementTransferBuffer_3D(mesh, r, nl)

      ! termination condition
      if (check_convergence) then
        if (present(r_max)) then
          rr_term = r_max ** 2
        else
          rr_term = ZERO
        end if
        if (present(r_red)) then
          rr_red = r_red ** 2
        else
          rr_red = ZERO
        end if
      else
        converged = .false.
      end if

      !$omp end master
      !$omp barrier

      ! coefficients ...........................................................

      ! element-averaged diffusivity
      if (present(nu_c)) then
        call SetArray(nu_avg, nu_c)
      else
        call TPO_Average(eop%w, nu_v, nu_avg)
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

      ! Schwarz iterations .....................................................

      do i = 1, i_max

        ! residual
        if (present(nu_c)) then
          call this % Residual(lambda, nu_c, f, bv, u, r)
        else
          call this % Residual(lambda, nu_v, f, bv, u, r)
        end if

        ! termination check
        if (check_convergence) then
          rr = ScalarProduct(r(:,:,:,:ne), r(:,:,:,:ne), mesh%comm_parts)
          !$omp master
          if (mesh%part == 0) then
            if (i == 1) then
              rr_term  = max(rr_term, rr * rr_red)
            end if
            converged = rr <= rr_term
          end if
          call XMPI_Bcast(converged, root=0, comm=mesh%comm_parts)
          !$omp end master
          !$omp barrier
        end if
        if (converged) exit

        ! Schwarz sweep
        select case(wp)
        case(RSP)
          call schwarz % RestrictResidual(mesh, buf_r, r, fs_sp)
          call TPO_Schwarz( schwarz % ops_sp % S      &
                          , schwarz % ops_sp % V      &
                          , schwarz % ops_sp % W      &
                          , schwarz % ops_sp % g      &
                          , schwarz % cfg             &
                          , lambda_sp                 &
                          , nu_sp                     &
                          , fs_sp                     &
                          , us_sp                     )
          call schwarz % MergeCorrections(mesh, buf_us, us_sp, u)
        case default
          call schwarz % RestrictResidual(mesh, buf_r, r, fs_dp)
          call TPO_Schwarz( schwarz % ops_dp % S      &
                          , schwarz % ops_dp % V      &
                          , schwarz % ops_dp % W      &
                          , schwarz % ops_dp % g      &
                          , schwarz % cfg             &
                          , lambda_dp                 &
                          , nu_dp                     &
                          , fs_dp                     &
                          , us_dp                     )
          call schwarz % MergeCorrections(mesh, buf_us, us_dp, u)
        end select

      end do

      if (present(ni)) ni = min(i, i_max)

      ! cleanup ................................................................

      !$omp master
      deallocate(nu_avg, buf_r, buf_us, r)
      select case(wp)
      case(RSP)
        deallocate(nu_sp, fs_sp, us_sp)
      case default
        deallocate(nu_dp, fs_dp, us_dp)
      end select
      !$omp end master

    end associate

  end subroutine Schwarz_Method_X

  !=============================================================================

end submodule MP_Schwarz_Method_X
