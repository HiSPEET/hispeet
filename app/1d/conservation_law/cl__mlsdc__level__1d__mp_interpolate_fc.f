!> summary:  Fine-to-coarse space-time solution interpolation
!> author:   Joerg Stiller
!> date:     2023/08/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_Interpolate_FC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse space-time interpolation of solution-like variables

  module subroutine Interpolate_FC(this, u_f, u_c)
    class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
    real(RNP), intent(in)    :: u_f(0:,:,:,0:,:) !< coarse solution variable
    real(RNP), intent(inout) :: u_c(0:,:,:,0:,:) !< fine solution variable

    real(RNP), allocatable :: u_i(:,:,:,:,:) ! intermediate interpolant
    real(RNP), allocatable :: u_t(:,:,:)     ! variable smoothed in time
    real(RNP), allocatable :: u_x(:,:)       ! variable smoothed in space
    real(RNP), allocatable :: jmp_u_t(:)     ! jump in time
    real(RNP)              :: jmp_u_x        ! jump in space

    integer :: po_c, po_f ! polynomial order in space
    integer :: ne_c, ne_f ! number of elements in space
    integer :: ns_c, ns_f ! number of subintervals in each time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: c, e, l, m, n
    integer :: i0, i1, i2, i3
    integer :: m0, m1, m2, m3

    associate( iop_x => this % iop_fc_x &
             , iop_t => this % iop_fc_t &
             , mask  => this % cl_operator % mask )

      ! initialization .........................................................

      ! fine dimensions
      po_f = ubound(u_f, 1)
      ne_f = ubound(u_f, 2)
      ns_f = ubound(u_f, 4)
      nt_f = ubound(u_f, 5)

      ! coarse dimensions
      po_c = ubound(u_c, 1)
      ne_c = ubound(u_c, 2)
      ns_c = ubound(u_c, 4)
      nt_c = ubound(u_c, 5)

      ! number of components
      nc = this % cl_problem % nc

      ! consistency of components
      is_consistent = ubound(u_f, 3) == nc .and. &
                      ubound(u_c, 3) == nc

      ! consistency of spatial dimensions
      is_consistent = is_consistent        .and. &
                      po_f == iop_x % po_f .and. &
                      po_c == iop_x % po_c

      ! consistency with spatial interpolation mode
      select case(iop_x % mode)
      case(1)
        is_consistent = is_consistent .and. ne_f == ne_c
      case(2)
        is_consistent = is_consistent .and. ne_f == ne_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Interpolate_FC'           &
                  , 'failed consistency check' &
                  , 'CL__MLSDC__Level__1D'     )
      end if

      ! workspace
      allocate(u_i(0:po_f,1:ne_f,nc,0:ns_c,1:nt_c), source = ZERO)

      ! temporal interpolation .................................................

      select case(iop_t % mode)

      case(0)
        u_i = u_f

      case(1)
        ! nt_f = nt_c
        do n = 1, nt_c
        do m = 0, ns_c
        do c = 1, nc
        do e = 1, ne_f
          if (mask(e)) then
            do l = 0, ns_f
              u_i(:,e,c,m,n) = u_i(:,e,c,m,n) + iop_t % A(m,l,1) * u_f(:,e,c,l,n)
            end do
          end if
        end do
        end do
        end do
        end do

      case(2)
        ! nt_f = 2 * nt_c

        ! index ranges
        m0 = 0                ! first coarse point interpolated from fine element 1
        m1 = ns_c / 2         ! last  coarse point interpolated from fine element 1
        m2 = m1 + mod(ns_c,2) ! first coarse point interpolated from fine element 2
        m3 = ns_c             ! last  coarse point interpolated from fine element 2

        ! workspace
        allocate(u_t(0:po_f,0:ns_f,2), jmp_u_t(0:po_f))

        ! interpolation
        do n = 1, nt_c
        do c = 1, nc
        do e = 1, ne_f
          if (mask(e)) then

            ! extract solution from involved fine elements
            do m = 0, ns_f
              u_t(:,m,1) = u_f(:,e,c,m,2*n - 1)
              u_t(:,m,2) = u_f(:,e,c,m,2*n    )
            end do

            ! optional smoothing, e.g. when using DG in time
            jmp_u_t = u_t(:,ns_f,1) - u_t(:,0,2)
            select case(iop_t % smoothing)
            case(1)
              ! linear
              do m = 0, ns_f
                u_t(:,m,1) = u_t(:,m,1) + jmp_u_t * iop_t % B(m,1)
                u_t(:,m,2) = u_t(:,m,2) + jmp_u_t * iop_t % B(m,2)
              end do
            case(2)
              ! averaging interface coefficients
              u_t(:,ns_f,1) = u_t(:,ns_f,1) - HALF * jmp_u_t
              u_t(:,   0,2) = u_t(:,   0,2) + HALF * jmp_u_t
            end select

            ! interpolation of smoothed variable
            do m = m0, m1
              do l = 0, ns_f
                u_i(:,e,c,m,n) = u_i(:,e,c,m,n) + iop_t % A(m,l,1) * u_t(:,l,1)
              end do
            end do
            do m = m2, m3
              do l = 0, ns_f
                u_i(:,e,c,m,n) = u_i(:,e,c,m,n) + iop_t % A(m-m2,l,2) * u_t(:,l,2)
              end do
            end do

          end if
        end do
        end do
        end do

      end select

      ! spatial interpolation ..................................................

      select case(iop_x % mode)

      case(0)
        u_c = u_i

      case(1)
        ! ne_f = ne_c
        do n = 1, nt_c
        do m = 0, ns_c
        do c = 1, nc
        do e = 1, ne_c
          if (mask(e)) then
            u_c(:,e,c,m,n) = matmul(iop_x % A(:,:,1), u_i(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do

      case(2)
        ! ne_f = 2 * ne_c

        ! index ranges
        i0 = 0                ! first coarse point interpolated from fine element 1
        i1 = po_c / 2         ! last  coarse point interpolated from fine element 1
        i2 = i1 + mod(po_c,2) ! first coarse point interpolated from fine element 2
        i3 = po_c             ! last  coarse point interpolated from fine element 2

        ! workspace
        allocate(u_x(0:po_f,2))

        ! interpolation
        do n = 1, nt_c
        do m = 0, ns_c
        do c = 1, nc
        do e = 1, ne_c
          if (mask(e)) then

            ! extract solution from involved fine elements
            u_x(:,1) = u_f(:,2*e - 1,c,m,n)
            u_x(:,2) = u_f(:,2*e    ,c,m,n)

            ! optional smoothing
            jmp_u_x = u_x(po_f,1) - u_x(0,2)
            select case(iop_x % smoothing)
            case(1)
              ! linear
              u_x(:,1) = u_x(:,1) + jmp_u_x * iop_x % B(:,1)
              u_x(:,2) = u_x(:,2) + jmp_u_x * iop_x % B(:,2)
            case(2)
              ! averaging interface coefficients
              u_x(po_f,1) = u_x(po_f,1) - HALF * jmp_u_x
              u_x(   0,2) = u_x(   0,2) + HALF * jmp_u_x
            end select

            ! interpolation of smoothed variable
            u_c(i0:i1,e,c,m,n) = u_c(i0:i1,e,c,m,n) &
                               + matmul(iop_x % A(:,:,1), u_x(:,1))
            u_c(i2:i3,e,c,m,n) = u_c(i2:i3,e,c,m,n) &
                               + matmul(iop_x % A(:,:,2), u_x(:,2))

          end if
        end do
        end do
        end do
        end do

      end select

    end associate

  end subroutine Interpolate_FC

  !=============================================================================

end submodule MP_Interpolate_FC
