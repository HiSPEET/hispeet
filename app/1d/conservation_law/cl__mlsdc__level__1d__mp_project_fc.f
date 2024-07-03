!> summary:  Fine-to-coarse space-time solution projection
!> author:   Joerg Stiller, Erik Pfister
!> date:     2023/08/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_Project_FC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse space-time projection of solution-like variables

  module subroutine Project_FC(this, u_f, u_c)
    class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
    real(RNP), intent(in)    :: u_f(0:,:,:,0:,:) !< fine solution variable
    real(RNP), intent(inout) :: u_c(0:,:,:,0:,:) !< coarse solution variable

    real(RNP), allocatable :: u_i(:,:,:,:,:) ! intermediate interpolant
    real(RNP), allocatable :: u_t(:,:,:)     ! variable smoothed in time
    real(RNP), allocatable :: u_x(:,:)       ! variable smoothed in space
    real(RNP), allocatable :: jmp_u_t(:)     ! jump in time
    real(RNP)              :: jmp_u_x        ! jump in space

    integer :: ps_c, ps_f ! polynomial degree in space
    integer :: ns_c, ns_f ! number of elements in space
    integer :: pt_c, pt_f ! polynomial degree in time
    integer :: mt_c, mt_f ! number of subintervals per time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: o_t        ! offset of time collocation points
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: c, e, l, m, n
    integer :: i0_1, i1_1, i0_2, i1_2
    integer :: m0_1, m1_1, m0_2, m1_2
    integer :: o_i1, o_i2

    associate( pop_x => this % pop_fc_x &
             , pop_t => this % pop_fc_t &
             , activity => this % cl_operator % activity )

      ! initialization .........................................................

      ! fine dimensions
      ps_f = ubound(u_f, 1)
      ns_f = ubound(u_f, 2)
      mt_f = ubound(u_f, 4)
      nt_f = ubound(u_f, 5)

      ! coarse dimensions
      ps_c = ubound(u_c, 1)
      ns_c = ubound(u_c, 2)
      mt_c = ubound(u_c, 4)
      nt_c = ubound(u_c, 5)

      ! polynomial degree in time and offset collocation points
      select case(this % cl_sdc % point_set)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        o_t   = 1
      case('E','L')
        ! equidistant or Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        o_t   = 0
      end select
      o_i1 = o_t

      ! number of components
      nc = this % cl_problem % nc

      ! consistency of components
      is_consistent = ubound(u_f, 3) == nc .and. &
                      ubound(u_c, 3) == nc

      ! consistency of spatial dimensions
      is_consistent = is_consistent        .and. &
                      ps_f == pop_x % po_f .and. &
                      ps_c == pop_x % po_c

      ! consistency with spatial projection mode
      select case(pop_x % mode)
      case(1)
        is_consistent = is_consistent .and. ns_f == ns_c
      case(2)
        is_consistent = is_consistent .and. ns_f == ns_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Project_FC'               &
                  , 'failed consistency check' &
                  , 'CL__MLSDC__Level__1D'     )
      end if

      ! workspace
      allocate(u_i(0:ps_f,1:ns_f,nc,0:mt_c,1:nt_c), source = ZERO)

      ! temporal projection ....................................................

      select case(pop_t % mode)

      case(0)
        u_i = u_f

      case(1)
        ! nt_f = nt_c
        do n = 1, nt_c

          ! inject boundary values (m = 0, mt)
          do c = 1, nc
          do e = 1, ns_f
            if (activity(e) > 0) then
              u_i(:,e,c, 0    ,n) = u_f(:,e,c, 0    ,n)
              u_i(:,e,c, mt_c ,n) = u_f(:,e,c, mt_f ,n)
            end if
          end do
          end do

          do m = 1, mt_c-1
            do c = 1, nc
            do e = 1, ns_f
              if (activity(e) > 0) then
                do l = o_t, mt_f
                  u_i(:,e,c,m,n) = u_i(:,e,c,m,n) &
                                 + pop_t % A(m-o_t,l-o_t,1) * u_f(:,e,c,l,n)
                end do
              end if
            end do
            end do
          end do
        end do

      case(2)
        ! nt_f = 2 * nt_c

        ! coarse point ranges
        select case(pop_t % method)
        case('P')
          ! L² projection: all fine points contribute to all coarse points
          m0_1 = 1
          m1_1 = mt_c
          m0_2 = 1
          m1_2 = mt_c
          o_i2 = m0_2 + o_t - 1
        case('I')
          ! interpolation: fine points contribute to left/right half only
          m0_1 = 1                        ! first point interpolated from element 1
          m1_1 = pt_c / 2 + o_t           ! last  ..
          m0_2 = m1_1 + mod(pt_c,2)       ! first point interpolated from element 2
          m1_2 = mt_c                     ! last  ..
          o_i2 = m0_2
        end select

        allocate(u_t(0:ps_f,0:mt_f,2), jmp_u_t(0:ps_f))

        ! projection
        do n = 1, nt_c

          ! left boundary values (m = 0)
          if (n == 1) then
            ! inject left boundary values interpolated from fine mesh
            do c = 1, nc
            do e = 1, ns_f
              if (activity(e) > 0) then
                u_i(:,e,c,0,1) = u_f(:,e,c,0,1)
              end if
            end do
            end do
          else
            ! adopt right boundary value from preceding step
            do c = 1, nc
            do e = 1, ns_f
              if (activity(e) > 0) then
                u_i(:,e,c,0,n) = u_i(:,e,c,mt_c,n-1)
              end if
            end do
            end do
          end if

          do c = 1, nc
          do e = 1, ns_f
            if (activity(e) > 0) then

              ! extract solution from involved fine elements
              do m = o_t, mt_f
                u_t(:,m,1) = u_f(:,e,c,m,2*n - 1)
                u_t(:,m,2) = u_f(:,e,c,m,2*n    )
              end do

              ! optional smoothing, e.g. when using DG in time
              jmp_u_t = u_t(:,mt_f,1) - u_t(:,0,2)
              select case(pop_t % smooth)
              case(1)
                ! linear
                do m = o_t, mt_f
                  u_t(:,m,1) = u_t(:,m,1) + jmp_u_t * pop_t % B(m,1)
                  u_t(:,m,2) = u_t(:,m,2) + jmp_u_t * pop_t % B(m,2)
                end do
              case(2)
                ! averaging interface coefficients
                u_t(:,mt_f,1) = u_t(:,mt_f,1) - HALF * jmp_u_t
                u_t(:,   0,2) = u_t(:,   0,2) + HALF * jmp_u_t
              end select

              ! projection of smoothed variable
              do m = m0_1, m1_1
                do l = o_t, mt_f
                  u_i(:,e,c,m,n) = u_i(:,e,c,m,n) &
                                 + pop_t % A(m-o_i1,l-o_t,1) * u_t(:,l,1)
                end do
              end do
              do m = m0_2, m1_2
                do l = o_t, mt_f
                  u_i(:,e,c,m,n) = u_i(:,e,c,m,n) &
                                 + pop_t % A(m-o_i2,l-o_t,2) * u_t(:,l,2)
                end do
              end do

            end if
          end do
          end do

        end do

      end select

      ! spatial projection .....................................................

      select case(pop_x % mode)

      case(0)
        u_c = u_i

      case(1)
        ! ns_f = ns_c
        do n = 1, nt_c
        do m = 0, mt_c
        do c = 1, nc
        do e = 1, ns_c
          if (activity(e) > 0) then
            u_c(:,e,c,m,n) = matmul(pop_x % A(:,:,1), u_i(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do

      case(2)
        ! ns_f = 2 * ns_c

        ! point ranges
        select case(pop_x % method)
        case('P')
          ! L² projection: all fine points contribute to all coarse points
          i0_1 = 0
          i1_1 = ps_c
          i0_2 = 0
          i1_2 = ps_c
        case('I')
          ! interpolation: fine points contribute to left/right half only
          i0_1 = 0                  ! first point interpolated from element 1
          i1_1 = ps_c / 2           ! last  ..
          i0_2 = i1_1 + mod(ps_c,2) ! first point interpolated from element 2
          i1_2 = ps_c               ! last  ..
        end select

        ! workspace
        allocate(u_x(0:ps_f,2))

        ! projection
        do n = 1, nt_c
        do m = 0, mt_c
        do c = 1, nc
        do e = 1, ns_c
          if (activity(e) > 0) then

            ! extract solution from involved fine elements
            u_x(:,1) = u_i(:,2*e - 1,c,m,n)
            u_x(:,2) = u_i(:,2*e    ,c,m,n)

            ! optional smoothing
            jmp_u_x = u_x(ps_f,1) - u_x(0,2)
            select case(pop_x % smooth)
            case(1)
              ! linear
              u_x(:,1) = u_x(:,1) + jmp_u_x * pop_x % B(:,1)
              u_x(:,2) = u_x(:,2) + jmp_u_x * pop_x % B(:,2)
            case(2)
              ! averaging interface coefficients
              u_x(ps_f,1) = u_x(ps_f,1) - HALF * jmp_u_x
              u_x(   0,2) = u_x(   0,2) + HALF * jmp_u_x
            end select

            ! projection of smoothed variable
            u_c(    :    ,e,c,m,n) = 0
            u_c(i0_1:i1_1,e,c,m,n) = u_c(i0_1:i1_1,e,c,m,n) &
                                   + matmul(pop_x % A(:,:,1), u_x(:,1))
            u_c(i0_2:i1_2,e,c,m,n) = u_c(i0_2:i1_2,e,c,m,n) &
                                   + matmul(pop_x % A(:,:,2), u_x(:,2))

          end if
        end do
        end do
        end do
        end do

      end select

    end associate

  end subroutine Project_FC

  !=============================================================================

end submodule MP_Project_FC
