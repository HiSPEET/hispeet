!> summary:  Coarse-to-fine space-time solution interpolation
!> author:   Joerg Stiller
!> date:     2023/08/25
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_Interpolate_CF
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Coarse-to-fine space-time interpolation of solution-like variables

  module subroutine Interpolate_CF(this, u_c, u_f, complete)
    class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
    real(RNP), intent(in)    :: u_c(0:,:,:,0:,:) !< coarse solution variable
    real(RNP), intent(inout) :: u_f(0:,:,:,0:,:) !< fine solution variable
    logical,   intent(in)    :: complete         !< T/F for all/refined elements

    real(RNP), allocatable, save :: u_i(:,:,:,:,:) ! intermediate interpolant
    logical,   allocatable, save :: eval_elem_c(:) ! switch for coarse elements
    logical,   allocatable, save :: eval_elem_f(:) ! switch for fine elements

    integer :: ps_c, ps_f ! polynomial degree in space
    integer :: ns_c, ns_f ! number of elements in space
    integer :: pt_c, pt_f ! polynomial degree in time
    integer :: mt_c, mt_f ! number of subintervals per time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: ot         ! offset of time collocation points
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: c, e, l, m, n, n1, n2

    associate( iop_x => this % iop_cf_x &
             , iop_t => this % iop_cf_t &
             , refinement => this % cl_operator % refinement )

      ! initialization .........................................................

      ! coarse dimensions
      ps_c = ubound(u_c, 1)
      ns_c = ubound(u_c, 2)
      mt_c = ubound(u_c, 4)
      nt_c = ubound(u_c, 5)

      ! fine dimensions
      ps_f = ubound(u_f, 1)
      ns_f = ubound(u_f, 2)
      mt_f = ubound(u_f, 4)
      nt_f = ubound(u_f, 5)

      ! polynomial degree in time and offset collocation points
      select case(this % cl_sdc % point_set)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        ot   = 1
      case('E','L')
        ! equidistant of Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        ot   = 0
      end select

      ! number of components
      nc = this % cl_problem % nc

      ! consistency of components
      is_consistent = ubound(u_c, 3) == nc .and. &
                      ubound(u_f, 3) == nc

      ! consistency of spatial dimensions
      is_consistent = is_consistent        .and. &
                      ps_c == iop_x % po_c .and. &
                      ps_f == iop_x % po_f

      ! consistency with spatial interpolation mode
      select case(iop_x % mode)
      case(1)
        is_consistent = is_consistent .and. ns_f == ns_c
      case(2)
        is_consistent = is_consistent .and. ns_f == ns_c * 2
      end select

      ! consistency of temporal dimensions
      is_consistent = is_consistent        .and. &
                      pt_c == iop_t % po_c .and. &
                      pt_f == iop_t % po_f

      ! consistency with spatial interpolation mode
      select case(iop_t % mode)
      case(1)
        is_consistent = is_consistent .and. nt_f == nt_c
      case(2)
        is_consistent = is_consistent .and. nt_f == nt_c * 2
      end select

      ! evaluate result of checks
      if (.not. is_consistent) then
        call Error( 'Interpolate_CF'           &
                  , 'failed consistency check' &
                  , 'CL__MLSDC__Level__1D'     )
      end if

      ! workspace
      allocate(u_i(0:ps_f,1:ns_f,nc,0:mt_c,1:nt_c))
      allocate(eval_elem_c(ns_c), source = complete)
      allocate(eval_elem_f(ns_f), source = complete)

      ! identify elements to be evaluated
      if (.not. complete) then
        do e = 1, ns_c
          if (refinement(e) < 0) cycle
          eval_elem_c(e      ) = .true.
          eval_elem_f(e*2 - 1) = .true.
          eval_elem_f(e*2    ) = .true.
        end do
      end if

      ! spatial interpolation ..................................................

      select case(iop_x % mode)

      case(0)

        u_i = u_c

      case(1)

        ! ns_f = ns_c
        do n = 1, nt_c
        do m = 0, mt_c
        do c = 1, nc
        do e = 1, ns_c
          if (eval_elem_c(e)) then
            u_i(:,e,c,m,n) = matmul(iop_x % A(:,:,1), u_c(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do

      case(2)

        ! ns_f = 2 * ns_c
        do n = 1, nt_c
        do m = 0, mt_c
        do c = 1, nc
        do e = 1, ns_c
          if (eval_elem_c(e)) then
            u_i(:,2*e-1,c,m,n) = matmul(iop_x % A(:,:,1), u_c(:,e,c,m,n))
            u_i(:,2*e  ,c,m,n) = matmul(iop_x % A(:,:,2), u_c(:,e,c,m,n))
          end if
        end do
        end do
        end do
        end do

      end select

      ! temporal interpolation .................................................

      select case(iop_t % mode)

      case(0)

        u_f = u_i

      case(1)

        ! nt_f = nt_c
        do n = 1, nt_c

          ! inject boundary values (m = 0, mt)
          do c = 1, nc
          do e = 1, ns_f
            if (eval_elem_f(e)) then
              u_f(:,e,c, 0    ,n) = u_i(:,e,c, 0    ,n)
              u_f(:,e,c, mt_f ,n) = u_i(:,e,c, mt_c ,n)
            end if
          end do
          end do

          ! interpolate interior values
          do m = 1, mt_f-1
            do c = 1, nc
            do e = 1, ns_f
              if (eval_elem_f(e)) then
                u_f(:,e,c,m,n) = 0
                do l = 0, pt_c
                  !! offset ot = 1 required with right-sided points !!
                  u_f(:,e,c,m,n) = u_f(:,e,c,m,n) &
                                 + iop_t % A(m-ot,l,1) * u_i(:,e,c,l+ot,n)
                end do
              end if
            end do
            end do
          end do
        end do

      case(2)

        ! nt_f = 2 * nt_c
        do n = 1, nt_c

          ! left boundary values (m = 0)
          if (n == 1) then
            ! inject left boundary values interpolated from coarse mesh
            do c = 1, nc
            do e = 1, ns_f
              if (eval_elem_f(e)) then
                u_f(:,e,c,0,1) = u_i(:,e,c,0,1)
              end if
            end do
            end do
          else
            ! adopt right boundary value from preceding step
            do c = 1, nc
            do e = 1, ns_f
              if (eval_elem_f(e)) then
                u_f(:,e,c,0,n) = u_f(:,e,c,mt_f,n-1)
              end if
            end do
            end do
          end if

          ! interior values
          do m = 1, mt_f
            do c = 1, nc
            do e = 1, ns_f
              if (eval_elem_f(e)) then
                n1 = 2 * n - 1
                n2 = 2 * n
                u_f(:,e,c,m,n1) = 0
                u_f(:,e,c,m,n2) = 0
                do l = 0, pt_c
                  !! offset ot = 1 required with right-sided points !!
                  u_f(:,e,c,m,n1) = u_f(:,e,c,m,n1) &
                                  + iop_t % A(m-ot,l,1) * u_i(:,e,c,l+ot,n)
                  u_f(:,e,c,m,n2) = u_f(:,e,c,m,n2) &
                                  + iop_t % A(m-ot,l,2) * u_i(:,e,c,l+ot,n)
                end do
              end if
            end do
            end do
          end do

        end do

      end select

      ! release workspace
      deallocate(u_i, eval_elem_c, eval_elem_f)

    end associate

  end subroutine Interpolate_CF

  !=============================================================================

end submodule MP_Interpolate_CF
