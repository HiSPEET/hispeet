!> summary:  Fine-to-coarse restriction of residual-like variables
!> author:   Joerg Stiller, Erik Pfister
!> date:     2023/09/07
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_Restrict_FC
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Fine-to-coarse restriction of residual-like variables

  module subroutine Restrict_FC(this, r_f, r_c)
    class(CL_MLSDC_Level_1D), intent(in) :: this !< coarse level
    real(RNP), intent(in)    :: r_f(0:,:,:,0:,:) !< fine residual variable
    real(RNP), intent(inout) :: r_c(0:,:,:,0:,:) !< coarse residual variable

    real(RNP), allocatable :: r_i(:,:,:,:,:) ! intermediate result

    integer :: ps_c, ps_f ! polynomial degree in space
    integer :: ns_c, ns_f ! number of elements in space
    integer :: pt_c, pt_f ! polynomial degree in time
    integer :: mt_c, mt_f ! number of subintervals per time step
    integer :: nt_c, nt_f ! number of time steps
    integer :: o_t        ! offset of time collocation points
    integer :: nc         ! number of components, must not change

    logical :: is_consistent
    integer :: c, e, l, m, n, n1, n2

    associate( iop_x => this % iop_cf_x &
             , iop_t => this % iop_cf_t &
             , refinement => this % cl_operator % refinement )

      ! initialization .........................................................

      ! coarse dimensions
      ps_c = ubound(r_c, 1)
      ns_c = ubound(r_c, 2)
      mt_c = ubound(r_c, 4)
      nt_c = ubound(r_c, 5)

      ! fine dimensions
      ps_f = ubound(r_f, 1)
      ns_f = ubound(r_f, 2)
      mt_f = ubound(r_f, 4)
      nt_f = ubound(r_f, 5)

      ! polynomial degree in time and offset collocation points
      select case(this % cl_sdc % nodes)
      case('RR')
        ! Radau-right
        pt_c = mt_c - 1
        pt_f = mt_f - 1
        o_t  = 1
      case('E','L')
        ! equidistant of Lobatto points
        pt_c = mt_c
        pt_f = mt_f
        o_t  = 0
      end select

      ! number of components
      nc = this % cl_problem % nc

      ! consistency of components
      is_consistent = ubound(r_c, 3) == nc .and. &
                      ubound(r_f, 3) == nc

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
        call Error( 'Restrict_FC'              &
                  , 'failed consistency check' &
                  , 'CL__MLSDC__Level__1D'     )
      end if

      ! workspace
      allocate(r_i(0:ps_c,1:ns_c,nc,0:mt_f,1:nt_f))

      ! spatial interpolation ..................................................

      select case(iop_x % mode)

      case(0)
        r_i = r_f

      case(1)
        ! ns_f = ns_c
        do n = 1, nt_f
        do m = 0, mt_f
        do c = 1, nc
        do e = 1, ns_c
          if (refinement(e) >= 0) then
            r_i(:,e,c,m,n) = matmul(r_f(:,e,c,m,n), iop_x % A(:,:,1))
          end if
        end do
        end do
        end do
        end do

      case(2)
        ! ns_f = 2 * ns_c
        do n = 1, nt_f
        do m = 0, mt_f
        do c = 1, nc
        do e = 1, ns_c
          if (refinement(e) >= 0) then
            r_i(:,e,c,m,n) = matmul(r_f(:,2*e-1,c,m,n), iop_x % A(:,:,1)) &
                           + matmul(r_f(:,2*e  ,c,m,n), iop_x % A(:,:,2))
          end if
        end do
        end do
        end do
        end do

      end select

      ! temporal interpolation .................................................

      select case(iop_t % mode)

      case(0)
        r_c = r_i

      case(1)
        ! nt_f = nt_c
        if (o_t == 1) then
          r_c(:,:,:,0,:) = 0
        end if
        do n = 1  , nt_c
        do m = o_t, mt_c
        do c = 1  , nc
        do e = 1  , ns_c
          if (refinement(e) >= 0) then
            r_c(:,e,c,m,n) = 0
            do l = o_t, mt_f
              r_c(:,e,c,m,n) = r_c(:,e,c,m,n) &
                             + iop_t % A(l-o_t,m-o_t,1) * r_i(:,e,c,l,n)
            end do
          end if
        end do
        end do
        end do
        end do

      case(2)
        ! nt_f = 2 * nt_c
        if (o_t == 1) then
          r_c(:,:,:,0,:) = 0
        end if
        do n = 1  , nt_c
        do m = o_t, mt_c
        do c = 1  , nc
        do e = 1  , ns_c
            if (refinement(e) >= 0) then
            n1 = 2 * n - 1
            n2 = 2 * n
            r_c(:,e,c,m,n) = 0
            do l = o_t, mt_f
              r_c(:,e,c,m,n) = r_c(:,e,c,m,n)                            &
                             + iop_t % A(l-o_t,m-o_t,1) * r_i(:,e,c,l,n1)  &
                             + iop_t % A(l-o_t,m-o_t,2) * r_i(:,e,c,l,n2)
            end do
          end if
        end do
        end do
        end do
        end do

      end select

    end associate

  end subroutine Restrict_FC

  !=============================================================================

end submodule MP_Restrict_FC
