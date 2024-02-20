!> summary:  Computation of the collocation residual for the time slice
!> author:   Joerg Stiller, Erik Pfister
!> date:     2023/09/04
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_GetResidual
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Multi-step collocation residual

  module subroutine GetResidual(this, dt, t_0, G, u, r)
    class(CL_MLSDC_Level_1D), intent(in) :: this
    real(RNP), intent(in)            :: dt             !< size of the time slice
    real(RNP), intent(in)            :: t_0            !< start time of the slice
    real(RNP), optional, intent(in)  :: G(0:,:,:,0:,:) !< FAS correction
    real(RNP), intent(in)            :: u(0:,:,:,0:,:) !< approximate solution
    real(RNP), intent(out)           :: r(0:,:,:,0:,:) !< residual

    real(RNP), allocatable, save :: t(:)       ! subinterval time nodes
    real(RNP), allocatable, save :: f(:)       ! RHS
    real(RNP), allocatable, save :: r_c(:,:,:) ! convection term
    real(RNP), allocatable, save :: r_d(:,:,:) ! diffusion term
    real(RNP), allocatable, save :: f_s(:,:,:) ! sources
    real(RNP), allocatable, save :: bv(:,:)    ! boundary values

    real(RNP) :: dt_step
    integer   :: nc, ns, nt, ps, pt, mt, ot
    integer   :: e, i, k, m, n

    associate( cl_problem  => this % cl_problem  &
             , cl_operator => this % cl_operator &
             , cl_sdc      => this % cl_sdc      )

      ! initialization .........................................................

      ! dimensions
      ps = ubound(u, 1)
      ns = ubound(u, 2)
      nc = ubound(u, 3)
      mt = ubound(u, 4)
      nt = ubound(u, 5)

      ! polynomial degree in time and offset collocation points
      select case(this % cl_sdc % point_set)
      case('RR')
        ! Radau-right
        pt = mt - 1
        ot = 1
      case('E','L')
        ! equidistant of Lobatto points
        pt = mt
        ot = 0
      end select

      ! work space
      allocate(t(0:mt))
      allocate(f(0:ps))
      allocate(r_c(0:ps,ns,nc))
      allocate(r_d, mold = r_c)
      allocate(f_s, mold = r_c)
      allocate(bv(nc,2))

      dt_step = dt / nt

      ! residual
      associate(Me => cl_operator % Me)
        r(:,:,:,0,:) = 0
        do n = 1, nt
        do m = 1, mt
        do k = 1, nc
        do e = 1, ns
          if(present(G)) then
            r(:,e,k,m,n) = Me * (u(:,e,k,m-1,n) - u(:,e,k,m,n)) &
                         + G(:,e,k,m,n)
          else
            r(:,e,k,m,n) = Me * (u(:,e,k,m-1,n) - u(:,e,k,m,n))
          end if
        end do
        end do
        end do
        end do
      end associate

      ! evaluation .............................................................

      do n = 1, nt

        t = cl_sdc % SubintervalPoints(t_0 + (n-1)*dt_step, dt_step)

        do i = 0, mt

          associate(ui => u(:,:,:,i,n))
            call cl_problem % GetBoundaryValues (t(i), bv)
            call cl_problem % GetConvectionTerm (cl_operator, bv, ui, r_c)
            call cl_problem % GetDiffusionTerm  (cl_operator, bv, ui, r_d)
            call cl_problem % GetSources        (cl_operator, t(i), ui, f_s)
          end associate

          associate(Me => cl_operator % Me, w_sub => cl_sdc % w_sub)
            do k = 1, nc
            do e = 1, ns
              f = r_c(:,e,k) + r_d(:,e,k) + Me * f_s(:,e,k)
              do m = 1, mt
                r(:,e,k,m,n) = r(:,e,k,m,n) + dt_step * w_sub(i,m) * f
              end do
            end do
            end do
          end associate

        end do
      end do

      deallocate(t, f, r_c, r_d, f_s, bv)

    end associate

  end subroutine GetResidual

  !=============================================================================

end submodule MP_GetResidual
