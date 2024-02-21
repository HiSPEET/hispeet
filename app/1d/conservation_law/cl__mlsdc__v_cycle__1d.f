!> summary:  MLSDC V-cycle for 1D conservation laws
!> author:   Erik Pfister, Joerg Stiller
!> date:     2023/11/30
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

module CL__MLSDC__V_Cycle__1D
  use Kind_Parameters
  use CL__MLSDC__1D
  use CL__MLSDC__Variable__1D
  implicit none
  private

  public :: CL_MLSDC_V_Cycle_1D

contains

  !-----------------------------------------------------------------------------
  !> Execution of one MLSDC V-cycle
  !> Variable names:
  !> R before a variable name indicates a restricted variable
  !> I bevore a variable name indicates an interpolated variable
  !> Index f indicates the finer grid of the iteration
  !> Index c indicates the coarser grid of the iteration

  subroutine CL_MLSDC_V_Cycle_1D(mlsdc, t, dt, n_s1, n_s2, n_coarse, n_cycle, u)
    class(CL_MLSDC_1D), intent(in) :: mlsdc
    real(RNP), intent(in) :: t        !< start time
    real(RNP), intent(in) :: dt       !< thickness of time slab
    integer,   intent(in) :: n_s1     !< number of pre-smoothing sweeps
    integer,   intent(in) :: n_s2     !< number of post-smoothing sweeps
    integer,   intent(in) :: n_coarse !< number of sweeps for coarse solution
    integer,   intent(in) :: n_cycle  !< number of cycles to perform
    class(CL_MLSDC_Variable_1D), intent(inout) :: u   !< approximate solution

    ! internal variables .......................................................

    type(CL_MLSDC_Variable_1D), allocatable, save :: g   ! FAS correction
    type(CL_MLSDC_Variable_1D), allocatable, save :: r   ! residual
    type(CL_MLSDC_Variable_1D), allocatable, save :: v   ! auxiliary

    integer                :: c, l, l_top
    real(RNP), allocatable :: r_old(:,:,:,:,:)
    logical                :: abort

    ! initialization ...........................................................

    ! identify top level
    do l_top = 1, size(mlsdc % level)
      if (mlsdc % level(l_top) % is_top) exit
    end do

    g   = CL_MLSDC_Variable_1D(mlsdc)
    r   = CL_MLSDC_Variable_1D(mlsdc)
    v   = CL_MLSDC_Variable_1D(mlsdc)

    ! hacky Abbruchkriterium
    abort     = .FALSE.

    do c = 1, n_cycle

      ! fine to coarse ...........................................................

      do l = l_top, 2, -1
        associate( r_f  => r % level(l  ) % val &
                 , Rr_f => r % level(l-1) % val &
                 , g_f  => g % level(l  ) % val &
                 , g_c  => g % level(l-1) % val &
                 , v_c  => v % level(l-1) % val &
                 , u_f  => u % level(l  ) % val )
    
          ! set FAS correction to zero for top level
          if (l == l_top) then
            g_f = 0.
          end if

          ! pre-smoothing
          call mlsdc % level(l) % ApplyCorrector(dt, t, g_f, u_f, n_s1)

          ! get residual on fine grid
          ! r_f = L(u_f)
          call mlsdc % level(l) % ApplyOperator(dt, t, u_f, r_f)
          ! r_f = f^ - L(u_f)
          call mlsdc % level(l) % GetResidual(G=g_f, v=r_f, r=r_f)

          ! restrict fine solution
          call mlsdc % level(l) % Project_FC(u_f, v_c)

          ! where regular refinement condition
          ! g_c = L(v_c)
          call mlsdc % level(l-1) % ApplyOperator(dt, t, v_c, g_c)
          ! g_c = 0 - L(v_c)
          call mlsdc % level(l-1) % GetResidual(v=g_c, r=g_c)

          ! restrict residual to coarse level
          call mlsdc % level(l-1) % Restrict_FC(r_f, Rr_f)

          ! compute FAS correction
          g_c = Rr_f - g_c

        end associate
      end do

      ! hacky Abbruchkriterium  -> auf Multilevel erweitern?
      associate( Rr_f => r % level(1) % val )
        if (allocated(r_old)) then
          if (maxval(Rr_f) > maxval(r_old)) then
            abort = .TRUE.
          endif
          deallocate(r_old)
        end if
        allocate(r_old, mold=Rr_f)
        r_old = Rr_f
      end associate

      ! coarse solution ..........................................................

      associate( g_c => g % level(1) % val &
               , u_c => u % level(1) % val )

        call mlsdc % level(1) % ApplyCorrector(dt, t, g_c, u_c, n_coarse)

      end associate

      ! coarse to fine ...........................................................

      do l = 2, l_top
        associate( g    => g % level(l  ) % val &
                 , v_cr => r % level(l-1) % val &
                 , Iv_c => r % level(l  ) % val &
                 , v_c  => v % level(l-1) % val &
                 , u_f  => u % level(l  ) % val &
                 , u_c  => u % level(l-1) % val )

          ! where regular refinement condition

          ! calculate correction v_cr
          v_cr = u_c - v_c

          ! interpolate to finer grid
          call mlsdc % level(l-1) % Interpolate_CF(v_cr, Iv_c, complete=.true.)
          ! hacky Abbruchkriterium
          ! turn off manually for tests
          abort = .false.
          if(.not. abort) then
            u_f = u_f + Iv_c
          end if

          if (l /= l_top .or. (l == l_top .and. c == n_cycle)) then
            call mlsdc % level(l) % ApplyCorrector(dt, t, g, u_f, n_s2)
          endif

        end associate
      end do
      
    ! finalization .............................................................

    end do

    deallocate(g, r, v)

  end subroutine CL_MLSDC_V_Cycle_1D

end module CL__MLSDC__V_Cycle__1D
