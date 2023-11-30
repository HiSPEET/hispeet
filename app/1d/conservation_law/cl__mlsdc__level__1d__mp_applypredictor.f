!> summary:  Application of the SDC predictor
!> author:   Joerg Stiller
!> date:     2023/09/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_ApplyPredictor
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of the SDC predictor

  module subroutine ApplyPredictor(this, dt, t_0, u)
    class(CL_MLSDC_Level_1D), intent(in) :: this
    real(RNP), intent(in)    :: dt             !< size of the time slice
    real(RNP), intent(in)    :: t_0            !< start time of the slice
    real(RNP), intent(inout) :: u(0:,:,:,0:,:) !< approximate solution

    real(RNP), allocatable :: t_sub(:), dt_sub(:)
    real(RNP) :: dt_step
    integer   :: m, n

    associate( cl_sdc => this % cl_sdc         &
             , n_sub  => this % cl_sdc % n_sub &
             , n_time => this % n_time         )

      ! initialization .........................................................

      allocate(t_sub(0:n_sub), dt_sub(1:n_sub))

      dt_step = dt / n_time
      ! time steps ..............................................................

      do n = 1, n_time
        t_sub  = cl_sdc % IntermediateTimes(t_0 + (n-1)*dt_step, dt_step)
        dt_sub = t_sub(1:n_sub) - t_sub(0:n_sub)

        ! sweep thru subintervals
        do m = 1, n_sub
          call cl_sdc % predictor % TimeStep( this % cl_problem     &
                                            , this % cl_operator    &
                                            , dt  = dt_sub(m)       &
                                            , t_0 = t_sub(m-1)      &
                                            , u_0 = u(:,:,:,m-1,n)  &
                                            , u   = u(:,:,:,m  ,n)  )
        end do

        ! set initial values for next step
        if (n < n_time) then
          u(:,:,:,0,n+1) = u(:,:,:,n_sub,n)
        end if

      end do

    end associate

  end subroutine ApplyPredictor

  !=============================================================================

end submodule MP_ApplyPredictor
