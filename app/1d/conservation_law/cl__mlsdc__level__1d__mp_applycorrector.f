!> summary:  Application of the SDC corrector
!> author:   Joerg Stiller
!> date:     2023/09/05
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!===============================================================================

submodule(CL__MLSDC__Level__1D) MP_ApplyCorrector
  implicit none

contains

  !-----------------------------------------------------------------------------
  !> Application of the SDC corrector

  module subroutine ApplyCorrector(this, dt, t_0, G, u, n_sweep)
    class(CL_MLSDC_Level_1D), intent(in) :: this
    real(RNP), intent(in)    :: dt             !< size of the time slice
    real(RNP), intent(in)    :: t_0            !< start time of the slice
    real(RNP), intent(in)    :: G(0:,:,:,0:,:) !< FAS defect correction
    real(RNP), intent(inout) :: u(0:,:,:,0:,:) !< approximate solution
    integer,   intent(in)    :: n_sweep        !< number of sweeps

    real(RNP), allocatable, save :: F        (:,:,:,:) ! [Fᵢ]ᵏ
    real(RNP), allocatable, save :: F_ex     (:,:,:,:) ! [F_exᵢ]ᵏ
    real(RNP), allocatable, save :: F_im     (:,:,:,:) ! [F_imᵢ]ᵏ
    real(RNP), allocatable, save :: F_ex_new (:,:,:,:) ! [F_exᵢ]ᵏ⁺¹
    real(RNP), allocatable, save :: F_im_new (:,:,:,:) ! [F_imᵢ]ᵏ⁺¹

    real(RNP), allocatable :: t_sub(:), dt_sub(:)
    real(RNP) :: dt_step
    integer   :: k, m, n

    if (n_sweep < 1) return

    associate( nc          => this % cl_problem  % nc       &
             , po          => this % cl_operator % eop % po &
             , ne          => this % cl_operator % ne       &
             , eop         => this % cl_operator % eop      &
             , n_sub       => this % cl_sdc % n_sub         &
             , n_time      => this % n_time                 &
             , cl_problem  => this % cl_problem             &
             , cl_operator => this % cl_operator            &
             , cl_sdc      => this % cl_sdc                 )

      ! initialization .........................................................

      allocate(F        (0:po,ne,nc,0:n_sub))
      allocate(F_ex     (0:po,ne,nc,0:n_sub))
      allocate(F_im     (0:po,ne,nc,0:n_sub))
      allocate(F_ex_new (0:po,ne,nc,0:n_sub))
      allocate(F_im_new (0:po,ne,nc,0:n_sub))

      allocate(t_sub(0:n_sub), dt_sub(1:n_sub))

      dt_step = dt / n_time

      ! time steps ..............................................................

      Steps: do n = 1, n_time

        t_sub  = cl_sdc % SubintervalPoints(t_0 + (n-1)*dt_step, dt_step)
        dt_sub = t_sub(1:n_sub) - t_sub(0:n_sub)

        ! prerequisites for SDC sweeps
        do m = 0, n_sub

          call cl_sdc % GetHighOrderRHS( cl_problem         &
                                       , cl_operator        &
                                       , t  = t_sub(m)      &
                                       , u  = u (:,:,:,m,n) &
                                       , F  = F (:,:,:,m)   )

          call cl_sdc % GetCorrectorRHS( cl_problem                         &
                                       , cl_operator                        &
                                       , t    = t_sub  (m)                  &
                                       , dt   = dt_sub (max(m,1))           &
                                       , u_0  = u      (:,:,:,max(m-1,0),n) &
                                       , u    = u      (:,:,:,m,n)          &
                                       , F_ex = F_ex   (:,:,:,m)            &
                                       , F_im = F_im   (:,:,:,m)            )

        end do

        ! initialization of new corrector RHS
        call SetArray(F_ex_new(:,:,:,0), F_ex(:,:,:,0), multi = .true.)
        call SetArray(F_im_new(:,:,:,0), F_im(:,:,:,0), multi = .true.)

        ! SDC sweeps
        Sweeps: do k = 1, n_sweep

          do m = 1, n_sub
            call cl_sdc % CorrectorStep( cl_problem              &
                                       , cl_operator             &
                                       , m        = m            &
                                       , t        = t_sub        &
                                       , u        = u(:,:,:,:,n) &
                                       , F        = F            &
                                       , F_ex     = F_ex         &
                                       , F_im     = F_im         &
                                       , F_ex_new = F_ex_new     &
                                       , F_im_new = F_im_new     &
                                       , G        = G(:,:,:,:,n) )
          end do

          if (k == n_sweep) exit

          do m = 1, n_sub

            call cl_sdc % GetHighOrderRHS( cl_problem         &
                                         , cl_operator        &
                                         , t  = t_sub(m)      &
                                         , u  = u (:,:,:,m,n) &
                                         , F  = F (:,:,:,m)   )

            call SetArray(F_ex(:,:,:,m), F_ex_new(:,:,:,m), multi = .true.)
            call SetArray(F_im(:,:,:,m), F_im_new(:,:,:,m), multi = .true.)

          end do

        end do Sweeps

        ! set initial values for next step
        if (n < n_time) then
          u(:,:,:,0,n+1) = u(:,:,:,n_sub,n)
        end if

      end do Steps

      ! clean-up ...............................................................

      deallocate(F, F_ex, F_im, F_ex_new, F_im_new)

    end associate

  end subroutine ApplyCorrector

  !=============================================================================

end submodule MP_ApplyCorrector
