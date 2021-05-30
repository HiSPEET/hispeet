! PROVISORIUM
program DG_Conservation_1D
  use Kind_Parameters
  use Constants
  use IP_Element_Operators__1D
  use Scalar_Problem_1D
  use Burgers_Problem_1D__Breaking_Wave
  implicit none

  class(ScalarProblem1D), allocatable :: problem
  type(IP_ElementOptions_1D) :: dg_opt
  real(RNP), allocatable    :: u(:,:,:)   ! discrete solution

  integer   :: po         =  15
  integer   :: po_cut_svv
  integer   :: ne         =  10
  real(RNP) :: penalty    =  2
  real(RNP) :: cfl        =  0.1
  real(RNP) :: t_end      =  0.35
  logical   :: svv        = .true.

  namelist /Discretization_Parameters/ po, po_cut_svv, ne, penalty, cfl,  &
                                       t_end, svv

  logical :: exists
  integer :: prm

  real(RNP) :: dt, t
  integer   :: i, k, ou

  print *, 'hello'

  po_cut_svv = dg_opt % po_cut_svv
  inquire(file='dg_conservation_1d.prm', exist=exists)
  if (exists) then
    open(newunit=prm, file='dg_conservation_1d.prm', action='READ')
    read(prm, nml=Discretization_Parameters)
    close(prm)
  end if

  dg_opt = IP_ElementOptions_1D( po         = po         &
                              , penalty    = penalty    &
                              , hybrid     = .true.     &
                              , svv        = svv        &
                              , po_cut_svv = po_cut_svv )

  allocate(BurgersProblem1D_BreakingWave :: problem)
  call problem % SetProblem('burgers_problem_1d__breaking_wave')
  call problem % SetSpaceDiscretization(dg_opt, ne)

  allocate(u(0:po, ne, problem%nc))
  u(0:,:,:) = problem % InitialValues()

  dt = min(cfl   *  problem%dx     /  po**2,                                  &
           cfl/2 * (problem%dx)**2 / (po**4 * (problem%nu_0r + problem%nu_0s)))

  k = 1
  t = 0
  do
    call TVDRK(problem, t, dt, u)

    if (10 * t >= k * t_end) then
      print '(2X,I3,"%")', 10*k
      k = k + 1
    end if
    if (t >= t_end) exit
  end do

  ! save results
  open(newunit=ou, file='dg_conservation_1d.dat')
  write(ou,'(A)') '# x, u'
  do k = 1, ne
  do i = 0, po
    write(ou,'(10(ES17.10,1X))') problem % x(i,k), u(i,k,:)
  end do
  end do
  close(ou)

contains

  !-----------------------------------------------------------------------------
  !> 3rd order TVD Runge-Kutta method

  subroutine TVDRK(problem, t, dt, u)
    class(ScalarProblem1D), intent(in)    :: problem
    real(RNP),                    intent(inout) :: t
    real(RNP),                    intent(in)    :: dt
    real(RNP),                    intent(inout) :: u(0:,:,:)

    real(RNP), allocatable, dimension(:,:,:), save :: f, u0, u1, u2, M_inv
    real(RNP) :: c0, c1, c2, ct
    integer   :: e, k

    ! bounds
    associate( po => problem % eop % po &
             , ne => problem % ne       &
             , nc => problem % nc       )

      ! workspace
      if (.not. allocated(f)) then
        allocate(u0, u1, u2, f, M_inv, mold = u)
        do k = 1, nc
        do e = 1, ne
          M_inv(:,e,k) = (2 / problem % dx) / problem % eop % w
        end do
        end do
      end if

      ! initial condition
      u0 = u

      ! step 1
      ct = dt
      f  = problem % RHS_Convection(t, u0)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u0)    ! does not matter with periodic BC
      u1 = u0 + ct * M_inv * f

      ! step 2
      c0 = 0.75_RNP
      c1 = 0.25_RNP
      ct = 0.25_RNP * dt
      f  = problem % RHS_Convection(t, u1)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u1)    ! does not matter with periodic BC
      u2 = c0 * u0 + c1 * u1 + ct * M_inv * f

      ! step 3
      c0 = THIRD
      c2 = TWO * THIRD
      ct = TWO * THIRD * dt
      f  = problem % RHS_Convection(t, u2)  & ! t is the wrong time here, but
         + problem % RHS_Diffusion (t, u2)    ! does not matter with periodic BC
      u  = c0 * u0 + c2 * u2 + ct * M_inv * f

      t  = t + dt

    end associate

  end subroutine TVDRK

end program DG_Conservation_1D
