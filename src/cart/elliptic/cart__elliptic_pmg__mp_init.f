!> summary:  Polynomial multigrid initialization routines
!> author:   Joerg Stiller
!> date:     2019/02/03
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Polynomial multigrid for use with elliptic solvers
!===============================================================================

submodule(CART__Elliptic_PMG) MP_Init
  implicit none

contains

!-------------------------------------------------------------------------------
!> PMG_Method3D initialization: IP, constant isotropic diffusivity

module subroutine Init__IP_CI( this, mesh, lambda, nu, bc,  &
                               ip_opt, schwarz_opt, pmg_opt )

  class(PMG_Method3D), intent(inout) :: this

  class(MeshPartition), target, intent(in) :: mesh        !< mesh partition
  real(RNP),                    intent(in) :: lambda      !< Helmholtz parameter
  real(RNP),                    intent(in) :: nu          !< diffusivity
  character,                    intent(in) :: bc(:)       !< boundary conditions
  class(IP_ElementOptions1D),   intent(in) :: ip_opt      !< IP/DG options
  class(SchwarzOptions3D),      intent(in) :: schwarz_opt !< Schwarz options
  class(PMG_Options3D),         intent(in) :: pmg_opt     !< PMG options

  class(IP_ElementOptions1D), allocatable :: ip_opt_l
  integer, allocatable :: po(:)
  integer :: l, l_top, ns1, ns2

  allocate(ip_opt_l, source = ip_opt)

  this % mesh => mesh
  call Assign_PMG_Options(pmg_opt, this)
  call CreatePolynomialLevels(pmg_opt, po)

  l_top = ubound(po,1)
  allocate(this % level(0:l_top))
  ns1 = pmg_opt % ns1
  ns2 = pmg_opt % ns2

  do l = l_top, 0, -1

    ip_opt_l % po = po(l)

    call this % level(l) % Init( bottom       =  l == 0       &
                               , top          =  l == l_top   &
                               , ns1          =  ns1          &
                               , ns2          =  ns2          &
                               , mesh         =  mesh         &
                               , lambda       =  lambda       &
                               , nu           =  nu           &
                               , bc           =  bc           &
                               , ip_opt       =  ip_opt_l     &
                               , schwarz_opt  =  schwarz_opt  &
                               )
    ns1 = ns1 * pmg_opt % mvs
    ns2 = ns2 * pmg_opt % mvs
  end do

  allocate(this % transfer(1:l_top))
  do l = 1, l_top
    call this % transfer(l) % New(po(l-1), po(l))
  end do

end subroutine Init__IP_CI

!-------------------------------------------------------------------------------
!> PMG_Method3D initialization: IP, variable isotropic diffusivity

module subroutine Init__IP_VI( this, mesh, lambda, nu, bc,  &
                               ip_opt, schwarz_opt, pmg_opt )

  class(PMG_Method3D), intent(inout) :: this

  class(MeshPartition), target, intent(in) :: mesh           !< mesh partition
  real(RNP),                    intent(in) :: lambda         !< Helmholtz param
  real(RNP),                    intent(in) :: nu(0:,0:,0:,:) !< diffusivity
  character,                    intent(in) :: bc(:)          !< boundary cond
  class(IP_ElementOptions1D),   intent(in) :: ip_opt         !< IP/DG opt
  class(SchwarzOptions3D),      intent(in) :: schwarz_opt    !< Schwarz opt
  class(PMG_Options3D),         intent(in) :: pmg_opt        !< PMG options

  class(IP_ElementOptions1D), allocatable :: ip_opt_l
  integer, allocatable :: po(:)
  integer :: l, l_top, ns1, ns2

  allocate(ip_opt_l, source = ip_opt)
!### CHECK
print *, '$ 01'
!### CHECK END

  this % mesh => mesh
  call Assign_PMG_Options(pmg_opt, this)
  call CreatePolynomialLevels(pmg_opt, po)
!### CHECK
print *, '$ 02'
!### CHECK END

  l_top = ubound(po,1)
  allocate(this % level(0:l_top))
  ns1 = pmg_opt % ns1
  ns2 = pmg_opt % ns2

  do l = l_top, 0, -1

!### CHECK
print *, '$ 03: l =',l
!### CHECK END
    ip_opt_l % po = po(l)
!### CHECK
print *, '$ 04: l =',l
!### CHECK END

    call this % level(l) % Init( bottom       =  l == 0       &
                               , top          =  l == l_top   &
                               , ns1          =  ns1          &
                               , ns2          =  ns2          &
                               , mesh         =  mesh         &
                               , lambda       =  lambda       &
                               , nu           =  nu           &
                               , bc           =  bc           &
                               , ip_opt       =  ip_opt_l     &
                               , schwarz_opt  =  schwarz_opt  &
                               )
    ns1 = ns1 * pmg_opt % mvs
    ns2 = ns2 * pmg_opt % mvs
!### CHECK
print *, '$ 05: l =',l
!### CHECK END
  end do

  allocate(this % transfer(1:l_top))
  do l = 1, l_top
    call this % transfer(l) % New(po(l-1), po(l))
!### CHECK
print *, '$ 06: l =',l
!### CHECK END
  end do

end subroutine Init__IP_VI

!-------------------------------------------------------------------------------
!> Assignment of options to PMG_Method3D object

subroutine Assign_PMG_Options(opt, pmg)
  class(PMG_Options3D),       intent(in)    :: opt
  class(PMG_Method3D), intent(inout) :: pmg

  pmg % i_max   = opt % i_max
  pmg % r_red   = opt % r_red
  pmg % r_max   = opt % r_max
  pmg % dr_min  = opt % dr_min

  pmg % solver  = opt % solver
  pmg % i0_max  = opt % i0_max
  pmg % r0_red  = opt % r0_red
  pmg % monitor = opt % monitor

end subroutine Assign_PMG_Options

!-------------------------------------------------------------------------------
!> Creates a set of ascending polynomial orders

subroutine CreatePolynomialLevels(opt, po)
  class(PMG_Options3D), intent(in)  :: opt
  integer, allocatable, intent(out) :: po(:) !< polynomial levels

  integer, allocatable :: q(:)
  integer :: lb, lt

  associate( po_top => opt % po_top, cr     => opt % cr,    &
             po_bot => opt % po_bot, cr_max => opt % cr_max )

    if (po_top > po_bot) then

      allocate(q(po_bot:po_top), source = po_top)
      lt = po_top
      lb = lt
      do
        lb = lb -1
        q(lb) = max( nint(q(lb+1)/cr), ceiling(q(lb+1)/cr_max), 1 )
        q(lb) = min( q(lb), q(lb+1) - 1 )
        if (q(lb) <= max(po_bot, 1)) then
          q(lb) = max(q(lb), po_bot)
          exit
        end if
      end do
      allocate(po(0:lt-lb), source=q(lb:lt))

    else if (po_top == po_bot) then

      allocate(po(0:0), source=po_top)

    else

      allocate(po(0:-1))

    end if

  end associate

end subroutine CreatePolynomialLevels

!===============================================================================

end submodule MP_Init
