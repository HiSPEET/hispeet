!> summary:  Implicit-explicit Runge-Kutta methods
!> author:   Susanne Stimpert, Joerg Stiller
!> date:     2017/02/06
!> license:  Institute of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany
!>
!>### Implicit-explicit Runge-Kutta schemes
!>
!> This module provides IMEX Runge-Kutta methods of the form
!>
!>     c | a_im     c | a_ex
!>     ––+––––--    ––+––––-
!>       | bᵀ         | bᵀ
!>
!> where `a_im` is the diagonally implicit part and `a_ex` the explicit part.
!> They all possess the first-same-as-last property, `c(1) = 0` and `c(ns) = 1`,
!> where `ns` is the number of stages. The first stage is always explicit, i.e.
!> `a_im(1,:) = 0` while, generally, `a_im(:,1) = 0`. The methods with 6 and 8
!> stages achieve stage-order 2.
!>
!> For accessing a particular method, an instance of the type `IMEX_RK_Method`
!> needs to be initialized using the type-bound procedure `New`.
!>
!> The implemented methods are described in
!>
!>   *  C.A. Kennedy, M.H. Carpenter, Appl Numer Math 44 (2003) 139–181, and
!>   *  D. Cavaglieri, T. Bewley, J Comput Phys 286 (2015) 172–193
!>
!> @note
!> In the references, coefficients are given as fractions. The numerators and
!> denominators of these fractions get too large for being represented by 8-byte
!> integers. Therefore, they were converted into a decimal form, which remains
!> precise with 16-byte reals and yields reasonable approximations, when only
!> 8-byte reals are available.
!> @endnote
!===============================================================================

module IMEX_Runge_Kutta_Method
  use, intrinsic :: ISO_Fortran_Env, only: OUTPUT_UNIT
  use Kind_Parameters, only: RNP, RHP
  use Constants,       only: ZERO
  use Execution_Control
  implicit none
  private

  public :: IMEX_RK_Method

  !-----------------------------------------------------------------------------
  !> Type for keeping the Butcher tableau of an IMEX Runge-Kutta method

  type IMEX_RK_Method
    character(len=80)      :: name    = ' ' !< name of RK method
    integer                :: n_stage = 0   !< number of stages
    integer                :: order   = 0   !< order of convergence
    real(RNP), allocatable :: a_im(:,:)     !< implicit RK matrix
    real(RNP), allocatable :: a_ex(:,:)     !< explicit RK matrix
    real(RNP), allocatable :: b(:)          !< RK weights
    real(RNP), allocatable :: c(:)          !< RK nodes
  contains
    procedure :: Init_IMEX_RK_Method
    procedure :: Show => Show_IMEX_RK_Method
  end type IMEX_RK_Method

  ! constructor
  interface IMEX_RK_Method
    module procedure New_IMEX_RK_Method
  end interface

contains

!-------------------------------------------------------------------------------
!> New IMEX_RK_Method

type(IMEX_RK_Method) function New_IMEX_RK_Method(ns, method) result(this)
  integer,           intent(in) :: ns     !< number of stages
  integer, optional, intent(in) :: method !< RK method [1]

  call Init_IMEX_RK_Method(this, ns, method)

end function New_IMEX_RK_Method

!-------------------------------------------------------------------------------
!> Initialize IMEX Butcher tableau
!>
!> The optional argument `method` allows to select from different RK methods
!> with the same number of stages `ns`.

subroutine Init_IMEX_RK_Method(this, ns, method)
  class(IMEX_RK_Method), intent(inout) :: this
  integer,               intent(in)    :: ns     !< number of stages
  integer,     optional, intent(in)    :: method !< RK scheme [1]

  integer :: method_

  if (present(method)) then
    method_ = method
  else
    method_ = 1
  end if

  ! set up components ..........................................................

  this % n_stage = min(8, max(2, ns))

  allocate( this % a_im ( this%n_stage, this%n_stage ), source = ZERO )
  allocate( this % a_ex ( this%n_stage, this%n_stage ), source = ZERO )
  allocate( this % b    ( this%n_stage )              , source = ZERO )
  allocate( this % c    ( this%n_stage )              , source = ZERO )

  ! select method ..............................................................

  select case(this % n_stage)

  case(2)

    this % name = 'Euler backward-forward'
    this % order = 1

    this % c(2) = 1
    this % b(2) = 1

    this % a_im(2,2) = 1
    this % a_ex(2,1) = 1

  case(3)

    this % name = 'IMEXRKCB2 (Cavaglieri & Bewley, JCP 286, 2015)'
    this % order = 2

    this % c(2) = real( 2._RHP / 5._RHP, RNP )
    this % c(3) = real( 1._RHP         , RNP )

    this % b(2) = real( 5._RHP / 6._RHP, RNP )
    this % b(3) = real( 1._RHP / 6._RHP, RNP )

    this % a_im(2,2) = real( 2._RHP / 5._RHP , RNP )
    this % a_im(3,:) = this % b

    this % a_ex(2,1) = real( 2._RHP / 5._RHP , RNP )
    this % a_ex(3,2) = real( 1._RHP          , RNP )

  case(4)

    select case(method_)

    case(1)

      this % name = 'IMEXRKCB3c (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 3

      this % c(2) = real( 337550982.9940_RHP / 452591907.6317_RHP , RNP )
      this % c(3) = real(  27277862.3835_RHP / 103945477.8728_RHP , RNP )
      this % c(4) = real(         1.0000_RHP                      , RNP )

      this % b(2) = real(  67348865.2607_RHP /     &
                          233403321.9546_RHP , RNP )
      this % b(3) = real(  49380121.9040_RHP /     &
                           85365302.6979_RHP , RNP )
      this % b(4) = real(  18481477.7513_RHP /     &
                          138966872.3319_RHP , RNP )

      this % a_im(2,2) = real(     337550982.99400000000_RHP /     &
                                   452591907.63170000000_RHP , RNP )
      this % a_im(3,2) = real( -117123838886.07531889907_RHP /     &
                                326945704956.02105556248_RHP , RNP )
      this % a_im(3,3) = real(      56613830.78810000000_RHP /     &
                                    91215372.11390000000_RHP , RNP )

      this % a_im(3,1) = this % b(1)
      this % a_im(4,:) = this % b

      this % a_ex(2,1) = this % c(2)
      this % a_ex(3,2) = this % c(3)
      this % a_ex(3,1) = this % b(1)
      this % a_ex(4,1) = this % b(1)
      this % a_ex(4,2) = this % b(2)
      this % a_ex(4,3) = real( 166054456.6939_RHP /233403321.9546_RHP , RNP )

   case(2)

      this % name = 'IMEXRKCB3e (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 3

      this % c(2) = real( 1._RHP / 3._RHP, RNP )
      this % c(3) = real( 1._RHP         , RNP )
      this % c(4) = real( 1._RHP         , RNP )

      this % b(2) = real( 3._RHP / 4._RHP, RNP )
      this % b(3) = real(-1._RHP / 4._RHP, RNP )
      this % b(4) = real( 1._RHP / 2._RHP, RNP )

      this % a_im(2,2) = real( 1._RHP / 3._RHP, RNP )
      this % a_im(3,2) = real( 1._RHP / 2._RHP, RNP )
      this % a_im(3,3) = real( 1._RHP / 2._RHP, RNP )
      this % a_im(4,:) = this % b

      this % a_ex(2,1) = real( 1._RHP / 3._RHP, RNP )
      this % a_ex(3,2) = real( 1._RHP         , RNP )
      this % a_ex(4,2) = real( 3._RHP / 4._RHP, RNP )
      this % a_ex(4,3) = real( 1._RHP / 4._RHP, RNP )

    case default
      call Error( 'Init_IMEX_RK_Method',            &
                  'requested method not available', &
                  'IMEX_Runge_Kutta_Method'         )
    end select

  case(6)

    this % name = 'IMEXRKCB4 (Cavaglieri & Bewley, JCP 286, 2015)'
    this % order = 4

    this % c(1) = ZERO
    this % c(2) = real( 1.0_RHP / 4.0_RHP , RNP )
    this % c(3) = real( 3.0_RHP / 4.0_RHP , RNP )
    this % c(4) = real( 3.0_RHP / 8.0_RHP , RNP )
    this % c(5) = real( 1.0_RHP / 2.0_RHP , RNP )
    this % c(6) = 1.0_RNP

    this % b(1) = real(  23.204908458700_RHP  / 137.713063006300_RHP, RNP )
    this % b(2) = real(   0.322009889509_RHP  /   2.243393849156_RHP, RNP )
    this % b(3) = real( -19.510967278700_RHP  / 123.316554581700_RHP, RNP )
    this % b(4) = real( -34.058241676100_RHP  /  70.541883231900_RHP, RNP )
    this % b(5) = real(  46.339607566100_RHP  /  40.997214447700_RHP, RNP )
    this % b(6) = real(  32.317794329400_RHP  / 162.664658063300_RHP, RNP )

    ! A_im . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    this % a_im(2,1) = real(  1.0000000000_RHP /  8.0000000000_RHP , RNP )
    this % a_im(3,1) = real( 21.6145252607_RHP / 96.1230882893_RHP , RNP )
    this % a_im(4,1) = this % b(1)
    this % a_im(5,1) = this % b(1)
    this % a_im(6,1) = this % b(1)

    this % a_im(2,2) = real(   1.0000000000_RHP /   8.0000000000_RHP , RNP )
    this % a_im(3,2) = real(  25.7479850128_RHP / 114.3310606989_RHP , RNP )
    this % a_im(4,2) = real( -38.1180097479_RHP / 127.6440792700_RHP , RNP )
    this % a_im(5,2) = this % b(2)
    this % a_im(6,2) = this % b(2)

    this % a_im(3,3) = real(   3.0481561667_RHP / 10.1628412017_RHP , RNP )
    this % a_im(4,3) = real(  -5.4660926949_RHP / 46.1115766612_RHP , RNP )
    this % a_im(5,3) = real( -10.0836174740_RHP / 86.1952129159_RHP , RNP )
    this % a_im(6,3) = this % b(3)

    this % a_im(4,4) = real(  34.4309628413_RHP /  55.2073727558_RHP , RNP )
    this % a_im(5,4) = real( -25.0423827953_RHP / 128.3875864443_RHP , RNP )
    this % a_im(6,4) = this % b(4)

    this % a_im(5,5) = 0.5_RNP
    this % a_im(6,5) = this % b(5)

    this % a_im(6,6) = this % b(6)

    ! A_ex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    this % a_ex(2,1) = real(  1.0000000000_RHP /   4.0000000000_RHP , RNP )
    this % a_ex(3,1) = real( 15.3985248130_RHP / 100.4999853329_RHP , RNP )
    this % a_ex(4,1) = this % b(1)
    this % a_ex(5,1) = this % b(1)
    this % a_ex(6,1) = this % b(1)

    this % a_ex(3,2) = real( 90.2825336800_RHP / 151.2825644809_RHP , RNP )
    this % a_ex(4,2) = real(  9.9316866929_RHP /  82.0744730663_RHP , RNP )
    this % a_ex(5,2) = this % b(2)
    this % a_ex(6,2) = this % b(2)

    this % a_ex(4,3) = real(  8.2888780751_RHP /  96.9573940619_RHP , RNP )
    this % a_ex(5,3) = real(  5.7501241309_RHP /  76.5040883867_RHP , RNP )
    this % a_ex(6,3) = this % b(3)

    this%a_ex(5,4) = real(    7.6345938311_RHP /  67.6824576433_RHP , RNP )
    this%a_ex(6,4) = real( -409.9309936455_RHP / 631.0162971841_RHP , RNP )

    this%a_ex(6,5) = real(  139.5992540491_RHP /  93.3264948679_RHP , RNP )

 case(8)

    this % name = 'ARK5(4)8L[2]SA (Kennedy & Carpenter, Appl Numer Math 44, 2003)'
    this % order = 5

    this % c(1) = ZERO
    this % c(2) = real(   4.1_RHP            /  10.0_RHP            , RNP )
    this % c(3) = real(  29.35347310677_RHP  / 112.92855782101_RHP  , RNP )
    this % c(4) = real(  14.26016391358_RHP  /  71.96633302097_RHP  , RNP )
    this % c(5) = real(   9.2_RHP            /  10.0_RHP            , RNP )
    this % c(6) = real(   2.4_RHP            /  10.0_RHP            , RNP )
    this % c(7) = real(   3.0_RHP            /   5.0_RHP            , RNP )
    this % c(8) = 1.0_RNP

    this % b(1) = real( -87.2700587467_RHP   / 913.3579230613_RHP   , RNP )
    this % b(2) = ZERO
    this % b(3) = ZERO
    this % b(4) = real( 223.482180632610_RHP /  95.558587375310_RHP , RNP )
    this % b(5) = real( -11.433695189920_RHP /  81.418160029310_RHP , RNP )
    this % b(6) = real( -39.379526789629_RHP /  19.018526304540_RHP , RNP )
    this % b(7) = real(  32.727382324388_RHP /  42.900044865799_RHP , RNP )
    this % b(8) = real(  4.1_RHP             /  20.0_RHP            , RNP )

    ! A_im . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    this % a_im(2,1) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(3,1) = real(   4.1_RHP            /   40.0_RHP            , RNP )
    this % a_im(4,1) = real(  68.37856364310_RHP  /  925.29203076860_RHP  , RNP )
    this % a_im(5,1) = real(  30.16520224154_RHP  /  100.81342136671_RHP  , RNP )
    this % a_im(6,1) = real(  21.88664790290_RHP  /  148.99783939110_RHP  , RNP )
    this % a_im(7,1) = real(  10.20004230633_RHP  /   57.15676835656_RHP  , RNP )
    this % a_im(8,1) = real( -87.27005874670_RHP  /  913.35792306130_RHP  , RNP )

    this % a_im(2,2) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(3,2) = real( -56.7603406766_RHP   / 1193.1857230679_RHP   , RNP )

    this % a_im(3,3) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(4,3) = real( -11.038504710300_RHP /  136.701519337300_RHP , RNP )
    this % a_im(5,3) = real(  30.586259806659_RHP /   12.414158314087_RHP , RNP )
    this % a_im(6,3) = real(  63.825689466800_RHP /  543.644631884100_RHP , RNP )
    this % a_im(7,3) = real(  25.762820946817_RHP /   25.263940353407_RHP , RNP )

    this % a_im(4,4) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(5,4) = real( -22.760509404356_RHP /   11.113319521817_RHP , RNP )
    this % a_im(6,4) = real( -11.797104745550_RHP /   53.211547248960_RHP , RNP )
    this % a_im(7,4) = real( -21.613759091450_RHP /   97.559073359090_RHP , RNP )
    this % a_im(8,4) = real( 223.482180632610_RHP /   95.558587375310_RHP , RNP )

    this % a_im(5,5) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(6,5) = real( -60.92811917200_RHP  / 8023.46106767100_RHP  , RNP )
    this % a_im(7,5) = real( -21.12173095930_RHP  /  584.68595025340_RHP  , RNP )
    this % a_im(8,5) = real( -11.43369518992_RHP  /   81.41816002931_RHP  , RNP )

    this % a_im(6,6) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(7,6) = real( -42.699250595730_RHP /   78.270590407490_RHP , RNP )
    this % a_im(8,6) = real( -39.379526789629_RHP /   19.018526304540_RHP , RNP )

    this % a_im(7,7) = real(   4.1_RHP            /   20.0_RHP            , RNP )
    this % a_im(8,7) = real(  32.727382324388_RHP /   42.900044865799_RHP , RNP )

    this % a_im(8,8) = real(   4.1_RHP            /   20.0_RHP            , RNP )

    ! A_ex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    this % a_ex(2,1) = real(    4.1_RHP            /  10.0_RHP            , RNP )
    this % a_ex(3,1) = real(   36.790274446400_RHP / 207.228047367700_RHP , RNP )
    this % a_ex(4,1) = real(   12.680235234080_RHP / 103.408227345210_RHP , RNP )
    this % a_ex(5,1) = real(  144.632819003510_RHP /  63.153537034770_RHP , RNP )
    this % a_ex(6,1) = real(   14.090043504691_RHP /  34.967701212078_RHP , RNP )
    this % a_ex(7,1) = real(   19.230459214898_RHP /  13.134317526959_RHP , RNP )
    this % a_ex(8,1) = real(  -19.977161125411_RHP /  11.928030595625_RHP , RNP )

    this % a_ex(3,2) = real(   67.762320755100_RHP / 822.414386656300_RHP , RNP )

    this % a_ex(4,3) = real(   10.299339394170_RHP / 136.365588504790_RHP , RNP )
    this % a_ex(5,3) = real(  661.144352112120_RHP /  58.794905890930_RHP , RNP )
    this % a_ex(6,3) = real(   15.191511035443_RHP /  11.219624916014_RHP , RNP )
    this % a_ex(7,3) = real(  212.753313583030_RHP /  29.424553649710_RHP , RNP )
    this % a_ex(8,3) = real( -407.959767960540_RHP /  63.849078235390_RHP , RNP )

    this % a_ex(5,4) = real( -540.531701528390_RHP /  42.847980215620_RHP , RNP )
    this % a_ex(6,4) = real(  -18.461159152457_RHP /  12.425892160975_RHP , RNP )
    this % a_ex(7,4) = real( -381.453459884190_RHP /  48.626203187230_RHP , RNP )
    this % a_ex(8,4) = real(  177.454434618887_RHP /  12.078138498510_RHP , RNP )

    this % a_ex(6,5) = real(  -28.166716381100_RHP / 901.161929587000_RHP , RNP )
    this % a_ex(7,5) = real(   -1.0_RHP            /   8.0_RHP            , RNP )
    this % a_ex(8,5) = real(   78.267220542500_RHP / 826.770190026100_RHP , RNP )

    this % a_ex(7,6) = real(   -1.0_RHP            /   8.0_RHP            , RNP )
    this % a_ex(8,6) = real( -695.630110598110_RHP /  96.465806942050_RHP , RNP )

    this % a_ex(8,7) = real(   73.566282105260_RHP /  49.421867764050_RHP , RNP )

  case default

    call Error( 'Init_IMEX_RK_Method',                      &
                'requested number of stages not supported', &
                'IMEX_Runge_Kutta_Method'                   )

  end select

end subroutine Init_IMEX_RK_Method

!-------------------------------------------------------------------------------
!> Deletes the given `IMEX_RK_Method` object.

subroutine Delete_IMEX_RK_Method(this)
  type(IMEX_RK_Method), intent(inout) :: this

  if (allocated( this % a_im )) deallocate( this % a_im )
  if (allocated( this % a_ex )) deallocate( this % a_ex )
  if (allocated( this % b    )) deallocate( this % b    )
  if (allocated( this % c    )) deallocate( this % c    )

end subroutine Delete_IMEX_RK_Method

!===============================================================================

subroutine Show_IMEX_RK_Method(this, unit)
  class(IMEX_RK_Method), intent(in) :: this
  integer,     optional, intent(in) :: unit  !< output unit

  character(len=*), parameter :: fmt_ca = '(F13.10," |",99F14.10)'
  character(len=*), parameter :: fmt_b =  '(13X,   " |",99F14.10)'
  integer :: i, io

  if (this % n_stage < 1) return

  if (present(unit)) then
    io = unit
  else
    io = OUTPUT_UNIT
  end if

  write(io,'(/,A,/)') 'IMEX Runge-Kutta method'
  write(io,'(2A,/)')  'name: ', trim(this % name)
  write(io,'(A,I0)')  'stages = ', this % n_stage
  write(io,'(A,I0)')  'order  = ', this % order

  write(io,'(/,A,/)') 'implicit part'
  do i = 1, this%n_stage
    write(io,fmt_ca) this % c(i), this % a_im(i,1:i)
  end do
  write(io,'(A)') repeat('-', 16 + 14*this%n_stage)
  write(io,fmt_b) this % b

  write(io,'(/,A,/)') 'explicit part'
  do i = 1, this%n_stage
    write(io,fmt_ca) this % c(i), this % a_ex(i,1:i-1)
  end do
  write(io,'(A)') repeat('-', 16 + 14*this%n_stage)
  write(io,fmt_b) this % b
  write(io,*)

end subroutine Show_IMEX_RK_Method

!===============================================================================

end module IMEX_Runge_Kutta_Method
