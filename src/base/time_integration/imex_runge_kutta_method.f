!------------------------------------------------------------------------------!
! This file is part of HiSPEET: High-order Spectral Element Techniques         !
!                                                                              !
! Copyright (C) 2026 by the HiSPEET authors and the                            !
! Chair of Fluid Mechanics, TU Dresden, 01062 Dresden, Germany                 !
!                                                                              !
! HiSPEET is free software: you can redistribute it and/or modify              !
! it under the terms of the GNU General Public License as published by         !
! the Free Software Foundation, either version 3 of the License, or            !
! (at your option) any later version.                                          !
!                                                                              !
! HiSPEET is distributed in the hope that it will be useful,                   !
! but WITHOUT ANY WARRANTY; without even the implied warranty of               !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                         !
! See the GNU General Public License for more details.                         !
!                                                                              !
! You should have received a copy of the GNU General Public License            !
! along with HiSPEET. If not, see <http://www.gnu.org/licenses/>.              !
!------------------------------------------------------------------------------!

!> summary:  Implicit-explicit Runge-Kutta methods
!> author:   Susanne Stimpert, Joerg Stiller
!> date:     2017/02/06
!>
!>### Implicit-explicit Runge-Kutta schemes
!>
!> This module provides IMEX Runge-Kutta methods of the form
!>
!>     c | a_im     c | a_ex
!>     ––|––––---   ––|––––--
!>       | b_imᵀ      | b_exᵀ
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
    real(RNP), allocatable :: c(:)          !< RK nodes
    real(RNP), allocatable :: b_im(:)       !< implicit RK weights
    real(RNP), allocatable :: b_ex(:)       !< explicit RK weights
    real(RNP), allocatable :: a_im(:,:)     !< implicit RK matrix
    real(RNP), allocatable :: a_ex(:,:)     !< explicit RK matrix
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
  real(RHP), allocatable :: c(:), b_im(:), b_ex(:), a_im(:,:), a_ex(:,:)

  if (present(method)) then
    method_ = method
  else
    method_ = 1
  end if

  ! set up components ..........................................................

  this % n_stage = min(8, max(2, ns))
  allocate( c    ( this%n_stage )              , source = 0.0_RHP )
  allocate( b_im ( this%n_stage )              , source = 0.0_RHP )
  allocate( b_ex ( this%n_stage )              , source = 0.0_RHP )
  allocate( a_im ( this%n_stage, this%n_stage ), source = 0.0_RHP )
  allocate( a_ex ( this%n_stage, this%n_stage ), source = 0.0_RHP )

  ! select method ..............................................................

  select case(this % n_stage)

  case(2)

    select case(method_)

    case(1)

      this % name  = 'Euler backward-forward'
      this % order = 1

      c(2)    = 1

      b_im(2) = 1
      b_ex(1) = 1

      a_im(2,2) = 1
      a_ex(2,1) = 1

    case(2)

      this % name  = 'IMEXRK22MP - IMEX mid-point rule'
      this % order = 2

      c(2)    = 1._RHP / 2._RHP

      b_im(2) = 1._RHP
      b_ex(:) = b_im

      a_im(2,2) = 1._RHP / 2._RHP
      a_ex(2,1) = 1._RHP / 2._RHP

    case(3)

      this % name  = 'IMEXRK22TR - IMEX trapezoidal rule'
      this % order = 2

      c(2) = 1

      b_im(1) = 1._RHP / 2._RHP
      b_im(2) = 1._RHP / 2._RHP

      b_ex(:) = b_im

      a_im(2,2) = 1
      a_ex(2,1) = 1

    case default

      call Error( 'Init_IMEX_RK_Method',            &
                  'requested method not available', &
                  'IMEX_Runge_Kutta_Method'         )

    end select

  case(3)

    select case(method_)

    case(1)

      this % name  = 'IMEXRKCB2 (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 2

      c(2) = 2._RHP / 5._RHP
      c(3) = 1._RHP

      b_im(2) = 5._RHP / 6._RHP
      b_im(3) = 1._RHP / 6._RHP

      b_ex(:) = b_im

      a_ex(2,1) = 2._RHP / 5._RHP
      a_ex(3,2) = 1._RHP

      a_im(2,2) = 2._RHP / 5._RHP
      a_im(3,:) = b_im

    case(2)

      this % name  = 'IMEXRK32TR - IMEX trapezoidal rule' ! Variant 2
      this % order = 2

      c(2) = 1._RHP
      c(3) = 1._RHP

      b_im(1) = 1._RHP / 2._RHP
      b_im(3) = 1._RHP / 2._RHP

      b_ex(1) = 1._RHP / 2._RHP
      b_ex(2) = 1._RHP / 2._RHP

      a_im(2,2) = 1._RHP
      a_im(3,:) = b_im

      a_ex(2,1) = 1._RHP
      a_ex(3,:) = b_ex

    case default

      call Error( 'Init_IMEX_RK_Method',            &
                  'requested method not available', &
                  'IMEX_Runge_Kutta_Method'         )

    end select

  case(4)

    select case(method_)

    case(1)

      this % name  = 'IMEXRKCB3c (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 3

      c(2) = 337550982.9940_RHP / 452591907.6317_RHP
      c(3) =  27277862.3835_RHP / 103945477.8728_RHP
      c(4) =         1.0000_RHP

      b_im(2) =  67348865.2607_RHP /     &
                233403321.9546_RHP
      b_im(3) =  49380121.9040_RHP /     &
                 85365302.6979_RHP
      b_im(4) =  18481477.7513_RHP /     &
                138966872.3319_RHP

      b_ex    = b_im

      a_im(2,2) =     337550982.99400000000_RHP /     &
                      452591907.63170000000_RHP
      a_im(3,2) = -117123838886.07531889907_RHP /     &
                   326945704956.02105556248_RHP
      a_im(3,3) =      56613830.78810000000_RHP /     &
                       91215372.11390000000_RHP

      a_im(3,1) = b_im(1)
      a_im(4,:) = b_im

      a_ex(2,1) = c(2)
      a_ex(3,2) = c(3)
      a_ex(3,1) = b_ex(1)
      a_ex(4,1) = b_ex(1)
      a_ex(4,2) = b_ex(2)
      a_ex(4,3) = 166054456.6939_RHP /233403321.9546_RHP

   case(2)

      this % name = 'IMEXRKCB3e (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 3

      c(2) =  1._RHP / 3._RHP
      c(3) =  1._RHP
      c(4) =  1._RHP

      b_im(2) =  3._RHP / 4._RHP
      b_im(3) = -1._RHP / 4._RHP
      b_im(4) =  1._RHP / 2._RHP

      b_ex = b_im

      a_im(2,2) = 1._RHP / 3._RHP
      a_im(3,2) = 1._RHP / 2._RHP
      a_im(3,3) = 1._RHP / 2._RHP
      a_im(4,:) = b_im

      a_ex(2,1) = 1._RHP / 3._RHP
      a_ex(3,2) = 1._RHP
      a_ex(4,2) = 3._RHP / 4._RHP
      a_ex(4,3) = 1._RHP / 4._RHP

   case(3)

      this % name = 'CN/RKW3 (Cavaglieri & Bewley, JCP 286, 2015)'
      this % order = 2

      c(2) =  8._RHP / 15._RHP
      c(3) =  2._RHP /  3._RHP
      c(4) =  1._RHP

      b_im(1) = 4._RHP / 15._RHP
      b_im(2) = 1._RHP /  3._RHP
      b_im(3) = 7._RHP / 30._RHP
      b_im(4) = 1._RHP /  6._RHP

      b_ex = b_im

      a_im(2,1) = 4._RHP / 15._RHP
      a_im(2,2) = 4._RHP / 15._RHP
      a_im(3,1) = 4._RHP / 15._RHP
      a_im(3,2) = 1._RHP /  3._RHP
      a_im(3,3) = 1._RHP / 15._RHP
      a_im(4,:) = b_im

      a_ex(2,1) = 8._RHP / 15._RHP
      a_ex(3,1) = 1._RHP /  4._RHP
      a_ex(3,2) = 5._RHP / 12._RHP
      a_ex(4,1) = 1._RHP /  4._RHP
      a_ex(4,3) = 3._RHP /  4._RHP

    case default

      call Error( 'Init_IMEX_RK_Method',            &
                  'requested method not available', &
                  'IMEX_Runge_Kutta_Method'         )

    end select

  case(5)

    select case(method_)

    case(1)

      this % name  = 'IMEXRK53 ARS443 (Ascher et al., ANM 25, 1997)'
      this % order = 3

      c(2) = 1._RHP / 2._RHP
      c(3) = 2._RHP / 3._RHP
      c(4) = 1._RHP / 2._RHP
      c(5) = 1._RHP

      b_im(2) =  3._RHP / 2._RHP
      b_im(3) = -3._RHP / 2._RHP
      b_im(4) =  1._RHP / 2._RHP
      b_im(5) =  1._RHP / 2._RHP

      b_ex(1) =  1._RHP / 4._RHP
      b_ex(2) =  7._RHP / 4._RHP
      b_ex(3) =  3._RHP / 4._RHP
      b_ex(4) = -7._RHP / 4._RHP

      a_im(2,2) =  1._RHP / 2._RHP
      a_im(3,2) =  1._RHP / 6._RHP
      a_im(3,3) =  1._RHP / 2._RHP
      a_im(4,2) = -1._RHP / 2._RHP
      a_im(4,3) =  1._RHP / 2._RHP
      a_im(4,4) =  1._RHP / 2._RHP
      a_im(5,:) =  b_im

      a_ex(2,1) =  1._RHP /  2._RHP
      a_ex(3,1) = 11._RHP / 18._RHP
      a_ex(3,2) =  1._RHP / 18._RHP
      a_ex(4,1) =  5._RHP /  6._RHP
      a_ex(4,2) = -5._RHP /  6._RHP
      a_ex(4,3) =  1._RHP /  2._RHP
      a_ex(5,:) =  b_ex

    case(2)

      this % name  = 'BHR(5,5,3) (Boscarino, ANM 59, 2009)'
      this % order = 3

      block

        real(RHP), parameter :: gamma = 0.435866521508482_RHP

        ! c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        ! Following private the communication with Sebastiano Boscarino
        ! c(3) = 2*gamma is set to avoid inconsistency with a_ex(3,1) + a_ex(3,2)

        c(2) = 2 * gamma
        c(3) = 2 * gamma

      ! c(3) =  902905985686._RHP &  ! value given
      !      / 1035759735069._RHP    ! in the paper

        c(4) = 2684624._RHP &
             / 1147171._RHP

        c(5) = 1

        ! b  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        b_im(1) =  48769850233674.0678603511_RHP &
                / 118115963692818.5920260208_RHP

        b_im(3) =  30298776308118.4622639300143137943089_RHP &
                / 153535994420329.3318639180129368156500_RHP

        b_im(4) = -10523592833510.0616072938218863_RHP &
                / 228255445206466.1756575727198000_RHP

        b_im(5) =  gamma

        b_ex(:) =  b_im(:)

        ! A_im . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        a_im(2,1) = gamma
        a_im(2,2) = gamma
        a_im(3,1) = gamma
        a_im(3,3) = gamma
        a_im(4,4) = gamma

      ! Following private the communication with Sebastiano Boscarino
      ! a_im(3,2) = 0 is set to avoid inconsistency with c(3)

      ! a_im(3,2) =             -31.733082319927313_RHP & ! value given
      !           / 455705377221960.889379854647102_RHP   ! in the paper

        a_im(4,1) = -30123785410849.22027361996761794919360516301377809610_RHP &
                  / 451233940565852.69977907753045030512597955897345819349_RHP

        a_im(4,2) =             -62.865589297807153294268_RHP &
                  / 102559673441610.672305587327019095047_RHP

        a_im(4,3) = 418769796920855.299603146267001414900945214277000_RHP &
                  / 212454360385257.708555954598099874818603217167139_RHP

        a_im(5,:) = b_im(:)

        ! A_ex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

        a_ex(2,1) = gamma * 2
        a_ex(3,1) = gamma
        a_ex(3,2) = gamma

        a_ex(4,1) = -475883375220285.986033264_RHP &
                 /  594112726933437.845704163_RHP

        a_ex(4,2) =  0

        a_ex(4,3) =  186623344982202.6827708736_RHP &
                 /   59411272693343.7845704163_RHP

        a_ex(5,1) =   62828845818073.169585635881686091391737610308247_RHP &
                 /  176112910684412.105319781630311686343715753056000_RHP

        a_ex(5,2) = -b_im(3)

        a_ex(5,3) =  262315887293043.739337088563996093207_RHP &
                 /  297427554730376.353252081786906492000_RHP

        a_ex(5,4) =  -9876182318941.76581438124717087_RHP &
                  / 238773376602029.69319526901856000_RHP

      end block

    case default

      call Error( 'Init_IMEX_RK_Method',            &
                  'requested method not available', &
                  'IMEX_Runge_Kutta_Method'         )

    end select

  case(6)

    this % name  = 'IMEXRKCB4 (Cavaglieri & Bewley, JCP 286, 2015)'
    this % order = 4

    ! c  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    c(1) = 0.0_RHP
    c(2) = 1.0_RHP / 4.0_RHP
    c(3) = 3.0_RHP / 4.0_RHP
    c(4) = 3.0_RHP / 8.0_RHP
    c(5) = 1.0_RHP / 2.0_RHP
    c(6) = 1.0_RHP

    ! b  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    b_im(1) =  23.204908458700_RHP  / 137.713063006300_RHP
    b_im(2) =   0.322009889509_RHP  /   2.243393849156_RHP
    b_im(3) = -19.510967278700_RHP  / 123.316554581700_RHP
    b_im(4) = -34.058241676100_RHP  /  70.541883231900_RHP
    b_im(5) =  46.339607566100_RHP  /  40.997214447700_RHP
    b_im(6) =  32.317794329400_RHP  / 162.664658063300_RHP

    b_ex = b_im

    ! A_im . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    a_im(2,1) =  1.0000000000_RHP /  8.0000000000_RHP
    a_im(3,1) = 21.6145252607_RHP / 96.1230882893_RHP
    a_im(4,1) = b_im(1)
    a_im(5,1) = b_im(1)
    a_im(6,1) = b_im(1)

    a_im(2,2) =   1.0000000000_RHP /   8.0000000000_RHP
    a_im(3,2) =  25.7479850128_RHP / 114.3310606989_RHP
    a_im(4,2) = -38.1180097479_RHP / 127.6440792700_RHP
    a_im(5,2) = b_im(2)
    a_im(6,2) = b_im(2)

    a_im(3,3) =   3.0481561667_RHP / 10.1628412017_RHP
    a_im(4,3) =  -5.4660926949_RHP / 46.1115766612_RHP
    a_im(5,3) = -10.0836174740_RHP / 86.1952129159_RHP
    a_im(6,3) = b_im(3)

    a_im(4,4) =  34.4309628413_RHP /  55.2073727558_RHP
    a_im(5,4) = -25.0423827953_RHP / 128.3875864443_RHP
    a_im(6,4) = b_im(4)

    a_im(5,5) = 0.5_RHP
    a_im(6,5) = b_im(5)

    a_im(6,6) = b_im(6)

    ! A_ex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    a_ex(2,1) =  1.0000000000_RHP /   4.0000000000_RHP
    a_ex(3,1) = 15.3985248130_RHP / 100.4999853329_RHP
    a_ex(4,1) = b_ex(1)
    a_ex(5,1) = b_ex(1)
    a_ex(6,1) = b_ex(1)

    a_ex(3,2) = 90.2825336800_RHP / 151.2825644809_RHP
    a_ex(4,2) =  9.9316866929_RHP /  82.0744730663_RHP
    a_ex(5,2) = b_ex(2)
    a_ex(6,2) = b_ex(2)

    a_ex(4,3) = 8.2888780751_RHP /  96.9573940619_RHP
    a_ex(5,3) = 5.7501241309_RHP /  76.5040883867_RHP
    a_ex(6,3) = b_ex(3)

    a_ex(5,4) =    7.6345938311_RHP /  67.6824576433_RHP
    a_ex(6,4) = -409.9309936455_RHP / 631.0162971841_RHP

    a_ex(6,5) =  139.5992540491_RHP /  93.3264948679_RHP

 case(8)

    this % name  = 'ARK5(4)8L[2]SA (Kennedy & Carpenter, ANM 44, 2003)'
    this % order = 5

    c(1) =  0.0_RHP
    c(2) =  4.1_RHP            /  10.0_RHP
    c(3) = 29.35347310677_RHP  / 112.92855782101_RHP
    c(4) = 14.26016391358_RHP  /  71.96633302097_RHP
    c(5) =  9.2_RHP            /  10.0_RHP
    c(6) =  2.4_RHP            /  10.0_RHP
    c(7) =  3.0_RHP            /   5.0_RHP
    c(8) =  1.0_RHP

    b_im(1) =  -87.2700587467_RHP   / 913.3579230613_RHP
    b_im(2) =    0.0_RHP
    b_im(3) =    0.0_RHP
    b_im(4) =  223.482180632610_RHP /  95.558587375310_RHP
    b_im(5) =  -11.433695189920_RHP /  81.418160029310_RHP
    b_im(6) =  -39.379526789629_RHP /  19.018526304540_RHP
    b_im(7) =   32.727382324388_RHP /  42.900044865799_RHP
    b_im(8) =    4.1_RHP            /  20.0_RHP

    b_ex    = b_im

    ! A_im . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    a_im(2,1) =   4.1_RHP            /   20.0_RHP
    a_im(3,1) =   4.1_RHP            /   40.0_RHP
    a_im(4,1) =  68.37856364310_RHP  /  925.29203076860_RHP
    a_im(5,1) =  30.16520224154_RHP  /  100.81342136671_RHP
    a_im(6,1) =  21.88664790290_RHP  /  148.99783939110_RHP
    a_im(7,1) =  10.20004230633_RHP  /   57.15676835656_RHP
    a_im(8,1) = -87.27005874670_RHP  /  913.35792306130_RHP

    a_im(2,2) =   4.1_RHP            /   20.0_RHP
    a_im(3,2) = -56.7603406766_RHP   / 1193.1857230679_RHP

    a_im(3,3) =   4.1_RHP            /   20.0_RHP
    a_im(4,3) = -11.038504710300_RHP /  136.701519337300_RHP
    a_im(5,3) =  30.586259806659_RHP /   12.414158314087_RHP
    a_im(6,3) =  63.825689466800_RHP /  543.644631884100_RHP
    a_im(7,3) =  25.762820946817_RHP /   25.263940353407_RHP

    a_im(4,4) =   4.1_RHP            /   20.0_RHP
    a_im(5,4) = -22.760509404356_RHP /   11.113319521817_RHP
    a_im(6,4) = -11.797104745550_RHP /   53.211547248960_RHP
    a_im(7,4) = -21.613759091450_RHP /   97.559073359090_RHP
    a_im(8,4) = 223.482180632610_RHP /   95.558587375310_RHP

    a_im(5,5) =   4.1_RHP            /   20.0_RHP
    a_im(6,5) = -60.92811917200_RHP  / 8023.46106767100_RHP
    a_im(7,5) = -21.12173095930_RHP  /  584.68595025340_RHP
    a_im(8,5) = -11.43369518992_RHP  /   81.41816002931_RHP

    a_im(6,6) =   4.1_RHP            /   20.0_RHP
    a_im(7,6) = -42.699250595730_RHP /   78.270590407490_RHP
    a_im(8,6) = -39.379526789629_RHP /   19.018526304540_RHP

    a_im(7,7) =   4.1_RHP            /   20.0_RHP
    a_im(8,7) =  32.727382324388_RHP /   42.900044865799_RHP

    a_im(8,8) =   4.1_RHP            /   20.0_RHP

    ! A_ex . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

    a_ex(2,1) =    4.1_RHP            /  10.0_RHP
    a_ex(3,1) =   36.790274446400_RHP / 207.228047367700_RHP
    a_ex(4,1) =   12.680235234080_RHP / 103.408227345210_RHP
    a_ex(5,1) =  144.632819003510_RHP /  63.153537034770_RHP
    a_ex(6,1) =   14.090043504691_RHP /  34.967701212078_RHP
    a_ex(7,1) =   19.230459214898_RHP /  13.134317526959_RHP
    a_ex(8,1) =  -19.977161125411_RHP /  11.928030595625_RHP

    a_ex(3,2) =   67.762320755100_RHP / 822.414386656300_RHP

    a_ex(4,3) =   10.299339394170_RHP / 136.365588504790_RHP
    a_ex(5,3) =  661.144352112120_RHP /  58.794905890930_RHP
    a_ex(6,3) =   15.191511035443_RHP /  11.219624916014_RHP
    a_ex(7,3) =  212.753313583030_RHP /  29.424553649710_RHP
    a_ex(8,3) = -407.959767960540_RHP /  63.849078235390_RHP

    a_ex(5,4) = -540.531701528390_RHP /  42.847980215620_RHP
    a_ex(6,4) =  -18.461159152457_RHP /  12.425892160975_RHP
    a_ex(7,4) = -381.453459884190_RHP /  48.626203187230_RHP
    a_ex(8,4) =  177.454434618887_RHP /  12.078138498510_RHP

    a_ex(6,5) =  -28.166716381100_RHP / 901.161929587000_RHP
    a_ex(7,5) =   -1.0_RHP            /   8.0_RHP
    a_ex(8,5) =   78.267220542500_RHP / 826.770190026100_RHP

    a_ex(7,6) =   -1.0_RHP            /   8.0_RHP
    a_ex(8,6) = -695.630110598110_RHP /  96.465806942050_RHP

    a_ex(8,7) =   73.566282105260_RHP /  49.421867764050_RHP

  case default

    call Error( 'Init_IMEX_RK_Method',                      &
                'requested number of stages not supported', &
                'IMEX_Runge_Kutta_Method'                   )

  end select

  allocate( this % c    ( this%n_stage )               )
  allocate( this % b_im ( this%n_stage )               )
  allocate( this % b_ex ( this%n_stage )               )
  allocate( this % a_im ( this%n_stage, this%n_stage ) )
  allocate( this % a_ex ( this%n_stage, this%n_stage ) )

  this % c    = real(c   , RNP)
  this % b_im = real(b_im, RNP)
  this % b_ex = real(b_ex, RNP)
  this % a_im = real(a_im, RNP)
  this % a_ex = real(a_ex, RNP)

end subroutine Init_IMEX_RK_Method

!===============================================================================

subroutine Show_IMEX_RK_Method(this, unit)
  class(IMEX_RK_Method), intent(in) :: this
  integer,     optional, intent(in) :: unit  !< output unit

  character(len=*), parameter :: fmt_ca = '(2X,F13.10," |",99F14.10)'
  character(len=*), parameter :: fmt_b =  '(2X,13X,   " |",99F14.10)'
  integer :: i, io

  if (this % n_stage < 1) return

  if (present(unit)) then
    io = unit
  else
    io = OUTPUT_UNIT
  end if

  write(io,'(/,2X,A)') 'IMEX Runge-Kutta method'
  write(io,'(2X,A,T15,A)')  'name:'  , trim(this % name)
  write(io,'(2X,A,T15,I0)') 'stages:', this % n_stage
  write(io,'(2X,A,T15,I0)') 'order:' , this % order

  write(io,'(/,2X,A,/)') 'implicit part'
  do i = 1, this%n_stage
    write(io,fmt_ca) this % c(i), this % a_im(i,1:i)
  end do
  write(io,'(2X,A)') repeat('-', 16 + 14*this%n_stage)
  write(io,fmt_b) this % b_im

  write(io,'(/,2X,A,/)') 'explicit part'
  do i = 1, this%n_stage
    write(io,fmt_ca) this % c(i), this % a_ex(i,1:i-1)
  end do
  write(io,'(2X,A)') repeat('-', 16 + 14*this%n_stage)
  write(io,fmt_b) this % b_ex
  write(io,*)

end subroutine Show_IMEX_RK_Method

!===============================================================================

end module IMEX_Runge_Kutta_Method
