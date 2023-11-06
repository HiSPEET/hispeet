module Smiling_Face
  use Kind_Parameters, only: RNP
  use Constants, only: ZERO, ONE, PI
  implicit none
  private

  public :: SmilingFace

  type SmilingFace
    real(RNP) :: r_e     = 0.2   !< radial position of eyes
    real(RNP) :: phi_e   = PI/3  !< half angle between eyes
    real(RNP) :: sigma_e = 0.04  !< standard deviation of eye potential
    real(RNP) :: h_m     = 0.10  !< height of mouth origin
    real(RNP) :: r_m     = 0.40  !< radius of mouth arc
    real(RNP) :: phi_m   = PI/5  !< half opening angle of mouth
    real(RNP) :: sigma_m = 0.02  !< standard deviation of mouth
  contains
    procedure :: Density
  end type SmilingFace

contains

  elemental real(RNP) function Density(this, x, y) result(f)
    class(SmilingFace), intent(in) :: this
    real(RNP), intent(in) :: x !< x-coordinate
    real(RNP), intent(in) :: y !< y-coordinate

    real(RNP) :: a, b, d, r, phi

    ! eyes .....................................................................

    associate( r_e     => this % r_e     &
             , phi_e   => this % phi_e   &
             , sigma_e => this % sigma_e )

      a = -ONE / (2 * sigma_e**2)

      f = ( exp( a * ( (x - r_e * sin(- phi_e))**2    &
                     + (y - r_e * cos(- phi_e))**2) ) &
          + exp( a * ( (x - r_e * sin(+ phi_e))**2    &
                     + (y - r_e * cos(+ phi_e))**2) ) )

    end associate

    ! mouth ....................................................................

    associate( h_m     => this % h_m     &
             , r_m     => this % r_m     &
             , phi_m   => this % phi_m   &
             , sigma_m => this % sigma_m )

      phi = atan(x, h_m - y)

      if (abs(phi) < phi_m) then
        r = sqrt(x**2 + (h_m - y)**2)    ! distance from mouth origin
        d = r - r_m                      ! distance from mouth arc
      else if (phi < ZERO) then
        a =  r_m * sin(-phi_m)           ! x-coordinate of left mouth corner
        b = -r_m * cos(-phi_m) + h_m     ! y-coordinate of left mouth corner
        d = sqrt((x-a)**2 + (y-b)**2)    ! distance from left mouth corner
      else
        a =  r_m * sin(phi_m)            ! x-coordinate of right mouth corner
        b = -r_m * cos(phi_m) + h_m      ! y-coordinate of right mouth corner
        d = sqrt((x-a)**2 + (y-b)**2)    ! distance from right mouth corner
      end if

      a = -ONE / (2 * sigma_m**2)

      f = f + exp(a * d**2)

    end associate

  end function Density

end module Smiling_Face
