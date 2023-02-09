module mod_copter

use, intrinsic :: iso_fortran_env ! real64
implicit none
integer, parameter :: wp = real64

type copter
    real(wp) :: R  ! Rotor radius (m)
    integer :: Nb  ! Number of blades
    real(wp) :: chord  ! chord (m)
    real(wp) :: Omega  ! angular velocity (rad/sec)
    real(wp) :: P_extra ! multiplicative factor for extra required power
end type

end module mod_copter