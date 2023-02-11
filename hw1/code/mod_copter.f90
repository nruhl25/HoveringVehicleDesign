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
    integer :: Nr      ! number of lifting rotors

    real(wp) :: sigma = 0.0  ! rotor solidity, yet to be defined
    real(wp) :: Ct = 0.0    ! Coefficient of thrust, yet to be defined
    real(wp) :: Ct_norm = 0.0    ! Coefficient of thrust normalized by solidity, yet to be defined
    real(wp) :: P_ideal = 0.0    ! Ideal power for the single lifting rotor
    real(wp) :: FM = 0.0    ! Figure of merit, yet to be defined
    real(wp) :: Cp = 0.0       ! Coefficient of power, yet to be defined
    real(wp) :: PL = 0.0    ! Power Loading, yet to be defined
    real(wp) :: DL = 0.0    ! Disk loading, yet to be defined
    real(wp) :: P_single = 0.0   ! Real power of each lifting rotor
    real(wp) :: P_tot = 0.0   ! Total power needed for hover
end type

end module mod_copter