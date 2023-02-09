program main
    use, intrinsic :: iso_fortran_env ! real64
    use mod_tools, only: linspace, plot2
    use mod_splines, only: spline
    use mod_copter, only: copter
    implicit none

    integer, parameter :: wp = real64
    real(wp), parameter :: pi = 4.*atan(1.)
    real(wp), parameter :: Tt = 500.    ! Total thrust required for hover

    ! data from problem assignment
    type(copter) :: trad   ! traditional main rotor copter
    type(copter) :: quad   ! quadcopter
    type(copter) :: tandem ! tandem rotor configuration
    type(copter) :: hex    ! hexacopter
    real(wp), dimension(13) :: CT_NORM, FM

    ! Variable declarations
    character(len=:), allocatable :: fn
    real(wp), dimension(100) :: FM_long, CT_NORM_long
    integer :: i
    type(spline) :: fm_sp
    real(wp) :: rho_ssl, rho_1    ! Density at SSL and 2000 ft

    ! Assign variables
    CT_NORM = [(0.0+i*0.1, i=0,12)]
    FM = (/0., 0.2, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59, 0.6, 0.6, 0.6, 0.6, 0.6/)

    trad = copter(2.1336, 6, 0.100584, 100.007, 1.20)
    quad = copter(0.61, 3, 1./6., 350., 1.15)
    tandem = copter(1.3716, 5, 0.0762, 155., 1.15)
    hex = copter(0.48768, 4, 0.0381, 350., 1.15)

    fm_sp = spline(CT_NORM, FM)

    !!!!!!!!!!! ANALYSIS BELOW !!!!!!!!!!!!!!!

    call analyze_copter(trad, rho_ssl)

    contains

    subroutine analyze_copter(heli, rho)
    type(copter), intent(in) :: heli   ! copter type
    real(wp) :: rho    ! density

    real(wp) :: A   ! Area of rotor disc

    A=pi*heli%R**2
    print *, A

    end subroutine analyze_copter

end program main

