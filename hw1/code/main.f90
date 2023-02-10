program main
    use, intrinsic :: iso_fortran_env ! real64
    use mod_tools, only: linspace, plot2
    use mod_splines, only: spline
    use mod_copter, only: copter
    implicit none

    integer, parameter :: wp = real64
    real(wp), parameter :: pi = 4.*atan(1.)
    real(wp), parameter :: Tt = 544.*9.81    ! Total thrust required for hover [N]

    ! data from problem assignment
    type(copter) :: trad   ! traditional main rotor copter
    type(copter) :: quad   ! quadcopter
    type(copter) :: tandem ! tandem rotor configuration
    type(copter) :: hex    ! hexacopter
    real(wp), dimension(13) :: CT_NORM, FM   ! data table in HW assignment
    ! real(wp), dimension(200) :: ct_norm_long, fm_long

    ! Variable declarations
    ! character(len=:), allocatable :: fn   ! target txt filename
    integer :: i
    type(spline) :: fm_sp         ! FM vs CT/sigma spline
    real(wp) :: rho_ssl, rho_1    ! Density at SSL and 2000 ft
    type(copter), dimension(4) :: copter_list
    type(character(len=4)), dimension(4) :: copter_names

    ! Assign variables
    CT_NORM = [(0.0+i*0.01, i=0,12)]
    FM = (/0., 0.2, 0.4, 0.5, 0.55, 0.57, 0.58, 0.59, 0.6, 0.6, 0.6, 0.6, 0.6/)
    rho_ssl = 1.22557 ! kg/m^3
    rho_1 = 1.15490 ! kg/m^3

    trad = copter(2.1336, 6, 0.100584, 100.007, 1.20, 1)
    quad = copter(0.61, 3, 1./6., 350., 1.15, 4)
    tandem = copter(1.3716, 5, 0.0762, 155., 1.15, 2)
    hex = copter(0.48768, 4, 0.0381, 350., 1.15, 6)
    copter_list = (/trad, quad, tandem, hex/)
    copter_names = (/"Trad", "Quad", "Tand", "Hexa"/)

    fm_sp = spline(CT_NORM, FM)

    !call linspace(ct_norm_long, minval(CT_NORM), maxval(CT_NORM), 200)
    !fm_long = fm_sp%value(ct_norm_long)

    !!!!!!!!!!! ANALYSIS BELOW !!!!!!!!!!!!!!!

    print*, "Ct/sigma   FM (at SSL)"
    do i = 1,4
        print*, copter_names(i), " results:"
        call analyze_copter(copter_list(i), rho_ssl, fm_sp)
    enddo

    contains

    subroutine analyze_copter(heli, rho, fm_sp)
        type(copter), intent(inout) :: heli   ! copter type object
        real(wp), intent(in) :: rho    ! density
        type(spline), intent(in) :: fm_sp

        real(wp) :: Ar, Tr, Ct, sigma, Ct_norm, FM, Cp, PL, DL

        ! Calculate Tr, Ct, Ct_norm, FM for single rotor
        Ar = pi*heli%R**2
        Tr = Tt/heli%Nr
        Ct = Tr/(0.5*Ar*rho*(heli%Omega*heli%R)**2)
        sigma = heli%Nb*heli%chord*heli%R/Ar
        Ct_norm = Ct/sigma
        FM = fm_sp%value(Ct_norm)

        ! Calculate PL, DL for single rotor
        Cp = Ct**(3./2.)/(sqrt(2.)*FM)
        PL = Ct/(heli%Omega*heli%R*Cp)   ! ADD IN EXTRA POWER HERE ??
        DL = Tr/Ar


        ! Add attributes to the copter object (add to zero)
        heli%sigma = heli%sigma + sigma
        heli%Ct = heli%Ct + Ct
        heli%Ct_norm = heli%Ct_norm + Ct_norm
        heli%FM = heli%FM + FM
        heli%Cp = heli%Cp + Cp
        heli%PL = heli%PL + PL
        heli%DL = heli%DL + DL

        print*, heli%Ct_norm, heli%FM, heli%PL, heli%DL

    end subroutine analyze_copter

end program main

