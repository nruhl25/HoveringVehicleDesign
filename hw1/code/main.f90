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

    trad = copter(2.1336, 6, 0.100584, 100.007, 1.15, 1)
    quad = copter(0.61, 3, 1./6., 350., 1.0, 4)
    tandem = copter(1.3716, 5, 0.0762, 155., 1.0, 2)
    hex = copter(0.48768, 4, 0.0381, 350., 1.0, 6)
    copter_list = (/trad, quad, tandem, hex/)
    copter_names = (/"Trad", "Quad", "Tand", "Hexa"/)

    fm_sp = spline(CT_NORM, FM)

    !call linspace(ct_norm_long, minval(CT_NORM), maxval(CT_NORM), 200)
    !fm_long = fm_sp%value(ct_norm_long)

    !!!!!!!!!!! ANALYSIS BELOW !!!!!!!!!!!!!!!

    ! print*, "# Results at Mean Sea Level (2000 ft)"
    ! print*, "# Ct, FM, Cp, DL [N/m^2], PL [N/W], P_ideal (rotor) [kW], P_single (rotor) [kW], P_total (system) [kW]"
    ! print*, "# ----------------------------------------------------------------------------------------------------"
    print*, "# Ct/sigma, FM, DL, PL"
    print *, "# -------------------"
    do i = 1,4
        call analyze_copter(copter_list(i), rho_ssl, fm_sp)
        ! print*, copter_list(i)%Ct, copter_list(i)%FM, copter_list(i)%Cp, copter_list(i)%DL, & 
        ! copter_list(i)%PL, copter_list(i)%P_ideal, copter_list(i)%P_single, copter_list(i)%P_tot
        print*, copter_list(i)%Ct_norm, copter_list(i)%FM, copter_list(i)%DL, copter_list(i)%PL
    enddo

    contains

    subroutine analyze_copter(heli, rho, fm_sp)
        type(copter), intent(inout) :: heli   ! copter type object
        real(wp), intent(in) :: rho    ! density
        type(spline), intent(in) :: fm_sp

        real(wp) :: Ar, Tr, Ct, sigma, Ct_norm, FM, Cp, PL, DL, P_ideal, P_single, P_tot, q_inf_v, Cp_ideal

        ! Calculate Tr, Ct, Ct_norm, FM for single rotor
        Ar = pi*heli%R**2
        Tr = Tt/heli%Nr
        Ct = Tr/(Ar*rho*(heli%Omega*heli%R)**2.)
        sigma = heli%Nb*heli%chord*heli%R/Ar
        Ct_norm = Ct/sigma
        FM = fm_sp%value(Ct_norm)

        ! Calculate PL, DL for single rotor (add 5% for inneficiency of power conversion)
        q_inf_v = rho*Ar*(heli%Omega*heli%R)**3.
        Cp_ideal = Ct**(3./2.)/sqrt(2.)
        Cp = Cp_ideal/FM
        PL = Ct/(heli%Omega*heli%R*Cp)  ! N/W
        DL = Tr/Ar  ! N/m^2
        P_ideal = Cp_ideal*q_inf_v/1000.
        P_single = Cp*q_inf_v/1000.   ! kW
        P_tot = (1./0.95)*heli%P_extra*(heli%Nr*P_single)   ! kW

        ! Add attributes to the copter object (add to zero)
        heli%sigma = heli%sigma + sigma
        heli%Ct = heli%Ct + Ct
        heli%Ct_norm = heli%Ct_norm + Ct_norm
        heli%FM = heli%FM + FM
        heli%Cp = heli%Cp + Cp
        heli%PL = heli%PL + PL
        heli%DL = heli%DL + DL
        heli%P_tot = heli%P_tot + P_tot
        heli%P_ideal = heli%P_ideal + P_ideal
        heli%P_single = heli%P_single + P_single

    end subroutine analyze_copter

end program main

