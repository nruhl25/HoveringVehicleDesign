module subs
implicit none

contains
real function f(lambda, C_T, mu, alpha)
    real, intent(in) :: lambda, C_T, mu, alpha
    f = lambda - mu*tan(alpha) - C_T/(2.*sqrt(mu**2+lambda**2))
end function f

real function fprime(lambda, C_T, mu, alpha)
    real, intent(in) :: lambda, C_T, mu, alpha
    fprime = 1 + (C_T/2.)*(mu**2+lambda**2)**(-3./2.)*lambda
end function fprime

end module subs

program main
    use, intrinsic :: iso_fortran_env ! real64
    use mod_tools, only: wp, pi, linspace, plot2, plotN
    use subs, only: f, fprime
    implicit none

    real, parameter :: T = 445.         ! [N]
    real, parameter :: rho = 1.22557    ! [kg/m^3]
    real, parameter :: vtip1 = 121.92   ! [m/s]
    real, parameter :: vtip2 = 160.02   ! [m/s]
    real, parameter :: R = 1.524/2.0    ! [m] blade radius
    real, parameter :: A = pi*R**2.      ! [m^2] rotor disc area
    integer, parameter :: N = 10        ! Number of of vtip to sample
    real, parameter :: tol = 0.0005     ! Tolerance for Newton-Raphson

    ! Declare variables
    real :: alpha      ! incidence angle (rad)
    real :: v_inf      ! [m/s] tunnel air speed
    real, dimension(N) :: vtip_list, mu, C_T, v_h
    real, dimension(N) :: lambda, lambda_i, lambda_h, P, P_h      ! OUTPUTS OF PROGRAM
    integer :: i

    ! Define variables
    alpha = 10.*pi/180.
    v_inf = 25.722222   ! 50 knots

    ! Depend on alpha, v_inf
    call linspace(vtip_list, vtip1, vtip2, N)
    mu = v_inf*cos(alpha)/vtip_list
    C_T = T/(rho*A*(vtip_list)**2.)
    lambda_h = sqrt(C_T/2.)
    v_h = lambda_h*vtip_list
    P_h = T*v_h

    ! Fill in lambda, lambda_i, P, P_h arrays
    do i=1,N
        lambda(i) = calc_lambda(alpha, v_inf, vtip_list(i))

        lambda_i(i) = lambda(i) - mu(i)*tan(alpha)
        P(i) = P_h(i)*((mu(i)/lambda_h(i))*tan(alpha) + (lambda_h(i)/sqrt(mu(i)**2+lambda(i)**2)))
        print*, mu(i)/lambda_h(i), lambda(i)/lambda_h(i), P(i)/P_h(i)
    end do

    contains

    ! This subroutine calculates the induced velocity coefficients
    real function calc_lambda(alpha, v_inf, vtip)

        real, intent(in) :: alpha, v_inf, vtip

        real :: C_T, v_h, mu, lambda, lambda_h
        real :: epsilon    ! error in the Newton-Raphson method
        real :: f_n, fprime_n, lambda_next
        integer :: num_iter

        mu = v_inf*cos(alpha)/vtip
        C_T = T/(rho*A*vtip**2.)
        lambda_h = sqrt(C_T/2.)
        v_h = lambda_h*vtip
        P_h = T*v_h

        ! Solve inflow ratio, lambda
        lambda_next = lambda_h       ! inital guess for Newton-Raphson
        epsilon = 1.
        num_iter = 0
        do while(epsilon>tol)
            lambda = lambda_next
            f_n = f(lambda, C_T, mu, alpha)
            fprime_n = fprime(lambda, C_T, mu, alpha)
            lambda_next = lambda - (f_n/fprime_n)
            num_iter = num_iter + 1
            epsilon = abs((lambda_next-lambda)/lambda_next)
        end do

        calc_lambda = lambda_next

    end function calc_lambda

end program main