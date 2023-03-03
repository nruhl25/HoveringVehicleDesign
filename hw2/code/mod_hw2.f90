module mod_hw2
use mod_tools, only: pi
implicit none
public :: f, fprime, calc_lambda

real, parameter :: T = 100.04      ! [lbf]
real, parameter :: rho = 0.002378  ! [slug/ft^3]
real, parameter :: vtip1 = 400.0   ! [ft/s]
real, parameter :: vtip2 = 525.0   ! [ft/s]
real, parameter :: R = 2.5         ! [ft] blade radius
real, parameter :: A = pi*R**2.    ! [ft^2] rotor disc area

contains
real function f(lambda, C_T, mu, alpha)
    real, intent(in) :: lambda, C_T, mu, alpha
    f = lambda - mu*tan(alpha) - C_T/(2.*sqrt(mu**2+lambda**2))
end function f

real function fprime(lambda, C_T, mu)
    real, intent(in) :: lambda, C_T, mu
    fprime = 1 + (C_T/2.)*(mu**2+lambda**2)**(-3./2.)*lambda
end function fprime

! This subroutine calculates the induced velocity coefficients
real function calc_lambda(alpha, v_inf, vtip)
    real, intent(in) :: alpha, v_inf, vtip
    real :: C_T, v_h, P_h, mu, lambda, lambda_h
    real :: epsilon    ! error in the Newton-Raphson method
    real :: tol = 0.0005    ! Tolerance for Newton-Raphson
    real :: f_n, fprime_n, lambda_next
    integer :: num_iter

    mu = v_inf*cos(alpha)/vtip

    ! Note for my future self: the values below could be accesed from the main program scope,
    ! but are re-defined here because that's not a great idea
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
        fprime_n = fprime(lambda, C_T, mu)
        lambda_next = lambda - (f_n/fprime_n)
        num_iter = num_iter + 1
        epsilon = abs((lambda_next-lambda)/lambda_next)
    end do

    calc_lambda = lambda_next

end function calc_lambda

end module mod_hw2