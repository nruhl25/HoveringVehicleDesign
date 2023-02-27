program main
    use, intrinsic :: iso_fortran_env ! real64
    use mod_tools, only: pi, linspace
    use mod_hw2, only: f, fprime
    implicit none

    real, parameter :: T = 100.04      ! [lbf]
    real, parameter :: rho = 0.002378  ! [slug/ft^3]
    real, parameter :: vtip1 = 525.0 ! 400.0   ! [ft/s]
    real, parameter :: vtip2 = 525.0   ! [ft/s]
    real, parameter :: R = 2.5         ! [ft] blade radius
    real, parameter :: A = pi*R**2.    ! [ft^2] rotor disc area
    integer, parameter :: N = 50       ! Number of of vtip to sample
    real, parameter :: tol = 0.0005    ! Tolerance for Newton-Raphson

    ! Declare variables
    real :: alpha      ! incidence angle (rad)
    real :: C_T, v_h, P_h, lambda_h
    real, dimension(N) :: v_inf, mu, lambda, lambda_i, P, v_i, v_net
    integer :: i

    ! Define variables
    alpha = 10.*pi/180.

    ! Arrays
    call linspace(v_inf, 0.0, 84.3905, N)
    mu = v_inf*cos(alpha)/vtip1

    ! Scalars dependenant on vtip
    C_T = T/(rho*A*(vtip1)**2.)
    lambda_h = sqrt(C_T/2.)
    v_h = lambda_h*vtip1
    P_h = T*v_h

    ! Fill in lambda, lambda_i, and P arrays
    !print*, 'alpha=-10 deg'
    !print*, '# v_inf, v_i, P'
    do i=1,N
        lambda(i) = calc_lambda(alpha, v_inf(i), vtip1)
        lambda_i(i) = lambda(i) - mu(i)*tan(alpha)
        P(i) = P_h*((mu(i)/lambda_h)*tan(alpha) + (lambda_h/sqrt(mu(i)**2+lambda(i)**2)))
        v_i(i) = vtip1*lambda_i(i)
        !print*, v_inf(i), vtip1*lambda_i(i), P(i)
    end do

    v_net = v_inf*cos((pi/2)+alpha) - v_i

    do i=1,N
        print*, v_inf(i), v_net(i)
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

        ! Note for my future self: the values below could be accesed from the main program scope,
        ! but are re-defined here
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