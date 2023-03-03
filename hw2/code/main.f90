program main
    use mod_tools, only: pi, linspace
    use mod_hw2, only: calc_lambda, T, rho, vtip1, vtip2, R, A
    implicit none

    ! Declare variables
    integer, parameter :: N = 50       ! Number of of vtip and v_inf to sample
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

    do i=1,N
        lambda(i) = calc_lambda(alpha, v_inf(i), vtip1)
        lambda_i(i) = lambda(i) - mu(i)*tan(alpha)
        P(i) = P_h*((mu(i)/lambda_h)*tan(alpha) + (lambda_h/sqrt(mu(i)**2+lambda(i)**2)))
        v_i(i) = vtip1*lambda_i(i)
        v_net(i) = lambda(i)*vtip1
        print*, v_inf(i), lambda(i)
    end do

end program main