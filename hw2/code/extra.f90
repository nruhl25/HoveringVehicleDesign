program extra
    use mod_tools, only: pi, linspace
    use mod_hw2, only: calc_lambda, T, rho, vtip1, vtip2, R, A
    implicit none

    ! In this program, I want to practice analyzing and plotting a function of two variables

    real, parameter :: alpha = 0.0   ! incidence angle
    integer, parameter :: N_vinf = 10
    integer, parameter :: N_vtip = 10

    real, dimension(N_vtip) :: vtip
    real, dimension(N_vinf) :: v_inf
    real, dimension(N_vinf, N_vtip) :: lambda
    integer :: i, j

    call linspace(v_inf, 0.0, 84.3905, N_vinf)
    call linspace(vtip, vtip1, vtip2, N_vtip)

    ! do i=1,N_vtip
    !     do j=1,N_vinf
    !         lambda(i,j) = calc_lambda(alpha, v_inf(j), vtip(i))
    !         print*, v_inf(j), lambda(i,j)
    !     enddo
    !     print*,''
    !     print*,''
    ! enddo
    print *, vtip

end program extra