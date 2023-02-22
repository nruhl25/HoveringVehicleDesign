module mod_tools
    use, intrinsic :: iso_fortran_env ! real64
    implicit none

    integer, parameter :: wp = real64
    real(wp), parameter :: pi = 4.*atan(1.)

    public :: linspace, plot2, plotN

    contains
    subroutine linspace(x, x_start, x_end, x_len)
        integer, intent(in) :: x_len
        real, intent(in) :: x_start, x_end
        real, dimension(1:x_len), intent(out) :: x
        real :: dx
        integer :: i

        dx = (x_end - x_start) / (x_len-1)
        x(1:x_len) = [(x_start + ((i-1)*dx), i=1,x_len)]
    end subroutine linspace

    subroutine plot2(x,y,filename)
        real(wp), intent(in), dimension(:) :: x, y
        character(len=:), allocatable, intent(in) :: filename
        integer :: size_x, size_y, i
        size_x = size(x)
        size_y = size(y)
        if (size_x /= size_y) then
            print *, "Array size mis-match"
        else
            open(unit=1, file=filename)
            do i = 1, size_x
                write(1,*) x(i), ' ', y(i)
            end do
        end if
    end subroutine plot2

    subroutine plotN(x, YY)
        real, dimension(:), intent(in) :: x
        real, dimension(:,:), intent(in) :: YY
        ! character(len=:), allocatable, intent(in) :: filename

        ! YY = [y1, y2,...] where yn are column vectors

        integer :: N            ! length of x and y arrays to plot
        integer :: Nyy          ! Number of y arrays to plot
        integer :: i, j
        real, dimension(:), allocatable :: yvals

        N = size(x)
        Nyy = size(YY, 2)
        allocate(yvals(Nyy))

        do i = 1,N
            yvals = [(YY(i,j), j=1,Nyy)]
            print*, x(i), yvals
        end do
        deallocate(yvals)
end subroutine plotN

end module mod_tools