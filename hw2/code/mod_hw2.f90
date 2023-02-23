module mod_hw2
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

end module mod_hw2