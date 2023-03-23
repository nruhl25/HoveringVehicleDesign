# Author: Nathaniel Ruhl
# MEAM 5460 HW3

import numpy as np

# Constants
rho0 = 0.002378 # slug/ft^3

def Theta(r, theta_75):
    '''Input collective as a function of non-dimensional radius and theta_75 input'''
    theta = 0.75*theta_75/r
    return theta

def Lambda(r, theta_75, rotor):
    '''Inflow ratio as a function of non-dimensional radius r, input collective at 75%R theta_75, and rotor object'''
    theta = Theta(r, theta_75)
    lmbda = (rotor.sigma*rotor.cl_slope/16)*(np.sqrt(1+(32/(rotor.sigma*rotor.cl_slope))*theta*r)-1)
    return lmbda

def Alpha(r, theta_75, rotor):
    '''Angle of attack as a function of non-dimensional radius r, input collective at 75%R theta_75, and rotor object'''
    return Theta(r, theta_75) - Phi(r, theta_75, rotor)

def Phi(r, theta_75, rotor):
    '''Inflow angle as a function of non-dimensional radius r, input collective at 75%R theta_75, and rotor object'''
    return Lambda(r, theta_75, rotor)/r

def dCT(r, theta_75, rotor):
    '''Local blade coefficient of thrust as a function of non-dimensional radius r, input collective at 75%R theta_75, and rotor object'''
    return 4*Lambda(r, theta_75, rotor)**2*r

def dCP(r, theta_75, rotor):
    '''Local blade coefficient of power as a function of non-dimensional radius r, input collective at 75%R theta_75, and rotor object'''
    return 4*Lambda(r, theta_75, rotor)**3*r

def calc_CT_CP(theta_75, rotor):
    '''Calculate total coefficients CT and CP for a given rotor object. This function makes used of the fact the Lambda is not a function of r, and uses r=1 to define the inflow (I originally made this function like this when it involved numerical integration)'''
    CT = 2*Lambda(1, theta_75, rotor)**2
    CP = 2*Lambda(1, theta_75, rotor)**3
    return CT, CP

def CP0(rotor):
    '''Profile power coefficient'''
    CP0 = (1/8)*rotor.sigma*rotor.cd0
    return CP0
