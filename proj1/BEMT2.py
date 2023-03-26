# Author: Nathaniel Ruhl
# MEAM 5460 Project 1

import matplotlib.pyplot as plt
import numpy as np

# Constants
rho0 = 0.002378 # slug/ft^3

# The functions below describe the linear variations of the rotor2 properties.
# They could be methods in the Rotor2 class, but I will keep them here to keep the Rotor2 class simple
#######################
def Cl_slope(r, rotor2):
    '''cl_slope as a function of r, rotor2 object'''
    cl_slope = rotor2.cl_slope_75 + rotor2.AV*(r-0.75)
    return cl_slope

def Theta(r, rotor2):
    '''Linear twist as a function of non-dimensional radius, rotor2 object'''
    theta = rotor2.theta_75 + rotor2.theta_tw*(r-0.75)
    return theta

def Chord(r, rotor2):
    '''Linear chord variation as a function of r
    TR = taper ratio (eg 0.5 for 2:1)
    chord_75: chord at 75%R'''
    chord_base = rotor2.chord_75/((0.75/rotor2.TR)+0.25)
    chord = chord_base*((1/rotor2.TR)-1)*r + chord_base
    return chord

def Sigma(r, rotor2):
        '''Rotor solidity as a function of r changes because of taper'''
        sigma = rotor2.Nb*Chord(r, rotor2)/(np.pi*rotor2.R)
        return sigma

def dCd(alpha, rotor2):
    '''Blade section coefficient of drag as a function of AoA (NACA0012 quadratic fit)'''
    return rotor2.d[0] + rotor2.d[1]*alpha + rotor2.d[2]*alpha**2
####################

###### BEMT functions ######

def Lambda(r, rotor2, F=1):
    '''Inflow ratio as a function of non-dimensional radius, rotor2 object, Prandtl tip loss factor'''
    theta = Theta(r, rotor2)
    sigma = Sigma(r, rotor2)
    cl_slope = Cl_slope(r, rotor2)
    lmbda = (sigma*cl_slope/(16*F))*(np.sqrt(1+((32*F)/(sigma*cl_slope))*theta*r)-1)
    return lmbda

def Alpha(r, rotor2, F=1):
    '''Angle of attack as a function of non-dimensional radius r, rotor2 object'''
    return Theta(r, rotor2) - Phi(r, rotor2, F)

def Phi(r, rotor2, F=1):
    '''Inflow angle as a function of non-dimensional radius r, input collective at 75%R theta_75'''
    return Lambda(r, rotor2, F)/r

def dCT(r, rotor2, F=1):
    '''Local blade coefficient of thrust as a function of non-dimensional radius r, rotor2 object, F: Prandtl tip loss factor'''
    return 4*F*Lambda(r, rotor2, F)**2*r

def dCPi(r, rotor2, F=1):
    '''Local blade coefficient of induced power as a function of non-dimensional radius r, rotor2 object, F: Prandtl tip loss factor (default: no tip losses)'''
    return 4*F*Lambda(r, rotor2, F)**3*r

def dCP0_GENERAL(r, rotor2, F=1):
    '''Differential expression for profile drag coefficient integral for a general airfoil (default: no tip losses)'''
    return (1/2)*rotor2.cd0*Sigma(r, rotor2)*r**3

def dCP0_NACA0012(r, rotor2, F=1):
    '''Differential expression for profile drag coefficient integral when NACA0012 is used (default: no tip losses)'''
    return (1/2)*dCd(Alpha(r, rotor2, F), rotor2)*Sigma(r, rotor2)*r**3

def Ff(r, lmbda, rotor2):
    '''Prandtl tip loss factor as a function of non-dimensional radius and inflow.
    Functioin used in the iterative method calcF(r, lmbda, rotors)'''
    f = (rotor2.Nb/2)*(1-r)/lmbda
    F = (2/np.pi)*np.arccos(np.exp(-f))
    return F

def calcF(r, rotor2):
    F = 1
    error = 100
    tol = 1e-6
    num_iter = 0
    while(error > tol):
        Flast = F
        lmbda = Lambda(r, rotor2, Flast)
        F = Ff(r, lmbda, rotor2)
        error = np.abs(F-Flast)
        num_iter += 1
    return F
