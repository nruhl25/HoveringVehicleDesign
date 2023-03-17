# Author: Nathaniel Ruhl
# MEAM 5460 HW3

import numpy as np

from gaussxw import gaussxwab

# Constants
rho0 = 0.002378 # slug/ft^3
Nb = 1

# Test rotor for hw3
c = 2.0 # Chord length, ft
vtip = 650 # ft/sec
R = 30 # ft
sigma = Nb*c/(np.pi*R)  # solidity
cl_slope = 2*np.pi

# Inflow ratio as a function of non-dimensional radius, r, and input collective at 75%R, theta_75
def Lambda(r, theta_75):
    theta = Theta(r, theta_75)
    lmbda = (sigma*cl_slope/16)*(np.sqrt(1+(32/(sigma*cl_slope))*theta*r)-1)
    return lmbda

def Theta(r, theta_75):
    theta = 0.75*theta_75/r
    return theta

# angle of attack as a function of non-dimensional radius
def Alpha(r, theta_75):
    return Theta(r, theta_75) - Phi(r, theta_75)

# Inflow angle as a function of non-dimensional radius
def Phi(r, theta_75):
    return Lambda(r, theta_75)/r

def dCT(r, theta_75):
    return 4*Lambda(r, theta_75)**2*r

def dCP(r, theta_75):
    return 4*Lambda(r, theta_75)**3*r

def calc_CT_CP(theta_75):
    # Calculate total coefficients CT and CP
    N = 10
    rs, w = gaussxwab(N, 0.0001, 1)
    CT = 0
    CP = 0
    for i in range(N):
        CT += w[i]*dCT(rs[i], theta_75)
        CP += w[i]*dCP(rs[i], theta_75)
    return CT, CP
