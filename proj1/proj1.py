# Author: Nathaniel Ruhl
# MEAM 5460 HW3

import matplotlib.pyplot as plt
import numpy as np

from gaussxw import gaussxwab
from Rotor2 import Rotor2

# Constants
rho0 = 0.002378 # slug/ft^3

#### These functions could be methods in the rotor2 class, but I don't want to make it too complicated so I'll keep them here
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
    INPUTS:
    TR = taper ratio (eg 2 for 2:1)
    c_75: chord at 75%R'''
    chord = (1/rotor2.TR)*(r-0.75)+rotor2.chord_75
    return chord

def Sigma(r, rotor2):
        '''Rotor solidity as a function of r changes because of taper'''
        sigma = rotor2.Nb*Chord(r, rotor2)/(np.pi*rotor2.R)
        return sigma
####

### BEMT functions ###

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

def dCP(r, rotor2, F=1):
    '''Local blade coefficient of power as a function of non-dimensional radius r, rotor2 object, F: Prandtl tip loss factor'''
    return 4*F*Lambda(r, rotor2, F)**3*r

def calc_CT_CP(rotor2):
    '''Calculate total coefficients CT and CP for a given rotor2 object'''
    N = 10
    rs, w = gaussxwab(N, 0, 1)
    CT = 0
    CP = 0
    for i in range(N):
        CT += w[i]*dCT(rs[i], rotor2)
        CP += w[i]*dCP(rs[i], rotor2)
    return CT, CP

def Ff(r, lmbda, rotor2):
    '''Prandtl tip loss factor as a function of non-dimensional radius and inflow.
    Functioin used in the iterative method calcF(r, lmbda, rotors)'''
    f = (rotor2.Nb/2)*(1-r)/lmbda
    F = (2/np.pi)*np.arccos(np.exp(-f))
    return F

def calcF(r, rotor2):
    F = 1
    error = 100
    num_iter = 0
    while(error > tol):
        Flast = F
        lmbda = Lambda(r, rotor2, Flast)
        F = Ff(r, lmbda, rotor2)
        error = np.abs(F-Flast)
        num_iter += 1
    print(f"F={F}")
    return F

# Testing

rotor2 = Rotor2()
rotor2.Nb=1
rotor2.theta_tw=0 # UNTWISTED BLADE!! (page 144)
tol = 1e-6
rs = np.linspace(0,0.99,100)
F_list = np.zeros(100)
lambda_list = np.zeros(100)
dCT_list = np.zeros(100)
for i in range(len(rs)):
    F = calcF(rs[i], rotor2)
    F_list[i] = F
    lambda_list[i] = Lambda(rs[i], rotor2, F)
    dCT_list[i] = dCT(rs[i], rotor2, F)

plt.figure(1)
plt.plot(rs, F_list)
plt.grid()
plt.xlim([0.7,1])
plt.xlabel('r')
plt.ylabel('F')

plt.figure(2)
plt.plot(rs, lambda_list, label="With Prandtl tip loss")
plt.plot(rs, Lambda(rs, rotor2, F=1), label="No tip loss")
plt.grid()
plt.xlabel('r')
plt.ylabel(r'$\lambda$')
plt.legend()

plt.figure(3)
plt.plot(rs, dCT_list, label="With Prandtl tip loss")
plt.plot(rs, dCT(rs, rotor2, F=1), label="No tip loss")
plt.grid()
plt.xlabel('r')
plt.ylabel(r'$c_T$')
plt.legend()
plt.show()

# rotor2 = Rotor2()

# theta75_list = np.deg2rad(np.array([3,6,9]))
# CT, CP = calc_CT_CP(rotor2)
# rs = np.linspace(0, 1, 100)  # Non-dimensional blade radius
# for theta75 in theta75_list:
#     rotor2.theta_75 = theta75
#     plt.figure(1)
#     plt.plot(rs, Lambda(rs, rotor2), label=fr"$\theta_{{75}}$={np.rad2deg(rotor2.theta_75):.1f}$^\circ$")

# plt.figure(1)
# plt.title("Inflow Distribution")
# plt.xlabel("Non-dimensional radial position, r")
# plt.ylabel(r"Local Inflow Ratio, $\lambda(r)$")
# plt.grid()
# plt.legend()
# plt.show()
