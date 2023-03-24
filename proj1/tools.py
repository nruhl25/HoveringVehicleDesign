# Varioius tools for proj1 that interfrac with BEMT2
import matplotlib.pyplot as plt
import numpy as np

from gaussxw import gaussxwab
# functions related to linear blade property variations
from BEMT2 import Cl_slope, Theta, Chord, Sigma
# BEMT functions
from BEMT2 import Lambda, Alpha, dCT, dCPi, dCP0_NACA0012, dCP0_GENERAL, calcF

def calc_rotor_profiles(rotor2, airfoil, use_F=False, rs=np.linspace(0.1, 0.99, 100)):
    '''This function is an interface to calculating performance profiless (default: no tip losses)'''
    # Choose method of calculating profile drag
    if airfoil == "GENERAL":
        dCP0_func = lambda r, rotor, Ff: dCP0_GENERAL(r, rotor, Ff)
    elif airfoil == "NACA0012":
        dCP0_func = lambda r, rotor, Ff: dCP0_NACA0012(r, rotor, Ff)
    else:
        dCP0_func = lambda x: x 
        raise RuntimeError(
            "User entered arguments are not valid. Must be either 'GENERAL' or 'NACA0012'")
    
    if use_F is False:
        F_list = np.ones(len(rs))
        lambda_list = Lambda(rs, rotor2, F=1)
        dCT_list = dCT(rs, rotor2, F=1)
        dCPi_list = dCPi(rs, rotor2, F=1)
        dCP0_list = dCP0_func(rs, rotor2, F=1)
    elif use_F is True:
        # Using Prandtl tip loss, cannot make use of vectorization
        F_list = np.zeros(len(rs))
        lambda_list = np.zeros(len(rs))
        dCT_list = np.zeros(len(rs))
        dCPi_list = np.zeros(len(rs))
        dCP0_list = np.zeros(len(rs))
        for i in range(len(rs)):
            F = calcF(rs[i], rotor2)
            F_list[i] = F
            lambda_list[i] = Lambda(rs[i], rotor2, F)
            dCT_list[i] = dCT(rs[i], rotor2, F)
            dCPi_list[i] = dCPi(rs[i], rotor2, F)
            dCP0_list[i] = dCP0_func(rs[i], rotor2, F)
    
    return lambda_list, dCT_list, dCPi_list, dCP0_list, F_list

def calc_CT_CPi_CP0(rotor2, airfoil="GENERAL", use_F=False):
    '''Calculate total coefficients CT, induced CP, CP0 profile for a given rotor2 object (default: no tip losses)'''
    if airfoil == "GENERAL":
        def dCP0_func(r, rotor, Ff): return dCP0_GENERAL(r, rotor, Ff)
    elif airfoil == "NACA0012":
        def dCP0_func(r, rotor, Ff): return dCP0_NACA0012(r, rotor, Ff)
    else:
        raise RuntimeError(
            "User entered arguments are not valid. Must be either None (default) or 'NACA0012'")

    N = 10
    rs, w = gaussxwab(N, 0, 1)
    CT = 0
    CPi = 0
    CP0 = 0
    for i in range(N):
        F = calcF(rs[i], rotor2)
        CT += w[i]*dCT(rs[i], rotor2, F)
        CPi += w[i]*dCPi(rs[i], rotor2, F)
        CP0 += dCP0_func(rs[i], rotor2, F)
    return CT, CPi, CP0
