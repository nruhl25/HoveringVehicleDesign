# Varioius tools for proj1
import matplotlib.pyplot as plt
import numpy as np
from BEMT2 import Cl_slope, Theta, Chord, Sigma

# functions related to linear blade property variations
from BEMT2 import Cl_slope, Theta, Chord, Sigma
from BEMT2 import Lambda, Alpha, Phi, dCT, dCPi, calc_CT_CPi, calcF  # BEMT functions

def plot_rotor2(rotor2):
    '''Function to plot rotor2 spec'''
    # Non-dimensional radius, shortened at ends to avoid discontinuities
    rs = np.linspace(0.1, 0.99, 100)

    plt.rcParams['figure.figsize'] = [14, 5]

    plt.subplot(1, 4, 1)
    plt.title("Blade Twist")
    plt.ylabel(r"$\theta(r)$")
    plt.xlabel("r")
    plt.plot(rs, np.rad2deg(Theta(rs, rotor2)))
    plt.grid()

    plt.subplot(1, 4, 2)
    plt.title("Lift-Curve Slope")
    plt.ylabel(r"$c_{l,\alpha}(r)$")
    plt.xlabel("r")
    plt.plot(rs, Cl_slope(rs, rotor2))
    plt.grid()

    plt.subplot(1, 4, 3)
    plt.title("Chord")
    plt.ylabel("c(r)")
    plt.xlabel("r")
    plt.plot(rs, Chord(rs, rotor2))
    plt.grid()

    plt.subplot(1, 4, 4)
    plt.title("Local Rotor Solidity")
    plt.ylabel(r"$\sigma(r)$")
    plt.xlabel("r")
    plt.plot(rs, Sigma(rs, rotor2))
    plt.grid()

    plt.tight_layout()
    plt.show()
    print(f"Other rotor2 properties: Nb={rotor2.Nb}, R={rotor2.R} ft, vtip = {rotor2.vtip} ft/sec, cd0={rotor2.cd0}")
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    return

def plot_prandtl_comparison(rotor2):
    rs = np.linspace(0.1, 0.99, 100)

    F_list = np.zeros(100)
    lambda_list = np.zeros(100)
    dCT_list = np.zeros(100)
    dCPi_list = np.zeros(100)
    for i in range(len(rs)):
        F = calcF(rs[i], rotor2)
        F_list[i] = F
        lambda_list[i] = Lambda(rs[i], rotor2, F)
        dCT_list[i] = dCT(rs[i], rotor2, F)
        dCPi_list[i] = dCPi(rs[i], rotor2, F)

    plt.figure()
    plt.title("Prandtl Tip Loss Function")
    plt.plot(rs, F_list)
    plt.grid()
    plt.xlim([0.7, 1])
    plt.xlabel('r')
    plt.ylabel('F')

    plt.figure()
    plt.title("Inflow Ratio Distribution")
    plt.plot(rs, lambda_list, label="With tip losses")
    plt.plot(rs, Lambda(rs, rotor2, F=1), label="Without tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$\lambda$')
    plt.legend()

    plt.figure()
    plt.title("Local Thrust Coefficient Distribution")
    plt.plot(rs, dCT_list, label="With Prandtl tip loss")
    plt.plot(rs, dCT(rs, rotor2, F=1), label="No tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$c_T$')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.title("Local Induced Power Coefficient Distribution")
    plt.plot(rs, dCPi_list, label="With Prandtl tip loss")
    plt.plot(rs, dCPi(rs, rotor2, F=1), label="No tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$c_{P,i}$')
    plt.legend()
    plt.show()
    return
