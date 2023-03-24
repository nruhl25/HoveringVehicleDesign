

# Various plotting tools for proj1 report
from IPython.display import HTML
import matplotlib.pyplot as plt
import numpy as np
from BEMT2 import Cl_slope, Theta, Chord, Sigma

# functions related to linear blade property variations
from BEMT2 import Cl_slope, Theta, Chord, Sigma
# BEMT functions
from BEMT2 import Lambda, Alpha, dCT, dCPi, dCP0_NACA0012, dCP0_GENERAL
import tools
from tools import calc_rotor_profiles, calc_CT_CPi_CP0

def plot_rotor2(rotor2):
    '''Function to plot rotor2 spec'''
    # Non-dimensional radius
    rs = np.linspace(0, 1, 100)

    plt.rcParams['figure.figsize'] = [14, 5]

    plt.subplot(1, 4, 1)
    plt.title("Blade Twist")
    plt.ylabel(r"$\theta(r)$ (deg)")
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
    plt.ylabel("c(r) (ft)")
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
    print(
        f"Other rotor2 properties: Nb={rotor2.Nb}, R={rotor2.R} ft, vtip = {rotor2.vtip} ft/sec, cd0={rotor2.cd0}")
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    return

def plot_prandtl_comparison(rotor2, airfoil, rs=np.linspace(0.1, 0.99, 100)):
    lambda_list, dCT_list, dCPi_list, dCP0_list, F_list = calc_rotor_profiles(
        rotor2, airfoil, use_F=True)
    
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
    plt.title("Section Thrust Coefficient Distribution")
    plt.plot(rs, dCT_list, label="With Prandtl tip loss")
    plt.plot(rs, dCT(rs, rotor2, F=1), label="No tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$c_T$')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.title("Section Induced Power Coefficient Distribution")
    plt.plot(rs, dCPi_list, label="With Prandtl tip loss")
    plt.plot(rs, dCPi(rs, rotor2, F=1), label="No tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$c_{P,i}$')
    plt.legend()
    plt.show()
    return


def plot_NACA0012_drag_effects(rotor2):
    rs = np.linspace(0.1, 0.99, 100)

    plt.rcParams['figure.figsize'] = [14, 5]
    plt.subplot(1,2,1)
    plt.title("AoA Distrution across rotor radius")
    plt.plot(rs, np.rad2deg(Alpha(rs, rotor2)))
    plt.ylabel(r"Section AoA, $\alpha$ (deg)")
    plt.xlabel("r")
    plt.grid()

    plt.subplot(1,2,2)
    plt.title("Section Power Coefficients")
    plt.plot(rs, dCP0_NACA0012(rs, rotor2), label="Profile (NACA0012)")
    plt.plot(rs, dCP0_GENERAL(rs, rotor2), label="Profile (GENERAL)")
    plt.plot(rs, dCPi(rs, rotor2), label="Induced Power")
    plt.ylabel("Section Power Coefficients across rotor radius")
    plt.xlabel("r")
    plt.legend()
    plt.grid()

    plt.show()
    plt.rcParams['figure.figsize'] = [14, 5]

    CT, CPi, CP0 = calc_CT_CPi_CP0(rotor2, airfoil="NACA0012")
    print(
        f"NACA0012 airfoils: Profile contribution to total power = {(CP0/(CPi+CP0)):.3f}")

    CT, CPi, CP0 = calc_CT_CPi_CP0(rotor2, airfoil="GENERAL")
    print(
        f"GENERAL airfoil: Profile contribution to total power = {(CP0/(CPi+CP0)):.3f}")
    return


def vary_propertyX_plot_profiles(rotor, propX, X_list, X_labels, airfoil="NACA0012", use_F=True):
    plt.rcParams['figure.figsize'] = [14, 5]

    rs = np.linspace(0.1, 0.99, 100)

    for i, xVal in enumerate(X_list):
        rotor.setProperty('_'+propX, xVal)

        lambda_list, dCT_list, dCPi_list, dCP0_list, F_list = tools.calc_rotor_profiles(
            rotor, airfoil=airfoil, use_F=use_F, rs=rs)
        dCP_list = dCPi_list + dCP0_list

        plt.subplot(1, 4, 1)
        plt.plot(rs, lambda_list,
                 label=X_labels[i])

        plt.subplot(1, 4, 2)
        plt.plot(rs, dCT_list,
                 label=X_labels[i])

        plt.subplot(1, 4, 3)
        plt.plot(rs, dCP_list,
                 label=X_labels[i])
        
        plt.subplot(1, 4, 4)
        plt.plot(rs, F_list,
                 label=X_labels[i])

    plt.subplot(1, 4, 1)
    plt.title("Inflow Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Local Inflow Ratio, $\lambda(r)$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 2)
    plt.title("Section Thrust Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Section thrust coefficient, $c_T$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 3)
    plt.title("Section Power Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Section power coefficient (total), $c_P$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 4)
    plt.title("Prandtl Tip Loss Function")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Tip Loss, $F$")
    plt.legend()
    plt.xlim([0.7,1.0])
    plt.grid()

    plt.tight_layout()

    plt.show()
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    return
