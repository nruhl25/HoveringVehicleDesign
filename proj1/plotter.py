

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

    plt.subplot(1, 3, 1)
    plt.title("Blade Twist")
    plt.ylabel(r"$\theta(r)$ (deg)")
    plt.xlabel("r")
    plt.plot(rs, np.rad2deg(Theta(rs, rotor2)))
    plt.grid()

    plt.subplot(1, 3, 2)
    plt.title("Lift-Curve Slope")
    plt.ylabel(r"$c_{l,\alpha}(r)$")
    plt.xlabel("r")
    plt.plot(rs, Cl_slope(rs, rotor2))
    plt.grid()

    plt.subplot(1, 3, 3)
    plt.title("Chord")
    plt.ylabel("c(r) (ft)")
    plt.xlabel("r")
    plt.plot(rs, Chord(rs, rotor2),
             label=rf"$\sigma_{{ave}}={Sigma(0, rotor2):.3f}$")
    plt.grid()
    plt.legend()

    plt.tight_layout()
    plt.show()
    print(
        f"Other rotor2 properties: Nb={rotor2.Nb}, R={rotor2.R} ft, vtip = {rotor2.vtip} ft/sec")
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
    plt.ylabel(r'$c_t$')
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.title("Section Induced Power Coefficient Distribution")
    plt.plot(rs, dCPi_list, label="With Prandtl tip loss")
    plt.plot(rs, dCPi(rs, rotor2, F=1), label="No tip loss")
    plt.grid()
    plt.xlabel('r')
    plt.ylabel(r'$c_{p,i}$')
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
        CT, CPi, CP0 = tools.calc_CT_CPi_CP0(rotor, airfoil=airfoil, use_F=use_F)
        CP = CPi + CP0
        p_ratio = CPi/CP0

        plt.subplot(1, 4, 1)
        plt.plot(rs, lambda_list,
                 label=X_labels[i])

        plt.subplot(1, 4, 2)
        plt.plot(rs, dCT_list,
                 label=X_labels[i]+fr", $C_T$={CT:.5f}")

        plt.subplot(1, 4, 3)
        plt.plot(rs, dCP_list,
                 label=X_labels[i]+fr", $C_P$={CP:.5f}")
        
        plt.subplot(1, 4, 4)
        plt.plot(rs, F_list,
                 label=X_labels[i])
        
        print(f"{X_labels[i]}: Power ratio: C_Pi/C_P0={p_ratio:.3f}, solidity = {Sigma(0,rotor):.3f}")


    plt.subplot(1, 4, 1)
    plt.title("Inflow Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Local Inflow Ratio, $\lambda(r)$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 2)
    plt.title("Section Thrust Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Section thrust coefficient, $c_t$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 3)
    plt.title("Section Power Distribution")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Section power coefficient (total), $c_p$")
    plt.legend()
    plt.grid()

    plt.subplot(1, 4, 4)
    plt.title("Prandtl Tip Loss Function")
    plt.xlabel("Non-dimensional radial position, r")
    plt.ylabel(r"Tip Loss, $F$")
    plt.legend()
    plt.xlim([0.8,1.0])
    plt.grid()

    plt.tight_layout()

    plt.show()
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    return

def compare_rotors(rotor_list):
    rs = np.linspace(0.01, 0.99, 100)


    CT_list = np.zeros(4)
    CPi_list = np.zeros(4)
    CP0_list = np.zeros(4)
    CP_list = np.zeros(4)

    plt.rcParams['figure.figsize'] = [14, 10]
    fig, ax = plt.subplots(2, 3)
    for i, rotor in enumerate(rotor_list):
        j = i+1
        lambda_list, dCT_list, dCPi_list, dCP0_list, F_list = tools.calc_rotor_profiles(
            rotor, airfoil="NACA0012", use_F=True, rs=rs)
        CT_list[i], CPi_list[i], CP0_list[i] = tools.calc_CT_CPi_CP0(
            rotor, airfoil="NACA0012", use_F=True)
        CP_list[i] = CPi_list[i] + CP0_list[i]

        ax[0, 0].plot(rs, lambda_list, label=f"Rotor {j}")

        ax[0, 1].plot(rs, dCT_list, label=fr"Rotor {j}: $C_T$={CT_list[i]:.6f}")

        ax[0, 2].plot(rs, F_list, label=fr"Rotor {j}")

        ax[1, 0].plot(
            rs, dCPi_list, label=fr"Rotor {j}: $C_{{P_i}}$={CPi_list[i]:.6f}")

        ax[1, 1].plot(
            rs, dCP0_list, label=fr"Rotor {j}: $C_{{P_0}}$={CP0_list[i]:.6f}")

        ax[1, 2].plot(rs, dCP0_list+dCPi_list,
                    label=fr"Rotor {j}: $C_{{P}}$={CP_list[i]:.6f}")


    ax[0, 0].set_title("Induced Inflow Ratio")
    ax[0, 0].set_ylabel(r"$\lambda$")
    ax[0, 0].set_xlabel("r")
    ax[0, 0].grid()
    ax[0, 0].legend()

    ax[0, 1].set_title("Section Thrust Coefficient")
    ax[0, 1].set_ylabel(r"$c_t$")
    ax[0, 1].set_xlabel("r")
    ax[0, 1].grid()
    ax[0, 1].legend()

    ax[0, 2].set_title("Prandtl Tip Loss")
    ax[0, 2].set_ylabel(r"$F$")
    ax[0, 2].set_xlabel("r")
    ax[0, 2].grid()
    ax[0, 2].legend()
    ax[0, 2].set_xlim([0.7,1.0])

    ax[1, 0].set_title("Section Coefficient of Induced Power")
    ax[1, 0].set_ylabel("$c_{p,i}$")
    ax[1, 0].set_xlabel("r")
    ax[1, 0].grid()
    ax[1, 0].legend()

    ax[1, 1].set_title("Section Coefficient of Profile Power")
    ax[1, 1].set_ylabel("$c_{p_0}$")
    ax[1, 1].set_xlabel("r")
    ax[1, 1].grid()
    ax[1, 1].legend()

    ax[1, 2].set_title("Total Section Coefficient of Power")
    ax[1, 2].set_ylabel("$c_{p}$")
    ax[1, 2].set_xlabel("r")
    ax[1, 2].grid()
    ax[1, 2].legend()

    plt.tight_layout()
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    return

def compare_aoas(rotor_list):
    rs = np.linspace(0.01, 0.99, 100)

    CT_list = np.zeros(4)
    CPi_list = np.zeros(4)
    CP0_list = np.zeros(4)
    CP_list = np.zeros(4)

    for i, rotor in enumerate(rotor_list):
        j = i+1
        lambda_list, dCT_list, dCPi_list, dCP0_list, F_list = tools.calc_rotor_profiles(
            rotor, airfoil="NACA0012", use_F=True, rs=rs)
        CT_list[i], CPi_list[i], CP0_list[i] = tools.calc_CT_CPi_CP0(
            rotor, airfoil="NACA0012", use_F=True)
        CP_list[i] = CPi_list[i] + CP0_list[i]

        plt.figure(1)
        plt.plot(rs, np.rad2deg(Alpha(
            rs, rotor, F=F_list)), label=f"Rotor {j}")


    plt.figure(1)
    plt.title("Section Angle of Attack along rotor blade")
    plt.xlabel("r")
    plt.ylabel("Section Angle of Attack")
    plt.grid()
    plt.legend()
    return

def convergence_test(rotor):
    N_list = np.linspace(100, 5000, 20, dtype=int)

    CT_list = []
    CPi_list = []
    CP0_list = []
    for N in N_list:
        CT, CPi, CP0 = tools.calc_CT_CPi_CP0(
            rotor, airfoil="GENERAL", use_F=True, N_int=N)
        CT_list.append(CT)
        CPi_list.append(CPi)
        CP0_list.append(CP0)

    plt.rcParams['figure.figsize'] = [14, 5]
    plt.suptitle("Convergence Test For Numerical Integration")

    plt.subplot(1, 3, 1)
    plt.plot(N_list, CT_list, label=r"$C_T$")
    plt.xlabel("Number of integration intervals")
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 2)
    plt.plot(N_list, CPi_list, label=r"$C_{P_i}$")
    plt.xlabel("Number of integration intervals")
    plt.grid()
    plt.legend()

    plt.subplot(1, 3, 3)
    plt.plot(N_list, CP0_list, label=r"$C_{P_0}$")
    plt.xlabel("Number of integration intervals")
    plt.grid()
    plt.legend()

    plt.tight_layout()
    return
