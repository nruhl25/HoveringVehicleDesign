# Author: Nathaniel Ruhl
# MEAM5460 Project 2

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Rotor properties
chord = 0.12192  # m
R = 1.524  # m
Nb = 3  # number of blades
A = np.pi*R**2  # m^2
sigma = Nb*chord/(np.pi*R)  # solidity
theta_tw = -8*np.pi/180  # rad
cl_alpha = 2*np.pi
vtip = 198.12  # m/s
gamma = 7  # enforce gamma
cd_0 = 0.01 # section coefficient of drag
rho = 1.225  # kg/m^3

# Desired output state
beta_1c = 0
beta_1s = 0

# initial guess of x for Newton iterations
n_dim = 5
x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.1])

# Rotor disc angle alpha is unique for given CX and CT
def calc_alpha(CX, CT, v_inf, nu_b):
    return np.arcsin((0.5*CX*v_inf**2)/(CT*vtip**2))

def calc_mu(v_inf, alpha):
    return v_inf*np.cos(alpha)/vtip

'''INPUT x = [beta_0, theta_0, theta_1c, theta_1s, lambda]'''
def f(x, CX, CT, v_inf, nu_b):
    # unpack input
    beta_0, theta_0, theta_1c, theta_1s, lmbda = x

    alpha = calc_alpha(CX,CT,v_inf,nu_b)
    mu = calc_mu(v_inf, alpha)

    # RHS c vector given in problem
    c = np.zeros(3)
    c[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2)*theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    c[1] = (1+(mu**2)/2)*theta_1c
    c[2] = (8/3)*mu*theta_0+2*mu*theta_tw+(1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    # Assemble the non-linear equations, evaluated given the inputs
    f = np.zeros(5)
    f[0] = beta_0*8*nu_b**2/gamma - c[0]
    f[1] = (4/3)*mu*beta_0 - c[1]
    f[2] = - c[2]
    f[3] = lmbda - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))
    f[4] = CT - (sigma*cl_alpha/2)*((theta_0/3)*(1+1.5*mu**2)+0.25*theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    return f

def fCQ(x, mu):
    # unpack input
    beta_0, theta_0, theta_1c, theta_1s, lmbda = x
    CQ = (sigma*cl_alpha/2)*((lmbda*theta_0/3)+(lmbda*theta_tw/4)-0.5*lmbda**2-(0.25*mu**2*beta_0**2)+cd_0*(1+mu**2)/(4*cl_alpha))
    return CQ

# User interface to return index for state variable (post-processing state)
def get_state_var_index(state_var):
    if state_var == "beta_0":
        return 0
    elif state_var == "theta_0":
        return 1
    elif state_var == "theta_1c":
        return 2
    elif state_var == "theta_1s":
        return 3
    elif state_var == "lambda":
        return 4
    elif state_var == "TP":
        return 5
    else:
        raise RuntimeWarning("Invalid user-entered argument")

# User interface to return plot ylabel for state variable
def get_state_ylabel(state_var):
    if state_var == "beta_0":
        return r"Coning Angle, $\beta_0$ (deg)"
    elif state_var == "theta_0":
        return r"Input Collective, $\theta_0$ (deg)"
    elif state_var == "theta_1c":
        return r"Longitudinal Input, $\theta_{1c}$ (deg)"
    elif state_var == "theta_1s":
        return r"Latitudianal Input, $\theta_{1s}$ (deg)"
    elif state_var == "lambda":
        return r"Inflow Ratio, $\lambda$"
    elif state_var == "TP":
        return r"Total Power, $P$ (hp)"
    else:
        raise RuntimeWarning("Invalid user-entered argument")
    
def mps2knots(v_mps):
    # converts velocity from mps to knots
    return v_mps*1.94384

def knots2mps(v_knots):
    return v_knots/1.94384

# Main function to calculate trim as a function of CX and CT
def calc_trims_array(CX_list, CT_list, v_inf, nu_b, x0=x0):
    trims = np.zeros((len(CX_list), len(CT_list), n_dim+1), dtype=float)
    for i in range(len(CX_list)):
        for j in range(len(CT_list)):
            alpha = calc_alpha(CX_list[i], CT_list[i], v_inf, nu_b)
            def func(x): return f(x, CX_list[i], CT_list[j], v_inf, nu_b)
            x_root = fsolve(func, x0)
            mu = calc_mu(v_inf, alpha)
            CP = fCQ(x_root, mu)
            P_tot = 0.0013410221*CP*rho*A*vtip**3  # horsepower
            x_root = np.append(x_root,P_tot)
            trims[i, j] = x_root
    return trims


# Function to plot a single state variable. Only use AFTER trim has been calculated and trim_tuple has been built accordingly
def plot_state_var(trim_tuple, state_var_string, label_type="CX"):
    CX_list, CT_list, v_inf, nu_b, trims = trim_tuple
    for i in range(len(CX_list)):
        state_index = get_state_var_index(state_var_string)
        state_ylabel = get_state_ylabel(state_var_string)
        state_array = trims[i, :, state_index]
        if (state_var_string != "lambda") or (state_var_string != "TP"):
            state_array = np.rad2deg(state_array)

        # Determing the plot label type
        if label_type == "CX":
            my_title = fr"Trim Analysis for $v_\infty$={mps2knots(v_inf):.1f} knots and $\nu_\beta$={nu_b:.3f}"
            my_label = fr"$C_X$={CX_list[i]:.3f}"
        elif label_type == "nu_b":
            my_title = fr"Trim Analysis for $v_\infty$={mps2knots(v_inf):.1f} knots and $C_X$={CX_list[0]:.3f}"
            my_label = fr"$\nu_\beta$={nu_b:.3f}"
        elif label_type == "v_inf":
            my_title = fr"Trim Analysis for $C_X$={CX_list[0]:.3f} and $\nu_\beta$={nu_b:.3f}"
            my_label = fr"$v_\infty$={mps2knots(v_inf):.1f} knots"
        
        plt.plot(CT_list/sigma, state_array, label=my_label)
    plt.grid()
    plt.legend()
    plt.xlabel(r"$C_T/\sigma$")
    plt.ylabel(state_ylabel)
    plt.title(my_title)
    return

# This is the main plotting function to plot the entire trim
def plot_rotor_trim(trim_tuple):
    CX_list, CT_list, v_inf, nu_b, trims = trim_tuple

    fig, ax = plt.subplots(2, 3, figsize=(18, 12))
    state_var_list = ["beta_0", "theta_0", "theta_1c", "theta_1s", "lambda", "TP"]
    r = [0,0,0,1,1,1]
    c = [0,1,2,0,1,2]
    for indx, state_var_string in enumerate(state_var_list):
        for i in range(len(CX_list)):
            state_index = get_state_var_index(state_var_string)
            state_ylabel = get_state_ylabel(state_var_string)
            state_array = trims[i, :, state_index]

            if (state_var_string != "lambda") or (state_var_string != "TP"):
                state_array = np.rad2deg(state_array)

            ax[r[indx],c[indx]].plot(CT_list/sigma, state_array,
                    label=fr"$C_X$={CX_list[i]:.3f}")
        ax[r[indx], c[indx]].legend()
        ax[r[indx], c[indx]].set_xlabel(r"$C_T/\sigma$")
        ax[r[indx], c[indx]].set_ylabel(state_ylabel)
        ax[r[indx], c[indx]].grid()
    
    ax[0,1].set_title(fr"Trim Analysis for $v_\infty$={mps2knots(v_inf): .2f} knots and $\nu_\beta$={nu_b: .3f}")
    plt.tight_layout()
    plt.show()
    return

