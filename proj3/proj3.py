# Author: Nathaniel Ruhl
# MEAM5460: Project 3

import numpy as np
from numpy import sin, cos
from scipy.optimize import fsolve

# Global variables
rho = 0.002378  # slug/ft^3
a = 2*np.pi  # cl_alpha
gamma = 8  # Locke number
W = 14000+6000  # lbf, wieght of helicopter (wet weight)
d_copter = 7.75  # ft, diameter of helicopter
A_copter = np.pi*(d_copter/2)**2

# Constants to choose
L_emp = 2  # ft
c_emp = 0.5  # ft
kappa = 1.15 # ideal power coefficient

# Main rotor
cd0 = 0.011  # section drag coefficient
c = 2  # ft chord
R = 27  # ft radius
A = np.pi*R**2  # ft^2 disc area
Omega = 27  # rad/sec
vtip = Omega*R  # ft/sec
Nb = 4
sigma = Nb*c/(np.pi*R)
theta_tw = -np.deg2rad(8)
A_blade = R*c

# Tail rotor
c_tr = 0.5  # ft chord
R_tr = 5.5  # ft
Omega_tr = 270  # rad/sec
vtip_tr = Omega_tr*R_tr  # ft/sec
Nb_tr = 4

################################################

def ff1(lmbda, mu, alpha, CT): 
    return lmbda - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))

def ff2(x, mu, lmbda, CT, nu_b):
    theta_0, theta_1c, theta_1s, beta_0 = x

    # RHS vector in the equations
    RHS = np.zeros(3)
    RHS[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2) * \
        theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    RHS[1] = (1+(mu**2)/2)*theta_1c
    RHS[2] = (8/3)*mu*theta_0+2*mu*theta_tw + \
        (1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    f2 = np.zeros(4)
    f2[0] = beta_0*8*nu_b**2/gamma - RHS[0]
    f2[1] = (4/3)*mu*beta_0 - RHS[1]
    f2[2] = RHS[2]
    f2[3] = CT - (sigma*a/2)*((theta_0/3)*(1+1.5*mu**2) +
                              0.25*theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    return f2

def ff3(x, xcg, ycg, hcg, xht, xtr, htr, Ttr, Lht, Y, Yf, Df, W, Mzr, alpha, phi):
    Mxf, Myf, Mzf = x
    f3 = np.zeros(3)
    f3[0] = Mxf - W*ycg + W*phi*hcg + Yf*hcg + Yf*phi*ycg + Ttr*(hcg - htr)
    f3[1] = Myf - W*cos(alpha)*xcg + W*sin(alpha)*hcg - \
        Df*cos(alpha)*hcg - Df*sin(alpha)*xcg - Lht*(xcg-xht)
    f3[2] = Mzr + Mzf + Ttr*(xtr-xcg)-Df*cos(alpha)*ycg - Yf*xcg
    return f3

def ffphi(x, T, W, Y, Lht, H, alpha):
    phi, Ttr = x
    Yf = W*phi
    f = np.zeros(2)
    f[0] = T*phi + Y + Yf + Ttr - Lht
    f[1] = T*cos(alpha) - Y*phi + H*sin(alpha) - \
        Lht - W - Ttr*(phi-psi_tr)
    return f


# Main equation for solving the trim
def solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b):

    # Estimate drag on airframe
    CD_sphere = 0.5
    Df = CD_sphere*A_copter*0.5*rho*v_inf**2

    # Assumption about thrust to start the iteration
    T = W
    CT = T/(rho*A*vtip**2)
    CQ = (sigma*cd0/8)+kappa*CT**(3./2.)/np.sqrt(2.)

    # Estimate rotor drag CH to start the iteration
    H = Nb*cd0*R*(rho*A_blade*v_inf**2)
    CH = H/(rho*A*vtip**2)

    # Step 1) Get CQ to converge by H changing
    for num_iter in range(50):
        CH_last = CH
        CT_last = CT
        CQ_last = CQ
        # A) Solve for alpha, Lht
        # Ty = W    # ASSUMPTION
        # Tx = np.sqrt(T**2-Ty**2)
        alpha = (Df+H)/T
        Lht = 2*np.pi*(psi_emp-alpha)*L_emp*(rho*v_inf**2)*L_emp*c_emp

        # B) Solve mu and lmbda
        mu = v_inf*cos(alpha)/vtip

        def f1(lmbda): return ff1(lmbda, mu, alpha, CT)
        sol = fsolve(f1, 0.1)
        lmbda = sol[0]

        # C) Solve theta_0, theta_1c, theta_1s, beta_0
        def f2(x): return ff2(x, mu, lmbda, CT, nu_b)
        sol2 = fsolve(f2, np.zeros(4))
        theta_0, theta_1c, theta_1s, beta_0 = sol2

        # D) Define CQ, Mzr, CY, Y
        CQ = (sigma*a/2)*((lmbda*theta_0/3)+(lmbda/4)*theta_tw-(1/2)
                        * lmbda**2-(1/4)*(beta_0*mu)**2+(cd0/(4*a))*(1+mu**2))
        CT_from_CQ = (np.sqrt(2)*(CQ-sigma*cd0/8)/kappa)**(2./3.)
        Mzr = -CQ*(rho*A*R*vtip**2)
        CY = (sigma*a/2)*(theta_0*(3/4)*mu*beta_0+theta_tw*(1/2)*mu*beta_0+theta_1c *
                        (1/4)*lmbda+theta_1s*(1/6)*beta_0*(1+3*mu**2)-(3/2)*mu*lmbda*beta_0)
        Y = CY*rho*A*vtip**2

        # F) Update H, T, Q
        CH = (sigma*a/2)*(0.5*theta_0*mu*lmbda+0.25*theta_tw*mu*lmbda -
                        (1/6)*theta_1c*beta_0+0.25*theta_1s*mu*lmbda+cd0*mu/(2*a))
        # CH = growth_rate*CH_new + (1-growth_rate)*CH
        H = CH*(rho*A*vtip**2)
        # T = CT*(rho*A*vtip**2)

        print(f"alpha={np.rad2deg(alpha)}")
        print(f"CH_last={CH_last}, CH={CH}")
        print(f"CQ_last={CQ_last}, CQ={CQ}")
        print(f"CT_last={CT_last}, CT={CT}")
        print("----------------")

        if (np.abs(CT-CT_from_CQ) and np.abs(CH-CH_last)) < 1e-9:
            print(f"num_iter={num_iter+1}")
            break

    # Step 2) Solve phi, Ttr
    def fphi(x): return ffphi(x, T, W, Y, Lht, H, alpha)
    sol = fsolve(fphi, np.array([0.01, 1]))
    phi, Ttr = sol
    Yf = W*phi

    # Step 3) Solve Moments Mxf, Myf, Mzf
    def f3(x): return ff3(x, xcg, ycg, hcg, xht, xtr,
                        htr, Ttr, Lht, Y, Yf, Df, W, Mzr, alpha, phi)
    sol = fsolve(f3, np.ones(3))
    Mxf, Myf, Mzf = sol
    print(f"--------------------")

    # Creat solution output form
    sol_dict = {}
    sol_dict['T'] = T
    sol_dict['alpha'] = alpha
    sol_dict["phi"] = phi
    sol_dict['beta_0'] = beta_0
    sol_dict['theta_0'] = theta_0
    sol_dict['theta_1c'] = theta_1c
    sol_dict['theta_1s'] = theta_1s

    return sol_dict

#### INPUTS TO SOLVER PROCESS #####
v_inf = 50 # ft/sec

xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = 0  # tail rotor angle
psi_emp = 0  # empenage angle
nu_b = 1

sol_dict = solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b)
print(sol_dict)

# xcg_list = np.linspace(-5,5,20)
# beta_0_list = []
# theta_0_list = []
# alpha_list = []
# phi_list = []
# for i in range(len(xcg_list)):
#     print(f"xcg = {xcg_list[i]}")
#     sol_dict = solve_trim_system(xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b)
#     beta_0_list.append(np.rad2deg(sol_dict["beta_0"]))
#     theta_0_list.append(np.rad2deg(sol_dict["theta_0"]))
#     alpha_list.append(np.rad2deg(sol_dict["alpha"]))
#     phi_list.append(np.rad2deg(sol_dict["phi"]))

# import matplotlib.pyplot as plt
# plt.plot(xcg_list, np.rad2deg(phi_list))
# plt.ylabel(r"$\phi$ (deg)")
# plt.xlabel(r"$x_{cg}$")
# plt.show()

