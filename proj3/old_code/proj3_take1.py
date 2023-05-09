# Author: Nathaniel Ruhl
# MEAM5460: Project 3

import numpy as np
from numpy import sin, cos
from scipy.optimize import fsolve

rho = 0.002378 # slug/ft^3
a = 2*np.pi # cl_alpha
gamma = 8 # Locke number
kappa = 1.15 # ideal power correction
W = 13650 # lbf, wieght of helicopter
d_copter = 7.75 # ft, diameter of helicopter
A_copter = np.pi*(d_copter/2)**2

nu_b = 1.0  # flapping frequency
v_inf = 238 # ft/sec

# Constants to choose
L_emp = 2 # ft
c_emp = 0.5 # ft

# Main rotor
cd0 = 0.011 # section drag coefficient
c = 2 # ft chord
R = 27 # ft radius
A = np.pi*R**2  # ft^2 disc area
Omega = 27 # rad/sec
vtip = Omega*R # ft/sec
Nb = 4
sigma = Nb*c/(np.pi*R)
theta_tw = -np.deg2rad(8)
A_blade = R*c

# Tail rotor
c_tr = 0.5 # ft chord
R_tr = 5.5 # ft
Omega_tr = 270 # rad/sec
vtip_tr = Omega_tr*R_tr # ft/sec
Nb_tr = 4

# Things to solve
CD_sphere = 0.5
Df = CD_sphere*A_copter*0.5*rho*v_inf**2
H = Nb*cd0*R*(rho*A_blade*v_inf**2)
CH = H/(rho*A*R*vtip**2)

# Thrust and disk-loading of UH-60
DL = 9.33 # lbf/ft^2
T = DL*A
CT = T/(rho*A*vtip**2)

def ff(x, xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp):
    theta_0, theta_1c, theta_1s, beta_0, lmbda, mu, alpha, phi, Lht, Ttr, Yf, Mxf, Myf, Mzr, Mzf, CY, CQ, CT = x
    print(f"theta_1c={theta_1c}")
    # Variations of input parameters
    Y = CY*(rho*A*R*vtip**2)
    T = CT*(rho*A*vtip**2)

    f = np.zeros(n_vars)
    aoa_tr = phi-psi_tr
    f[0] = Df + H*cos(alpha) - T*sin(alpha)
    f[1] = T*sin(phi) + Y*cos(phi) + Yf + Ttr*cos(aoa_tr) - Lht*sin(phi)
    f[2] = T*cos(alpha)*cos(phi) - Y*sin(phi) + H*sin(alpha)-Lht-W-Ttr*sin(aoa_tr)
    f[3] = Mxf - W*cos(phi)*ycg + W*sin(phi)*hcg + Yf*cos(phi)*hcg + Yf*sin(phi)*ycg + Ttr*(hcg - htr)*cos(aoa_tr)
    f[4] = Myf - W*cos(alpha)*xcg + W*sin(alpha)*hcg - Df*cos(alpha)*hcg - Df*sin(alpha)*xcg - Lht*(xcg-xht)
    f[5] = Mzr + Mzf + Ttr*(xtr-xcg)*cos(aoa_tr)-Df*cos(alpha)*ycg - Yf*cos(phi)*xcg

    f[6] = CT - (sigma*a/2)*((theta_0/3)*(1+1.5*mu**2) +
                                    0.25*theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    f[7] = CH - (sigma*a/2)*(0.5*theta_0*mu*lmbda+0.25*theta_tw*mu*lmbda-(1/6)*theta_1c*beta_0+0.25*theta_1s*mu*lmbda+cd0*mu/(2*a))
    f[8] = CY - (sigma*a/2)*(theta_0*(3/4)*mu*beta_0+theta_tw*(1/2)*mu*beta_0+theta_1c*(1/4)*lmbda+theta_1s*(1/6)*beta_0*(1+3*mu**2)-(3/2)*mu*lmbda*beta_0)
    f[9] = CQ - (sigma*a/2)*((lmbda*theta_0/3)+(lmbda/4)*theta_tw-(1/2)*lmbda**2-(1/4)*(beta_0*mu)**2+(cd0/(4*a))*(1+mu**2))
    f[10] = lmbda - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))

    # Equations 12-14: 
    # RHS vector in the equations
    RHS = np.zeros(3)
    RHS[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2) * \
        theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    RHS[1] = (1+(mu**2)/2)*theta_1c
    RHS[2] = (8/3)*mu*theta_0+2*mu*theta_tw+(1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    f[11] = beta_0*8*nu_b**2/gamma - RHS[0]
    f[12] = (4/3)*mu*beta_0 - RHS[1]
    f[13] = RHS[2]

    f[14] = CQ + Mzr/(rho*A*R*vtip**2)
    f[15] = mu - v_inf*cos(alpha)/vtip
    f[16] = Lht - 2*np.pi*(psi_emp-alpha)*L_emp*(rho*v_inf**2)*L_emp*c_emp
    f[17] = CQ - ((sigma*cd0/8)+(kappa*CT**(3./2.)/np.sqrt(2)))
    return f

n_vars=18
x0 = 0.1*np.ones(n_vars) # initial guess

# INPUTS TO NON-LINEAR SOLVER
xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = 0  # tail rotor angle
psi_emp = 0  # empenage angle

def f(x): return ff(x, xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp) # in-line
sol = fsolve(f, x0, full_output=True)
root, infodict, ier, mesg = sol
#print(f(root))

print(f(0.1*np.ones(n_vars)))

# xcg_list = np.arange(-0.5, 0.5, 0.05)
# sol_array = np.zeros((n_vars,len(xcg_list)))
# for i, xcg in enumerate(xcg_list):
#     def f(x): return ff(x, xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp) # in-line 
#     sol = fsolve(f, x0, full_output=True)
#     root, infodict, ier, mesg = sol
#     print(f"xcg={xcg}:")
#     print(f(root))
#     print(f"number of f(x) evals = {infodict['nfev']}")
#     print(mesg)
#     print("---")
#     sol_array[:,i] = root

# import matplotlib.pyplot as plt
# plt.plot(xcg_list, np.rad2deg(sol_array[3,:]))
# plt.ylabel(r"$\beta_0$ (deg)")
# plt.xlabel(r"$x_{cg}$ (ft)")
# plt.show()

