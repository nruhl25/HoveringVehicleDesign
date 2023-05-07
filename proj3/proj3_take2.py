# Author: Nathaniel Ruhl
# MEAM5460: Project 3

import numpy as np
from numpy import sin, cos
from scipy.optimize import fsolve

rho = 0.002378  # slug/ft^3
a = 2*np.pi  # cl_alpha
gamma = 8  # Locke number
kappa = 1.15  # ideal power correction
W = 14000+6000  # lbf, wieght of helicopter (wet weight)
d_copter = 7.75  # ft, diameter of helicopter
A_copter = np.pi*(d_copter/2)**2

nu_b = 1.0  # flapping frequency
v_inf = 205  # ft/sec

# Constants to choose
L_emp = 2  # ft
c_emp = 0.5  # ft

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

# Things to solve
CD_sphere = 0.5
Df = CD_sphere*A_copter*0.5*rho*v_inf**2
H = Nb*cd0*R*(rho*A_blade*v_inf**2)
CH = H/(rho*A*R*vtip**2)

# Thrust and disk-loading of UH-60
DL = 9.33  # lbf/ft^2
T = DL*A
CT = T/(rho*A*vtip**2)

#### INPUTS TO SOLVER PROCESS #####
xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = 0  # tail rotor angle
psi_emp = 0  # empenage angle

################################################

#1) Solve for alpha, Lht
Ty = W
Tx = np.sqrt(T**2-Ty**2)
alpha = np.arccos((Tx-Df)/H)
Lht = 2*np.pi*(psi_emp-alpha)*L_emp*(rho*v_inf**2)*L_emp*c_emp

# 2) Solve mu and lmbda
mu = v_inf*cos(alpha)/vtip

def f(lmbda): return lmbda - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))

sol = fsolve(f,0.1)
lmbda = sol[0]

#3) Solve theta_0, theta_1c, theta_1s, beta_0
def f1(x):
    theta_0, theta_1c, theta_1s, beta_0 = x

    # RHS vector in the equations
    RHS = np.zeros(3)
    RHS[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2) * \
        theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    RHS[1] = (1+(mu**2)/2)*theta_1c
    RHS[2] = (8/3)*mu*theta_0+2*mu*theta_tw+(1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    f1 = np.zeros(4)
    f1[0] = beta_0*8*nu_b**2/gamma - RHS[0]
    f1[1] = (4/3)*mu*beta_0 - RHS[1]
    f1[2] = RHS[2]
    f1[3] = CT - (sigma*a/2)*((theta_0/3)*(1+1.5*mu**2) +
                             0.25*theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    return f1

sol1 = fsolve(f1,np.zeros(4))
theta_0, theta_1c, theta_1s, beta_0 = sol1

# 4) Define CQ
CQ = (sigma*a/2)*((lmbda*theta_0/3)+(lmbda/4)*theta_tw-(1/2)
                  * lmbda**2-(1/4)*(beta_0*mu)**2+(cd0/(4*a))*(1+mu**2))

# 5) Define Mzr, CY, Y
Mzr = -CQ*(rho*A*R*vtip**2)
CY = (sigma*a/2)*(theta_0*(3/4)*mu*beta_0+theta_tw*(1/2)*mu*beta_0+theta_1c *
             (1/4)*lmbda+theta_1s*(1/6)*beta_0*(1+3*mu**2)-(3/2)*mu*lmbda*beta_0)
Y = CY*rho*A*vtip**2

#6) Solve 2 eqn system for Ttr and phi
Yf = Y
Ttr = 5000
phi = -(Y+Yf+Ttr)/(T-Lht)
for i in range(1000):
    phi_last = phi
    Ttr_last = Ttr
    phi = (Y+Yf+Ttr)/(T-Lht)
    Ttr = (T*cos(alpha)-Y*phi+H*sin(alpha)-Lht-W)/(phi-psi_tr)
    if (np.abs(Ttr - Ttr_last)<1):
        print(f"Ttr={Ttr}, Ttr_last={Ttr_last}")
        break

def f3(x):
    Mxf, Myf, Mzf = x
    f3 = np.zeros(3)
    f3[0] = Mxf - W*ycg + W*phi*hcg + Yf*hcg + Yf*phi*ycg + Ttr*(hcg - htr)
    f3[1] = Myf - W*cos(alpha)*xcg + W*sin(alpha)*hcg - Df*cos(alpha)*hcg - Df*sin(alpha)*xcg - Lht*(xcg-xht)
    f3[2] = Mzr + Mzf + Ttr*(xtr-xcg)-Df*cos(alpha)*ycg - Yf*xcg
    return f3

sol = fsolve(f3,np.ones(3))