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
kappa = 1.15

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
A_tr = np.pi*R_tr**2
Omega_tr = 270  # rad/sec
vtip_tr = Omega_tr*R_tr  # ft/sec
Nb_tr = 4
sigma_tr = Nb_tr*c_tr/(np.pi*R_tr)

################################################


def ff1(lmbda, mu, alpha, CT):
    return lmbda - mu*np.tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))


def ff2(x, xcg, ycg, hcg, xht, xtr, htr, Lht, Df, W, alpha, lmbda, mu, CT, CH):
    theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, phi, Mxf, Myf, Mzf = x
    # Define helpful relationships using primary unknowns
    CQ = (sigma*a/2)*((lmbda*theta_0/3)+(lmbda/4)*theta_tw-(1/2)*lmbda**2-(1/8)*(beta_1c**2+beta_1s)**2-0.5*(0.5*beta_0**2+(3./8.)*beta_1c**2+(1/8)*beta_1s**2)+(cd0/(4*a))*(1+mu**2)-0.5*mu*lmbda*beta_1c-mu*theta_0*beta_1s/3)
    CH_rhs = (sigma*a/2)*(theta_0*(-beta_1c/3+0.5*mu*lmbda)+theta_tw*(-beta_1c/4+mu*lmbda/4) -
                      (1/6)*theta_1c*beta_0+theta_1s*(-0.25*beta_1c+0.25*mu*lmbda)-0.75*mu*lmbda*beta_1c+(1/6)*beta_0*beta_1s+0.25*mu*(beta_0**2+beta_1c**2)+cd0*mu/(2*a))
    CY = (sigma*a/2)*(theta_0*(0.75*mu*beta_0+beta_1s*(1+1.5*mu**2)/3)+theta_tw*(0.5*mu*beta_0+0.25*beta_1s*(1+mu**2))+theta_1c*((1/4)*lmbda+(1/4)*mu*beta_1c)+theta_1s*((1/6)*beta_0*(1+3*mu**2)+0.5*mu*beta_1s)-(3/4)*lmbda*beta_1s+beta_0*beta_1c*((1/6)-mu**2)-(3/2)*mu*lmbda*beta_0-0.25*mu*beta_1c*beta_1s)
    Y = CY*(rho*A*vtip**2)
    H = CH*(rho*A*vtip**2)  # enforced
    T = CT*(rho*A*vtip**2)   # enforced
    CQ_tr = 0.05*CQ
    Ttr = CQ_tr*rho*A_tr*vtip_tr**2
    Yf = W*phi
    CMx = (sigma*a/(2*gamma))*(nu_b**2-1)*beta_1s
    CMy = (sigma*a/(2*gamma))*(nu_b**2-1)*beta_1c
    Mxr = CMx*(rho*A*R*vtip**2)
    Myr = CMy*(rho*A*R*vtip**2)
    Mzr = -CQ*(rho*A*R*vtip**2)

    f2 = np.zeros(10)
    # RHS vector in the equations
    RHS = np.zeros(3)
    MAT = np.zeros((3,3))
    RHS[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2) * \
        theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    RHS[1] = (1+(mu**2)/2)*theta_1c
    RHS[2] = (8/3)*mu*theta_0+2*mu*theta_tw + \
        (1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    MAT[0,0] = (8/gamma)*nu_b**2
    MAT[1,0] = 4*mu/3
    MAT[1,1] = 8*(nu_b**2-1)/gamma
    MAT[1,2] = 1+0.5*mu**2
    MAT[2,1] = -(1-0.5*mu**2)
    MAT[2, 2] = 8*(nu_b**2-1)/gamma

    # Assemble 10 equations
    f2[0] = T*phi + Y + Yf + Ttr - Lht*phi
    f2[1] = T*cos(alpha) - Y*phi + H*sin(alpha) - Lht - W - Ttr*(phi-psi_tr)
    f2[2] = Mxr + Mxf - W*ycg + W*phi*hcg + \
        Yf*hcg + Yf*phi*ycg + Ttr*(hcg - htr)
    f2[3] = Myr + Myf - W*cos(alpha)*xcg + W*sin(alpha)*hcg - \
        Df*cos(alpha)*hcg - Df*sin(alpha)*xcg - Lht*(xcg-xht)
    f2[4] = Mzr + Mzf + Ttr*(xtr-xcg)-Df*cos(alpha)*ycg - Yf*xcg
    f2[5] = CT - (sigma*a/2)*((theta_0/3)*(1+1.5*mu**2) + 0.25 *
                              theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    f2[6] = CH - CH_rhs
    f2[7] = MAT[0, 0]*beta_0+MAT[0, 1]*beta_1c+MAT[0, 2]*beta_1s - RHS[0]
    f2[8] = MAT[1, 0]*beta_0+MAT[1, 1]*beta_1c+MAT[1, 2]*beta_1s - RHS[1]
    f2[9] = MAT[2, 0]*beta_0+MAT[2, 1]*beta_1c+MAT[2, 2]*beta_1s - RHS[2]
    return f2

# Main function for solving the trim
def solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b):

    # Estimate drag on airframe
    CD_sphere = 0.5
    Df = CD_sphere*A_copter*0.5*rho*v_inf**2

    # Assumption about thrust to start the iteration
    T = W
    CT = T/(rho*A*vtip**2)

    # Estimate rotor drag CH (get it from no flap solution)
    H = Nb*cd0*R*(rho*A_blade*v_inf**2)
    CH = H/(rho*A*vtip**2)

    # Step 1) Get CQ to converge by H changing
    # A) Solve for alpha, Lht
    alpha = (Df+H)/T
    Lht = 2*np.pi*(psi_emp-alpha)*L_emp*(rho*v_inf**2)*L_emp*c_emp

    # B) Solve mu and lmbda
    mu = v_inf*cos(alpha)/vtip
    def f1(lmbda): return ff1(lmbda, mu, alpha, CT)
    sol = fsolve(f1, 0.1)
    lmbda = sol[0]

    # Step 2) Solve everything else
    def f2(x): return ff2(x, xcg, ycg, hcg, xht, xtr, htr, Lht, Df, W, alpha, lmbda, mu, CT, CH)
    root = fsolve(f2, np.zeros(10))
    theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, phi, Mxf, Myf, Mzf = root
    print(theta_0)

    # Create solution output form
    sol_dict = {}

    return sol_dict

#### INPUTS TO SOLVER PROCESS #####
v_inf = 205 # ft/sec

xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = np.deg2rad(6)  # tail rotor angle
psi_emp = 0  # empenage angle
nu_b = 1
sol_dict = solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b)
