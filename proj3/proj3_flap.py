# Author: Nathaniel Ruhl
# MEAM5460: Project 3

import numpy as np
from numpy import sin, cos, tan
from scipy.optimize import fsolve

# Global variables
rho = 0.002378  # slug/ft^3
a = 2*np.pi  # cl_alpha
gamma = 8  # Locke number
W = 14000+6000  # lbf, wieght of helicopter (wet weight)
d_copter = 7.75  # ft, diameter of helicopter
A_copter = np.pi*(d_copter/2)**2

# Constants to choose
c_emp = 0.5  # ft
# L_emp baseline is 2 ft... it will be solved in the loop

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
kappa = 1.15

# Tail rotor
c_tr = 0.5  # ft chord
R_tr = 5.5  # ft
A_tr = np.pi*R_tr**2
Omega_tr = 270  # rad/sec
vtip_tr = Omega_tr*R_tr  # ft/sec
Nb_tr = 4
sigma_tr = Nb_tr*c_tr/(np.pi*R_tr)
theta_tw_tr = 0
kappa_tr = 1.25

################################################

# Equation given in problem set


def fCT(theta_0, theta_1s, mu, lmbda):
    CT = (sigma*a/2)*((theta_0/3)*(1+1.5*mu**2) + 0.25 *
                      theta_tw*(1+mu**2)+0.5*mu*theta_1s-0.5*lmbda)
    return CT


def fCY(theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, mu, lmbda):
    CY = (sigma*a/2)*(theta_0*(0.75*mu*beta_0+beta_1s*(1+1.5*mu**2)/3)+theta_tw*(0.5*mu*beta_0+0.25*beta_1s*(1+mu**2))+theta_1c*((1/4)*lmbda+(1/4)*mu*beta_1c) +
                      theta_1s*((1/6)*beta_0*(1+3*mu**2)+0.5*mu*beta_1s)-(3/4)*lmbda*beta_1s+beta_0*beta_1c*((1/6)-mu**2)-(3/2)*mu*lmbda*beta_0-0.25*mu*beta_1c*beta_1s)
    return CY

# Equation given in problem set


def fCH(theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, mu, lmbda):
    CH = (sigma*a/2)*(theta_0*(-beta_1c/3+0.5*mu*lmbda)+theta_tw*(-beta_1c/4+mu*lmbda/4) -
                      (1/6)*theta_1c*beta_0+theta_1s*(-0.25*beta_1c+0.25*mu*lmbda)-0.75*mu*lmbda*beta_1c+(1/6)*beta_0*beta_1s+0.25*mu*(beta_0**2+beta_1c**2)+cd0*mu/(2*a))
    return CH

# Equation given in problem set


def fCQ(theta_0, beta_0, beta_1c, beta_1s, mu, lmbda):
    CQ = (sigma*a/2)*((lmbda*theta_0/3)+(lmbda/4)*theta_tw-(1/2)*lmbda**2-(1/8)*(beta_1c**2+beta_1s**2)**2-0.5 *
                      (0.5*beta_0**2+(3./8.)*beta_1c**2+(1/8)*beta_1s**2)+(cd0/(4*a))*(1+mu**2)-0.5*mu*lmbda*beta_1c-mu*theta_0*beta_1s/3)
    return CQ


def fCMx(nu_b, beta_1s):
    CMx = (sigma*a/(2*gamma))*(nu_b**2-1)*beta_1s
    return CMx


def fCMy(nu_b, beta_1c):
    CMy = (sigma*a/(2*gamma))*(nu_b**2-1)*beta_1c
    return CMy


def ff1(lmbda, mu, alpha, CT):
    return lmbda - mu*tan(alpha) - CT/(2*np.sqrt(mu**2+lmbda**2))

# Tail rotor trim function


def ff_trim_tr(x, CT_tr, vtip_tr, psi_tr, v_inf):
    alpha, theta_0 = x
    # Pre-processing of inputs
    mu = v_inf*cos(alpha+psi_tr)/vtip_tr
    def f1(lmbda): return ff1(lmbda, mu, alpha+psi_tr, CT_tr)
    sol = fsolve(f1, 0.1)
    lmbda = sol[0]

    # System of equations
    f = np.zeros(2)
    f[0] = (1+mu**2)*theta_0+((4/5)+(2/3)*mu**2)*theta_tw_tr - (4./3.)*lmbda
    f[1] = (8/3)*mu*theta_0+2*mu*theta_tw-2*mu*lmbda
    return f


def ffphi(x, T, W, Lht, H, alpha, psi_tr):
    phi, Ttr = x
    Y = -Ttr
    Yf = W*phi
    f = np.zeros(2)
    f[0] = T*phi + Y + Yf + Ttr - Lht*phi
    f[1] = T*cos(alpha) - Y*phi + H*sin(alpha) - \
        Lht - W - Ttr*(phi-psi_tr)
    return f


def ff2(x, xcg, ycg, hcg, xht, xtr, htr, Lht, Df, W, alpha, lmbda, mu, CT, CQ, CH, phi, nu_b):
    theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, Mxf, Myf, Mzf = x
    # Make assumption that phi is zero
    Yf = W*phi
    # Define helpful relationships using primary unknowns
    CQ_rhs = fCQ(theta_0, beta_0, beta_1c, beta_1s, mu, lmbda)
    CH_rhs = fCH(theta_0, theta_1c, theta_1s,
                 beta_0, beta_1c, beta_1s, mu, lmbda)
    CT_rhs = fCT(theta_0, theta_1s, mu, lmbda)

    # Determine tail rotor thrust
    CQ_tr = 0.05*CQ*A*R*vtip**2/(A_tr*R_tr*vtip_tr**2)
    CT_tr = ((np.sqrt(2)/kappa)*(CQ_tr-sigma_tr*cd0/8))**(2./3.)
    Ttr = CT_tr*rho*A_tr*vtip_tr**2

    CMx = fCMx(nu_b, beta_1s)
    CMy = fCMy(nu_b, beta_1c)
    Mxr = CMx*(rho*A*R*vtip**2)
    Myr = CMy*(rho*A*R*vtip**2)
    Mzr = -CQ*(rho*A*R*vtip**2)

    f2 = np.zeros(9)
    # RHS vector in the equations
    RHS = np.zeros(3)
    MAT = np.zeros((3, 3))
    RHS[0] = (1+mu**2)*theta_0 + ((4/5)+(2/3)*mu**2) * \
        theta_tw + (4/3)*mu*theta_1s - (4/3)*lmbda
    RHS[1] = (1+(mu**2)/2)*theta_1c
    RHS[2] = (8/3)*mu*theta_0+2*mu*theta_tw + \
        (1+(3/2)*mu**2)*theta_1s-2*mu*lmbda

    MAT[0, 0] = (8/gamma)*nu_b**2
    MAT[1, 0] = 4*mu/3
    MAT[1, 1] = 8*(nu_b**2-1)/gamma
    MAT[1, 2] = 1+0.5*mu**2
    MAT[2, 1] = -(1-0.5*mu**2)
    MAT[2, 2] = 8*(nu_b**2-1)/gamma

    # Assemble equations
    f2[0] = Mxr + Mxf - W*ycg + W*phi*hcg + \
        Yf*hcg + Yf*phi*ycg + Ttr*(hcg - htr)
    f2[1] = Myr + Myf - W*cos(alpha)*xcg + W*sin(alpha)*hcg - \
        Df*cos(alpha)*hcg - Df*sin(alpha)*xcg - Lht*(xcg-xht)
    f2[2] = Mzr + Mzf + Ttr*(xtr-xcg)-Df*cos(alpha)*ycg - Yf*xcg
    f2[3] = CT - CT_rhs
    f2[4] = CQ - CQ_rhs
    f2[5] = CH - CH_rhs
    f2[6] = MAT[0, 0]*beta_0+MAT[0, 1]*beta_1c+MAT[0, 2]*beta_1s - RHS[0]
    f2[7] = MAT[1, 0]*beta_0+MAT[1, 1]*beta_1c+MAT[1, 2]*beta_1s - RHS[1]
    f2[8] = MAT[2, 0]*beta_0+MAT[2, 1]*beta_1c+MAT[2, 2]*beta_1s - RHS[2]
    return f2


def ffphi(x, T, W, Y, Lht, H, alpha, psi_tr):
    phi, Ttr = x
    Yf = W*phi
    f = np.zeros(2)
    f[0] = T*phi + Y + Yf + Ttr - Lht*phi
    f[1] = T*cos(alpha) - Y*phi + H*sin(alpha) - \
        Lht - W - Ttr*(phi-psi_tr)
    return f

# Main function for solving the helicopter trim


def solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp):
    # Estimate drag on airframe
    CD_sphere = 0.5
    Df = CD_sphere*A_copter*0.5*rho*v_inf**2

    # Assumption about thrust and power of helicopter and tail rotor
    T = W
    CT = T/(rho*A*vtip**2)
    CQ = sigma*cd0/8+(kappa/np.sqrt(2))*CT**(3./2.)
    Mzr = -CQ*(rho*A*R*vtip**2)
    P_mr = CQ*rho*A*R*vtip**2

    # Define rotor drag CH
    H = Nb*cd0*R*(rho*A_blade*v_inf**2)
    CH = H/(rho*A*vtip**2)

    # Consumes 5% of main rotor power
    # CQ_tr = 0.05*CQ*A*R*vtip**2/(A_tr*R_tr*vtip_tr**2)
    # CT_tr = ((np.sqrt(2)/kappa_tr)*(CQ_tr-sigma_tr*cd0/8))**(2./3.)
    # Ttr = CT_tr*rho*A_tr*vtip_tr**2
    # P_tr = CQ_tr*rho*A_tr*R_tr*vtip_tr**2

    Y = 100  # initial guess
    for num_iter in range(20):
        Y_last = Y

        # Step 1)
        # A) Solve for alpha, L_emp, Lht
        alpha = (Df+H)/T
        Lht = -2*np.pi*(psi_emp-alpha)*L_emp*(rho*v_inf**2)*L_emp*c_emp

        # B) Solve mu and lmbda
        mu = v_inf*cos(alpha)/vtip
        def f1(lmbda): return ff1(lmbda, mu, alpha, CT)
        sol = fsolve(f1, 0.1)
        lmbda = sol[0]

        # Step 2) Solve for phi, Ttr
        def fphi(x): return ffphi(x, T, W, Y, Lht, H, alpha, psi_tr)
        sol = fsolve(fphi, np.array([0, 1]))
        phi, Ttr = sol
        CT_tr = Ttr/(rho*A_tr*vtip_tr**2)
        Yf = W*phi

        # Step 3) Solve everything else
        def f2(x): return ff2(x, xcg, ycg, hcg, xht, xtr, htr,
                              Lht, Df, W, alpha, lmbda, mu, CT, CQ, CH, phi, nu_b)
        x0 = np.zeros(9)
        x0[0] = 1
        x0[6] = 0
        x0[7] = x0[8] = 100
        root = fsolve(f2, x0)  # beta_1c has the hardest time converging
        theta_0, theta_1c, theta_1s, beta_0, beta_1c, beta_1s, Mxf, Myf, Mzf = root
        if not any(np.isclose(f2(root), np.zeros(9), atol=1e-4)):
            print(f"CHECK CONVERGENCE!")

        CY = fCY(theta_0, theta_1c, theta_1s,
                 beta_0, beta_1c, beta_1s, mu, lmbda)
        Y = CY*rho*A*vtip**2
        if np.abs(Y-Y_last) < 1e-9:
            break

    # Step 5) Solve tail rotor trim
    def f_tr(x): return ff_trim_tr(x, CT_tr, vtip_tr, psi_tr, v_inf)
    soln = fsolve(f_tr, 0.1*np.ones(2))
    alpha_tr, theta_0_tr = soln

    # Create solution output form
    CP_tr = sigma_tr*cd0/8 + (kappa_tr/np.sqrt(2))*CT_tr**(3./2.)
    P_tr = CP_tr*rho*A_tr*R_tr*vtip_tr**2
    sol_dict = {}
    sol_dict['CT'] = CT
    sol_dict['CT_tr'] = CT_tr
    sol_dict['P_tr'] = P_tr*0.0018433993923729736  # HP
    sol_dict['P_mr'] = P_mr*0.0018433993923729736  # HP
    sol_dict['P'] = (P_mr+P_tr)*0.0018433993923729736  # HP
    sol_dict['CH'] = CH
    sol_dict['CY'] = CY
    sol_dict['Yf'] = Yf
    sol_dict['Lht'] = Lht
    sol_dict['mu'] = mu
    sol_dict['lmbda'] = lmbda
    sol_dict['alpha'] = np.rad2deg(alpha)
    sol_dict["phi"] = np.rad2deg(phi)
    sol_dict['theta_0'] = np.rad2deg(theta_0)
    sol_dict['theta_1c'] = np.rad2deg(theta_1c)
    sol_dict['theta_1s'] = np.rad2deg(theta_1s)
    sol_dict['beta_0'] = np.rad2deg(beta_0)
    sol_dict['beta_1c'] = np.rad2deg(beta_1c)
    sol_dict['beta_1s'] = np.rad2deg(beta_1s)
    sol_dict['theta_0_tr'] = np.rad2deg(theta_0_tr)
    sol_dict['alpha_tr'] = np.rad2deg(alpha_tr)
    sol_dict['Mxf'] = Mxf
    sol_dict['Myf'] = Myf
    sol_dict['Mzf'] = Mzf
    sol_dict['Mxr'] = fCMx(nu_b, beta_1s)
    sol_dict['Myr'] = fCMy(nu_b, beta_1c)
    sol_dict['Mzr'] = Mzr

    return sol_dict
