import numpy as np
import matplotlib.pyplot as plt

from proj3_flap import solve_trim_system

#### INPUTS TO SOLVER PROCESS #####
v_inf = 60  # ft/sec
nu_b = 1
L_emp = 2 # ft
xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = 0  # tail rotor angle
psi_emp = 0  # empenage angle

key_list_fig1 = ["theta_0", "theta_1c", "theta_1s", "beta_0", "beta_1c", "beta_1s"]
key_label_list_fig1 = [r"Main Rotor Collective, $\theta_0$ (deg)", r"Longitudinal Cyclic Input, $\theta_{1c}$ (deg)", r"Lateral Cyclic Input, $\theta_{1s}$ (deg)",
                       r"Coning Angle, $\beta_0$ (deg)", r"Longitudinal Cyclic Flapping $\beta_{1c}$ (deg)", r"Lateral Cyclic Flapping, $\beta_{1s}$ (deg)"]

key_list_fig2 = ["theta_0_tr", "alpha", "phi", "P"]
key_label_list_fig2 = [r"$\theta_{0_{tr}}$", r"$\alpha$ (deg)", r"$\phi$ (deg)", r"$P$ (Hp)"]
