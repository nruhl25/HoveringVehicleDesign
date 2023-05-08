import numpy as np
import matplotlib.pyplot as plt

from proj3_no_flap import solve_trim_system

#### INPUTS TO SOLVER PROCESS #####
v_inf = 100  # ft/sec
nu_b = 1

xcg = 0
ycg = 0
hcg = -3
xht = 15
xtr = 45
htr = 3
psi_tr = 0  # tail rotor angle
psi_emp = 0  # empenage angle

psi_tr_list = np.deg2rad(np.linspace(-5,5,20))
test_list = []
test_key = 'alpha'
for i in range(len(psi_tr_list)):
    sol_dict = solve_trim_system(xcg, ycg, hcg, xht, xtr, htr, psi_tr_list[i], psi_emp, v_inf, nu_b)
    test_list.append(sol_dict[test_key])

plt.plot(np.rad2deg(psi_tr_list), test_list)
plt.ylabel(f"{test_key}")
plt.xlabel(r"$\psi_{tr}$")
plt.show()
