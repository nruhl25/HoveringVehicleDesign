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


def vary_L_emp():
    L_emp_list = np.arange(1, 3, 0.1)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    plt.title(
        "Control Inputs and Blade Flapping When Varying Non-dimensional Flap Frequency")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(L_emp_list))
        for i in range(len(L_emp_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp_list[i])
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 3, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(L_emp_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$\nu_\beta$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.title(
        "Control Inputs and Vehicale State When Varying Empenage Size")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig2, key_label_list_fig2)):
        test_list = np.zeros(len(L_emp_list))
        for i in range(len(L_emp_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp_list[i])
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(L_emp_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$L_{emp}$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.rcParams['figure.figsize'] = [6, 4]
    P_ratio_list = []
    for i in range(len(L_emp_list)):
        sol_dict = solve_trim_system(
            xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp_list[i])
        P_ratio_list.append(sol_dict['P_tr']/sol_dict['P_mr'])

    plt.plot(L_emp_list, P_ratio_list)
    plt.ylabel(r"$P_{tr}/P_{mr}$")
    plt.xlabel(r"$L_{emp}$ (ft)")
    plt.title("Power of tail rotor vs Power of main rotor")
    plt.show()
    return
