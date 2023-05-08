import matplotlib.pyplot as plt
import numpy as np

from proj3_flap import solve_trim_system
# Import default configuration
from sensitivity_analysis import v_inf, nu_b, xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, key_list_fig1, key_list_fig2, key_label_list_fig1, key_label_list_fig2

def vary_xcg():
    xcg_list = np.arange(-4, 7, 0.1)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(xcg_list))
        for i in range(len(xcg_list)):
            sol_dict = solve_trim_system(
                xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 3, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(xcg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$x_{cg}$ (ft)")
    plt.suptitle(
        "Control Inputs and Blade Flapping When Varying Longitudinal Center of Gravity")
    plt.tight_layout()

    plt.figure()
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig2, key_label_list_fig2)):
        test_list = np.zeros(len(xcg_list))
        for i in range(len(xcg_list)):
            sol_dict = solve_trim_system(
                xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(xcg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$x_{cg}$ (ft)")
    plt.suptitle(
        "Control Inputs and Vehicale State When Varying Longitudinal Center of Gravity")
    plt.tight_layout()

    plt.show()
    return
