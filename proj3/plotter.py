import matplotlib.pyplot as plt
import numpy as np

from proj3_flap import solve_trim_system  # version that doesn't assume phi=0
# Import default configuration
from sensitivity_analysis import v_inf, nu_b, xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, L_emp, key_list_fig1, key_list_fig2, key_label_list_fig1, key_label_list_fig2

def vary_xcg():
    xcg_list = np.arange(-2.2, 6, 0.1)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(xcg_list))
        for i in range(len(xcg_list)):
            sol_dict = solve_trim_system(
                xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
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
                xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(xcg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$x_{cg}$ (ft)")
    plt.suptitle(
        "Control Inputs and Vehicale State When Varying Longitudinal Center of Gravity")
    plt.tight_layout()

    plt.figure()
    plt.rcParams['figure.figsize'] = [6, 4]
    P_ratio_list = []
    for i in range(len(xcg_list)):
        sol_dict = solve_trim_system(
            xcg_list[i], ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
        P_ratio_list.append(sol_dict['P_tr']/sol_dict['P_mr'])
    
    plt.plot(xcg_list, P_ratio_list)
    plt.ylabel(r"$P_{tr}/P_{mr}$")
    plt.xlabel(r"$x_{cg}$ (ft)")
    plt.title("Power of tail rotor vs Power of main rotor")

    plt.show()
    return


def vary_ycg():
    ycg_list = np.arange(-5, 2, 0.1)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    plt.title("Control Inputs and Blade Flapping When Varying Lateral Center of Gravity")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(ycg_list))
        for i in range(len(ycg_list)):
            sol_dict = solve_trim_system(
                xcg, ycg_list[i], hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 3, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(ycg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$y_{cg}$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.title("Control Inputs and Vehicale State When Varying Lateral Center of Gravity")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig2, key_label_list_fig2)):
        test_list = np.zeros(len(ycg_list))
        for i in range(len(ycg_list)):
            sol_dict = solve_trim_system(
                xcg, ycg_list[i], hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(ycg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$y_{cg}$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.rcParams['figure.figsize'] = [6, 4]
    P_ratio_list = []
    for i in range(len(ycg_list)):
        sol_dict = solve_trim_system(
            xcg, ycg_list[i], hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
        P_ratio_list.append(sol_dict['P_tr']/sol_dict['P_mr'])

    plt.plot(ycg_list, P_ratio_list)
    plt.ylabel(r"$P_{tr}/P_{mr}$")
    plt.xlabel(r"$y_{cg}$ (ft)")
    plt.title("Power of tail rotor vs Power of main rotor")

    plt.show()
    return


def vary_hcg():
    hcg_list = np.arange(-8, 0, 0.1)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    plt.title("Control Inputs and Blade Flapping When Vertical Center of Gravity")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(hcg_list))
        for i in range(len(hcg_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg_list[i], xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 3, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(hcg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$h_{cg}$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.title("Control Inputs and Vehicle State When Varying Vertical Center of Gravity")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig2, key_label_list_fig2)):
        test_list = np.zeros(len(hcg_list))
        for i in range(len(hcg_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg_list[i], xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(hcg_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$h_{cg}$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.rcParams['figure.figsize'] = [6, 4]
    P_ratio_list = []
    for i in range(len(hcg_list)):
        sol_dict = solve_trim_system(
            xcg, ycg, hcg_list[i], xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b, L_emp)
        P_ratio_list.append(sol_dict['P_tr']/sol_dict['P_mr'])

    plt.plot(hcg_list, P_ratio_list)
    plt.ylabel(r"$P_{tr}/P_{mr}$")
    plt.xlabel(r"$h_{cg}$ (ft)")
    plt.title("Power of tail rotor vs Power of main rotor")

    plt.show()
    return


def vary_L_emp():
    L_emp_list = np.arange(2, 10, 0.5)
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


def vary_nu_beta():
    nu_b_list = np.arange(1, 1.15, 0.01)
    plt.rcParams['figure.figsize'] = [12, 8]
    plt.figure()
    plt.title("Control Inputs and Blade Flapping When Varying Non-dimensional Flap Frequency")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig1, key_label_list_fig1)):
        test_list = np.zeros(len(nu_b_list))
        for i in range(len(nu_b_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b_list[i], L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 3, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(nu_b_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$\nu_\beta$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.title("Control Inputs and Vehicale State When Varying Non-dimensional Flap Frequency")
    for indx, (test_key, key_label) in enumerate(zip(key_list_fig2, key_label_list_fig2)):
        test_list = np.zeros(len(nu_b_list))
        for i in range(len(nu_b_list)):
            sol_dict = solve_trim_system(
                xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b_list[i], L_emp)
            test_list[i] = sol_dict[test_key]
        plt.subplot(2, 2, indx+1)
        plt.ticklabel_format(style='plain')
        plt.plot(nu_b_list, test_list)
        plt.ylabel(f"{key_label}")
        plt.xlabel(r"$\nu_\beta$ (ft)")
    plt.tight_layout()

    plt.figure()
    plt.rcParams['figure.figsize'] = [6, 4]
    P_ratio_list = []
    for i in range(len(nu_b_list)):
        sol_dict = solve_trim_system(
            xcg, ycg, hcg, xht, xtr, htr, psi_tr, psi_emp, v_inf, nu_b_list[i], L_emp)
        P_ratio_list.append(sol_dict['P_tr']/sol_dict['P_mr'])

    plt.plot(nu_b_list, P_ratio_list)
    plt.ylabel(r"$P_{tr}/P_{mr}$")
    plt.xlabel(r"$\nu_\beta$ (ft)")
    plt.title("Power of tail rotor vs Power of main rotor")

    plt.show()
    return
