import numpy as np
import matplotlib.pyplot as plt

from proj2 import sigma, plot_state_var, calc_trims_array, knots2mps

def vary_flight_velocity(state_to_plot, v_inf_list, CX=0.05, nu_b=1):
    plt.figure()
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    # (v_inf must be in m/s)
    # Typical values of CT and CX when performing the trim
    N_CT = 50
    CT_norm_list = np.linspace(0.005, 0.01, N_CT)
    CT_list = CT_norm_list*sigma

    CX_fake_list = np.array([CX])

    for v_inf in v_inf_list:
        trims = calc_trims_array(CX_fake_list, CT_list, v_inf, nu_b)
        trim_tuple = (CX_fake_list, CT_list, v_inf, nu_b,
                    trims)  # info relevant to trim
        plot_state_var(trim_tuple, state_to_plot, label_type="v_inf")
    return

def vary_flap_frequency(state_to_plot, nu_b_list,v_inf=0.0, CX=0.05):
    plt.figure()
    plt.rcParams['figure.figsize'] = [6.4, 4.8]   # default
    # (v_inf must be in m/s)
    # Typical values of CT and CX when performing the trim
    N_CT = 50
    CT_norm_list = np.linspace(0.005, 0.01, N_CT)
    CT_list = CT_norm_list*sigma
    CX_fake_list = np.array([CX])

    for nu_b in nu_b_list:
        trims = calc_trims_array(CX_fake_list, CT_list, v_inf, nu_b)
        trim_tuple = (CX_fake_list, CT_list, v_inf, nu_b,
                    trims)  # info relevant to trim

        plot_state_var(trim_tuple, state_to_plot, label_type="nu_b")
    return
