import numpy as np
import matplotlib.pyplot as plt

from proj2 import sigma, plot_state_var, calc_trims_array, knots2mps, calc_alpha, mps2knots

def vary_flight_velocity(state_to_plot, v_inf_list, CX=0.05, nu_b=1):
    # (v_inf must be in m/s)
    # Typical values of CT and CX when performing the trim
    N_CT = 50
    CT_norm_list = np.linspace(0.005, 0.01, N_CT)
    CT_list = CT_norm_list*sigma

    CX_fake_list = np.array([CX])

    for indx, v_inf in enumerate(v_inf_list):
        trims = calc_trims_array(CX_fake_list, CT_list, v_inf, nu_b)
        trim_tuple = (CX_fake_list, CT_list, v_inf, nu_b,
                    trims)  # info relevant to trim
        plot_state_var(trim_tuple, state_to_plot, label_type="v_inf")
    return

def vary_flap_frequency(state_to_plot, nu_b_list,v_inf=0.0, CX=0.05):
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

# I like the below workflow for plotting subfigures. Keep this in mind in the future.
def plot_trim_vary_velocity(v_inf_list, CX=0.05, nu_b=1.0):
    plt.rcParams['figure.figsize'] = [18, 12]
    plt.figure()
    state_var_list = ["beta_0", "theta_0",
                        "theta_1c", "theta_1s", "lambda", "TP"]
    for indx, state_var in enumerate(state_var_list):
        plt.subplot(2, 3, indx+1)
        vary_flight_velocity(state_var, v_inf_list, CX, nu_b)
    plt.tight_layout()
    plt.show()
    return


def plot_trim_vary_flap_frequency(nu_b_list, v_inf=knots2mps(25), CX=0.05):
    plt.rcParams['figure.figsize'] = [18, 12]
    plt.figure()
    state_var_list = ["beta_0", "theta_0",
                      "theta_1c", "theta_1s", "lambda", "TP"]
    for indx, state_var in enumerate(state_var_list):
        plt.subplot(2, 3, indx+1)
        vary_flap_frequency(state_var, nu_b_list, v_inf, CX)
    plt.tight_layout()
    plt.show()
    return


def plot_rotor_disc_tilt(CX_list, CT_list, v_inf, nu_b):
    plt.rcParams['figure.figsize'] = [6.4, 4.8]
    for i in range(len(CX_list)):
        plt.figure(1)
        alpha = calc_alpha(CX_list[i], CT_list, v_inf, nu_b)
        plt.plot(CT_list/sigma, np.rad2deg(alpha), label=fr"$C_X$={CX_list[i]:.3f}")
    plt.grid()
    plt.title(fr"Required rotor rotor shaft incidence angle to achieve desired forward thrust ($v_\infty$={mps2knots(v_inf):.2f}, $\nu_b=${nu_b:.3f})")
    plt.xlabel(r"$C_T/\sigma$")
    plt.ylabel(r"Rotor Shaft Incidence Angle, $\alpha$ (deg)")
    plt.legend()
    return
