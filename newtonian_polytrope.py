#!/usr/bin/env python3
"""
Newtonian Neutron Star Model with Non-Relativistic Polytrope EOS

This script models neutron stars using Newtonian gravity and a
non-relativistic polytropic equation of state.

It integrates the mass and pressure profiles using the Runge–Kutta method
and plots the final radius and mass vs central pressure.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as con

# --- Physical constants (SI units) ---
G = 6.67e-11        # gravitational constant
M_SUN = 1.989e30    # solar mass

# --- EOS and polytrope parameters ---
K_NR = 2.98e-25     # polytropic constant
A_POLY = 5 / 3      # polytropic index

# --- Numerical settings ---
PRESSURE_DROP_FACTOR = 1e-3
R_MAX = 1.0e6       # m
STEP_SIZE = 10      # m
R_START = 1         # m


def f1(r, P, M):
    """
    First ODE: dP/dr for the Newtonian polytrope.
    """
    return -(G * M * P**(1 / A_POLY)) / \
           (r**2 * K_NR**(1 / A_POLY) * con.c**2)


def f2(r, P, M):
    """
    Second ODE: dM/dr for the Newtonian polytrope.
    """
    rho = (P / K_NR)**(1 / A_POLY) / (con.c**2)
    return 4 * con.pi * r**2 * rho


def runge_kutta(P0, M0, r0, r_max, h):
    """
    Integrate the Newtonian polytrope equations using 4th-order Runge–Kutta.

    Parameters
    ----------
    P0 : float
        Central pressure (Pa).
    M0 : float
        Initial mass at r0 (kg).
    r0 : float
        Starting radius (m).
    r_max : float
        Maximum integration radius (m).
    h : float
        Step size (m).

    Returns
    -------
    tuple
        (r_values, P_values, M_values) arrays up to surface.
    """
    r_values = np.arange(r0, r_max, h)
    P_values = np.zeros_like(r_values)
    M_values = np.zeros_like(r_values)

    P_values[0] = P0
    M_values[0] = M0

    for i in range(1, len(r_values)):
        r = r_values[i - 1]
        P = P_values[i - 1]
        M = M_values[i - 1]

        k1_P = h * f1(r, P, M)
        k1_M = h * f2(r, P, M)

        k2_P = h * f1(r + h / 2, P + k1_P / 2, M + k1_M / 2)
        k2_M = h * f2(r + h / 2, P + k1_P / 2, M + k1_M / 2)

        k3_P = h * f1(r + h / 2, P + k2_P / 2, M + k2_M / 2)
        k3_M = h * f2(r + h / 2, P + k2_P / 2, M + k2_M / 2)

        k4_P = h * f1(r + h, P + k3_P, M + k3_M)
        k4_M = h * f2(r + h, P + k3_P, M + k3_M)

        dP = (k1_P + 2 * k2_P + 2 * k3_P + k4_P) / 6
        dM = (k1_M + 2 * k2_M + 2 * k3_M + k4_M) / 6

        P_values[i] = P + dP
        M_values[i] = M + dM

        # Stop if pressure falls below threshold or numerical failure
        if (P_values[i] < PRESSURE_DROP_FACTOR * P0 or
                np.isnan(P_values[i]) or
                np.isnan(M_values[i])):
            return r_values[:i], P_values[:i], M_values[:i]

    return r_values, P_values, M_values


def run_simulation():
    """
    Run the Newtonian neutron star simulation for a range of central pressures.

    Returns
    -------
    tuple
        (central_pressures, final_radii, final_masses)
    """
    central_pressures = np.logspace(29, 41, 50)
    final_radii = []
    final_masses = []

    for P0 in central_pressures:
        rho0 = (P0 / K_NR)**(1 / A_POLY) / (con.c**2)
        M0 = (4 / 3) * con.pi * R_START**3 * rho0

        r_vals, _, M_vals = runge_kutta(
            P0, M0=M0, r0=R_START, r_max=R_MAX, h=STEP_SIZE
        )

        final_radii.append(r_vals[-1] / 1000)   # km
        final_masses.append(M_vals[-1] / M_SUN)  # solar masses

    return central_pressures, final_radii, final_masses


def plot_results(central_pressures, final_radii, final_masses):
    """
    Plot neutron star mass and radius vs central pressure.
    """
    fig, ax1 = plt.subplots()
    ax1.set_xscale("log")
    ax1.plot(
        central_pressures, final_radii,
        color='black', linestyle='--', label='Radius'
    )
    ax1.set_xlabel('Central Pressure (Pa)')
    ax1.set_ylabel('Radius (km)', color='black')
    ax1.tick_params(axis='y', labelcolor='black')

    ax2 = ax1.twinx()
    ax2.plot(
        central_pressures, final_masses,
        color='black', linestyle='-', label='Mass'
    )
    ax2.set_ylabel('Mass (Solar Masses)', color='black')
    ax2.tick_params(axis='y', labelcolor='black')

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(
        lines + lines2, labels + labels2,
        loc='lower right', bbox_to_anchor=(1, 0.1)
    )

    fig.tight_layout()
    plt.show()


def main():
    """
    Main entry point: run simulation and plot results.
    """
    cp, radii, masses = run_simulation()
    plot_results(cp, radii, masses)


if __name__ == "__main__":
    main()
