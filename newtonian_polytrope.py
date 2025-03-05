#!/usr/bin/env python3
"""
Newtonian polytrope: simulation loop for multiple central pressures.
"""
import numpy as np
import scipy.constants as con

# Constants
G = 6.67e-11
K_NR = 2.98e-25
A_POLY = 5/3
R_START = 1
R_MAX = 1.0e6
STEP_SIZE = 10

def f1(r, P, M):
    return -(G * M * P**(1 / A_POLY)) / (r**2 * K_NR**(1 / A_POLY) * con.c**2)

def f2(r, P, M):
    rho = (P / K_NR)**(1 / A_POLY) / con.c**2
    return 4 * np.pi * r**2 * rho

def runge_kutta(P0, M0, r0, r_max, h):
    r_values = np.arange(r0, r_max, h)
    P_values = np.zeros_like(r_values)
    M_values = np.zeros_like(r_values)
    P_values[0] = P0
    M_values[0] = M0

    for i in range(1, len(r_values)):
        r, P, M = r_values[i-1], P_values[i-1], M_values[i-1]
        k1_P = h*f1(r,P,M)
        k1_M = h*f2(r,P,M)
        k2_P = h*f1(r+h/2,P+k1_P/2,M+k1_M/2)
        k2_M = h*f2(r+h/2,P+k1_P/2,M+k1_M/2)
        k3_P = h*f1(r+h/2,P+k2_P/2,M+k2_M/2)
        k3_M = h*f2(r+h/2,P+k2_P/2,M+k2_M/2)
        k4_P = h*f1(r+h,P+k3_P,M+k3_M)
        k4_M = h*f2(r+h,P+k3_P,M+k3_M)
        P_values[i] = P + (k1_P + 2*k2_P + 2*k3_P + k4_P)/6
        M_values[i] = M + (k1_M + 2*k2_M + 2*k3_M + k4_M)/6

    return r_values, P_values, M_values

def run_simulation():
    central_pressures = np.logspace(29, 41, 20)
    final_radii, final_masses = [], []

    for P0 in central_pressures:
        rho0 = (P0 / K_NR)**(1/A_POLY) / con.c**2
        M0 = (4/3) * np.pi * R_START**3 * rho0
        r_vals, _, M_vals = runge_kutta(P0, M0, R_START, R_MAX, STEP_SIZE)
        final_radii.append(r_vals[-1]/1000)
        final_masses.append(M_vals[-1]/1.989e30)

    return central_pressures, final_radii, final_masses

if __name__ == "__main__":
    cp, radii, masses = run_simulation()
    for p, r, m in zip(cp, radii, masses):
        print(f"P={p:.2e}, R={r:.2f} km, M={m:.2f} Msun")
