#!/usr/bin/env python3
"""
Initial Newtonian polytrope script.
Basic ODEs and Euler integration for one central pressure.
"""
import numpy as np
import scipy.constants as con

# Physical constants
G = 6.67e-11
K_NR = 2.98e-25
A_POLY = 5/3

def f1(r, P, M):
    return -(G * M * P**(1 / A_POLY)) / (r**2 * K_NR**(1 / A_POLY) * con.c**2)

def f2(r, P, M):
    rho = (P / K_NR)**(1 / A_POLY) / con.c**2
    return 4 * np.pi * r**2 * rho

def integrate_euler(P0, M0, r0, r_max, h):
    r_values = np.arange(r0, r_max, h)
    P_values = np.zeros_like(r_values)
    M_values = np.zeros_like(r_values)
    P_values[0] = P0
    M_values[0] = M0

    for i in range(1, len(r_values)):
        P_values[i] = P_values[i-1] + h * f1(r_values[i-1], P_values[i-1], M_values[i-1])
        M_values[i] = M_values[i-1] + h * f2(r_values[i-1], P_values[i-1], M_values[i-1])

    return r_values, P_values, M_values

if __name__ == "__main__":
    r, P, M = integrate_euler(1e34, 0, 1, 1e6, 100)
    print("Final mass:", M[-1])
