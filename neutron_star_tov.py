#!/usr/bin/env python3
"""
TOV Neutron Star Model - Run simulation and save data
"""

import numpy as np
from scipy.integrate import solve_ivp

G = 6.67430e-11
C = 2.998e8
M_N = 1.675e-27
HBAR = 1.055e-34
EPSILON_0 = (M_N**4 * C**5)/(np.pi**2 * HBAR**3)
SURFACE_PRESSURE = 1e10
TOL = 1e-5
M_SUN = 1.989e30

def epsilon_x(x):
    return EPSILON_0 / 8 * ((2*x**3 + x)*np.sqrt(1+x**2) - np.arcsinh(x))

def pressure_x(x):
    return EPSILON_0 / 24 * ((2*x**3 - 3*x)*np.sqrt(1+x**2) + 3*np.arcsinh(x))

def bisection_method(p):
    return 1e-5  # simplified placeholder

def tov_equations(r, y):
    p, m = y
    if p <= SURFACE_PRESSURE:
        return [0,0]
    x_n = bisection_method(p)
    eps = epsilon_x(x_n)
    dpdr = - (G * eps * m / (C**2 * r**2)) * (1 + p/eps)
    dmdr = 4 * np.pi * r**2 * eps / C**2
    return [dpdr, dmdr]

def surface_event(r, y):
    return y[0] - SURFACE_PRESSURE
surface_event.terminal = True

def run_simulation():
    central_pressures = np.logspace(28,40,50)
    final_masses, final_radii = [], []
    for p_central in central_pressures:
        sol = solve_ivp(
            tov_equations,
            (1e-6,2e5),
            [p_central,1e-6],
            method='RK45',
            events=surface_event,
            rtol=TOL,
            atol=TOL
        )
        final_masses.append(sol.y[1][-1]/M_SUN)
        final_radii.append(sol.t[-1]/1e3)
    with open('neutron_star_data.txt','w') as f:
        f.write("Pressure\tMass\tRadius\n")
        for p,m,r in zip(central_pressures,final_masses,final_radii):
            f.write(f"{p}\t{m}\t{r}\n")
    print("Data saved.")

def main():
    run_simulation()

if __name__ == "__main__":
    main()
