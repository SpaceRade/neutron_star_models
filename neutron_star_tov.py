#!/usr/bin/env python3
"""
TOV Neutron Star Model - TOV equations added
"""

import numpy as np

G = 6.67430e-11
C = 2.998e8
M_N = 1.675e-27
HBAR = 1.055e-34
EPSILON_0 = (M_N**4 * C**5) / (np.pi**2 * HBAR**3)
SURFACE_PRESSURE = 1e10
TOL = 1e-5

def epsilon_x(x):
    return EPSILON_0 / 8 * ((2*x**3 + x)*np.sqrt(1+x**2) - np.arcsinh(x))

def pressure_x(x):
    return EPSILON_0 / 24 * ((2*x**3 - 3*x)*np.sqrt(1+x**2) + 3*np.arcsinh(x))

def bisection_method(p):
    return 1e-5  # placeholder

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
    pass

def main():
    run_simulation()

if __name__ == "__main__":
    main()
