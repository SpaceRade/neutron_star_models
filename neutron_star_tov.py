#!/usr/bin/env python3
"""
TOV Neutron Star Model - Added bisection method for EOS inversion
"""

import numpy as np

M_N = 1.675e-27
C = 2.998e8
HBAR = 1.055e-34
EPSILON_0 = (M_N**4 * C**5) / (np.pi**2 * HBAR**3)
TOL = 1e-5
MIN_X = 1e-5
X_LOW = 1e-8
X_HIGH = 1e30
MAX_ITER = 1000

def epsilon_x(x):
    return EPSILON_0 / 8 * ((2*x**3 + x)*np.sqrt(1+x**2) - np.arcsinh(x))

def pressure_x(x):
    return EPSILON_0 / 24 * ((2*x**3 - 3*x)*np.sqrt(1+x**2) + 3*np.arcsinh(x))

def bisection_method(p, tol=TOL, max_iter=MAX_ITER):
    """Solve for x given pressure p."""
    if p <= 0: return MIN_X
    x_low, x_high = X_LOW, X_HIGH
    p_low, p_high = pressure_x(x_low), pressure_x(x_high)
    for _ in range(max_iter):
        x_mid = (x_low + x_high)/2
        p_mid = pressure_x(x_mid)
        if abs(p_mid - p) < tol:
            return x_mid
        if (p_mid - p)*(p_low - p) < 0:
            x_high = x_mid
        else:
            x_low, p_low = x_mid, p_mid
    return x_mid

def run_simulation():
    pass

def main():
    run_simulation()

if __name__ == "__main__":
    main()
