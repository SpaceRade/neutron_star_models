#!/usr/bin/env python3
"""
TOV Neutron Star Model - EOS functions added
"""

import numpy as np

# Physical constants
M_N = 1.675e-27
C = 2.998e8
HBAR = 1.055e-34
EPSILON_0 = (M_N**4 * C**5) / (np.pi**2 * HBAR**3)

def epsilon_x(x):
    """Energy density as a function of dimensionless momentum x."""
    return EPSILON_0 / 8 * ((2*x**3 + x)*np.sqrt(1+x**2) - np.arcsinh(x))

def pressure_x(x):
    """Pressure as a function of dimensionless momentum x."""
    return EPSILON_0 / 24 * ((2*x**3 - 3*x)*np.sqrt(1+x**2) + 3*np.arcsinh(x))

def run_simulation():
    pass

def main():
    run_simulation()

if __name__ == "__main__":
    main()
