# Neutron Star Stability Models

This project models neutron stars under two physical frameworks to determine their **mass–radius relationship** and the **maximum stable configuration** before collapse:

1. **Newtonian Model with Non-Relativistic Polytropic EOS**  
2. **General Relativistic Model using TOV Equations with a Dynamic EOS**  

By varying the central pressure over many orders of magnitude, each model predicts the star’s radius and total mass, enabling stability comparisons between classical and relativistic treatments.

---

## Features
- ✅ **Newtonian Polytrope Model** — baseline stability estimates using classical gravity and a non-relativistic polytropic equation of state.
- ✅ **TOV (General Relativity) Model** — more realistic calculation of neutron star structure using the Tolman–Oppenheimer–Volkoff equations.
- ✅ **Dynamic EOS Support** in the TOV model for more physically accurate results.
- ✅ **Automated Parameter Sweep** — evaluates a range of central pressures in each model.
- ✅ **Visualization** — generates mass–radius and radius–pressure plots for comparison.

---

## Requirements
- **Python 3.8+**
- **Libraries:**
  - `numpy`
  - `scipy`
  - `matplotlib`

Install dependencies:
```bash
pip install numpy scipy matplotlib

