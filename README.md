![Build](https://github.com/username/repo/actions/workflows/main.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

# FeynPy

FeynPy is a Python library for symbolic manipulation of Feynman amplitudes in Quantum Field Theory. It provides tools for spinors, polarization vectors, gamma matrices, and momenta, along with example notebooks for QED, Standard Model, and EFT processes.

---

## Project Structure

```
.
├── feynpy                  # Core library
│   │
│   ├── __init__.py
│   ├── gamma_matrices.py   # Gamma matrices and related operations
│   ├── misc.py             # Generic helper functions
│   ├── momenta.py          # Momentum handling and manipulations
│   ├── polarization_vectors.py # Polarization vectors (transverse and longitudinal)
│   └── spinors.py          # Fermion spinors
│
├── processes               # Example notebooks for physical processes
│   │
│   ├── EFT
│   │   └── HHZZ.ipynb      
│   ├── QED
│   │   ├── e+e->aa.ipynb
│   │   └── e+e->mu+mu-.ipynb
│   └── SM
│       ├── ee~->zz.ipynb
│       └── Untitled-1.nb
│
└── pyproject.toml          # Project configuration and dependencies
```


## Installation

It is recommended to use a Python 3.11+ virtual environment.

```
git clone https://github.com/Francesco-Montagno/feynpy.git
cd feynpy
pip install -e .
```

This installs the library in editable mode, allowing local modifications.

---

## Basic Usage
In this section, we outline the main functions and definitions provided by FeynPy.
### Functions and Definitions
```python
from feynpy import * 
#----------- Momenta -----------#
MomUP(En, p, θ, φ) # Returns a 4-momentum with upper indices (En, p*sp.sin(θ)*sp.cos(φ), p*sp.sin(θ)*sp.sin(φ), p*sp.cos(θ))
MomDOWN(En, p, θ, φ) # Returns a 4-momentum with lower indices (En, -p*sp.sin(θ)*sp.cos(φ), -p*sp.sin(θ)*sp.sin(φ), -p*sp.cos(θ))

#------------- Misc ------------#
metric(mu, nu) # Returns the metric tensor g_{mu nu} (most-minus convention, g = diag(1, -1, -1, -1))
lower_index(p) # Returns the 4-vector p with lowered index
upper_index(p) # Returns the 4-vector p with upper index
minkowski_dot(a, b) # Returns the Minkowski dot product a.T * g * b

#------- Polarization Vectors -------#
EpsilonUP(λ, En, p, θ, φ) # Returns the polarization vector with upper indices for helicity λ = +1, -1, 0
 

#------- Spinors -------#
u_spinor(h, En, p, θ, φ) # Returns the Dirac spinor u with helicity h = +1/2, -1/2
v_spinor(h, En, p, θ, φ) # Returns the Dirac spinor v with helicity h = +1/2, -1/2

#------- Gamma Matrices -------#
GammaUP(mu) # Returns the gamma matrix with upper index mu
GammaDOWN(mu) # Returns the gamma matrix with lower index mu
Gamma5 # Is the gamma^5 matrix
Gamma6 # Is (1 + gamma^5)/2
Gamma7 # Is (1 - gamma^5)/2
Slash(p) # Returns the slashed momentum p
Bar(psi) # Returns psi^dagger gamma^0
```
---

### Example: $ e^+ e^- \to \gamma \gamma $ Amplitude Calculation
This Jupyter notebook snippet demonstrates how to compute the amplitude for the process \( e^+ e^- \to \gamma \gamma \) using FeynPy. The amplitude expression is obtained using the `FeynCalc` package in `Mathematica` as a reference.

```python
from feynpy import *
print("FeynPy imports OK ✅")

import sympy as sp

# Define symbols
θ, φ = sp.symbols('theta phi', real=True)
p = sp.symbols('p', real=True, positive=True)

# Momenta
p1 = MomUP(p, p, 0, 0)          # e-
p2 = MomUP(p, p, sp.pi, 0)      # e+
p3 = MomUP(p, p, θ, 0)          # a
p4 = MomUP(p, p, sp.pi+ θ, 0)   # a
display(p1, p2, p3, p4)

# Define amplitude function for e+ e- -> γ γ
def Amplitude(h1, h2, λ3, λ4):
    u1 = u_spinor(h1, p, p, 0, 0)
    v2 = v_spinor(h2, p, p, sp.pi, 0)
    eps1 = EpsilonUP(λ3, p, p, θ, 0)
    eps2 = EpsilonUP(λ4, p, p, sp.pi+ θ, 0)

    e = sp.symbols('e', real=True)

    t = minkowski_dot(p1 - p3, p1 - p3)
    u = minkowski_dot(p1 - p4, p1 - p4)

    M = sp.ZeroMatrix(1, 1)

    M = -e**2 * (
        -2 * (minkowski_dot(p2, eps1) * Bar(v2) * Slash(eps2) * Gamma6 * u1 / u)
        -2 * (minkowski_dot(p2, eps1) * Bar(v2) * Slash(eps2) * Gamma7 * u1 / u)
        -2 * (minkowski_dot(p2, eps2) * Bar(v2) * Slash(eps1) * Gamma6 * u1 / t)
        -2 * (minkowski_dot(p2, eps2) * Bar(v2) * Slash(eps1) * Gamma7 * u1 / t)
        + (Bar(v2) * Slash(eps1) * Slash(p3) * Slash(eps2) * Gamma6 * u1 / u)
        + (Bar(v2) * Slash(eps1) * Slash(p3) * Slash(eps2) * Gamma7 * u1 / u)
        + (Bar(v2) * Slash(eps2) * Slash(p4) * Slash(eps1) * Gamma6 * u1 / t)
        + (Bar(v2) * Slash(eps2) * Slash(p4) * Slash(eps1) * Gamma7 * u1 / t)
    )

    return sp.simplify(M)

# Example helicities
h1 = +1/2
h2 = -1/2
λ1 = -1
λ2 = +1

M = Amplitude(h1, h2, λ1, λ2)

# Helicity display helper
def hel_to_pm(h):
    return '+' if h == +1/2 or h == +1 else '-' if h == -1/2 or h == -1 else '0'

for h1 in [+1/2, -1/2]:
    for h2 in [+1/2, -1/2]:
        for λ1 in [+1, -1]:
            for λ2 in [+1, -1]:
                M = Amplitude(h1, h2, λ1, λ2)
                print(f"Helicities: {hel_to_pm(h1)}, {hel_to_pm(h2)}, {hel_to_pm(λ1)}, {hel_to_pm(λ2)}")
                display(M.simplify())
                print("--------------------------------\n")
```

---

## Example Notebooks
In the `processes` directory, you can find example Jupyter notebooks demonstrating various physical processes:
- **QED**
  - `e+e->aa.ipynb` 
  - `e+e->mu+mu-.ipynb` 

- **Standard Model**
  - `ee~->zz.ipynb`

- **EFT**
  - `HHZZ.ipynb` — example of the EFT corrections to $ q\bar q \to ZZ $ coming from the EFT operator $ H^\dagger H W_{\mu\nu}^i W^{\mu\nu\,, i}$.

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub. You can also send me an email at `fmontagno@ifae.es`

---

## License

Distributed under the MIT License. See LICENSE for details.

