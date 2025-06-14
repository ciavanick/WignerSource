# WignerSource

**WignerSource** is a C++ module for computing and handling Wigner functions using the ROOT framework. It provides utilities for simulating quantum mechanical systems and visualizing phase space distributions.

---
## Wigner Function Theory and Implementation

This software models quantum coalescence using the Wigner function, which provides a phase-space distribution consistent with quantum mechanics, capturing both position and momentum information.

### Wigner Function from a Gaussian Source

The source is modeled as a 3D Gaussian wave packet modulated by a plane wave:

    φ(r) = A * exp( -r² / (8 * R₀²) ) * exp( i * k·r / ħ )

The corresponding Wigner function is:

    W(r, p) = (1 / (2πħ)³) * exp( -r² / (4 * R₀²) ) * exp( -4 * R₀² * (k + p/ħ)² )

This function is Gaussian in both r and p, centered around the classical values.

To obtain the radial momentum distribution, the function is integrated over angles assuming the momentum vector k is aligned along the z-axis:

    W(r, p) ∝ exp( -r² / (4 * R₀²) - 4 * R₀² * (k² + p²) / ħ² )
             * sinh( 8 * R₀² * k * p / ħ² ) / (8 * R₀² * k * p / ħ²)

---

### Energy Terms and Wigner-Weighted Integrals

This software evaluates key observables by integrating the Wigner source function weighted by classical energy terms.

#### Kinetic Energy

    T(p) = p² / (2 * μ)      where μ = pm[3]

The kinetic contribution is:

    W_K = ∬ T(p) * W_source(r, p) * J(r, p) dr dp

#### Potential Energy (Square-Well)

The potential is modeled as a square well of depth V₀ and radius R = pm[4]:

    V(r) = -V₀   if r < R
           0     otherwise

The potential energy contribution is:

    W_V = ∬ V(r) * W_source(r, p) * J(r, p) dr dp

#### Total Hamiltonian

The total energy is:

    H(r, p) = T(p) + V(r)

The Wigner-weighted total energy is:

    W_H = ∬ H(r, p) * W_source(r, p) * J(r, p) dr dp

These scalar integrals quantify the energy distribution of the quantum source in phase space.

---

### Deuteron Coalescence Probability

The coalescence probability represents the likelihood that a proton-neutron pair forms a deuteron. It is computed as the overlap of the source and deuteron Wigner functions:

    P_coal = ∬ W_source(r, p) * W_deuteron(r, p) * J(r, p) dr dp

Where:
- W_source(r, p) is the Gaussian source Wigner function
- W_deuteron(r, p) is provided via a ROOT file containing a precomputed 2D histogram from numerical integration of the deuteron wave function
- J(r, p) is the full spherical Jacobian including angular corrections:

    J(r, p) = (r * p)² * (4π)² * [ (1 - exp(-2α)) / (2α) ]

With:

    α = (8 * R₀² * k * p) / ħ²

This formulation incorporates both classical volume factors and quantum angular integration.

## Requirements

- [ROOT Framework](https://root.cern/)
  - Ensure `root-config` is available in your environment
- C++17 or higher
- cmake

---

## Compilation and Installation

Clone the repository:

```bash
git clone https://github.com/ciavanick/WignerSource.git
cd WignerSource
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=..
cmake --build .
cmake --install .
cd ..
```

## How to Run

### Workflow

The software is structured so that each source file is individually compiled and linked into shared libraries. These libraries, along with the ROOT dictionary, are then wrapped into a custom executable called `wigneroot`. This executable acts as a ROOT interpreter with all project-specific classes preloaded, allowing users to run macros directly without manual library loading or setup.

### Usage

```bash
source wignerenv.sh
```

```bash
wigneroot -l wignertest2.cpp
```
## Docker

```bash
docker build -t wignerutils .
docker run --rm -it wignerutils 
```

If you have an arm achitecture:

```bash
docker build --platform linux/amd64 -t wignerutils .
docker run --platform linux/amd64 --rm -it wignerutils 
```
