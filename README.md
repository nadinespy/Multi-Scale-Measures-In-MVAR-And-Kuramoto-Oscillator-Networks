# Multi-Scale Measures of Emergence and Complexity in Autoregressive Networks and Kuramoto Oscillators

This repository contains code for calculating multi-scale measures of complexity and emergence in Multivariate Autoregressive (MVAR) networks and Kuramoto oscillators.

## Approach

We use a hybrid approach—combining numerical and analytical methods—to compute various measures across both MVAR networks and simulated Kuramoto oscillator time-series:

**Emergence measures (selection):**
- Dynamical Dependence
- Whole-Parts Emergence
- Co-Information
- Emergence Capacity

**Complexity measures:**
- Various versions of Integrated Information

The emergence measures include both Shannon information-based approaches (Dynamical Dependence, Whole-Parts Emergence, Co-Information) and measures based on Partial and Integrated Information Decomposition (Emergence Capacity).

## Current Status

**This work is ongoing and the code is not yet sufficiently documented or tested for replication and wider use.** The project is led by Nadine Spychala in collaboration with Lionel Barnett.

## Usage

### Main Scripts

**MVAR Networks:**
- Main script: `mec_mvar_rand_cross_coup_rmi.m`
- Plotting: `plot_mec_mvar_rand_cross_coup_rmi.m`

**Kuramoto Oscillators:**
- Main script: `mec_km_rand_coup_ratio_phase_lag_rmi.m`
- Plotting: `plot_mec_km_rand_coup_ratio_phase_lag_rmi.m`

Both scripts use functions from the `functions` directory.

### Requirements

- MATLAB
- [kuramoto](https://github.com/lcbarnett/kuramoto) (efficient C implementation of the standard Kuramoto coupled oscillators system)
- [kvar](https://github.com/lcbarnett/kvar) (Vector autoregressive (VAR) analysis for Kuramoto oscillator networks)
- Extended and heavily adapted version of [PhiID](https://github.com/pmediano/PhiID) (Integrated Information Decomposition)
- [MVGC](https://github.com/lcbarnett/MVGC2) (Multivariate Granger Causality MATLAB software suite)
- [ssdi](https://github.com/lcbarnett/ssdi) (An implementation of Dynamical Independence computation and optimisation for linear state-space systems)

### Setup

You'll need to configure local-specific directories in:
- `mec_mvar_rand_cross_coup_rmi_setup.m` (for MVAR networks)
- `mec_km_rand_coup_ratio_phase_lag_rmi_setup.m` (for Kuramoto oscillators)

## Publication

A publication titled "Making Sense of Measures of Emergence and Complexity in Autoregressive Networks and Kuramoto Oscillators" is in preparation. Corresponding updates to this repository will follow.

