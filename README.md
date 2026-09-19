# NSBFS_ETDRK-G

This MATLAB package provides a two-dimensional implementation of the third-order exponential time differencing Runge--Kutta gauge (ETDRK3-G) scheme for the hydrodynamics-coupled binary fluid--surfactant (BFS) phase-field model.

The included example simulates two surfactant-laden droplets in simple shear flow on a rectangular domain using a MAC staggered grid and FFT/DCT/DST-based fast solvers.

## Associated manuscript

This code accompanies the manuscript:

**Mass-Conserving Exponential Time Differencing Gauge Method for the Hydrodynamics-Coupled Binary Fluid-Surfactant Phase-Field Model**  
Jiayi Duan, Lili Ju, Xiao Li, and Zhonghua Qiao

The manuscript has not yet been published.

## Main files

- `run.m` — sets the basic numerical parameters and runs the example.
- `NSBFS_ETDRK3G.m` — main ETDRK3-G solver.
- `draw_figure.m` — plots the energy, extrema of the phase variables, and saved phase-field snapshots.
- `F1.m`, `F2.m`, `F3u.m`, `F3v.m` — nonlinear terms in the phase-field and momentum equations.
- `mu_phi.m`, `mu_rho.m`, `get_energy.m` — chemical potentials and discrete energy.
- `varphi*.m` — functions used in the ETDRK3 time discretization.
- `dst*byfft.m`, `idst*byfft.m`, `Leigenv2d_*.m` — fast-transform and eigenvalue routines.
- The remaining `grad*`, `aver*`, and `Laplace*` files provide finite-difference operators on the staggered grid.

## Requirements

- MATLAB with `fft`, `dct`, and `idct` available.

## How to run

1. Place all `.m` files in the same folder and set this folder as the MATLAB working directory.
2. Run

```matlab
run
```

The default parameters in `run.m` are

```matlab
N   = 128;
T   = 0.9;
dt  = 0.001;
dtn = 0.1;
```

Here, the spatial grid is `Nx = 2*N` and `Ny = N`; `T` is the final time, `dt` is the time-step size, and `dtn` is the interval used to save snapshots.

The main solver can also be called directly as

```matlab
[time,u,v,phi,rho] = NSBFS_ETDRK3G(N,T,dt,dtn);
```

## Output

The simulation saves several MATLAB data files in the current folder, including

- `energy.mat` — discrete energy history;
- `max_min.mat` — maximum and minimum values of `phi` and `rho`;
- `phase_seperation.mat` — saved phase-field snapshots;
- `velocity.mat` — saved velocity snapshots;
- `pressure_T.mat` — pressure at the final time;
- `time.mat` — total elapsed wall-clock time.

After the simulation, `draw_figure.m` plots the energy evolution, extrema of `phi` and `rho`, and the saved `phi` snapshots.

## Notes

The physical and numerical parameters for the example are specified near the beginning of `NSBFS_ETDRK3G.m` and can be modified there. The supplied code is intended as a reproducible reference implementation for the ETDRK3-G scheme used in the accompanying manuscript.
