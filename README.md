# Optimized Schwarz Waveform Relaxation for the Damped Wave Equation

This repository contains the MATLAB implementation and reproducibility scripts for the research presented in the conference paper **"Optimized Schwarz Waveform Relaxation for the Damped Wave Equation"**.

The code implements an Overlapping Schwarz Waveform Relaxation (SWR) method with optimized transmission conditions to solve the one-dimensional viscoelastic-telegrapher damped wave equation.

## 📄 Reference Paper

This code supports the findings detailed in:

> **Optimized Schwarz Waveform Relaxation for the Damped Wave Equation**
> *Gerardo Cicalese, Gabriele Ciaramella, Ilario Mazzieri, and Martin J. Gander*
> Proceedings of the 29th International Conference on Domain Decomposition Methods (DD29).

## 🧮 Mathematical Model

The solver addresses the general damped wave equation on a domain $\Omega=(0,L)$ for time $t \ge 0$. The governing equation includes both **telegrapher damping** ($\gamma$) and **viscoelastic damping** ($\nu$):

$$
\partial_{tt}u + \gamma\partial_{t}u = c^{2}\partial_{xx}u + \nu\partial_{t}\partial_{xx}u + f
$$

Where:
* $c$: Wave speed.
* $\gamma$: Telegrapher damping coefficient (viscous).
* $\nu$: Viscoelastic damping coefficient (structural).

The problem is closed with homogeneous Dirichlet boundary conditions $u(0,t)=u(L,t)=0$.

## ⚙️ Algorithm: Schwarz Waveform Relaxation

In general, the domain $\Omega$ is decomposed into $N$ overlapping subdomains $\Omega_j$ that cover the original interval $(0, L)$. The algorithm iterates in parallel over these subdomains, exchanging information at the boundaries.

For the mathematical analysis and core testing presented in the paper, the setup often focuses on the two-subdomain case ($\Omega_1 = (0, b)$ and $\Omega_2 = (a, L)$ with overlap $\delta = b - a$), but the method scales to multiple subdomains.

The iteration uses **Robin transmission conditions** parameterized by $\Lambda = p\partial_t + q$ to optimize convergence. For a general subdomain $\Omega_j$ with neighbors $\Omega_{j-1}$ and $\Omega_{j+1}$:

* **Left Boundary:** Matches transmission condition $(\partial_x - \Lambda)u_{j-1}$ from the left neighbor.
* **Right Boundary:** Matches transmission condition $(\partial_x + \Lambda)u_{j+1}$ from the right neighbor.

### Parameter Optimization

The repository includes scripts to optimize the transmission parameters $(p, q)$ using spectral analysis. The optimization minimizes the convergence factor $\rho$ over relevant frequencies $[\omega_{min}, \omega_{max}]$.

Two optimization strategies are implemented:
1.  **$L_\infty$ Optimization:** Minimizes the maximum convergence factor $\hat{\rho}_{\infty}$.
2.  **$L_2$ Optimization:** Minimizes the root-mean-square (RMS) of the convergence factor $\hat{\rho}_{2}$.

## 📂 Repository Structure

The repository is organized into root-level execution scripts, a utility library for core logic, and a figures folder for results.

### 1. Root Directory: Simulation & Analysis Drivers
These scripts are the primary entry points for running simulations and reproducing paper results.

**Configuration & Interactive Tools**
* `get_sim_params.m`: Configuration file for setting physical parameters ($c, \gamma, \nu$) and grid discretization ($\Delta x, \Delta t$).
* `interactive_rho_all_params.m`: Interactive GUI tool to explore how changes in transmission parameters and damping coefficients affect the theoretical convergence factor $\rho$.
* `compare_surfaces.m`: Visualization utility that compares theoretical global contraction factor surfaces against experimental results.

**Solvers & Testers**
* `run_swr_tester.m`: Main wrapper script to run the SWR simulation (calls `utils/run_swr.m`).
* `run_fdtd_tester.m`: Runs the reference Finite Difference Time Domain (FDTD) solver (calls `utils/run_fdtd.m`).
* `run_fdtd_tester_varying_gamma.m`: Runs FDTD simulations while systematically varying the viscous damping coefficient $\gamma$.
* `run_fdtd_tester_varying_nu.m`: Runs FDTD simulations while systematically varying the viscoelastic damping coefficient $\nu$.
* `run_fdtd_conv_test.m`: Verifies the numerical convergence order of the FDTD scheme.
* `weak_scalability_test.m`: Tests the algorithm's scalability by analyzing performance as the number of subdomains and problem size increase.

**Optimization & Plotting**
* `optimal_parameters_vs_damping_L2.m`: Computes optimal $(p, q)$ parameters minimizing the $L_2$ norm of the convergence factor.
* `optimal_parameters_vs_damping_Linf.m`: Computes optimal parameters minimizing the $L_\infty$ norm.
* `SWR_error_vs_iteration_L2.m`: Simulates and plots error convergence behavior ($L_2$ norm) over iterations.
* `SWR_error_vs_iteration_Linf.m`: Simulates and plots error convergence behavior ($L_\infty$ norm) over iterations.

### 2. Utils Directory (`utils/`)
Contains the core backend functions and mathematical definitions.

**Core Solvers**
* `run_swr.m`: Implementation of the Schwarz Waveform Relaxation loop with optimized Robin transmission conditions.
* `run_swr_dirichlet.m`: Implementation of SWR using classical Dirichlet transmission conditions.
* `run_fdtd.m`: Core implementation of the FDTD solver for a single domain.

**Mathematical definitions**
* `rho.m`: Computes the frequency-dependent convergence factor $\rho(\omega)$ for the 2 (analytical) or N (numerical) subdomains case.
* `matrix_G.m`: Constructs the convergence matrix used in the analysis for the N subdomains case.
* `analytic_solution_single_mode.m`: Computes the exact analytical solution for single-mode verification.

**Optimization Objectives**
* `obj_L2.m`: Calculates the global convergence factor using the $L_2$ norm (RMS).
* `obj_Linf.m`: Calculates the global convergence factor using the $L_{\infty}$ norm (Maximum).
* `obj_exp.m`: Computes the experimentally determined global convergence factor from simulation data.

**Plotting Helpers**
* `contraction_surface.m`: Helper to generate data for theoretical contraction factor surfaces.
* `swr_error_surface.m`: Helper to generate data for experimental error surfaces.

## 📊 Key Results

The implementation demonstrates the distinct behaviors of the two damping types identified in the paper:

1.  **Viscoelastic Damping ($\nu > 0, \gamma=0$):**
    * The method exhibits approximately linear convergence.
    * Optimization strategies successfully identify parameters that drastically improve performance.

2.  **Telegrapher Damping ($\gamma > 0, \nu=0$):**
    * Convergence is sublinear and resembles the behavior of the undamped wave equation.
    * Optimized parameters do not significantly alter the underlying convergence mechanism, which is limited by finite wave propagation speed.

## 💻 Usage

To reproduce the convergence results (e.g., Figure 2 in the paper):

1.  Clone the repository:
    ```bash
    git clone [https://github.com/mazamin7/damped_wave_equation_SWR.git](https://github.com/mazamin7/damped_wave_equation_SWR.git)
    cd damped_wave_equation_SWR
    ```

2.  Open MATLAB.

3.  Set the physical parameters in `get_sim_params.m` (e.g., set $\nu > 0$ for viscoelastic tests).

4.  Compare the SWR solution against the FDTD reference using the plotting scripts:
    ```matlab
    SWR_error_vs_iteration_Linf
    SWR_error_vs_iteration_L2
    ```

## 🔗 Citation

If you use this code or the methods described, please cite the associated paper:

```bibtex
@inproceedings{Cicalese2026SWR,
  author    = {Gerardo Cicalese and Gabriele Ciaramella and Ilario Mazzieri and Martin J. Gander},
  title     = {Optimized Schwarz Waveform Relaxation for the Damped Wave Equation},
  booktitle = {Proceedings of the 29th International Conference on Domain Decomposition Methods},
  year      = {2026}
}
