# Fractional Burgers' Equation Solver

This repository contains a high-performance **Python solver** for the **time-fractional Burgers' equation**, a non-linear PDE used to model **shock waves**, **memory effects**, and **anomalous diffusion**.

It uses a hybrid **L1 fractional derivative + Crank–Nicolson + predictor–corrector** numerical scheme, achieving **unconditional stability** and **second-order spatial accuracy**.

## Governing Equation

The solver computes numerical solutions to the time-fractional Burgers' equation:

$$
\frac{\partial^\alpha u}{\partial t^\alpha} + au\frac{\partial u}{\partial x} - c\frac{\partial^2 u}{\partial x^2} = 0, \quad (0 < \alpha \le 1)
$$

where:

- \(u(x,t)\): velocity field  
- $\alpha$: fractional order of the time derivative (Caputo sense)  
- $a$: advection coefficient  
- $c$: viscosity/diffusion coefficient  

---

## Features

- **Numerical method**
  - L1 approximation for the Caputo fractional time derivative  
  - Crank–Nicolson scheme for the diffusion term  
  - Predictor–corrector iteration for handling the non-linear advection term  

- **Efficiency**
  - Reduces the semi-discrete system to a tridiagonal linear system at each time step  
  - Solved efficiently with `scipy.linalg.solve_banded` (\(O(N)\) complexity per step)  

- **Flexibility**
  - User-defined simulation parameters:
    - $x_{\max}$, $t_{\max}$, $\Delta x$, $\Delta t\$
  - User-defined equation coefficients:
    - $\alpha$, $a$, $c$

- **Output**
  - Generates a tabulated solution u(x,t)  
  - Saves the full solution matrix to a CSV file  

---

## Installation

You need Python and the following libraries:

- `numpy`  
- `scipy`  
- `pandas`  

Install them directly:

```bash
pip install numpy scipy pandas
```

Or, if you have a `requirements.txt`:

```bash
pip install -r requirements.txt
```

---

## Usage

### 1. Configure Initial Conditions (optional)

By default, the solver uses the sinusoidal initial condition  
$u(x, O) = \sin(\pi x)$.

To simulate a different scenario (for example, a step function to model a shock), modify the `initial_condition` function in your Python script:

```python
def initial_condition(x):
    # Example: step function
    return np.where(x < 0.5, 1.0, 0.0)
```

### 2. Run the Solver

From a terminal in the project directory, run:

```bash
python fractional_burgers_solver.py
```

Replace `fractional_burgers_solver.py` with the actual name of your solver script if it is different.

### 3. Enter Simulation Parameters

When prompted, provide the simulation settings. For example:

```text
--- 1. INPUTS ---
Enter the max spatial value for x (e.g., 1.0): 1.0
Enter the max time value for t (e.g., 1.0): 1.0
Enter the space step (dx, e.g., 0.1): 0.02
Enter the time step (dt, e.g., 0.1): 0.02
Enter alpha (0 < alpha <= 1, e.g., 0.5): 0.9
Enter coefficient a (e.g., 1.0): 1.0
Enter coefficient c (e.g., 1.0): 0.1
```

---

## Output

- **CSV file**  
  A file such as `fractional_burgers_solution1.csv` is written to the current directory.  
  This file contains the full solution matrix (rows = time levels, columns = spatial nodes) and can be:
  - Opened in Excel or other spreadsheet software  
  - Used for plotting and post-processing in Python, MATLAB, etc.  

---

## Theoretical Validation

The numerical method has been validated against exact analytical solutions, including the classical Burgers' equation obtained for $\alpha$ (via the Cole–Hopf transformation).

The scheme exhibits a temporal convergence rate of $O(\Delta t^{2-alpha})$,, which is consistent with the expected accuracy of L1-based fractional time-stepping methods and confirms that the solver captures the **memory effects** inherent to fractional-order dynamics.

---
