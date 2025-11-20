# Fractional Burgers' Equation Solver

This repository contains a high-performance **Python solver** for the **Time-Fractional Burgers' Equation**, a non-linear PDE used to model **shock waves**, **memory effects**, and **anomalous diffusion**.

It uses a hybrid **L1 fractional derivative + Crank–Nicolson + Predictor–Corrector** numerical scheme, achieving **unconditional stability** and **second-order spatial accuracy**.

## Governing Equation

The solver computes solutions to:

$$
\frac{\partial^\alpha u}{\partial t^\alpha} 
+ a\,u\,\frac{\partial u}{\partial x} 
- c\,\frac{\partial^2 u}{\partial x^2} 
= 0,\qquad 0 < \alpha \le 1
$$

where:

- \( u(x,t) \): velocity  
- \( \alpha \): fractional order (Caputo)  
- \( a \): advection coefficient  
- \( c \): viscosity coefficient  

# Features

- **L1 Caputo fractional derivative**
- **Crank–Nicolson diffusion**
- **Predictor–Corrector for nonlinearity**
- **Scipy tridiagonal solver (`solve_banded`)**
- User-defined:
  - \( x_{\max}, t_{\max}, \Delta x, \Delta t \)
  - \( \alpha, a, c \)
- Outputs:
  - Console summary
  - CSV file of full solution matrix

# Installation (Dependencies Included Here)

Dependencies:

```
numpy
scipy
pandas
```

Install:

```bash
pip install numpy scipy pandas
```

OR using a `requirements.txt`:

```bash
pip install -r requirements.txt
```

# Usage

## Step 1: Initial Condition

Default:

\[
u(x,0)=\sin(\pi x)
\]

To change (example: step function):

```python
def initial_condition(x):
    return np.where(x < 0.5, 1.0, 0.0)
```

## Step 2: Run

```bash
python fractional_burgers_solver.py
```

Example input:

```
--- 1. INPUTS ---
1.0
1.0
0.02
0.02
0.9
1.0
0.1
```

# Output

- Console table of \(u(x,t)\)
- CSV file:

```
fractional_burgers_solution1.csv
```

(rows = time, columns = space)

# Theoretical Validation

- Verified against Cole–Hopf analytical solution when \( \alpha = 1 \)
- Temporal accuracy of:

\[
O(\Delta t^{2-\alpha})
\]

# License

MIT (recommended)
