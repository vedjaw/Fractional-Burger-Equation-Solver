Fractional Burgers' Equation Solver

This repository contains a high-performance Python solver for the Time-Fractional Burgers' Equation, a non-linear partial differential equation (PDE) used to model shock waves and anomalous diffusion.

The solver utilizes a hybrid L1/Crank-Nicolson Predictor-Corrector scheme, ensuring unconditional stability and second-order spatial accuracy.

The Governing Equation

The code solves the following Time-Fractional Burgers' Equation:

$$\frac{\partial^\alpha u}{\partial t^\alpha} + a u \frac{\partial u}{\partial x} - c \frac{\partial^2 u}{\partial x^2} = 0, \quad (0 < \alpha \le 1)$$

Where:

$u(x,t)$: The velocity field.

$\alpha$: Fractional order of the time derivative (Caputo sense).

$a$: Advection coefficient.

$c$: Diffusion (viscosity) coefficient.

Features

Numerical Method: L1 approximation for the fractional time derivative + Crank-Nicolson for spatial discretization + Predictor-Corrector for non-linearity.

Efficiency: Reduces the problem to a tridiagonal linear system, solved efficiently ($O(N)$ complexity) using scipy.linalg.solve_banded.

Flexibility: Allows user-defined simulation parameters ($x_{max}, t_{max}, dx, dt$) and equation coefficients ($\alpha, a, c$).

Output: Generates a tabulated solution of $u(x,t)$ saved to a CSV file.

Installation

Clone the repository (or download the files to a local directory).

Install dependencies: Ensure you have Python installed. Then, install the required libraries using the provided requirements.txt file:

pip install -r requirements.txt


Usage

Configure Initial Conditions (Optional):
By default, the solver uses a sinusoidal initial condition: $u(x, 0) = \sin(\pi x)$.
If you wish to simulate a different scenario (e.g., a step function for shock waves), open the python script and modify the initial_condition function:

# Example: Change to a step function
def initial_condition(x):
    return np.where(x < 0.5, 1.0, 0.0) 


Run the Solver:
Execute the Python script from your terminal:

python fractional_burgers_solver.py


(Replace fractional_burgers_solver.py with the actual name of your python file).

Enter Simulation Parameters:
The program will prompt you for inputs. Example for a standard test case:

--- 1. INPUTS ---
Enter the max spatial value for x (e.g., 1.0): 1.0
Enter the max time value for t (e.g., 1.0): 1.0
Enter the space step (dx, e.g., 0.1): 0.02
Enter the time step (dt, e.g., 0.1): 0.02
Enter alpha (0 < alpha <= 1, e.g., 0.5): 0.9
Enter coefficient a (e.g., 1.0): 1.0
Enter coefficient c (e.g., 1.0): 0.1


Output

Console: The program will print a summary table of the computed $u(x,t)$ values to the console.

CSV File: A file named fractional_burgers_solution1.csv will be saved in the current directory. This file contains the full solution matrix, which can be opened in Excel or used for plotting.

Theoretical Validation

The method has been validated against exact analytical solutions (e.g., Cole-Hopf transformation for $\alpha=1$). The solver exhibits a temporal convergence rate of $O(\Delta t^{2-\alpha})$, confirming the accurate modeling of memory effects.
