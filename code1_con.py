import numpy as np
import pandas as pd
from scipy.special import gamma
from scipy.linalg import solve_banded 

# Helper function for input
def get_float_input(prompt):
    """Gets a valid float from the user."""
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Invalid input. Please enter a number.")

print("1. INPUTS")

# Simulation parameters
x_max = get_float_input("Enter the max spatial value for x (e.g., 1.0): ")
t_max = get_float_input("Enter the max time value for t (e.g., 1.0): ")
dx = get_float_input("Enter the space step (dx, e.g., 0.1): ")
dt = get_float_input("Enter the time step (dt, e.g., 0.1): ")

# Equation parameters
alpha = get_float_input("Enter alpha (0 < alpha <= 1, e.g., 0.5): ")
a = get_float_input("Enter coefficient a (e.g., 1.0): ")
c = get_float_input("Enter coefficient c (e.g., 1.0): ")

# SCHEME IMPLEMENTATION

print("\n 2. APPLYING NUMERICAL SCHEME ")

# A: Setup the domain and solution array
x = np.linspace(0, x_max, int(x_max / dx) + 1)
t = np.linspace(0, t_max, int(t_max / dt) + 1)
N = len(x)  # Number of space points
M = len(t)  # Number of time points

# 2D solution array (M rows, N columns)
U = np.zeros((M, N))

# B: Apply Initial and Boundary Conditions
def initial_condition(x):
    # Example: u(x, 0) = sin(πx)
    return np.sin(np.pi * x)

def boundary_left(t):
    # Example: u(0, t) = 0
    return 0 * t # Returns a vector of zeros

def boundary_right(t):
    # Example: u(x_max, t) = 1.0
    return 1.0 + 0 * t # Returns a vector of ones

U[0, :] = initial_condition(x)
U[:, 0] = boundary_left(t)
U[:, -1] = boundary_right(t)

print(f"Grid setup: {M} time steps, {N} space points.")

# C: Pre-calculate constants for efficiency
K = (dt**(-alpha)) / gamma(2 - alpha)
# Pre-calculate all b_k weights (memory kernel)
k_vec = np.arange(M)
b = (k_vec + 1)**(1 - alpha) - k_vec**(1 - alpha)

# Number of internal points to solve for
N_internal = N - 2

# D: The Main Time-Stepping Loop
for n in range(0, M - 1): # Loop from t_0 to t_{M-1} to find t_{n+1}
    
    # We are solving for the unknown U[n+1, 1:-1]
    
    # 1. Build the Right-Hand-Side (RHS) vector R_j^n
    
    # Get current u_j^n values (internal points)
    u_n_j = U[n, 1:-1]
    
    # H_j^n (History Term)
    history_sum = np.zeros(N_internal)
    for k in range(1, n + 1):
        u_diff = U[n + 1 - k, 1:-1] - U[n - k, 1:-1]
        history_sum += b[k] * u_diff
    
    H_j = K * u_n_j - K * history_sum

    # S_j^n (Explicit Spatial Term)
    u_n_jm1 = U[n, 0:-2] # j-1
    u_n_jp1 = U[n, 2:]   # j+1
    
    u_n_x = (u_n_jp1 - u_n_jm1) / (2 * dx)
    u_n_xx = (u_n_jp1 - 2 * u_n_j + u_n_jm1) / (dx**2)
    
    S_j = 0.5 * (-a * u_n_j * u_n_x + c * u_n_xx)
    
    # Total RHS
    R_j = H_j + S_j

    # 2. Build the Left-Hand-Side (LHS) Tridiagonal Matrix
    
    # A_j, B_j, C_j (as vectors for all internal j)
    # We use u_n_j (the predictor) to linearize
    A_j = 0.5 * (-a * u_n_j / (2 * dx) - c / (dx**2)) # Coeff for u_{j-1}^{n+1}
    B_j = K - 0.5 * (c * (-2 / (dx**2)))              # Coeff for u_j^{n+1}
    C_j = 0.5 * (a * u_n_j / (2 * dx) - c / (dx**2))  # Coeff for u_{j+1}^{n+1}

    # 3. Account for Boundaries in RHS
    R_j[0] -= A_j[0] * U[n + 1, 0]   # Move A_1 * u_0^{n+1} to RHS
    R_j[-1] -= C_j[-1] * U[n + 1, -1] # Move C_{N-2} * u_{N-1}^{n+1} to RHS

    # 4. Create the banded matrix for solve_banded
    A_banded = np.zeros((3, N_internal))
    A_banded[0, 1:] = C_j[:-1]  # Upper diagonal (C_1, C_2, ...)
    A_banded[1, :] = B_j        # Main diagonal (B_1, B_2, ...)
    A_banded[2, :-1] = A_j[1:]  # Lower diagonal (A_2, A_3, ...)

    # 5. Solve the tridiagonal system
    try:
        solution = solve_banded((1, 1), A_banded, R_j)
        # Store the solution in our main array
        U[n + 1, 1:-1] = solution
    except np.linalg.LinAlgError:
        print(f"Error: Linear system solve failed at time step {n}.")
        break
    
print("Scheme implementation complete.")

# 3. OUTPUT THE RESULTS

print("\n--- 3. FINAL RESULTS TABLE ---")

# Creating a pandas DataFrame for output
# Rounding the values for cleaner display
df = pd.DataFrame(U, index=np.round(t, 4), columns=np.round(x, 4))
df.index.name = "Time (t)"
df.columns.name = "Space (x)"


print(df)