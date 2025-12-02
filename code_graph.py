import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.linalg import solve_banded

# --- 1. INPUT FUNCTIONS ---
def get_float_input(prompt):
    """Safely gets a float input from the user."""
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Invalid input. Please enter a number.")

def get_int_input(prompt):
    """Safely gets an integer input from the user."""
    while True:
        try:
            return int(input(prompt))
        except ValueError:
            print("Invalid input. Please enter an integer.")

# --- 2. USER INPUTS ---
print("1. INPUTS")
x_max = get_float_input("Enter the max spatial value for x : ")
t_max = get_float_input("Enter the max time value for t : ")
dx = get_float_input("Enter the space step (dx): ")
dt = get_float_input("Enter the time step (dt): ")
alpha = get_float_input("Enter alpha (0 < alpha <= 1): ")
a = get_float_input("Enter coefficient a : ")
c = get_float_input("Enter coefficient c : ")

# --- 3. GRID SETUP ---
print("\n2. APPLYING NUMERICAL SCHEME ")

# Calculate grid points (ensure x_max/dx and t_max/dt are integers for consistent spacing)
N = int(round(x_max / dx)) + 1
M = int(round(t_max / dt)) + 1
if abs(x_max / dx - (N - 1)) > 1e-9 or abs(t_max / dt - (M - 1)) > 1e-9:
     print("Warning: Max values may not be perfectly divisible by step sizes.")

x = np.linspace(0, x_max, N)
t = np.linspace(0, t_max, M)
U = np.zeros((M, N))

# --- 4. INITIAL AND BOUNDARY CONDITIONS ---
def initial_condition(x):
    # Using the initial condition u(x,0) = x from the presentation's example
    return x 

def boundary_left(t):
    # Using the boundary condition u(0,t) = 0 from the presentation's example
    return 0 * t

def boundary_right(t):
    # Using the boundary condition u(1,t) = 1 from the presentation's example
    return 1.0 + 0 * t

U[0, :] = initial_condition(x)
U[:, 0] = boundary_left(t)
U[:, -1] = boundary_right(t)

print(f"Grid setup: {M} time steps, {N} space points.")

# --- 5. L1 FORMULA CONSTANTS ---
K = (dt**(-alpha)) / gamma(2 - alpha)
k_vec = np.arange(M)
# b_k = (k+1)^(1-alpha) - k^(1-alpha)
b = (k_vec + 1)**(1 - alpha) - k_vec**(1 - alpha)

N_internal = N - 2 # Number of unknown points to solve for

# --- 6. TIME-STEPPING LOOP (TRIDIAGONAL SYSTEM SOLVE) ---
for n in range(0, M - 1):
    
    # 6.1. CALCULATE HISTORY TERM (H_j)
    u_n_j = U[n, 1:-1] # u_j^n (interior points)
    
    # Summation for the memory term: sum_{k=1}^{n} b_k (u_j^{n+1-k} - u_j^{n-k})
    history_sum = np.zeros(N_internal)
    for k in range(1, n + 1):
        # Indexing for U: time index n+1-k goes from n down to 1.
        # time index n-k goes from n-1 down to 0.
        u_diff = U[n + 1 - k, 1:-1] - U[n - k, 1:-1]
        history_sum += b[k] * u_diff
    
    # H_j = K * u_j^n - K * history_sum (The history term for the RHS)
    H_j = K * u_n_j - K * history_sum

    # 6.2. CALCULATE SOURCE TERM (S_j)
    u_n_jm1 = U[n, 0:-2] # u_{j-1}^n
    u_n_jp1 = U[n, 2:]   # u_{j+1}^n
    
    # Spatial derivatives at time n
    u_n_x = (u_n_jp1 - u_n_jm1) / (2 * dx) # Central difference for u_x
    u_n_xx = (u_n_jp1 - 2 * u_n_j + u_n_jm1) / (dx**2) # Central difference for u_xx
    
    # S_j = 0.5 * (-a * u_j^n * u_x^n + c * u_xx^n)
    S_j = 0.5 * (-a * u_n_j * u_n_x + c * u_n_xx)
    
    R_j = H_j + S_j # The complete RHS vector

    # 6.3. CALCULATE TRIDIAGONAL MATRIX COEFFICIENTS (LHS)
    # The coefficients are based on u_j^n (the predictor step)
    
    # A_j = -(a*u_j^n)/(4*dx) - c/(2*dx^2) 
    A_j = 0.5 * (-a * u_n_j / (2 * dx) - c / (dx**2))
    
    # B_j = K + c/dx^2 (K is from L1, c/dx^2 is from the C-N implicit u_xx term)
    B_j = K + c / (dx**2)
    
    # C_j = (a*u_j^n)/(4*dx) - c/(2*dx^2)
    C_j = 0.5 * (a * u_n_j / (2 * dx) - c / (dx**2))

    # 6.4. HANDLE BOUNDARY CONDITIONS ON RHS
    # For j=1 (first internal point): A_1 * u_0^{n+1} is moved to the RHS.
    R_j[0] -= A_j[0] * U[n + 1, 0]
    
    # For j=N-2 (last internal point): C_{N-2} * u_{N-1}^{n+1} is moved to the RHS.
    R_j[-1] -= C_j[-1] * U[n + 1, -1]

    # 6.5. CONSTRUCT BANDED MATRIX FOR SOLVE_BANDED
    # solve_banded expects a (3, N_internal) matrix for a tridiagonal system (1 lower, 1 upper).
    # Row 0: Upper diagonal (C_j)
    # Row 1: Main diagonal (B_j)
    # Row 2: Lower diagonal (A_j)
    
    A_banded = np.zeros((3, N_internal))
    
    # Upper Diagonal (C_j) - C_j[0] is not used, so shift C_j by 1
    A_banded[0, 1:] = C_j[:-1] 
    
    # Main Diagonal (B_j)
    A_banded[1, :] = B_j 
    
    # Lower Diagonal (A_j) - A_j[-1] is not used, so shift A_j by 1
    A_banded[2, :-1] = A_j[1:] 

    # 6.6. SOLVE THE SYSTEM
    try:
        # solve_banded((l, u), a_banded, b) where l=1, u=1
        solution = solve_banded((1, 1), A_banded, R_j)
        U[n + 1, 1:-1] = solution
    except np.linalg.LinAlgError:
        print(f"Error: Linear system solve failed at time step {n}.")
        break
    
print("Scheme implementation complete.")

# --- 7. PLOTTING THE RESULTS ---
print("\n3. PLOTTING RESULTS ")

# Define which time steps to plot
plot_indices = [0]
# Automatically select a few intermediate time steps and the final time step
intermediate_steps = np.linspace(1, M - 2, 3, dtype=int)
plot_indices.extend(intermediate_steps)
if M > 1:
    plot_indices.append(M - 1)

# Ensure indices are unique and within bounds
plot_indices = sorted(list(set(i for i in plot_indices if i >= 0 and i < M)))

plt.figure(figsize=(10, 6))
colors = plt.cm.viridis(np.linspace(0, 1, len(plot_indices)))

for i, n in enumerate(plot_indices):
    t_val = t[n]
    # Label each line with the time t
    plt.plot(x, U[n, :], label=f'$t = {t_val:.4f}$', color=colors[i], linewidth=2, marker='o', markersize=4)

plt.title(f'Solution $u(x,t)$ for Fractional Burgers\' Equation $(\\alpha={alpha}, a={a}, c={c})$')
plt.xlabel('Spatial Position $x$')
plt.ylabel('Solution $u$')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(title='Time Steps')
plt.tight_layout()
plt.show()

print("\nPlot displayed successfully.")