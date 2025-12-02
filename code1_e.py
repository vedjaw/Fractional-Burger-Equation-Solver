import numpy as np
import pandas as pd
from scipy.special import gamma
from scipy.linalg import solve_banded

def get_float_input(prompt):
    while True:
        try:
            return float(input(prompt))
        except ValueError:
            print("Invalid input. Please enter a number.")

print("1. INPUTS")
x_max = get_float_input("Enter the max spatial value for x : ")
t_max = get_float_input("Enter the max time value for t : ")
dx = get_float_input("Enter the space step : ")
dt = get_float_input("Enter the time step : ")
alpha = get_float_input("Enter alpha (0 < alpha <= 1): ")
a = get_float_input("Enter coefficient a : ")
c = get_float_input("Enter coefficient c : ")

print("\n 2. APPLYING NUMERICAL SCHEME ")

x = np.linspace(0, x_max, int(x_max / dx) + 1)
t = np.linspace(0, t_max, int(t_max / dt) + 1)
N = len(x)
M = len(t)

U = np.zeros((M, N))

def initial_condition(x):
    return x

def boundary_left(t):
    return 0 * t

def boundary_right(t):
    return 1.0 + 0 * t

U[0, :] = initial_condition(x)
U[:, 0] = boundary_left(t)
U[:, -1] = boundary_right(t)

print(f"Grid setup: {M} time steps, {N} space points.")

K = (dt**(-alpha)) / gamma(2 - alpha)
k_vec = np.arange(M)
b = (k_vec + 1)**(1 - alpha) - k_vec**(1 - alpha)

N_internal = N - 2

for n in range(0, M - 1):
    
    u_n_j = U[n, 1:-1]
    
    history_sum = np.zeros(N_internal)
    for k in range(1, n + 1):
        u_diff = U[n + 1 - k, 1:-1] - U[n - k, 1:-1]
        history_sum += b[k] * u_diff
    
    H_j = K * u_n_j - K * history_sum

    u_n_jm1 = U[n, 0:-2]
    u_n_jp1 = U[n, 2:]
    
    u_n_x = (u_n_jp1 - u_n_jm1) / (2 * dx)
    u_n_xx = (u_n_jp1 - 2 * u_n_j + u_n_jm1) / (dx**2)
    
    S_j = 0.5 * (-a * u_n_j * u_n_x + c * u_n_xx)
    
    R_j = H_j + S_j

    A_j = 0.5 * (-a * u_n_j / (2 * dx) - c / (dx**2))
    B_j = K - 0.5 * (c * (-2 / (dx**2)))
    C_j = 0.5 * (a * u_n_j / (2 * dx) - c / (dx**2))

    R_j[0] -= A_j[0] * U[n + 1, 0]
    R_j[-1] -= C_j[-1] * U[n + 1, -1]

    A_banded = np.zeros((3, N_internal))
    A_banded[0, 1:] = C_j[:-1]
    A_banded[1, :] = B_j
    A_banded[2, :-1] = A_j[1:]

    try:
        solution = solve_banded((1, 1), A_banded, R_j)
        U[n + 1, 1:-1] = solution
    except np.linalg.LinAlgError:
        print(f"Error: Linear system solve failed at time step {n}.")
        break
    
print("Scheme implementation complete.")

print("\n 3. FINAL RESULTS TABLE ")

df = pd.DataFrame(U, index=np.round(t, 4), columns=np.round(x, 4))
df.index.name = "Time (t)"
df.columns.name = "Space (x)"

print(df)

try:
    df.to_csv("fractional_burgers_solution1.csv")
    print("\nResults successfully saved to 'fractional_burgers_solution.csv'")
except PermissionError:
    print("\nError: Could not save file. Check if it's open in another program.")