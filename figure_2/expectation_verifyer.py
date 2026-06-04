import numpy as np

k = 200
m = 0.2
d = 200
N = 1000

def calculate_expected_times(k, m, d, N):
    alpha = 2*m*(1-m)
    # beta = (1 - m)**2 + 2*m*(1-m)*(d-2)/(d-1) + m**2*(d-2)/(d-1)
    beta = 1

    Q = np.array([[(1-1/k)*(1-1/N)*(1-alpha), (1-1/k) * alpha],
                [(1/d)*(1-1/N), (1-1/d)]])

    I = np.eye(2)
    N_matrix = np.linalg.inv(I - Q)

    # Calculate expected absorption times
    expected_times = N_matrix @ np.ones((2, 1))

    T_S = expected_times[0, 0]
    T_B = expected_times[1, 0]
    T_T = T_S/d + T_B*(1-1/d)

    return T_S, T_B, T_T

def analytical_expectations(k, m, d, N):
    alpha = 2*m*(1-m)
    T_0 = (N*(alpha*d*k - alpha*d + k))/(N + k - 1)
    T_1 = (N*alpha*d*k - N*alpha*d + N*d + N*k - alpha*d*k + alpha*d + d*k - d - k)/(N+k-1)
    T_T = (N*d*(alpha*d*(k-1) + d + k - 1) + (d-1)*((1-alpha)*d*(k-1) - k ))/(d*(N + k - 1))

    return T_0, T_1, T_T

T_S, T_B, T_T_matrix = calculate_expected_times(k, m, d, N)
T_0, T_1, T_T = analytical_expectations(k, m, d, N)

print(f"T_0 residual: {T_S - T_0}")
print(f"T_1 residual: {T_B - T_1}")
print(f"T_T residual: {T_T_matrix - T_T}")
