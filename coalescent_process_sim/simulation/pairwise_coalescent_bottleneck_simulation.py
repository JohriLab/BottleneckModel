import random
import pandas as pd
import numpy as np

def run_simulation(k, m):
    # Initialize the population
    generation = 0
    same_deme = random.random() < 1/d

    if same_deme:
        if random.random() < 1/N:
            return generation

    alpha = 2*m*(1-m)
    beta = (1-m)**2 + 2*m*(1-m)*(d-2)/(d-1) + m**2 * (d-2)/(d-1)

    # Run backwards in time
    while True:
        generation += 1
        if same_deme and random.random() < alpha:
            same_deme = False
        if not same_deme and random.random() < 1/d:
            same_deme = True
        if not same_deme and random.random() < 1 - beta:
            same_deme = True

        if same_deme:
            if random.random() < 1/k:
                return generation
            if random.random() < 1/N:
                return generation


m_values = [0.0001, 0.001, 0.01, 0.1, 0.4, 0.8]
k_values = [1, 2, 5, 10]
d = 100
N = 1e3
u = 1e-6

output_file = "analysis/pi_data_plotting.csv"
output_df = pd.DataFrame(columns=["k", "m", "pi"])

for m in m_values:
    for k in k_values:
        results = pd.DataFrame([{"k": k, "m": m, "T_mrca": run_simulation(k, m)} for _ in range(1000000)])
        # For each T_mrca, draw a pi value from the Poisson distribution with lambda = T_mrca * u * 2
        results["pi"] = np.random.poisson(lam=results["T_mrca"] * u * 2)
        mean_pi = results["pi"].mean()
        se_pi = results["pi"].std() / np.sqrt(1000000)
        output_df = pd.concat([output_df, pd.DataFrame([{"k": k, "m": m, "pi": mean_pi, "se_pi": se_pi}])], ignore_index=True)
        print("Completed m = ", m, "k = ", k)

output_df.to_csv(output_file, index=False)