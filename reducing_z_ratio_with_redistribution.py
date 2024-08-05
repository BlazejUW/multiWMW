import numpy as np
from scipy.stats import multivariate_normal, t, gamma
from scipy.spatial.distance import pdist, squareform
import pandas as pd
from multiprocessing import Pool
import time
import itertools

def dist_matrix(data):
    distances = squareform(pdist(data, metric='euclidean'))
    return distances

def rank_distance_upgraded(X, Y, Z):
    n = X.shape[0]
    m = Y.shape[0]
    z_size = Z.shape[0]
    ZN = Z

    DZ = dist_matrix(ZN)
    DZX = dist_matrix(np.vstack([X, ZN]))
    DZY = dist_matrix(np.vstack([Y, ZN]))

    DX = DZX[:n, n:(n + z_size)]
    DY = DZY[:m, m:(m + z_size)]

    Tz = np.zeros(z_size)

    for i in range(z_size):
        WX = 0
        WY = 0
        sorted_DZ = np.sort(DZ[i, :])

        for k in range(n):
            WX += np.sum(sorted_DZ <= DX[k, i])

        for l in range(m):
            WY += np.sum(sorted_DZ <= DY[l, i])

        Tz[i] = (WX / n - WY / m) ** 2 / (m + n) ** 2

    stat = np.mean(Tz)
    return stat

def bootstrap_iteration(i, combined, n, m, Z):
    np.random.seed(i)
    indices = np.random.choice(n + m, n + m, replace=True)
    X_bootstrap = combined[indices[:n], :]
    Y_bootstrap = combined[indices[n:], :]
    return rank_distance_upgraded(X_bootstrap, Y_bootstrap, Z)

def bootstrap_test(X, Y, Z, B=500):
    n = X.shape[0]
    m = Y.shape[0]
    combined = np.vstack([X, Y])
    T_obs = rank_distance_upgraded(X, Y, Z)
    T_bootstrap = np.zeros(B)

    with Pool() as pool:
        T_bootstrap = pool.starmap(bootstrap_iteration, [(i, combined, n, m, Z) for i in range(B)])

    p_value = np.mean(np.array(T_bootstrap) >= T_obs)
    return p_value

def generate_data(total_size_X, total_size_Y, d, distribution_X='normal', distribution_Y='normal'):
    if distribution_X == 'normal':
        data_X = multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=total_size_X)
    elif distribution_X == 't':
        data_X = t.rvs(df=5, size=(total_size_X, d))
    elif distribution_X == 'gamma':
        data_X = gamma.rvs(a=2, size=(total_size_X, d))
    elif distribution_X == 'uniform':
        data_X = np.random.uniform(0, 1, size=(total_size_X, d))
    else:
        raise ValueError("Unknown distribution")

    if distribution_Y == 'normal':
        data_Y = multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=total_size_Y)
    elif distribution_Y == 't':
        data_Y = t.rvs(df=5, size=(total_size_Y, d))
    elif distribution_Y == 'gamma':
        data_Y = gamma.rvs(a=2, size=(total_size_Y, d))
    elif distribution_Y == 'uniform':
        data_Y = np.random.uniform(0, 1, size=(total_size_Y, d))
    else:
        raise ValueError("Unknown distribution")

    return data_X, data_Y

def generate_scenarios(initial_X, initial_Y):
    proportions = [1.0, 0.75, 0.6, 0.5, 0.37, 0.25]
    scenarios = []
    total_size = initial_X + initial_Y

    for proportion in proportions:
        z_size = int(proportion * total_size)
        remaining = 2 * total_size - z_size
        x_size = int(remaining / 2)
        y_size = remaining - x_size
        scenarios.append((x_size, y_size, z_size, f"{x_size}/{y_size}/{z_size}"))

    return scenarios

def run_scenario(data_X, data_Y, n, m, z_size, scenario_name, B=500):
    Z_X = data_X[n:n + z_size // 2]
    Z_Y = data_Y[m:m + z_size // 2]
    Z = np.vstack([Z_X, Z_Y])

    start_time = time.time()
    p_value = bootstrap_test(data_X[:n], data_Y[:m], Z, B)
    end_time = time.time()
    time_taken = end_time - start_time
    print(f"{scenario_name} p-value: {p_value} Time taken: {time_taken}")
    return {'p_value': p_value, 'time_taken': time_taken, 'scenario_name': scenario_name}

def run_experiments(initial_X, initial_Y, d=10, distribution_X='normal', distribution_Y='normal', B=500):
    total_size_X = initial_X * 2
    total_size_Y = initial_Y * 2
    while True:
        data_X, data_Y = generate_data(total_size_X, total_size_Y, d, distribution_X, distribution_Y)
        results = []

        scenarios = generate_scenarios(initial_X, initial_Y)

        first_scenario = run_scenario(data_X, data_Y, scenarios[0][0], scenarios[0][1], scenarios[0][2], scenarios[0][3], B)
        if first_scenario['p_value'] >= 0.1 or distribution_X != distribution_Y:
            results.append(first_scenario)
            for scenario in scenarios[1:]:
                result = run_scenario(data_X, data_Y, scenario[0], scenario[1], scenario[2], scenario[3], B)
                results.append(result)
            break
        else:
            print(f"Low p-value ({first_scenario['p_value']}) for scenario {scenarios[0][3]} with distributions {distribution_X} and {distribution_Y}, repeating the experiment...")

    results_df = pd.DataFrame(results)
    return results_df

def run_and_save_all_experiments(initial_X=300, initial_Y=300, d_values=[10], distribution_types=None, B=500):
    if distribution_types is None:
        distribution_types = ['normal', 't', 'gamma', 'uniform']
    
    distribution_pairs = list(itertools.combinations_with_replacement(distribution_types, 2))

    all_results = []

    for pair in distribution_pairs:
        for d in d_values:
            print(f"Running experiments for distribution pair: {pair} and dimension: {d}")
            results = run_experiments(initial_X, initial_Y, d=d, distribution_X=pair[0], distribution_Y=pair[1], B=B)
            results['Distribution_X'] = pair[0]
            results['Distribution_Y'] = pair[1]
            results['Dimension'] = d
            all_results.append(results)
    
    final_results_df = pd.concat(all_results, ignore_index=True)
    final_results_df.to_csv(f"combined_results_{initial_X}.csv", index=False)


# Uruchomienie eksperymentów dla wszystkich par rozkładów
sample_sizes = [10, 50, 100, 200, 300, 500]
distribution_types = ['normal', 't', 'gamma', 'uniform']

for size in sample_sizes:
    print(f"Running experiments for sample size: {size}")
    run_and_save_all_experiments(initial_X=size, initial_Y=size, d_values=[20], distribution_types=distribution_types, B=500)
