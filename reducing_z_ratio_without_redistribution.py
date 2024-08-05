import numpy as np
from scipy.stats import multivariate_normal, t, gamma, uniform
from scipy.spatial import distance_matrix
from joblib import Parallel, delayed
from sklearn.utils import resample


def rank_distance_upgraded(X, Y, Z):
    n = X.shape[0]  # Number of rows in X
    m = Y.shape[0]  # Number of rows in Y
    z_size = Z.shape[0]
    
    ZN = Z
    DZ = distance_matrix(ZN, ZN)

    DZX = distance_matrix(np.vstack([X, ZN]), np.vstack([X, ZN]))
    DX = DZX[:n, n:(n + z_size)]

    DZY = distance_matrix(np.vstack([Y, ZN]), np.vstack([Y, ZN]))
    DY = DZY[:m, m:(m + z_size)]

    Tz = np.zeros(z_size)

    # Calculate Wx and Wy
    for i in range(z_size):
        WX = np.sum(DZ[i, :] <= DX[:, i].reshape(-1, 1), axis=0)
        WY = np.sum(DZ[i, :] <= DY[:, i].reshape(-1, 1), axis=0)

        Tz[i] = np.sum(((WX / n - WY / m) ** 2)) * ((m + n) ** -2)
    
    # Calculate the final statistic
    stat = np.mean(Tz)
    return stat

# Function for bootstrap test
def bootstrap_test(X, Y, Z, B=500):
    combined = np.vstack((X, Y))
    n, m = X.shape[0], Y.shape[0]
    T_obs = rank_distance_upgraded(X, Y, Z)
    T_bootstrap = np.zeros(B)
    
    for i in range(B):
        indices = np.random.choice(n + m, n + m, replace=True)
        X_bootstrap = combined[indices[:n], :]
        Y_bootstrap = combined[indices[n:], :]
        T_bootstrap[i] = rank_distance_upgraded(X_bootstrap, Y_bootstrap, Z)
    
    p_value = np.mean(T_bootstrap >= T_obs)
    return p_value

# Function to generate samples
def generate_samples(n, m, d, distribution='normal', different_distributions=False):
    if distribution == 'normal':
        X = multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=n)
        if not different_distributions:
            Y = multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=m)
        else:
            Y = multivariate_normal.rvs(mean=np.ones(d), cov=np.eye(d), size=m)
    elif distribution == 't':
        X = t.rvs(df=5, size=(n, d))
        if not different_distributions:
            Y = t.rvs(df=5, size=(m, d))
        else:
            Y = t.rvs(df=10, size=(m, d))
    elif distribution == 'gamma':
        X = gamma.rvs(a=2, scale=1, size=(n, d))
        if not different_distributions:
            Y = gamma.rvs(a=2, scale=1, size=(m, d))
        else:
            Y = gamma.rvs(a=5, scale=1, size=(m, d))
    elif distribution == 'uniform':
        X = uniform.rvs(loc=0, scale=1, size=(n, d))
        if not different_distributions:
            Y = uniform.rvs(loc=0, scale=1, size=(m, d))
        else:
            Y = uniform.rvs(loc=0.5, scale=1, size=(m, d))
    else:
        raise ValueError("Unknown distribution")
    return X, Y

# Function to generate Z sample
def generate_Z(n, m, d, distribution='normal', different_distributions=False):
    size = n + m
    if distribution == 'normal':
        Z = multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=size)
    elif distribution == 't':
        Z = t.rvs(df=5, size=(size, d))
    elif distribution == 'gamma':
        Z = gamma.rvs(a=2, scale=1, size=(size, d))
    elif distribution == 'uniform':
        Z = uniform.rvs(loc=0, scale=1, size=(size, d))
    else:
        raise ValueError("Unknown distribution")
    return Z

# Function to generate joint Z sample for different distributions
def generate_joint_Z(n, m, d, distribution='normal', different_distributions=False):
    size = n + m
    if distribution == 'normal':
        Z = np.vstack((multivariate_normal.rvs(mean=np.zeros(d), cov=np.eye(d), size=size//2),
                       multivariate_normal.rvs(mean=np.ones(d), cov=np.eye(d), size=size - size//2)))
    elif distribution == 't':
        Z = np.vstack((t.rvs(df=5, size=(size//2, d)),
                       t.rvs(df=10, size=(size - size//2, d))))
    elif distribution == 'gamma':
        Z = np.vstack((gamma.rvs(a=2, scale=1, size=(size//2, d)),
                       gamma.rvs(a=5, scale=1, size=(size - size//2, d))))
    elif distribution == 'uniform':
        Z = np.vstack((uniform.rvs(loc=0, scale=1, size=(size//2, d)),
                       uniform.rvs(loc=0.5, scale=1, size=(size - size//2, d))))
    else:
        raise ValueError("Unknown distribution")
    return Z

# Function to run an experiment
def run_experiment(n, m, d, distribution, Z_ratios, different_distributions=False, B=500):
    results = []
    X, Y = generate_samples(n, m, d, distribution, different_distributions)
    
    if different_distributions:
        Z_full = generate_joint_Z(n, m, d, distribution, different_distributions)
    else:
        Z_full = generate_Z(n, m, d, distribution, different_distributions)
    
    for ratio in Z_ratios:
        Z_sample = Z_full[:int((n + m) * ratio), :]
        p_value = bootstrap_test(X, Y, Z_sample, B)
        results.append({
            'n': n, 'm': m, 'd': d, 'distribution': distribution,
            'ratio': ratio, 'different_distributions': different_distributions,
            'p_value': p_value
        })
    
    return results

# Function to run the experiments with different settings
def run(x):
    if x > 30:
        Z_ratios = [1.0, 0.75, 0.60, 0.50, 0.37, 0.25, 0.10, 0.05, 0.01]
    else:
        Z_ratios = [1.0, 0.75, 0.60, 0.50, 0.37, 0.25, 0.10]
    n_values = [x]
    m_values = [x]
    d_values = [10, 40]
    distributions = ['normal', 't', 'gamma', 'uniform']
    
    results = []
    for n in n_values:
        for m in m_values:
            for d in d_values:
                for distribution in distributions:
                    for different_distributions in [False, True]:
                        result = run_experiment(n, m, d, distribution, Z_ratios, different_distributions)
                        results.extend(result)
    
    output_file = f"nowe_wyniki/{x}_results.csv"
    import pandas as pd
    final_results = pd.DataFrame(results)
    final_results.to_csv(output_file, index=False)

# Running the experiments
run(10)
run(30)
run(100)
run(300)
run(500)
run(700)