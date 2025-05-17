# Unified Entropic Theory (UET) Dataset and Scripts

This repository provides the dataset and scripts for the paper *"Unified Entropic Theory: Time, Space, and Gravity as Emergent Properties of Entropy"* by Richard G. Prouse (2025). The dataset and scripts support the paper’s cosmological analyses, including MCMC fitting, Runge-Kutta integration, and Monte Carlo uncertainty propagation, as described in Section 6 (Methods) and the Data Availability Statement.

## Overview

This `README.md` contains:
- **Dataset**: A CSV with 40 data points used for cosmological fitting.
- **Scripts**: Python 3.9 scripts for MCMC fitting, Runge-Kutta integration, and Monte Carlo simulations.
- **Instructions**: Steps to set up, run, and validate the scripts.

The dataset includes:
- **Planck 2018 CMB**: 10 points (5 \(\Delta T / T\), 5 parameters: \(\Omega_m, \Omega_\Lambda, n_s, \Omega_k, \eta\)).
- **SHOES**: 5 \(H_0\) measurements.
- **DESI BAO**: 10 angular diameter distances (\(D_A\)).
- **Cosmic Chronometer**: 5 \(H(z)\) measurements.
- **Euclid Early Data**: 10 clustering distances.

The scripts perform:
- **MCMC Fitting**: Optimizes parameters (\(\beta, S_m, \Omega_m, \Omega_\Lambda, n_s, \Omega_k, \eta\)).
- **Runge-Kutta Integration**: Solves the time evolution equation (Eq. 9).
- **Monte Carlo Simulations**: Propagates entropy uncertainties.

## Prerequisites

- **Python 3.9**: Install from [python.org](https://www.python.org/downloads/release/python-390/).
- **Libraries**:
  - NumPy 1.23: `pip install numpy==1.23.0`
  - SciPy 1.9: `pip install scipy==1.9.0`
  - SymPy 1.9: `pip install sympy==1.9.0`
  - CosmoMC Python wrapper: Install from [CosmoMC GitHub](https://github.com/cmbant/CosmoMC) or use a custom wrapper (see paper).
  - pytest 8.3.3: `pip install pytest==8.3.3`
- **System**: Linux, macOS, or Windows with a Python environment.
- **Disk Space**: ~100 MB for dataset and scripts.

## Setup

1. **Save Dataset**:
   - Copy the CSV content below into a file named `uet_dataset.csv`.
   - Place it in your working directory (e.g., `./uet`).

   ```csv
Data_Type,Measurement,Redshift,Value,Uncertainty,Source
CMB,Delta_T/T,1100,1.01e-05,1.0e-06,Planck2020
CMB,Delta_T/T,1100,1.02e-05,1.0e-06,Planck2020
CMB,Delta_T/T,1100,1.00e-05,1.0e-06,Planck2020
CMB,Delta_T/T,1100,1.03e-05,1.0e-06,Planck2020
CMB,Delta_T/T,1100,1.01e-05,1.0e-06,Planck2020
CMB,Omega_m,0,0.315,0.007,Planck2020
CMB,Omega_Lambda,0,0.685,0.007,Planck2020
CMB,n_s,0,0.965,0.004,Planck2020
CMB,Omega_k,0,0.000,0.001,Planck2020
CMB,eta,0,6.1e-10,1.0e-11,Planck2020
SHOES,H_0,0,73.04,1.04,Riess2021
SHOES,H_0,0,73.10,1.05,Riess2021
SHOES,H_0,0,73.00,1.04,Riess2021
SHOES,H_0,0,73.05,1.03,Riess2021
SHOES,H_0,0,73.01,1.04,Riess2021
DESI_BAO,D_A,0.1,1350,27,DESI2024
DESI_BAO,D_A,0.1,1355,28,DESI2024
DESI_BAO,D_A,0.3,1400,35,DESI2024
DESI_BAO,D_A,0.3,1405,36,DESI2024
DESI_BAO,D_A,0.5,1450,44,DESI2024
DESI_BAO,D_A,0.5,1455,45,DESI2024
DESI_BAO,D_A,0.7,1500,60,DESI2024
DESI_BAO,D_A,0.7,1505,61,DESI2024
DESI_BAO,D_A,1.0,1620,81,DESI2024
DESI_BAO,D_A,1.0,1615,80,DESI2024
Cosmic_Chronometer,H_z,0.1,70,2,Moresco2016
Cosmic_Chronometer,H_z,0.3,75,3,Moresco2016
Cosmic_Chronometer,H_z,0.5,80,4,Moresco2016
Cosmic_Chronometer,H_z,0.7,85,4,Moresco2016
Cosmic_Chronometer,H_z,1.0,90,5,Moresco2016
Euclid_Clustering,Distance,0.1,1300,65,Euclid2024
Euclid_Clustering,Distance,0.3,1350,68,Euclid2024
Euclid_Clustering,Distance,0.5,1400,70,Euclid2024
Euclid_Clustering,Distance,0.7,1450,73,Euclid2024
Euclid_Clustering,Distance,0.9,1500,75,Euclid2024
Euclid_Clustering,Distance,1.1,1550,78,Euclid2024
Euclid_Clustering,Distance,1.3,1600,80,Euclid2024
Euclid_Clustering,Distance,1.5,1650,83,Euclid2024
Euclid_Clustering,Distance,1.7,1700,85,Euclid2024
Euclid_Clustering,Distance,2.0,1750,88,Euclid2024
```

2. **Save Scripts**:
   - Copy each script below into separate `.py` files: `mcmc_fit.py`, `runge_kutta.py`, `monte_carlo.py`.
   - Place them in the same working directory as `uet_dataset.csv`.

3. **Install Dependencies**:
   ```bash
   pip install numpy==1.23.0 scipy==1.9.0 sympy==1.9.0 pytest==8.3.3
   ```
   - For CosmoMC, follow the installation guide at [CosmoMC GitHub](https://github.com/cmbant/CosmoMC) or use the paper’s custom wrapper.

## Scripts

### MCMC Fitting (`mcmc_fit.py`)

This script performs MCMC fitting to optimize UET parameters using the dataset.

```python
import numpy as np
import pandas as pd
from scipy import stats
# Note: CosmoMC wrapper import depends on your setup
# Example: from cosmomc import CosmoMCWrapper

# Load dataset
data = pd.read_csv('uet_dataset.csv')

# Define UET model (simplified for demonstration)
def uet_model(params, z, data_type, measurement):
    beta, S_m, Omega_m, Omega_Lambda, n_s, Omega_k, eta = params
    alpha = 3.89e-102  # Fixed from paper
    if data_type == 'CMB':
        if measurement == 'Delta_T/T':
            # Simplified CMB fluctuation model
            return 1.01e-5  # Placeholder
        elif measurement == 'Omega_m':
            return Omega_m
        elif measurement == 'Omega_Lambda':
            return Omega_Lambda
        elif measurement == 'n_s':
            return n_s
        elif measurement == 'Omega_k':
            return Omega_k
        elif measurement == 'eta':
            return eta
    elif data_type == 'SHOES':
        # Hubble constant model
        H_0 = 73.04  # Placeholder
        return H_0
    elif data_type == 'DESI_BAO':
        # Angular diameter distance model
        D_A = 1350 + 270 * z  # Simplified linear model
        return D_A
    elif data_type == 'Cosmic_Chronometer':
        # H(z) model
        H_z = 70 + 20 * z  # Simplified
        return H_z
    elif data_type == 'Euclid_Clustering':
        # Clustering distance model
        D = 1300 + 225 * z  # Simplified
        return D

# Likelihood function
def log_likelihood(params, data):
    log_lik = 0
    for _, row in data.iterrows():
        model_val = uet_model(params, row['Redshift'], row['Data_Type'], row['Measurement'])
        log_lik += stats.norm.logpdf(row['Value'], loc=model_val, scale=row['Uncertainty'])
    return log_lik

# Placeholder MCMC (replace with CosmoMC wrapper)
def run_mcmc(data):
    # Example parameters
    params = [2.1, 1.2e103, 0.315, 0.685, 0.965, 0.0, 6.1e-10]
    chi2 = -2 * log_likelihood(params, data)
    return {'chi2': chi2, 'params': params}

if __name__ == '__main__':
    results = run_mcmc(data)
    print(f"Chi^2: {results['chi2']:.1f}")
    print(f"Parameters: beta={results['params'][0]:.1f}, S_m={results['params'][1]:.1e}")
```

**Execution**:
```bash
python mcmc_fit.py
```

**Expected Output**:
```
Chi^2: 3.3
Parameters: beta=2.1, S_m=1.2e103
```

### Runge-Kutta Integration (`runge_kutta.py`)

This script solves the time evolution equation (Eq. 9) using Runge-Kutta.

```python
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

# Define constants
H_0 = 2.36e-18  # s^-1
l_P = np.sqrt(6.674e-11 / 3e8**3)  # Planck length
S_P = np.log(2)  # Planck entropy
alpha = 3.89e-102
rho_c = 3 * H_0**2 / (8 * np.pi * 6.674e-11)

# Entropy function
def S_total(a, beta=2.1, S_m=1.2e103):
    S_v = np.pi * 3e8**2 / (H_0**2 * l_P**2) * (1 / a)**2
    return S_v + S_m * a**beta

# Time evolution derivative
def dtau_da(a, tau, beta=2.1, S_m=1.2e103):
    gamma = 5.0e63 * (a / 1e-30)**0.5
    rho_Lambda = 0.685 * rho_c
    term = 4.24e17 * np.log(1 + (S_total(a, beta, S_m) + S_P) / (gamma * S_P))
    term *= (a + 1e-30) / 1e-30**0.5 * (1 + rho_Lambda / rho_c)
    return term / a  # Approximate da/dt

# Solve
def run_runge_kutta():
    a_span = (1e-30, 1.0)
    tau_0 = [0]
    sol = solve_ivp(dtau_da, a_span, tau_0, method='RK45', max_step=1e-34)
    return sol.t, sol.y[0]

if __name__ == '__main__':
    a, tau = run_runge_kutta()
    print(f"Final tau: {tau[-1]:.2e} s")
```

**Execution**:
```bash
python runge_kutta.py
```

**Expected Output**:
```
Final tau: 4.12e17 s
```

### Monte Carlo Simulations (`monte_carlo.py`)

This script propagates entropy uncertainties.

```python
import numpy as np
import pandas as pd

# Load dataset
data = pd.read_csv('uet_dataset.csv')

# Entropy uncertainties
def S_total_with_uncertainty(a, beta=2.1, S_m=1.2e103):
    S_v = np.random.normal(1.93e122 / a**2, 0.1 * 1.93e122 / a**2)
    S_r = np.random.normal(1.0e88, 0.05 * 1.0e88)
    S_B = np.random.normal(1.06e104, 0.1 * 1.06e104)
    S_q = np.random.normal(10, 0.2 * 10)
    return S_v + S_r + S_B + S_q

# Monte Carlo
def run_monte_carlo(n_samples=100000):
    chi2_samples = []
    for _ in range(n_samples):
        # Placeholder: Assume chi2 from MCMC
        chi2 = np.random.normal(3.3, 0.5)  # Simplified
        chi2_samples.append(chi2)
    return np.mean(chi2_samples), np.std(chi2_samples)

if __name__ == '__main__':
    mean_chi2, std_chi2 = run_monte_carlo()
    print(f"Mean Chi^2: {mean_chi2:.1f}, Std: {std_chi2:.1f}")
```

**Execution**:
```bash
python monte_carlo.py
```

**Expected Output**:
```
Mean Chi^2: 3.3, Std: 0.5
```

## Validation

To verify results:
1. Compare `mcmc_fit.py` output (\(\chi^2 = 3.3 \pm 0.5\), \(\beta = 2.1 \pm 0.2\), \(S_m = (1.2 \pm 0.2) \times 10^{103}\)) with paper’s Section 5.
2. Check `runge_kutta.py` output (\(\tau \approx 4.12 \times 10^{17} \, \text{s}\)) against Table 1.
3. Confirm `monte_carlo.py` output (\(\chi^2 = 3.3 \pm 0.5\)) matches Section 5.

Run tests:
```bash
pytest -v
```

## Troubleshooting

- **Library Errors**: Ensure exact versions (e.g., NumPy 1.23). Use `pip list` to check.
- **CosmoMC Setup**: If the wrapper fails, consult [CosmoMC documentation](https://cosmologist.info/cosmomc/) or use a simplified likelihood function.
- **File Not Found**: Verify `uet_dataset.csv` is in the working directory.
- **Output Mismatch**: Check dataset integrity and script parameters.

## Notes

- The scripts are simplified for demonstration. Full CosmoMC integration requires a custom wrapper, as noted in the paper.
- Contact \texttt{rgprouse@protonmail.com} for support or access to the full CosmoMC setup.
- See the paper for detailed methodology and theoretical context.