import osqp
import os
import sys
import numpy as np
import scipy.io as spio
import gurobipy as gp
import time

# Global configuration constants
PATHS           = None
EPS_ABS         = 1e-6
EPS_REL         = 1e-6
EPS_PRIM_INF    = 1e-9
EPS_DUAL_INF    = 1e-9

def determine_prob_date(path):

    # Define problem data
    name = os.path.basename(path)[:-4]

    if path.lower().endswith('.mat'):
        P, q, A, l, u, n, m = determine_prob_date_mat(path)
    else:
        print(f"Error the file is not a .mat file")
        sys.exit(1)
    assert A.shape == (m, n)
    
    # Infinity constant is 1e20
    A[A > +9e19] = +np.inf
    l[l > +9e19] = +np.inf
    u[u > +9e19] = +np.inf
    A[A < -9e19] = -np.inf
    l[l < -9e19] = -np.inf
    u[u < -9e19] = -np.inf
        
    return P, q, A, l, u, n, m, name

def determine_prob_date_mat(path):
    try:
        mat_dict = spio.loadmat(path)
    except Exception as e:
        print(f"Error loading problem data: {e}")
        sys.exit(1)
        
    P = mat_dict["P"].astype(float).tocsc()
    q = mat_dict["q"].T.flatten().astype(float)
    A = mat_dict["A"].astype(float).tocsc()
    l = mat_dict["l"].T.flatten().astype(float)
    u = mat_dict["u"].T.flatten().astype(float)
    
    n = mat_dict["n"].T.flatten().astype(int)[0]
    m = mat_dict["m"].T.flatten().astype(int)[0]
    
    return P, q, A, l, u, n, m


def run_osqp():
    
    prob = osqp.OSQP()

    # Problem Data
    P, q, A, l, u, n, m, name = determine_prob_date(PATHS)
            
    prob.setup(P, q, A, l, u, verbose=True,
               eps_abs=EPS_ABS, eps_rel=EPS_REL, eps_prim_inf=EPS_PRIM_INF, eps_dual_inf=EPS_DUAL_INF,
               time_limit=300, max_iter=2000000000, linsys_solver='mkl pardiso')
    
    res = prob.solve(raise_error=False)
    
    print('Problem name:', name, flush=True)
    print('Run time:', res.info.run_time, flush=True)
    print('Status:', res.info.status, flush=True)
    
    return res


if __name__ == '__main__':
    
    PATHS = '/scratch/gpfs/lb3825/Mittelman_mat/L1_sixm250obs_Mittelman.mat'
    # PATHS = 'L1_sixm1000obs_Mittelman.mat'
    
    results = run_osqp()