import osqp
import os
import sys
import numpy as np
from scipy import sparse
import scipy.io as spio
import scipy.sparse as spa
import csv
import pandas as pd
import random
import time
import argparse

def generate_random_qp(n, m, p_scale=0.1, p_rank=None, seed=42):
    """
    Generate a random QP problem of the form:
        min 0.5 x^T P x + q^T x
        s.t. l <= A x <= u

    Parameters
    ----------
    n : int
        Number of variables
    m : int
        Number of constraints
    p_scale : float
        Scaling factor for P matrix
    p_rank : int or None
        Rank of P matrix (if None, full rank is used)
    seed : int
        Random seed for reproducibility

    Returns
    -------
    P : scipy.sparse.csc_matrix
        Quadratic cost matrix (n x n), positive semidefinite
    q : numpy.ndarray
        Linear cost vector (n)
    A : scipy.sparse.csc_matrix
        Constraint matrix (m x n)
    l : numpy.ndarray
        Lower bound vector (m)
    u : numpy.ndarray
        Upper bound vector (m)
    """
    # Set random seed for reproducibility
    np.random.seed(seed)
    random.seed(seed)

    # Generate random data for the problem
    z = np.random.randn(m)
    y = np.maximum(z, 0)  # Element-wise maximum between z and 0
    s = y - z

    P = np.random.randn(n, n)
    P = p_scale * P.T @ P

    # Make P low rank if p_rank is set
    if p_rank is not None:
        eigs, V = np.linalg.eig(P)
        eigs[p_rank:] = 0
        P = (V * eigs) @ V.T

    P = sparse.csc_matrix(P)

    # Make problem slightly more numerically challenging:
    A = np.random.randn(m, n)
    U, S, V = np.linalg.svd(A, full_matrices=False)
    S = S**2
    S /= np.max(S)
    A = (U * S) @ V
    A = sparse.csc_matrix(A)

    x = np.random.randn(n)
    q = -A.T @ y - P @ x
    u = A.dot(x) + s
    
    u_mag = np.linalg.norm(u)

    u /= u_mag
    # u = np.inf * np.ones(len(u))
    # q /= np.linalg.norm(q)
    # l = -np.inf * np.ones(len(u))
    l = u - u_mag * np.abs(np.random.randn(u.shape[0]))
    # l = np.zeros(len(u))

    return P, q, A, l, u

def halpern_combinations():
    """
    Generates all of the possible combinations of paramters for 
    osqp algorithm with restart_type = "halpern". All combinations
    are stored in an array as dictionaries
    """
    # Define array of dictionaries
    param_dicts = []
    
    # Start timer
    start = time.time()
    
    # Parameters
    # alpha_vals = np.linspace(1, 2, 11)[:-1]
    alpha_vals = [1., 1.4, 1.6, 1.8]
    # rho_is_vec_vals = [0, 1]
    # Omit 0
    adaptive_rho_tolerance_vals = np.linspace(0, 10, 5)[1:]
    rho_custom_condition_vals = [0, 1]
    # Omit 0
    rho_custom_tolerance_vals = np.linspace(0, 10, 5)[1:]
    restart_type = "halpern"
    adapt_rho_on_restart_vals = [0, 1]
    beta_vals = np.linspace(0, 1, 11)[1:-1]
    ini_rest_len_vals = np.arange(10, 50, 5)
    halpern_scheme_vals = ["none", "adaptive", "adaptive only before init_rest_len"]
    adaptive_rest_vals = [0, 1]
    # Omit 0
    restart_necessary_vals = np.linspace(0, 1, 5)[1:]
    restart_artificial_vals = np.linspace(0, 1, 5)[1:]
    halpern_step_first_inner_iter_vals = [0, 1]
    
    for alpha in alpha_vals:
        for adaptive_rho_tolerance_greater in adaptive_rho_tolerance_vals:
            for adaptive_rho_tolerance_less in adaptive_rho_tolerance_vals:
                for adapt_rho_on_restart in adapt_rho_on_restart_vals:
                    for beta in beta_vals:
                        for ini_rest_len in ini_rest_len_vals:
                            for halpern_scheme in halpern_scheme_vals:
                                for halpern_step_first_inner_iter in halpern_step_first_inner_iter_vals:
                                    for rho_custom_condition in rho_custom_condition_vals:
                                        if rho_custom_condition == 0:
                                            for adaptive_rest in adaptive_rest_vals:
                                                if adaptive_rest == 0:
                                                    param_dicts.append({
                                                                "alpha": alpha,
                                                                "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                "beta": beta,
                                                                "ini_rest_len": ini_rest_len,
                                                                "halpern_scheme": halpern_scheme,
                                                                "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                "rho_custom_condition": rho_custom_condition,
                                                                "adaptive_rest": adaptive_rest,
                                                                "restart_type": restart_type
                                                            })
                                                else:
                                                    for restart_necessary in restart_necessary_vals:
                                                        for restart_artificial in restart_artificial_vals:
                                                            param_dicts.append({
                                                                "alpha": alpha,
                                                                "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                "beta": beta,
                                                                "ini_rest_len": ini_rest_len,
                                                                "halpern_scheme": halpern_scheme,
                                                                "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                "rho_custom_condition": rho_custom_condition,
                                                                "adaptive_rest": adaptive_rest,
                                                                "restart_necessary": restart_necessary,
                                                                "restart_artificial": restart_artificial,
                                                                "restart_type": restart_type
                                                            })
                                        else:
                                            for rho_custom_tolerance in rho_custom_tolerance_vals:
                                                for adaptive_rest in adaptive_rest_vals:
                                                    if adaptive_rest == 0:
                                                        param_dicts.append({
                                                                    "alpha": alpha,
                                                                    "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                    "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                    "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                    "beta": beta,
                                                                    "ini_rest_len": ini_rest_len,
                                                                    "halpern_scheme": halpern_scheme,
                                                                    "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                    "rho_custom_condition": rho_custom_condition,
                                                                    "rho_custom_tolerance": rho_custom_tolerance,
                                                                    "adaptive_rest": adaptive_rest,
                                                                    "restart_type": restart_type
                                                                })
                                                    else:
                                                        for restart_necessary in restart_necessary_vals:
                                                            for restart_artificial in restart_artificial_vals:
                                                                param_dicts.append({
                                                                    "alpha": alpha,
                                                                    "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                    "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                    "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                    "beta": beta,
                                                                    "ini_rest_len": ini_rest_len,
                                                                    "halpern_scheme": halpern_scheme,
                                                                    "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                    "rho_custom_condition": rho_custom_condition,
                                                                    "rho_custom_tolerance": rho_custom_tolerance,
                                                                    "adaptive_rest": adaptive_rest,
                                                                    "restart_necessary": restart_necessary,
                                                                    "restart_artificial": restart_artificial,
                                                                    "restart_type": restart_type
                                                                })
    end = time.time()
    print(f"Time it took to generate all of the combinations for halpern is: {end - start} seconds")
    print(f"Generated {len(param_dicts)} parameter dictionaries")
    
    return param_dicts


def reflected_halpern_combinations():
    """
    Generates all of the possible combinations of paramters for 
    osqp algorithm with restart_type = "reflected halpern". All combinations
    are stored in an array as dictionaries. This requires over
    10 Gb of memory
    """
    # Define array of dictionaries
    param_dicts = []
    
    # Start timer
    start = time.time()
    
    # Parameters
    # alpha_vals = np.linspace(1, 2, 11)[:-1]
    alpha_vals = [1., 1.4, 1.6, 1.8]
    # rho_is_vec_vals = [0, 1]
    # Omit 0
    adaptive_rho_tolerance_vals = np.linspace(0, 10, 5)[1:]
    rho_custom_condition_vals = [0, 1]
    # Omit 0
    rho_custom_tolerance_vals = np.linspace(0, 10, 5)[1:]
    restart_type = "reflected halpern"
    adapt_rho_on_restart_vals = [0, 1]
    beta_vals = np.linspace(0, 1, 11)[1:-1]
    ini_rest_len_vals = np.arange(10, 50, 5)
    halpern_scheme_vals = ["none", "adaptive", "adaptive only before init_rest_len"]
    # alpha_adjustment_reflected_halpern_vals = [0, 1]
    alpha_adjustment_reflected_halpern_vals = [1]
    adaptive_rest_vals = [0, 1]
    # Omit 0
    restart_necessary_vals = np.linspace(0, 1, 5)[1:]
    restart_artificial_vals = np.linspace(0, 1, 5)[1:]
    halpern_step_first_inner_iter_vals = [0, 1]
    


    for alpha in alpha_vals:
        for adaptive_rho_tolerance_greater in adaptive_rho_tolerance_vals:
            for adaptive_rho_tolerance_less in adaptive_rho_tolerance_vals:
                for adapt_rho_on_restart in adapt_rho_on_restart_vals:
                    for beta in beta_vals:
                        for ini_rest_len in ini_rest_len_vals:
                            for halpern_scheme in halpern_scheme_vals:
                                for halpern_step_first_inner_iter in halpern_step_first_inner_iter_vals:
                                    for lambd in np.linspace(0.0, (2.0 / alpha) - 1.0, 5):
                                        for alpha_adjustment_reflected_halpern in alpha_adjustment_reflected_halpern_vals:
                                            for rho_custom_condition in rho_custom_condition_vals:
                                                if rho_custom_condition == 0:
                                                    for adaptive_rest in adaptive_rest_vals:
                                                        if adaptive_rest == 0:
                                                            param_dicts.append({
                                                                        "alpha": alpha,
                                                                        "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                        "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                        "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                        "beta": beta,
                                                                        "ini_rest_len": ini_rest_len,
                                                                        "halpern_scheme": halpern_scheme,
                                                                        "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                        "rho_custom_condition": rho_custom_condition,
                                                                        "adaptive_rest": adaptive_rest,
                                                                        "restart_type": restart_type,
                                                                        "alpha_adjustment_reflected_halpern": alpha_adjustment_reflected_halpern,
                                                                        "lambd": lambd
                                                                    })
                                                        else:
                                                            for restart_necessary in restart_necessary_vals:
                                                                for restart_artificial in restart_artificial_vals:
                                                                    param_dicts.append({
                                                                        "alpha": alpha,
                                                                        "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                        "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                        "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                        "beta": beta,
                                                                        "ini_rest_len": ini_rest_len,
                                                                        "halpern_scheme": halpern_scheme,
                                                                        "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                        "rho_custom_condition": rho_custom_condition,
                                                                        "adaptive_rest": adaptive_rest,
                                                                        "restart_necessary": restart_necessary,
                                                                        "restart_artificial": restart_artificial,
                                                                        "restart_type": restart_type,
                                                                        "alpha_adjustment_reflected_halpern": alpha_adjustment_reflected_halpern,
                                                                        "lambd": lambd
                                                                    })
                                                else:
                                                    for rho_custom_tolerance in rho_custom_tolerance_vals:
                                                        for adaptive_rest in adaptive_rest_vals:
                                                            if adaptive_rest == 0:
                                                                param_dicts.append({
                                                                            "alpha": alpha,
                                                                            "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                            "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                            "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                            "beta": beta,
                                                                            "ini_rest_len": ini_rest_len,
                                                                            "halpern_scheme": halpern_scheme,
                                                                            "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                            "rho_custom_condition": rho_custom_condition,
                                                                            "rho_custom_tolerance": rho_custom_tolerance,
                                                                            "adaptive_rest": adaptive_rest,
                                                                            "restart_type": restart_type,
                                                                            "alpha_adjustment_reflected_halpern": alpha_adjustment_reflected_halpern,
                                                                            "lambd": lambd
                                                                        })
                                                            else:
                                                                for restart_necessary in restart_necessary_vals:
                                                                    for restart_artificial in restart_artificial_vals:
                                                                        param_dicts.append({
                                                                            "alpha": alpha,
                                                                            "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                            "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                            "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                            "beta": beta,
                                                                            "ini_rest_len": ini_rest_len,
                                                                            "halpern_scheme": halpern_scheme,
                                                                            "halpern_step_first_inner_iter": halpern_step_first_inner_iter,
                                                                            "rho_custom_condition": rho_custom_condition,
                                                                            "rho_custom_tolerance": rho_custom_tolerance,
                                                                            "adaptive_rest": adaptive_rest,
                                                                            "restart_necessary": restart_necessary,
                                                                            "restart_artificial": restart_artificial,
                                                                            "restart_type": restart_type,
                                                                            "alpha_adjustment_reflected_halpern": alpha_adjustment_reflected_halpern,
                                                                            "lambd": lambd
                                                                        })
                                                                        
    end = time.time()
    print(f"Time it took to generate all of the combinations for reflected halpern is: {end - start} seconds")
    print(f"Generated {len(param_dicts)} parameter dictionaries")
    
    return param_dicts



def averaged_combinations():
    """
    Generates all of the possible combinations of paramters for 
    osqp algorithm with restart_type = "averaged". All combinations
    are stored in an array as dictionaries. This requires over
    10 Gb of memory
    """
    # Define array of dictionaries
    param_dicts = []
    
    # Start timer
    start = time.time()
    
    # Parameters
    # alpha_vals = np.linspace(1, 2, 11)[:-1]
    alpha_vals = [1., 1.4, 1.6, 1.8]
    # rho_is_vec_vals = [0, 1]
    # Omit 0
    adaptive_rho_tolerance_vals = np.linspace(0, 10, 5)[1:]
    rho_custom_condition_vals = [0, 1]
    # Omit 0
    rho_custom_tolerance_vals = np.linspace(0, 10, 5)[1:]
    restart_type = "averaged"
    adapt_rho_on_restart_vals = [0, 1]
    custom_average_rest_vals = [0, 1]
    beta_vals = np.linspace(0, 1, 11)[1:-1]
    ini_rest_len_vals = np.arange(10, 50, 5)
    xi_vals = np.linspace(0, 2, 11)[1:]
    vector_rho_in_averaged_KKT_vals = [0, 1]
    # halpern_scheme_vals = ["none", "adaptive", "adaptive only before init_rest_len"]
    # # alpha_adjustment_reflected_halpern_vals = [0, 1]
    # alpha_adjustment_reflected_halpern_vals = [1]
    adaptive_rest_vals = [0, 1]
    # Omit 0
    restart_necessary_vals = np.linspace(0, 1, 5)[1:]
    restart_artificial_vals = np.linspace(0, 1, 5)[1:]
    # halpern_step_first_inner_iter_vals = [0, 1]
    
    
    
    for alpha in alpha_vals:
        for adaptive_rho_tolerance_greater in adaptive_rho_tolerance_vals:
            for adaptive_rho_tolerance_less in adaptive_rho_tolerance_vals:
                for adapt_rho_on_restart in adapt_rho_on_restart_vals:
                    for custom_average_rest in custom_average_rest_vals:
                        for beta in beta_vals:
                            for ini_rest_len in ini_rest_len_vals:
                                for xi in xi_vals:
                                    for vector_rho_in_averaged_KKT in vector_rho_in_averaged_KKT_vals:
                                        for rho_custom_condition in rho_custom_condition_vals:
                                            if rho_custom_condition == 0:
                                                for adaptive_rest in adaptive_rest_vals:
                                                    if adaptive_rest == 0:
                                                        param_dicts.append({
                                                                    "alpha": alpha,
                                                                    "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                    "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                    "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                    "custom_average_rest": custom_average_rest,
                                                                    "beta": beta,
                                                                    "ini_rest_len": ini_rest_len,
                                                                    "xi": xi,
                                                                    "rho_custom_condition": rho_custom_condition,
                                                                    "adaptive_rest": adaptive_rest,
                                                                    "restart_type": restart_type,
                                                                    "vector_rho_in_averaged_KKT": vector_rho_in_averaged_KKT
                                                                })
                                                    else:
                                                        for restart_necessary in restart_necessary_vals:
                                                            for restart_artificial in restart_artificial_vals:
                                                                param_dicts.append({
                                                                    "alpha": alpha,
                                                                    "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                    "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                    "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                    "custom_average_rest": custom_average_rest,
                                                                    "beta": beta,
                                                                    "ini_rest_len": ini_rest_len,
                                                                    "xi": xi,
                                                                    "rho_custom_condition": rho_custom_condition,
                                                                    "adaptive_rest": adaptive_rest,
                                                                    "restart_necessary": restart_necessary,
                                                                    "restart_artificial": restart_artificial,
                                                                    "restart_type": restart_type,
                                                                    "vector_rho_in_averaged_KKT": vector_rho_in_averaged_KKT
                                                                })
                                            else:
                                                for rho_custom_tolerance in rho_custom_tolerance_vals:
                                                    for adaptive_rest in adaptive_rest_vals:
                                                        if adaptive_rest == 0:
                                                            param_dicts.append({
                                                                        "alpha": alpha,
                                                                        "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                        "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                        "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                        "custom_average_rest": custom_average_rest,
                                                                        "beta": beta,
                                                                        "ini_rest_len": ini_rest_len,
                                                                        "xi": xi,
                                                                        "vector_rho_in_averaged_KKT": vector_rho_in_averaged_KKT,
                                                                        "rho_custom_condition": rho_custom_condition,
                                                                        "rho_custom_tolerance": rho_custom_tolerance,
                                                                        "adaptive_rest": adaptive_rest,
                                                                        "restart_type": restart_type
                                                                    })
                                                        else:
                                                            for restart_necessary in restart_necessary_vals:
                                                                for restart_artificial in restart_artificial_vals:
                                                                    param_dicts.append({
                                                                        "alpha": alpha,
                                                                        "adaptive_rho_tolerance_greater": adaptive_rho_tolerance_greater,
                                                                        "adaptive_rho_tolerance_less": adaptive_rho_tolerance_less,
                                                                        "adapt_rho_on_restart": adapt_rho_on_restart,
                                                                        "custom_average_rest": custom_average_rest,
                                                                        "beta": beta,
                                                                        "ini_rest_len": ini_rest_len,
                                                                        "xi": xi,
                                                                        "vector_rho_in_averaged_KKT": vector_rho_in_averaged_KKT,
                                                                        "rho_custom_condition": rho_custom_condition,
                                                                        "rho_custom_tolerance": rho_custom_tolerance,
                                                                        "adaptive_rest": adaptive_rest,
                                                                        "restart_necessary": restart_necessary,
                                                                        "restart_artificial": restart_artificial,
                                                                        "restart_type": restart_type
                                                                    })



    end = time.time()
    print(f"Time it took to generate all of the combinations for averaged is: {end - start} seconds")
    print(f"Generated {len(param_dicts)} parameter dictionaries")
    
    return param_dicts



if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print("script_dir =", script_dir)
    
    # Check for command line argument
    # if len(sys.argv) != 2 and sys.argv[1] != "random":
    #     print("Usages:")
    #     print("Existing Problem: py maros_meszaros_benchmarks.py <problem_name>")
    #     print("Random Problem: py maros_meszaros_benchmarks.py random")
    #     print("Example: py maros_meszaros_benchmarks.py HS21")
    #     print("Available problems are .mat files in qpbenchmark/maros_meszaros_qpbenchmark/data/")
    #     sys.exit(1)
    
    parser = argparse.ArgumentParser(description="Run Maros-Meszaros benchmarks.")
    parser.add_argument("problem_name", type=str, help="Name of the problem (e.g., HS21 or 'random').")
    parser.add_argument("--n", type=int, default=100, help="Number of variables for random QP generation.")
    parser.add_argument("--seed", type=int, default=1, help="Random seed for reproducibility.")
    parser.add_argument("--p_scale", type=float, default=0.01, help="Scaling factor for P matrix.")
    parser.add_argument("--p_rank", type=int, default=25, help="Rank of P matrix (None for full rank).")
    parser.add_argument("--plot", type=int, default=1, help="Parameter to pass to the osqp solver in order for it to temporary store each iteration for ploting purposes")
    parser.add_argument("--save_stat", type=int, default=0, help="Save the setup time, solve time, run time, primal residual, dual residual, duality gap, and the number of restarts. This is a binary either 0 or 1.")
    parser.add_argument("--save_all", type=int, default=0, help="Store all of the iterations for each problem ran (For usecase where problem_name = 'ALL'). This is a binary either 0 or 1.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Access arguments
    problem_name = args.problem_name.upper()
    n = args.n
    m = int(n * 1.5)
    seed = args.seed
    p_scale = args.p_scale
    p_rank = args.p_rank
    plot = args.plot
    save_stats = args.save_stat
    save_all = args.save_all
    
    if (save_stats != 0) and (save_stats != 1):
        print(f"save_stats must be 0 or 1")
        sys.exit(1)
        
    if (save_all != 0) and (save_all != 1):
        print(f"save_all must be 0 or 1")
        sys.exit(1)
    
    
    print(f"problem_name {problem_name}")
    if problem_name == "RANDOM":
        print(f"n: {n}, m: {m}, seed: {seed}, p_scale: {p_scale}, p_rank: {p_rank}")
    
    # Construct path to the qpbenchmark data directory
    # Navigate from osqp-python/examples to qpbenchmark/maros_meszaros_qpbenchmark/data
    qpbenchmark_data_dir = os.path.join(script_dir, "..", "..", "qpbenchmark", "maros_meszaros_qpbenchmark", "data")
    qpbenchmark_data_dir = os.path.abspath(qpbenchmark_data_dir)
    
    if problem_name != "RANDOM":
        # Add .mat extension if not provided
        if not problem_name.endswith('.mat'):
            problem_name += '.mat'
            
        if (problem_name == "ALL.mat"):
            files = [f for f in os.listdir(qpbenchmark_data_dir) if f.endswith('.mat')]
            paths = [os.path.join(qpbenchmark_data_dir, f) for f in files]
        
        else:
            problem_path = os.path.join(qpbenchmark_data_dir, problem_name)
            
            # Check if the problem file exists
            if not os.path.exists(problem_path):
                print(f"Error: Problem file '{problem_path}' not found.")
                print(f"Available problems in '{qpbenchmark_data_dir}':")
                try:
                    files = [f for f in os.listdir(qpbenchmark_data_dir) if f.endswith('.mat')]
                    for f in sorted(files):
                        print(f"  {f[:-4]}")  # Remove .mat extension for display
                except FileNotFoundError:
                    print(f"  Directory not found: {qpbenchmark_data_dir}")
                sys.exit(1)
            
            paths = [problem_path]
    
        print(f"Running problem: {problem_name}")
    else:
        # Paths is not important for randomly generated problems
        paths = ["random"]
        print("Running randomly generated problem")

    for path in paths:
        if problem_name != "RANDOM":
            # Define problem data
            name = os.path.basename(path)[:-4]
            print(f"name: {name}")
            print(f"Solving problem: {name}")

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
            assert A.shape == (m, n)
            
            # Infinity constant is 1e20
            A[A > +9e19] = +np.inf
            l[l > +9e19] = +np.inf
            u[u > +9e19] = +np.inf
            A[A < -9e19] = -np.inf
            l[l < -9e19] = -np.inf
            u[u < -9e19] = -np.inf
        else:
            name = "Random Problem"
            print(f"name: {name}")
            print(f"Solving problem: {name}")
            P, q, A, l, u = generate_random_qp(
                n,
                m,
                p_scale=p_scale,
                p_rank=p_rank,
                seed=seed,
            )
            
        
        # Create an OSQP object
        prob = osqp.OSQP()
        
        # Generate the set of all possible combinations
        #       This requires a lot of free memory on the system 
        #       (my pc can't load all of these combinations, roughly 60 million comb.)
        
        halpern_combination_array = halpern_combinations()
        # reflected_halpern_combination_array = reflected_halpern_combinations()
        # averaged_combination_array = averaged_combinations()
          
        all_param_dict = {
            "restart_type": "-1",
            "alpha": -1,
            "rho_is_vec": -1,
            "adaptive_rho_tolerance_greater": -1,
            "adaptive_rho_tolerance_less": -1,
            "rho_custom_condition": -1,
            "rho_custom_tolerance": -1,
            "adapt_rho_on_restart": -1,
            "custom_average_rest": -1,
            "beta": -1,
            "ini_rest_len": -1,
            "lambd": -1,
            "xi": -1,
            "vector_rho_in_averaged_KKT": -1,
            "halpern_scheme": "-1",
            "alpha_adjustment_reflected_halpern": -1,
            "adaptive_rest": -1,
            "restart_necessary": -1,
            "restart_artificial": -1,
            "halpern_step_first_inner_iter": -1
        }
        

        # Setup workspace
        """
        All of the parameters
                                                                
        Restart Types:
        none
        halpern
        reflected halpern
        averaged
        
        Halpern Schemes:
        none
        adaptive
        adaptive only before init_rest_len
        
        Parameter values:
        alpha                               (0, 2)
        beta                                (0, 1]
        lambd                               [0, 1]
        restart_necessary                   (0, 1]
        restart_artificial                  (0, 1]
        ini_rest_len                        any # > 0
        # adaptive_rho_tolerance            any # > 1
        adaptive_rho_tolerance_greater      any # > 0
        adaptive_rho_tolerance_less         any # > 0
        adaptive_rest                       0 or 1, either off or on
        alpha_adjustment_reflected_halpern  0 or 1, either off or on
        rho_custom_condition                0 or 1, either off or on
        rho_custom_tolerance                any # > 0
        custom_average_rest                 0 or 1, either off or on
        vector_rho_in_averaged_KKT          0 or 1, either off or on
        adapt_rho_on_restart                0 or 1, either off or on
        xi                                  any # > 0
        halpern_step_first_inner_iter       0 or 1, either off or on
        """
        
        # prob.setup(P, q, A, l, u, alpha=1.6, verbose=True, 
        #            beta=0.1, lambd=0.2, restart_necessary=0.3, 
        #            restart_artificial=0.4, ini_rest_len=5,
        #            restart_type="none", halpern_scheme="none",
        #            adaptive_rho_tolerance_greater=0.6,
        #            adaptive_rho_tolerance_less=0.7, adaptive_rest=1,
        #            alpha_adjustment_reflected_halpern=0, rho_custom_condition=1,
        #            custom_average_rest=0, vector_rho_in_averaged_KKT=0,
        #            adapt_rho_on_restart=0, rho_custom_tolerance=0.8, xi=1.0
        #            scaling=0)
        
        alpha_val=1.6
        
        # prob.setup(P, q, A, l, u, plot=1, alpha=alpha_val, lambd=max(min((2. / alpha_val) - 1., 1), 0),
        #            restart_type="reflected halpern", alpha_adjustment_reflected_halpern=1, adapt_rho_on_restart=1,
        #            halpern_scheme="none", beta=0.7, max_iter=250000, rho_is_vec=1, rho_custom_condition=0,
        #            rho_custom_tolerance=0.5, adaptive_rho_tolerance_greater=5, adaptive_rho_tolerance_less=5,
        #            eps_abs=1e-6, eps_rel=1e-6)
                                            
        
        current_set = {
            'alpha': alpha_val,
            'rho_is_vec': 1,
            'adaptive_rho_tolerance_greater': 5,
            'adaptive_rho_tolerance_less': 5,
            'rho_custom_condition': 0,
            'rho_custom_tolerance': 0.05,
            'restart_type': "halpern",
            'adapt_rho_on_restart': 1,
            'custom_average_rest': 1,
            'beta': 0.8,
            'ini_rest_len': 5,
            'xi': 0.2,
            'adaptive_rest': 1,
            'restart_necessary': 0.9,
            'restart_artificial': 0.7,
            'vector_rho_in_averaged_KKT': 1
        }
        
        # for i in range(len(halpern_combination_array)):
        for i in range(1):
            current_comb = halpern_combination_array[i]
            # prob.setup(P, q, A, l, u, plot=plot, scaling=0, eps_abs=1e-6, eps_rel=0, max_iter=25000, **current_comb)
            prob.setup(P, q, A, l, u, plot=plot, scaling=0, eps_abs=1e-6, eps_rel=0, max_iter=25000, **current_set)
            
            
            # Settings can be changed using .update_settings()
            # prob.update_settings(polishing=1)

            # Solve problem
            res = prob.solve(raise_error=False)
            
            print('Setup time:', res.info.setup_time)
            print('Solve time:', res.info.solve_time)
            print('Run time:', res.info.run_time)
            print('Primal Residual:', res.info.prim_res)
            print('Dual Residual:', res.info.dual_res)
            print('Duality Gap:', res.info.duality_gap)
            print('Number of restarts:', res.info.restart)
            
            if save_stats:
                stat_file = "../../osqp/plot/stats.csv"
                
                # Solution details
                stats_data = {
                    "problem name": problem_name,
                    "setup time": res.info.setup_time,
                    "solve time": res.info.solve_time,
                    "run time": res.info.run_time,
                    "primal residual": res.info.prim_res,
                    "dual residual": res.info.dual_res,
                    "duality gap": res.info.duality_gap,
                    "restart": res.info.restart
                }
                
                # Adding paramters
                halpern_comb_keys = list(halpern_combination_array[i].keys())
                all_param_dict_keys = list(all_param_dict.keys())
                
                for key in all_param_dict_keys:
                    stats_data[key] = halpern_combination_array[i][key] if key in halpern_comb_keys else all_param_dict[key]
                
                # Writing to .csv file and adding header if needed
                file_exists = os.path.exists(stat_file)
                file_is_empty = os.path.getsize(stat_file) == 0 if file_exists else True
                with open(stat_file, mode='a', newline='', encoding='utf8') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=stats_data.keys())
                    if not file_exists or file_is_empty:
                        writer.writeheader()
                    writer.writerow(stats_data)
                
                
                
            # Transfer the iterates values for the current problem from the .csv file generated by 
            #   the osqp sovler to a master .csv file that stores the iterates for all problems
            if (problem_name == "ALL.mat") and (save_all == 1):
                residuals_file = "../../osqp/plot/residuals.csv"
                # master_file = "../../osqp/plot/osqp_reflected_halpern_residuals.csv"
                master_file = "../../osqp/plot/osqp_default_residuals.csv"
                
                df = pd.read_csv(residuals_file, header=0)
                df.insert(1, 'prob_name', name)
                
                if os.path.exists(master_file):
                    df.to_csv(master_file, mode='a', header=True, index=False, float_format="%.3e")
                    # df.to_csv(master_file, mode='a', header=True, index=False)