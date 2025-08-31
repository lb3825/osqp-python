import osqp
import os
import sys
import numpy as np
import numpy.linalg as la
from scipy import sparse
import scipy.io as spio
import scipy.sparse as spa
import csv
import pandas as pd
import random
import time
import argparse
import multiprocessing
from functools import partial
import smtplib
from email.message import EmailMessage
import optuna
from optuna.storages import JournalStorage
from optuna.storages.journal import JournalFileBackend
import traceback
import gurobipy as gp
from gurobipy import GRB
import cvxpy as cp
# import scs

# Global configuration constants
START_TIMER = None
# Set to 3 hours
JOB_TIME_LIMIT = 3 * 60 * 60
TIME_THRESHOLD = 2500
PATHS = None
PROBLEM_NAME = None
N = None
M = None
P_SCALE = None
P_RANK = None
SEED = None
PLOT = None
SAVE_STATS = None
SAVE_ALL = None
ALL_PARAM_DICT = None
CLUSTER = None
PATH_LESS_200 = [
    "AUG2D.mat", "AUG2DC.mat", "AUG3D.mat", "AUG3DC.mat",
    "AUG3DCQP.mat", "AUG3DQP.mat","CVXQP2_M.mat", "CVXQP2_S.mat",
    "CVXQP3_S.mat", "DPKLO1.mat", "DUAL1.mat", "DUAL2.mat",
    "DUAL3.mat", "DUAL4.mat", "DUALC8.mat", "GENHS28.mat",
    "GOULDQP3.mat", "HS21.mat", "HS35.mat", "HS35MOD.mat",
    "HS51.mat", "HS52.mat", "HS53.mat", "HS76.mat", "LASER.mat",
    "LOTSCHD.mat", "PRIMAL1.mat", "PRIMAL2.mat", "PRIMAL3.mat",
    "PRIMAL4.mat", "QAFIRO.mat", "QPTEST.mat", "QSC205.mat",
    "STCQP1.mat", "STCQP2.mat", "TAME.mat"
]
EPS_ABS         = 1e-6
EPS_REL         = 1e-6
EPS_PRIM_INF    = 1e-6
EPS_DUAL_INF    = 1e-6


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

def geom_mean(t, shift=10.0):
    """Compute the shifted geometric mean using formula from
    http://plato.asu.edu/ftp/shgeom.html

    NB. Use logarithms to avoid numeric overflows
    """
    return np.exp(np.sum(np.log(np.maximum(1, t + shift))/len(t))) - shift

def send_email(subject, body):
    # Write an email message when the work is finished    
    msg = EmailMessage()
    msg['Subject'] = subject
    msg['From'] = 'emailsender062@gmail.com'
    msg['To'] = 'lb3825@princeton.edu'
    msg.set_content(body)

    # Gmail SMTP server setup
    smtp_server = 'smtp.gmail.com'
    smtp_port = 587

    # Sending the email
    with smtplib.SMTP(smtp_server, smtp_port) as server:
        server.starttls()
        server.login('emailsender062@gmail.com', 'qull hdcu ivlb hubu')
        server.send_message(msg)
        
def accuracy_ratios(prob_name, x, y):
    """
    Computes the current primal residual, dual residual, and duality gap
    relative to the values neeeded for optimality
    Returns:
    (||primal residual||_inf / eps_pri)
    (||dual residual||_inf / eps_dua)
    (dua_gap / eps_dua_gap)
    """
    
    eps_abs = EPS_ABS
    eps_rel = EPS_ABS
    
    # Get problem path
    prob_path = get_problem_path(prob_name)
    
    # Get problem data
    P, q, A, l, u, _, _, _ = determine_prob_date(prob_path)
    
    # Primal feasibility
    Ax = A.dot(x)
    
    eps_pri = eps_abs + eps_rel * la.norm(Ax, np.inf)
    pri_res = np.minimum(Ax - l, 0) + np.maximum(Ax - u, 0)
    
    # Dual feasibility
    Px = P.dot(x)
    Aty = A.T.dot(y)
    eps_dua = eps_abs + eps_rel * np.max(
        [la.norm(Px, np.inf), la.norm(q, np.inf), la.norm(Aty, np.inf)]
    )
    dua_res = Px + q + Aty
    
    # Duality gap
    u_inf = np.isinf(u)
    l_inf = np.isinf(l)
    
    u_notinf = u[~u_inf]
    l_notinf = l[~l_inf]
    y_notinf_u = y[~u_inf]
    y_notinf_l = y[~l_inf]
    
    y_plus = np.maximum(y_notinf_u, 0)
    y_minus = np.minimum(y_notinf_l, 0)
    
    supp_func = u_notinf.dot(y_plus) + l_notinf.dot(y_minus)
    xPx = x.dot(Px)
    qx = q.dot(x)
    dua_gap = np.abs(xPx + qx + supp_func)
    eps_dua_gap = eps_abs + eps_rel * np.max(
        [np.abs(xPx), np.abs(qx), np.abs(supp_func)]
    )
    
    primal_ratio = la.norm(pri_res, np.inf) / eps_pri
    dua_ratio = la.norm(dua_res, np.inf) / eps_dua
    duality_gap_ratio = dua_gap / eps_dua_gap
    
    return primal_ratio, dua_ratio, duality_gap_ratio
    

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
    # alpha_vals = [1., 1.4, 1.6, 1.8]
    alpha_vals = [1., 1.6]
    # rho_is_vec_vals = [0, 1]
    # Omit 0
    # adaptive_rho_tolerance_vals = np.linspace(0, 10, 5)[1:]
    adaptive_rho_tolerance_vals = np.linspace(0, 10, 3)[1:]
    rho_custom_condition_vals = [0, 1]
    # Omit 0
    # rho_custom_tolerance_vals = np.linspace(0, 10, 5)[1:]
    rho_custom_tolerance_vals = np.linspace(0, 10, 3)[1:]
    restart_type = "halpern"
    adapt_rho_on_restart_vals = [0, 1]
    # beta_vals = np.linspace(0, 1, 11)[1:-1]
    beta_vals = np.linspace(0, 1, 6)[1:-1]
    # ini_rest_len_vals = np.arange(10, 50, 5)
    ini_rest_len_vals = np.arange(10, 50, 3)
    halpern_scheme_vals = ["none", "adaptive", "adaptive only before init_rest_len"]
    adaptive_rest_vals = [0, 1]
    # Omit 0
    # restart_necessary_vals = np.linspace(0, 1, 5)[1:]
    restart_necessary_vals = np.linspace(0, 1, 3)[1:]
    # restart_artificial_vals = np.linspace(0, 1, 5)[1:]
    restart_artificial_vals = np.linspace(0, 1, 3)[1:]
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


# def determine_prob_date(problem_name, path, n, m, p_scale, p_rank, seed):
def determine_prob_date(path):
    if PROBLEM_NAME != "RANDOM":
        # Define problem data
        name = os.path.basename(path)[:-4]
        # print(f"name: {name}")
        # print(f"Solving problem: {name}")

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
        # print(f"name: {name}")
        print(f"Solving problem: {name}")
        P, q, A, l, u = generate_random_qp(
            n,
            m,
            # p_scale=p_scale,
            # p_rank=p_rank,
            # seed=seed,
            p_scale=P_SCALE,
            p_rank=P_RANK,
            seed=SEED,
        )
        
    return P, q, A, l, u, n, m, name

def get_problem_path(problem_name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    qpbenchmark_data_dir = os.path.join(script_dir, "..", "..", "maros_meszaros_qpbenchmark", "data")
    qpbenchmark_data_dir = os.path.abspath(qpbenchmark_data_dir)
    
    if not problem_name.upper().endswith('.mat'):
        problem_name = problem_name.upper() + '.mat'
    else:
        problem_name = problem_name.upper()
        
    problem_path = os.path.join(qpbenchmark_data_dir, problem_name)
    
    if not os.path.exists(problem_path):
        try:
            files = [f for f in os.listdir(qpbenchmark_data_dir) if f.endswith('.mat')]
            available = sorted([f[:-4] for f in files])  # Remove .mat extension
            raise FileNotFoundError(
                f"Problem '{problem_name[:-4]}' not found. "
                f"Available problems: {available}"
            )
        except OSError:
            raise FileNotFoundError(f"Problem directory not found: {qpbenchmark_data_dir}")
        
    return problem_path


def problem_solver(
    # paths, problem_name, n, m, p_scale, p_rank, seed, plot, save_stats, 
    # save_all, combination_array, all_param_dict, current_comb, task_id
    current_comb, task_id
):
    problem_time = time.time()
    prob = osqp.OSQP()
    
    # Prepare lists for each metric
    setup_times = []
    solve_times = []
    run_times = []
    statuses = []
    prim_res = []
    dual_res = []
    duality_gaps = []
    restarts = []
    iterations = []
    task_ids = []
    prim_res_ratios = []
    dua_res_ratios = []
    duality_gap_res_ratios = []
    integral_sums = []
    
    # for path in paths:
    for path in PATHS:
        # Problem Data
        # P, q, A, l, u, n, m, name = determine_prob_date(problem_name, path, n, m, p_scale, p_rank, seed)
        P, q, A, l, u, n, m, name = determine_prob_date(path)
            
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
        
        # prob.setup(P, q, A, l, u, plot=plot, scaling=0, eps_abs=1e-6, eps_rel=0, max_iter=25000, time_limit=1e3, **current_comb)
        # prob.setup(P, q, A, l, u, plot=plot, verbose=False, time_limit=1e3, **current_comb)
        prob.setup(P, q, A, l, u, plot=PLOT, verbose=False, eps_abs=EPS_ABS, eps_rel=EPS_REL, eps_prim_inf=EPS_PRIM_INF, eps_dual_inf=EPS_DUAL_INF, time_limit=1e3, max_iter=25000, **current_comb)
        # prob.setup(P, q, A, l, u, verbose=True, eps_abs=1e-6, eps_rel=1e-6, time_limit=1e3, pid_controller=1, pid_controller_log=1, max_iter=25000)
        prob.setup(P, q, A, l, u, plot=PLOT, eps_abs=EPS_ABS, eps_rel=EPS_ABS, max_iter=100, **current_set)
        
        
        # # Using CVXPY
        # x = cp.Variable(P.shape[0])
        # obj = cp.Minimize(0.5 * cp.quad_form(x, P) + q.T @ x)
        # constraints = [l <= A @ x, A @ x <= u]
        # prob = cp.Problem(obj, constraints)
        # prob.solve(solver=cp.GUROBI, FeasibilityTol=1e-6, OptimalityTol=1e-6)
        
        
        # # Using Gurobi
        # model = gp.Model("problems")
        # x = model.addMVar(P.shape[0], lb=-GRB.INFINITY, name="x")
        # obj = x @ P @ x + q @ x
        # model.setObjective(obj, GRB.MINIMIZE)
        # model.addConstr(A @ x >= l, name="lower_bounds")
        # model.addConstr(A @ x <= u, name="upper_bounds")
        # model.setParam("FeasibilityTol", 1e-7)
        # model.setParam("OptimalityTol", 1e-9)
        # model.setParam("TimeLimit", 60)
        # model.optimize()
        
        
        # Settings can be changed using .update_settings()
        # prob.update_settings(polishing=1)

        # Solve problem
        res = prob.solve(raise_error=False)
        
        # print('Setup time:', res.info.setup_time)
        # print('Solve time:', res.info.solve_time)
        print('Problem name:', name, flush=True)
        print('Run time:', res.info.run_time, flush=True)
        print('Status:', res.info.status, flush=True)
        # print('Primal Residual:', res.info.prim_res)
        # print('Dual Residual:', res.info.dual_res)
        # print('Duality Gap:', res.info.duality_gap)
        # print('Number of restarts:', res.info.restart)
        
        # Compute the primal, dual, and duality gap ratios
        primal_ratio, dua_ratio, duality_gap_ratio = accuracy_ratios(name, res.x, res.y)
        
        # Collect results
        setup_times.append(res.info.setup_time)
        solve_times.append(res.info.solve_time)
        run_times.append(res.info.run_time)
        statuses.append(res.info.status)
        prim_res.append(res.info.prim_res)
        dual_res.append(res.info.dual_res)
        duality_gaps.append(res.info.duality_gap)
        restarts.append(res.info.restart)
        iterations.append(res.info.iter)
        task_ids.append(task_id)
        prim_res_ratios.append(primal_ratio)
        dua_res_ratios.append(dua_ratio)
        duality_gap_res_ratios.append(duality_gap_ratio)
        integral_sums.append(res.info.total_integral)
        
        # Save result statistics
        # if save_stats:
        if SAVE_STATS:
            # saving_stats_csv(name, res, combination_array, all_param_dict, task_id)
            saving_stats_csv(name, res, current_comb, task_id, primal_ratio, dua_ratio, duality_gap_ratio)
            # saving_stats_csv(name, prob, current_comb, task_id)
            # saving_stats_csv(name, model, current_comb, task_id)
            # stat_file = "../../osqp/plot/stats.csv"
            
            # # Solution details
            # stats_data = {
            #     "problem name": name,
            #     "setup time": res.info.setup_time,
            #     "solve time": res.info.solve_time,
            #     "run time": res.info.run_time,
            #     "primal residual": res.info.prim_res,
            #     "dual residual": res.info.dual_res,
            #     "duality gap": res.info.duality_gap,
            #     "restart": res.info.restart,
            #     "iterations": res.info.iter
            # }
            
            # # Adding paramters
            # comb_keys = list(combination_array[i].keys())
            # all_param_dict_keys = list(all_param_dict.keys())
            
            # for key in all_param_dict_keys:
            #     stats_data[key] = combination_array[i][key] if key in comb_keys else all_param_dict[key]
            
            # # Writing to .csv file and adding header if needed
            # file_exists = os.path.exists(stat_file)
            # file_is_empty = os.path.getsize(stat_file) == 0 if file_exists else True
            # with open(stat_file, mode='a', newline='', encoding='utf8') as csvfile:
            #     writer = csv.DictWriter(csvfile, fieldnames=stats_data.keys())
            #     if not file_exists or file_is_empty:
            #         writer.writeheader()
            #     writer.writerow(stats_data)
            
            
            
        # Transfer the iterates values for the current problem from the .csv file generated by 
        #   the osqp sovler to a master .csv file that stores the iterates for all problems
        # if (problem_name == "ALL.mat") and (save_all == 1):
        if (PROBLEM_NAME == "ALL.mat") and (SAVE_ALL == 1):
            saving_iter_csv(name)
            # residuals_file = "../../osqp/plot/residuals.csv"
            # # master_file = "../../osqp/plot/osqp_reflected_halpern_residuals.csv"
            # master_file = "../../osqp/plot/osqp_default_residuals.csv"
            
            # df = pd.read_csv(residuals_file, header=0)
            # df.insert(1, 'prob_name', name)
            
            # if os.path.exists(master_file):
            #     df.to_csv(master_file, mode='a', header=True, index=False, float_format="%.3e")
            #     # df.to_csv(master_file, mode='a', header=True, index=False)
            
        
        # Struggles to solve simple problems
        if (((res.info.iter >= 24900) or (res.info.status == 'unsolved')) and (path in PATH_LESS_200)):
            print("pruned as iter >= 24900 or status == unsolved for an easy problem", flush=True)
            raise optuna.TrialPruned()
        
        # The run takes too long (over 5 times the pure OSQP implementation)
        if time.time() - problem_time >= TIME_THRESHOLD:
            print(f"pruned as the run took too long. It took {time.time() - problem_time} seconds", flush=True)
            raise optuna.TrialPruned()
        
        # Tells that a problem is not convex (solver fail)
        if res.info.status == 'problem non convex':
            print(f"pruned as status == prolem non convex", flush=True)
            raise optuna.TrialPruned()


    # Convert lists to numpy arrays
    results = {
        'setup_time': np.array(setup_times),
        'solve_time': np.array(solve_times),
        'run_time': np.array(run_times),
        'status': np.array(statuses, dtype=object),
        'prim_res': np.array(prim_res),
        'dual_res': np.array(dual_res),
        'duality_gap': np.array(duality_gaps),
        'restart': np.array(restarts),
        'iterations': np.array(iterations),
        'task_id': np.array(task_ids),
        'prim_res_ratios': np.array(prim_res_ratios),
        'dua_res_ratios': np.array(dua_res_ratios),
        'duality_gap_res_ratios': np.array(duality_gap_res_ratios),
        'integral_sums': np.array(integral_sums)
    }
    
    return results


# def saving_stats_csv(name, res, combination_array, all_param_dict, task_id):
def saving_stats_csv(name, res, current_comb, task_id, primal_ratio, dua_ratio, duality_gap_ratio):
    if CLUSTER == "della_stellato":
        stat_file = "../../osqp/plot/stats.csv"
    elif CLUSTER == "della":
        stat_file = "../../osqp/plot/stats_della.csv"
    
    # stat_file = "../../osqp/plot/stats_orig_osqp.csv"
            
    # Solution details
    stats_data = {
        "problem name": name,
        "setup time": res.info.setup_time,
        "solve time": res.info.solve_time,
        "run time": res.info.run_time,
        "status": res.info.status,
        "primal residual": res.info.prim_res,
        "primal ratio": primal_ratio,
        "dual residual": res.info.dual_res,
        "dual ratio": dua_ratio,
        "duality gap": res.info.duality_gap,
        "duality gap ratio": duality_gap_ratio,
        "integral_sum": res.info.total_integral,
        "restart": res.info.restart,
        "iterations": res.info.iter,
        "task_id": task_id,
        
        # "run time": res.solver_stats.solve_time,
        # "status": res.status == gp.GRB.OPTIMAL
    }
    
    # Adding paramters
    comb_keys = list(current_comb.keys())
    # all_param_dict_keys = list(all_param_dict.keys())
    all_param_dict_keys = list(ALL_PARAM_DICT.keys())
    
    for key in all_param_dict_keys:
        # stats_data[key] = combination_array[task_id][key] if key in comb_keys else all_param_dict[key]
        stats_data[key] = current_comb[key] if key in comb_keys else ALL_PARAM_DICT[key]
    
    # Writing to .csv file and adding header if needed
    file_exists = os.path.exists(stat_file)
    file_is_empty = os.path.getsize(stat_file) == 0 if file_exists else True
    with open(stat_file, mode='a', newline='', encoding='utf8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=stats_data.keys())
        if not file_exists or file_is_empty:
            writer.writeheader()
        writer.writerow(stats_data)
        
    return


def saving_iter_csv(name):
    residuals_file = "../../osqp/plot/residuals.csv"
    # master_file = "../../osqp/plot/osqp_reflected_halpern_residuals.csv"
    master_file = "../../osqp/plot/osqp_default_residuals.csv"
    
    df = pd.read_csv(residuals_file, header=0)
    df.insert(1, 'prob_name', name)
    
    if os.path.exists(master_file):
        df.to_csv(master_file, mode='a', header=True, index=False, float_format="%.3e")
        # df.to_csv(master_file, mode='a', header=True, index=False)
        
    return

def dispatch(
    # task_id, combination_array, paths, problem_name, n, m, p_scale, 
    # p_rank, seed, plot, save_stats, save_all, all_param_dict
    task_id, combination_array
):
    current_comb = combination_array[task_id]
    
    return problem_solver(
        # paths, problem_name, n, m, p_scale, p_rank, seed, plot, save_stats, 
        # save_all, combination_array, all_param_dict, current_comb, task_id
        current_comb, task_id
    )
    
    
def run_locally(
    # combination_array, paths, problem_name, n, m, p_scale, 
    # p_rank, seed, plot, save_stats, save_all, all_param_dict
    combination_array
):
    # num_cpus = int(int(multiprocessing.cpu_count()) / 2.0)
    num_cpus = int(multiprocessing.cpu_count())
    num_cpus = 30
    print(f"Using {num_cpus} CPUs for {len(combination_array)} total tasks")
    # print(f"task_id = {task_id}")
    
    task_func = partial(
        # dispatch, combination_array=combination_array, paths=paths,
        # problem_name=problem_name, n=n, m=m, p_scale=p_scale,
        # p_rank=p_rank, seed=seed, plot=plot, save_stats=save_stats,
        # save_all=save_all, all_param_dict=all_param_dict
        dispatch, combination_array=combination_array
    )
    
    with multiprocessing.Pool(processes=num_cpus) as pool:
        # results = pool.map(task_func, range(len(combination_array)))
        results = pool.map(task_func, range(60))
    
    return results


def objective(trial):
    # Optuna does not allow for bounds of the form (0, 1) or (0, 1]
    small_num = 1e-4
    
    # All of the combinations set to default values
    current_comb = {}
    
    # Optuna determines to what values to change the parameters
    current_comb["alpha"] = trial.suggest_float("alpha", 0 + small_num, 2 - small_num)
    current_comb["adaptive_rho_tolerance_greater"] = trial.suggest_float("adaptive_rho_tolerance_greater", 0 + small_num, 100)
    current_comb["adaptive_rho_tolerance_less"] = trial.suggest_float("adaptive_rho_tolerance_less", 0 + small_num, 100)
    current_comb["pid_controller"] = trial.suggest_int("pid_controller", 0, 1)
    if current_comb["pid_controller"] == 1:
        current_comb["KP"] = trial.suggest_float("KP", 0, 100)
        current_comb["KI"] = trial.suggest_float("KI", 0, 100)
        
    current_comb["adapt_rho_on_restart"] = 1
    current_comb["beta"] = trial.suggest_float("beta", 0 + small_num, 1 - small_num)
    current_comb["ini_rest_len"] = trial.suggest_int("ini_rest_len", 25, 100, step=25)
    current_comb["adaptive_rest"] = 1
    current_comb["restart_necessary"] = trial.suggest_float("restart_necessary", 0 + small_num, 1)
    current_comb["restart_artificial"] = trial.suggest_float("restart_artificial", 0 + small_num, 1)
        # current_comb["negate_K"] = trial.suggest_int("negate_K", 0, 1)
        # if current_comb["negate_K"]:
        #     current_comb["KP"] = -current_comb["KP"]
        #     current_comb["KI"] = -current_comb["KI"]
        # current_comb["pid_controller_sqrt"] = trial.suggest_int("pid_controller_sqrt", 0, 1)
        # current_comb["pid_controller_sqrt_mult"] = trial.suggest_int("pid_controller_sqrt_mult", 0, 1)
        # current_comb["pid_controller_sqrt_mult_2"] = trial.suggest_int("pid_controller_sqrt_mult_2", 0, 1)
        # current_comb["pid_controller_log"] = trial.suggest_int("pid_controller_log", 0, 1)
    
    # current_comb["rho_custom_condition"] = trial.suggest_int("rho_custom_condition", 0, 1)
    # if current_comb["rho_custom_condition"] == 1:
    #     current_comb["rho_custom_tolerance"] = trial.suggest_float("rho_custom_tolerance", 0 + small_num, 100)
        
    # current_comb["restart_type"] = trial.suggest_categorical("restart_type", ["none", "halpern", "reflected halpern", "averaged"])
    # current_comb["restart_type"] = trial.suggest_categorical("restart_type", ["halpern", "reflected halpern", "averaged"])
    current_comb["restart_type"] = trial.suggest_categorical("restart_type", ["reflected halpern", "averaged"])
    # if current_comb["restart_type"] == "halpern":
    #     # current_comb["adapt_rho_on_restart"] = trial.suggest_int("adapt_rho_on_restart", 0, 1) 
    #     current_comb["adapt_rho_on_restart"] = 1
    #     current_comb["beta"] = trial.suggest_float("beta", 0 + small_num, 1 - small_num)
    #     # current_comb["ini_rest_len"] = trial.suggest_int("ini_rest_len", 1, 100)
    #     current_comb["ini_rest_len"] = trial.suggest_int("ini_rest_len", 25, 100, step=25)
    #     # current_comb["halpern_scheme"] = trial.suggest_categorical("halpern_scheme", ["none", "adaptive", "adaptive only before init_rest_len"])
    #     # current_comb["halpern_scheme"] = trial.suggest_categorical("halpern_scheme", ["none", "adaptive"])
    #     current_comb["halpern_scheme"] = trial.suggest_categorical("halpern_scheme", ["none"])
        
    #     # current_comb["adaptive_rest"] = trial.suggest_int("adaptive_rest", 0, 1)
    #     # if current_comb["adaptive_rest"] == 1:
    #     #     current_comb["restart_necessary"] = trial.suggest_float("restart_necessary", 0 + small_num, 1)
    #     #     current_comb["restart_artificial"] = trial.suggest_float("restart_artificial", 0 + small_num, 1)
    #     current_comb["adaptive_rest"] = 1
    #     current_comb["restart_necessary"] = trial.suggest_float("restart_necessary", 0 + small_num, 1)
    #     current_comb["restart_artificial"] = trial.suggest_float("restart_artificial", 0 + small_num, 1)
    #     current_comb["halpern_step_first_inner_iter"] = trial.suggest_int("halpern_step_first_inner_iter", 0, 1)
    # elif current_comb["restart_type"] == "reflected halpern":
    if current_comb["restart_type"] == "reflected halpern":
        # current_comb["adapt_rho_on_restart"] = trial.suggest_int("adapt_rho_on_restart", 0, 1)
        # current_comb["halpern_scheme"] = trial.suggest_categorical("halpern_scheme", ["none", "adaptive", "adaptive only before init_rest_len"])
        current_comb["halpern_scheme"] = trial.suggest_categorical("halpern_scheme", ["none"])
        current_comb["lambd"] = trial.suggest_float("lambd", 0, max(0, min((2.0 / current_comb["alpha"]) - 1.0, 1)))
        # if current_comb["alpha"] > 1:
        #     current_comb["alpha_adjustment_reflected_halpern"] = 1
        # else:
        #     current_comb["alpha_adjustment_reflected_halpern"] = trial.suggest_int("alpha_adjustment_reflected_halpern", 0, 1)
        current_comb["alpha_adjustment_reflected_halpern"] = 1
        # current_comb["adaptive_rest"] = trial.suggest_int("adaptive_rest", 0, 1)
        # if current_comb["adaptive_rest"] == 1:
        #     current_comb["restart_necessary"] = trial.suggest_float("restart_necessary", 0 + small_num, 1)
        #     current_comb["restart_artificial"] = trial.suggest_float("restart_artificial", 0 + small_num, 1)
        current_comb["halpern_step_first_inner_iter"] = trial.suggest_int("halpern_step_first_inner_iter", 0, 1)
    elif current_comb["restart_type"] == "averaged":
        # current_comb["adapt_rho_on_restart"] = trial.suggest_int("adapt_rho_on_restart", 0, 1)
        # current_comb["custom_average_rest"] = trial.suggest_int("custom_average_rest", 0, 1)
        current_comb["custom_average_rest"] = 1
        current_comb["xi"] = trial.suggest_float("xi", 0 + small_num, 100)
        current_comb["vector_rho_in_averaged_KKT"] = trial.suggest_int("vector_rho_in_averaged_KKT", 0, 1)
        # current_comb["adaptive_rest"] = trial.suggest_int("adaptive_rest", 0, 1)
        # if current_comb["adaptive_rest"] == 1:
        #     current_comb["restart_necessary"] = trial.suggest_float("restart_necessary", 0 + small_num, 1)
        #     current_comb["restart_artificial"] = trial.suggest_float("restart_artificial", 0 + small_num, 1)
    
    # Feeding the parameters to solve all of the desired problems
    res = problem_solver(
        # paths, problem_name, n, m, p_scale, p_rank, seed, plot, save_stats, 
        # save_all, combination_array, all_param_dict, current_comb, task_id
        current_comb, trial.number
    )
    
    # Geometric mean of solve time
    conditions = [
        res['status'] == 'solved',
        res['status'] == 'problem non convex',
        (res['status'] == 'maximum iterations reached') | (res['status'] == 'solved inaccurate')
    ]
    default = 1e6
    
    # For pure OSQP total run time is 543.947886066
    
    # solve_time_geom = geom_mean(res['solve_time'])
    # For pure OSQP it is 2.234657228416994
    # solve_time_geom = geom_mean(np.select(conditions, [res['solve_time'], res['solve_time'] + 1e10, res['solve_time'] + 1e3], default=default))
    # For pure OSQP it is 
    # The prim_res and dual_res are typically small compared to 1e4, so I want to scale them up a little. Dividing by 1e-6 (the eps_abs/tolerance)
    #       will work but the resulting geom_mean will be extremely large for all solvers as some prim_res and dual_res are 100 or even 1,000 even for
    #       the pure OSQP solver. This will then create an over emthasis on solvign a small subset of the hard problems rather than the overall solver
    #       efficiency. For this reason we divide by 1e-3 instead of 1e-6, this still puts an emthasis on solving those hard problems but evens out the
    #       emthasis on those problems. The min(np.abs(...), 1e10) is purely for numerical stability purposes.
    scaled_prim_dual = np.minimum(np.abs((res['prim_res'] + res['dual_res']) / (1e-3)), 1e10)
    # For pure OSQP solve_time_geom = 255.32010499705535
    solve_time_geom = geom_mean(np.select(conditions, [res['solve_time'], 
                                                       1e15, 
                                                       1e4 + res['solve_time'] + scaled_prim_dual + np.abs(res['duality_gap'])],
                                          default=default))
    
    # Geometric mean of integral sum
    integral_geom = geom_mean(np.select(conditions, [res['integral_sums'], 
                                                       1e10, 
                                                       1e2 + res['integral_sums']],
                                          default=default))
    
    # Geometric mean of primal residual
    # prim_res_geom = geom_mean(res['prim_res'])
    prim_res_geom = geom_mean(np.select(conditions, [res['prim_res'], res['prim_res'] + 1e10, res['prim_res'] + 1e3], default=default))
    
    # Geometric mean of dual residual
    # dual_res_geom = geom_mean(res['dual_res'])
    dual_res_geom = geom_mean(np.select(conditions, [res['dual_res'], res['dual_res'] + 1e10, res['dual_res'] + 1e3], default=default))
    
    # Geometric mean of duality gap
    # duality_gap_geom = geom_mean(np.abs(res['duality_gap']))
    duality_gap_geom = geom_mean(np.select(conditions, [np.abs(res['duality_gap']), np.abs(res['duality_gap']) + 1e10, np.abs(res['duality_gap']) + 1e3], default=default))
    
    # Geometric mean of iteration count
    # iterations_geom = geom_mean(res['iterations'])
    iterations_geom = geom_mean(np.select(conditions, [res['iterations'], res['iterations'] + 1e10, res['iterations'] + 1e3], default=default))
    
    # return solve_time_geom
    return integral_geom


def run_optimization(_):
    if CLUSTER == "della_stellato":
        # journal_log_file_path = "./journal.log"
        db_path = "./optuna_study.db"
    elif CLUSTER == "della":
        # journal_log_file_path = "./journal_della.log"
        db_path = "./optuna_study_della.db"
    study = optuna.create_study(
        study_name="journal_storage_multiprocess",
        # storage=JournalStorage(JournalFileBackend(file_path=journal_log_file_path)),
        storage=f"sqlite:///{db_path}",
        load_if_exists=True, # Useful for multi-process or multi-node optimization.
        sampler=optuna.samplers.CmaEsSampler()
    )
    study.optimize(objective, n_trials=1)


def run_optuna(use_multiprocessing=True):
    # num_cpus = int(int(multiprocessing.cpu_count()) / 2.0)
    # num_cpus = max(1, int(multiprocessing.cpu_count()) - 2)
    if CLUSTER == "della_stellato":
        # num_cpus = 34
        num_cpus = 1
    elif CLUSTER == "della":
        num_cpus = 1
    # print(f"Using {num_cpus} CPUs for {len(combination_array)} total tasks")
    
    if use_multiprocessing:
        with multiprocessing.Pool(processes=num_cpus) as pool:
            if CLUSTER == "della_stellato":
                # results = pool.map(run_optimization, range(250))
                # results = pool.map(run_optimization, range(40))
                results = run_optimization(1)
            elif CLUSTER == "della":
                results = pool.map(run_optimization, range(512))
                # results = pool.map(run_optimization, range(1))
            
            # # Dynamically allocate jobs to ensure CPUs dont sit idly
            # elif CLUSTER == "della":
            #     results = []
            #     active_tasks = []
                
            #     for i in range(num_cpus):
            #         active_tasks.append(pool.apply_async(run_optimization, (i,)))
                    
            #     while active_tasks:
            #         ready_task = next((task for task in active_tasks if task.ready()), None)
                    
            #         if ready_task:
            #             try:
            #                 result = ready_task.get()
            #                 results.append(result)
            #             except Exception as e:
            #                 print(f"A task failed with an exception: {e}")
            #             active_tasks.remove(ready_task)
                        
            #             elapsed_time = time.time() - START_TIMER
                        
            #             if elapsed_time < JOB_TIME_LIMIT - TIME_THRESHOLD:
            #                 new_task_id = len(result) + len(active_tasks)
            #                 active_tasks.append(pool.apply_async(run_optimization, (new_task_id,)))
            #             else:
            #                 print("Not enough time to start a new taks, waiting till all existing tasks finish")
                            
            #         else:
            #             # Just wait
            #             time.sleep(1)
                        
            #     # Save the number of tasks completed
            #     stats_file = "../../osqp/plot/stats_orig_osqp_old.csv"
            #     results_len = len(results)
            #     file_exists = os.path.exists(stats_file)
            #     with open(stats_file, mode='a', newline='', encoding='utf8') as csvfile:
            #         writer = csv.writer(csvfile)
            #         if not file_exists:
            #             writer.writerow(['results_len'])
            #         writer.writerow([results_len])
    else:
        results = run_optimization(1)
    
    return results



if __name__ == '__main__':
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        print("script_dir =", script_dir, flush=True)
        
        # Check for command line argument
        # if len(sys.argv) != 2 and sys.argv[1] != "random":
        #     print("Usages:")
        #     print("Existing Problem: py maros_meszaros_benchmarks.py <problem_name>")
        #     print("Random Problem: py maros_meszaros_benchmarks.py random")
        #     print("Example: py maros_meszaros_benchmarks.py HS21")
        #     print("Available problems are .mat files in qpbenchmark/maros_meszaros_qpbenchmark/data/")
        #     sys.exit(1)

        array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID")
        # print(f"Current SLURM array task ID: {array_task_id}")
        
        parser = argparse.ArgumentParser(description="Run Maros-Meszaros benchmarks.")
        parser.add_argument("problem_name", type=str, help="Name of the problem (e.g., HS21 or 'random').")
        parser.add_argument("--n", type=int, default=100, help="Number of variables for random QP generation.")
        parser.add_argument("--seed", type=int, default=1, help="Random seed for reproducibility.")
        parser.add_argument("--p_scale", type=float, default=0.01, help="Scaling factor for P matrix.")
        parser.add_argument("--p_rank", type=int, default=25, help="Rank of P matrix (None for full rank).")
        parser.add_argument("--plot", type=int, default=1, help="Parameter to pass to the osqp solver in order for it to temporary store each iteration for ploting purposes")
        parser.add_argument("--save_stat", type=int, default=0, help="Save the setup time, solve time, run time, primal residual, dual residual, duality gap, and the number of restarts. This is a binary either 0 or 1.")
        parser.add_argument("--save_all", type=int, default=0, help="Store all of the iterations for each problem ran (For usecase where problem_name = 'ALL'). This is a binary either 0 or 1.")
        parser.add_argument("--restart_comb", type=str, default="none", help="Restart framework used for the combination generation (e.g. 'none', 'halpern', 'reflected_halpern', or 'averaged').")
        parser.add_argument("--comb_count", type=int, default=1, help="Number of combinations to try (-1 means all).")
        parser.add_argument("--cluster", type=str, help="Name of the cluster used, either 'della' or 'della_stellato'")
        
        # Parse arguments
        args = parser.parse_args()
        
        # Access arguments
        # problem_name = args.problem_name.upper()
        PROBLEM_NAME = args.problem_name.upper()
        # n = 
        N = args.n
        # m = int(n * 1.5)
        M = int(N * 1.5)
        # seed = args.seed
        SEED = args.seed
        # p_scale = args.p_scale
        P_SCALE= args.p_scale
        # p_rank = args.p_rank
        P_RANK = args.p_rank
        # plot = args.plot
        PLOT = args.plot
        # save_stats = args.save_stat
        SAVE_STATS = args.save_stat
        # save_all = args.save_all
        SAVE_ALL = args.save_all
        restart_comb = args.restart_comb
        comb_count = args.comb_count
        CLUSTER = args.cluster
        
        # if (save_stats != 0) and (save_stats != 1):
        if (SAVE_STATS != 0) and (SAVE_STATS != 1):
            print(f"save_stats must be 0 or 1", flush=True)
            sys.exit(1)
            
        # if (save_all != 0) and (save_all != 1):
        if (SAVE_ALL != 0) and (SAVE_ALL != 1):
            print(f"save_all must be 0 or 1", flush=True)
            sys.exit(1)
            
        if (restart_comb != "none") and (restart_comb != "halpern") and (restart_comb != "reflected_halpern") and (restart_comb != "averaged"):
            print(f"restart_comb must be either 'none', 'halpern', 'reflected_halpern', or 'averaged' ", flush=True)
            sys.exit(1)
        
        if (comb_count <= 0) and (comb_count != -1):
            print(f"comb_count must be a positive integer or -1", flush=True)
            sys.exit(1)
            
        if (CLUSTER != "della") and (CLUSTER != "della_stellato"):
            print(f"cluster must be either 'della' or 'della_stellato'", flush=True)
            sys.exit(1)
        
        
        # print(f"problem_name {problem_name}")
        # print(f"problem_name {PROBLEM_NAME}")
        # if problem_name == "RANDOM":
        if PROBLEM_NAME == "RANDOM":
            # print(f"n: {n}, m: {m}, seed: {seed}, p_scale: {p_scale}, p_rank: {p_rank}")
            print(f"n: {N}, m: {M}, seed: {SEED}, p_scale: {P_SCALE}, p_rank: {P_RANK}", flush=True)
        
        # Construct path to the qpbenchmark data directory
        # Navigate from osqp-python/examples to qpbenchmark/maros_meszaros_qpbenchmark/data
        qpbenchmark_data_dir = os.path.join(script_dir, "..", "..", "maros_meszaros_qpbenchmark", "data")
        qpbenchmark_data_dir = os.path.abspath(qpbenchmark_data_dir)
        
        # if problem_name != "RANDOM":
        if PROBLEM_NAME != "RANDOM":
            # Add .mat extension if not provided
            # if not problem_name.endswith('.mat'):
            if not PROBLEM_NAME.endswith('.mat'):
                # problem_name += '.mat'
                PROBLEM_NAME += '.mat'
                
            # if (problem_name == "ALL.mat"):
            if (PROBLEM_NAME == "ALL.mat"):
                exclude_probs = ['CONT-300.mat', 'CONT-201.mat', 'CONT-200.mat', 'BOYD2.mat']
                
                files = [f for f in os.listdir(qpbenchmark_data_dir) if (f.endswith('.mat')) and (f not in exclude_probs)]
                files.sort()
                
                # paths = [os.path.join(qpbenchmark_data_dir, f) for f in files]
                PATHS = [os.path.join(qpbenchmark_data_dir, f) for f in files]
            
            else:
                # problem_path = os.path.join(qpbenchmark_data_dir, problem_name)
                problem_path = os.path.join(qpbenchmark_data_dir, PROBLEM_NAME)
                
                # Check if the problem file exists
                if not os.path.exists(problem_path):
                    print(f"Error: Problem file '{problem_path}' not found.", flush=True)
                    print(f"Available problems in '{qpbenchmark_data_dir}':", flush=True)
                    try:
                        files = [f for f in os.listdir(qpbenchmark_data_dir) if f.endswith('.mat')]
                        for f in sorted(files):
                            print(f"  {f[:-4]}", flush=True)  # Remove .mat extension for display
                    except FileNotFoundError:
                        print(f"  Directory not found: {qpbenchmark_data_dir}", flush=True)
                    sys.exit(1)
                
                # paths = [problem_path]
                PATHS = [problem_path]
        
            # print(f"Running problem: {problem_name}")
            # print(f"Running problem: {PROBLEM_NAME}")
        else:
            # Paths is not important for randomly generated problems
            # paths = ["random"]
            PATHS = ["random"]
            # print("Running randomly generated problem")

        
        # Generate the set of all possible combinations
        #       This requires a lot of free memory on the system 
        #       (my pc can't load all of these combinations, roughly 60 million comb.)

        if (restart_comb == "halpern"):
            # Combinations 4700160 (if use 5 choices per hyperparemeter)
            # Combinations 80640 (if use 3 choices per hyperparemeter)
            combination_array = halpern_combinations()
        elif (restart_comb == "reflected_halpern"):
            # Combinations 23500800
            combination_array = reflected_halpern_combinations()
        elif (restart_comb == "averaged"):
            # Combinations 31334400
            combination_array = averaged_combinations()
        elif (restart_comb == "none"):
            combination_array = [{}]
            
        # all_param_dict = {
        ALL_PARAM_DICT = {
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
            "halpern_step_first_inner_iter": -1,
            "pid_controller": -1,
            "KP": -1,
            "KI": -1,
            "negate_K": -1
        }

        # Create an OSQP object
        prob = osqp.OSQP()
        
        if (comb_count > len(combination_array)):
            comb_count = len(combination_array)
            print(f"comb_count is set to the len(combination_array) [{len(combination_array)}]", flush=True)
        
        if (comb_count == -1):
            comb_count = len(combination_array)


        # Time the execution process
        START_TIMER = time.time()


        # Run it in parallel
        if array_task_id is not None:
            # SLURM Mode: Run only one task
            # task_id = int(array_task_id)
            # results = dispatch(
            #     # task_id, combination_array, paths, problem_name, n, m, p_scale, 
            #     # p_rank, seed, plot, save_stats, save_all, all_param_dict
            #     task_id, combination_array
            # )
            results = run_optuna(use_multiprocessing=False)
        else:
            # # Multithreading
            # results = run_locally(
            #     # combination_array[:comb_count], paths, problem_name, n, m, p_scale, 
            #     # p_rank, seed, plot, save_stats, save_all, all_param_dict
            #     combination_array[:comb_count]
            # )
            
            # Optuna (Multi-process)
            results = run_optuna()
        
        
            
        # for i in range(comb_count):
        #     current_comb = combination_array[i]
        
        #     problem_solver(
        #         paths, problem_name, n, m, p_scale, p_rank, seed, plot, save_stats, 
        #         save_all, combination_array, all_param_dict, current_comb, i
        #     )
        
        if CLUSTER == "della_stellato":
            end = time.time()
            print(f"To solve all of the instance it took {end - START_TIMER} seconds")
            print(f"To solve all of the instance it took {(end - START_TIMER) / 60.0} minutes")
            print(f"To solve all of the instance it took {((end - START_TIMER) / 60.0) / 60.0} hours")
            
            
            
            # Write an email message when the work is finished    
            msg = EmailMessage()
            msg['Subject'] = 'Python Jobs Report'
            msg['From'] = 'emailsender062@gmail.com'
            msg['To'] = 'lb3825@princeton.edu'
            msg.set_content(f'Job Finished !!!! It took {(end - START_TIMER) / 60.0} minutes')

            # Gmail SMTP server setup
            smtp_server = 'smtp.gmail.com'
            smtp_port = 587

            # Sending the email
            with smtplib.SMTP(smtp_server, smtp_port) as server:
                server.starttls()
                server.login('emailsender062@gmail.com', 'qull hdcu ivlb hubu')
                server.send_message(msg)

            print("Email sent successfully.")
    # except Exception as e:
    except BaseException as e:
        if CLUSTER == "della_stellato":
            tb = traceback.format_exc()
            send_email('Python Jobs Report (ERROR)', f'Job failed with error:\n{tb}')
            
            print("Error email sent.")
        raise