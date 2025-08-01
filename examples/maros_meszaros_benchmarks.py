import osqp
import os
import sys
import numpy as np
from scipy import sparse
import scipy.io as spio
import scipy.sparse as spa

if __name__ == '__main__':
    script_dir = os.path.dirname(os.path.abspath(__file__))
    print("script_dir =", script_dir)
    
    # Check for command line argument
    if len(sys.argv) != 2:
        print("Usage: py maros_meszaros_benchmarks.py <problem_name>")
        print("Example: py maros_meszaros_benchmarks.py HS21")
        print("Available problems are .mat files in qpbenchmark/maros_meszaros_qpbenchmark/data/")
        sys.exit(1)
    
    problem_name = sys.argv[1].upper()
    
    # Construct path to the qpbenchmark data directory
    # Navigate from osqp-python/examples to qpbenchmark/maros_meszaros_qpbenchmark/data
    qpbenchmark_data_dir = os.path.join(script_dir, "..", "..", "qpbenchmark", "maros_meszaros_qpbenchmark", "data")
    qpbenchmark_data_dir = os.path.abspath(qpbenchmark_data_dir)
    
    # Add .mat extension if not provided
    if not problem_name.endswith('.mat'):
        problem_name += '.mat'
    
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
    print(f"Problem path: {problem_path}")
    print()
    
    for path in paths:
        # Define problem data
        name = os.path.basename(path)[:-4]
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
        
        # Create an OSQP object
        prob = osqp.OSQP()

        # Setup workspace and change alpha parameter
        prob.setup(P, q, A, l, u, alpha=1.0)
        
        # # Set the problem name for CSV logging
        # prob.set_problem_name(problem_name)

        # Settings can be changed using .update_settings()
        prob.update_settings(polishing=1)

        # Solve problem
        res = prob.solve(raise_error=False)
        
        print('Setup time:', res.info.setup_time)
        print('Solve time:', res.info.solve_time)

        # Additional Information
        # if res.info.status_val == osqp.SolverStatus.OSQP_SOLVED:
        #     print('Setup time:', res.info.setup_time)
        #     print('Solve time:', res.info.solve_time)
            # print('Optimal solution x (first 10 elements):')
            # print(res.x[:min(10, len(res.x))])
            # if len(res.x) > 10:
            #     print(f"... (total {len(res.x)} variables)")
        # else:
        #     print(f"Problem not solved successfully. Status: {res.info.status}")
        #     if hasattr(res.info, 'obj_val') and res.info.obj_val is not None:
        #         print('Best objective value found:', res.info.obj_val)
