# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the Python wrapper for OSQP (Operator Splitting Quadratic Program solver), a numerical optimization package. The package uses scikit-build-core with CMake to build C/C++ extensions via pybind11.

## Build System Architecture

### Multi-Backend System
OSQP supports multiple algebra backends with different performance characteristics:
- **builtin**: Default CPU backend, always available (built into main package)
- **mkl**: Intel MKL-accelerated backend (separate package: `osqp-mkl`)
- **cuda**: CUDA GPU backend (separate package: `osqp-cu12`)

The backend selection happens at runtime via:
1. `OSQP_ALGEBRA_BACKEND` environment variable (if set)
2. Auto-detection in priority order: cuda → mkl → builtin (see `src/osqp/interface.py:14-18`)

Each backend compiles to a separate extension module:
- `ext_builtin` (in osqp package)
- `osqp_mkl` (in osqp-mkl package)
- `osqp_cuda` (in osqp-cu12 package)

### Build Configuration
The root `CMakeLists.txt` fetches OSQP from https://github.com/lb3825/osqp.git (branch: `b/diagonal_step_sizes`), not the official repo. Custom memory/printing routines are injected via `cmake/printing.h` and `cmake/memory.h`.

The `bindings.cpp.in` template is configured at build time to create `src/bindings.cpp`, which is the pybind11 module source.

## Development Commands

### Installing from Source
The recommended approach is to use `uv` for environment management and installation. Avoid editable installs (`-e`) as they can break other parts of the build system.

```bash
# Recommended: use uv
uv venv                          # Create virtual environment
source .venv/bin/activate        # Activate it
uv pip install .                 # Install main package with builtin backend
uv pip install .[dev]            # Install with development dependencies

# Alternative: standard pip
pip install .                    # Install main package with builtin backend
pip install .[dev]               # Install with development dependencies
```

### Building Backend Packages
```bash
# Build MKL backend
cd backend/mkl && pip install .

# Build CUDA backend
cd backend/cuda && pip install .
```

### Testing
```bash
# Run all tests
pytest src/osqp/tests

# Run single test file
pytest src/osqp/tests/basic_test.py

# Run with specific algebra backend
OSQP_TEST_ALGEBRA_INCLUDE="builtin" pytest src/osqp/tests

# Skip specific backends
OSQP_TEST_ALGEBRA_SKIP="cuda" pytest src/osqp/tests
```

Test parametrization is configured in `src/osqp/tests/conftest.py`, which generates test variations for different algebra backends and solver types (direct/indirect). Environment variables `OSQP_TEST_ALGEBRA_INCLUDE` and `OSQP_TEST_ALGEBRA_SKIP` control which backends are tested.

### Linting and Formatting
```bash
pre-commit run --all-files      # Run all pre-commit hooks
flake8                          # Run flake8 linter
blue --line-length=120 .        # Run blue formatter
```

Code style: uses `blue` formatter (120 char line length) and `absolufy-imports` for absolute imports.

### Building Wheels
```bash
# Build source distribution
python -m build --sdist

# Build wheel for current platform
python -m build --wheel
```

CI uses cibuildwheel (config in `cibuildwheel.toml`) to build wheels for multiple platforms.

## Key Architecture Points

### Extension Module Loading
The `interface.py` module contains the algebra backend loading logic:
- `_ALGEBRAS` tuple defines priority order (cuda > mkl > builtin)
- `_ALGEBRA_MODULES` maps algebra names to importable module names
- `algebra_available()` checks if a backend can be imported
- `default_algebra()` selects backend based on env var or availability

### Python Interface
- Main class: `OSQP` in `src/osqp/interface.py` (algebra-agnostic)
- Convenience classes: `osqp.builtin.OSQP`, `osqp.mkl.OSQP`, `osqp.cuda.OSQP` (pre-set algebra)
- All convenience classes inherit from the main `OSQP` class, just passing the `algebra` kwarg

### Code Generation
The package supports generating standalone C code from a problem. Codegen templates are in `src/osqp/codegen/pywrapper/` and controlled by the `OSQP_CODEGEN` environment variable during build.

### CMake Build Variables
Key CMake defines (set in `pyproject.toml` or backend-specific configs):
- `OSQP_ALGEBRA_BACKEND`: builtin/mkl/cuda
- `OSQP_EXT_MODULE_NAME`: Name of the compiled extension module
- `OSQP_ENABLE_INTERRUPT`: Enable interrupt handling (default ON)
- `OSQP_CODEGEN`: Enable code generation support
- `OSQP_USE_LONG`: Use long integers (currently ON in CMakeLists.txt:7)
- `OSQP_CUSTOM_PRINTING/MEMORY`: Paths to custom header files

## Important File Locations

- `src/osqp/interface.py`: Main Python interface, backend selection logic
- `src/bindings.cpp.in`: Pybind11 binding template (becomes src/bindings.cpp)
- `CMakeLists.txt`: Root CMake configuration
- `src/osqp/tests/conftest.py`: Pytest parametrization for multi-backend testing
- `backend/{cuda,mkl}/pyproject.toml`: Backend-specific build configs
- `pyproject.toml`: Main package configuration

## Notes on Current Branch

Current branch: `feature/pr8-compat` (main branch for PRs: `master`)
Recent commits indicate work on updating the Python wrapper for PR #8 (diagonal step sizes), CUDA backends, and GitHub workflow fixes.
