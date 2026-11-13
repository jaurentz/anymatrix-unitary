# Test strategy and numerical guarantees

This repository contains **only MATLAB source code**, plus documentation and
JSON metadata. The tests run in the repository checkout and after installation
as an Anymatrix remote group. No Python or Statistics Toolbox dependency is
required.

## Local invocation

```matlab
addpath(pwd)
results = test_run;
```

The command fails with an assertion when any verification fails. Tests use
`matlab.unittest` function-based test cases and are independently discoverable
under `tests/`. For CI, `.github/workflows/matlab-tests.yml` runs on pushes,
pull requests and manually dispatched workflows.

## What is tested

| Test class | Independent or diagnostic reference |
|---|---|
| `test_jacobi_measure` | Weight normalization, explicit Verblunsky formulas, forward/inverse Geronimus, classical specializations and moments |
| `test_independent_stieltjes` | Golub–Welsch quadrature followed by an independent Stieltjes recurrence, then inverse Geronimus |
| `test_core_representations` | Unitarity, sparse structure, spectral moments and distinguished-vector quadrature moments across four forms |
| `test_generators` | Existing cycle, CMV, quantum walk, Haar and public interface regression |
| `test_input_validation` | Invalid orders, coefficients, boundaries, and default construction |
| `test_new_families` | Analytical spectra, moments, limiting cases, determinant/symmetry and structural identities of added families |
| `test_opuc_additions` | Bernstein–Szegő density/moments and circular-beta coefficient sampling law |
| `test_opuc_extended` | Generalized circular Jacobi and Rogers–Szegő weights, finite Bernstein–Szegő, constant/periodic real-complex Verblunsky sequences, unitarity and finite closure |
| `test_floquet_families` | Floquet unitarity, exact limits, spin-layer isospectrality, and Kitaev particle-hole symmetry |
| `test_floquet_extended` | Rotor exact phases, embedded two-site gates, local dual-unitarity, XXZ Yang-Baxter, and global magnetization sectors |
| `test_signal_processing_extended` | DFT fourth-power identity, DST I–IV, Walsh sequency, Slepian concentration and graph eigenspaces; only useful fast implicit actions |
| `test_signal_processing` | Exact identities, spectra, autocorrelations, multi-RHS forward/adjoint, implicit matrix materialization and large-order actions |
| `test_catalog` | Catalog/documentation completeness and consistency of Anymatrix property declarations |
| `test_math_markup` | GitHub display-math fence consistency and guards against bare dollar delimiters / unsupported notation |

## Accuracy versus identity

Unitary and orthogonal identities are checked numerically in double precision;
nontrivial eigenvalue reference tests often compare traces of powers against
independently specified known spectra. Random generators are checked for
seed reproducibility and restoration of the global MATLAB RNG state.

Passing finite-dimensional numerical tests does **not** prove ensemble limits,
large-order numerical stability, or mathematical theorems beyond those already
established in the cited literature. The supplied tests cover representative
orders and parameter settings, not every floating-point regime. Further tests
for extreme Jacobi parameters, large phase wrapping and near-singular polar
problems remain desirable.

## Anymatrix installation detail

Anymatrix clones the remote group inside `unitary/private`. MATLAB private
functions are not generally visible to tests in a nested `tests/` directory.
Our `test_run.m` therefore copies the MATLAB sources, metadata, docs and
regression suite to a temporary ordinary directory for the integration run,
executes the **same tests**, then restores the caller's path and removes the
temporary copy. No alternate test implementation is maintained.
