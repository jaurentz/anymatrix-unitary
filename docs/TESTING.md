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
| `test_structural_spectral_identities` | Analytical spectra, moments, limiting cases, determinant/symmetry and structural identities of added families |
| `test_opuc_additions` | Bernstein–Szegő density/moments and circular-beta coefficient sampling law |
| `test_opuc_extended` | Generalized circular Jacobi and Rogers–Szegő weights, finite Bernstein–Szegő, constant/periodic real-complex Verblunsky sequences, unitarity and finite closure |
| `test_floquet_families` | Floquet unitarity, exact limits, spin-layer isospectrality, and Kitaev particle-hole symmetry |
| `test_floquet_extended` | Rotor exact phases, embedded two-site gates, local dual-unitarity, XXZ Yang-Baxter, and global magnetization sectors |
| `test_quantum_expansion` | Dimer exact spectra, Jaynes–Cummings invariant blocks, Chalker–Coddington network and pi-pairing |
| `test_quantum_matrix_free` | Explicit/implicit multi-RHS forward and adjoint actions; sparse tridiagonal EXPM, MATLAB EXPMV and Chebyshev fallback; large bosonic, kicked-top, Jaynes and network tests |
| `test_spectral_examples` | Runs the two example scripts without figures in a full checkout; for Anymatrix installations without `examples/`, verifies equivalent generator invariants. In both contexts checks period-two CMV transfer discriminant, Dirac support/density, and tight-binding Bessel moment |
| `test_signal_processing_extended` | DFT fourth-power identity, DST I–IV, Walsh sequency, Slepian concentration and graph eigenspaces; only useful fast implicit actions |
| `test_signal_processing` | Exact identities, spectra, autocorrelations, multi-RHS forward/adjoint, implicit matrix materialization and large-order actions |
| `test_eigenvalue_problems` | Strict pencil Gram identities, numerical quotients, nonlinear secular determinants, exact spectra, block-CMV consistency and input validation |
| `test_polynomial_rotation_transforms` | Gauss quadrature identities, Wigner exact characters, spherical blocks and FFT-sinc adjoints |
| `test_catalog` | Catalog/documentation completeness, Anymatrix property declarations, and cross-checks for the Dirac operator API and catalog flag |
| `test_math_markup` | GitHub math-fence checks, guards against unsupported macros, malformed inline delimiters and nonstandard little-o complexity claims |

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

Wavelet-packet tests cover db2/db3/db4 and custom QMF filters, full and
incomplete binary trees, n=2 periodic wraparound, adjoint multi-RHS
operators, exact reconstruction, and preservation of the legacy DWT API.

The finite Christoffel OPUC helper tests positive reweighting,
preservation of support eigenphases versus changed cyclic weights,
exact support loss and dimension reduction, conjugation covariance,
operator/adjoint actions and annihilation guards. It delegates to
`inverse_opuc` rather than registering a redundant Anymatrix class.

The `test_integrable_spectral_flows.m` suite verifies the discrete finite CMV
positive-measure steps (#105), the Jacobi–Toda QR/Cayley isospectral evolution
(#104), and the unitary Hessenberg positive-diagonal QR/Schur flow (#103).
It separately tests fixed trace eigenphases, evolving e1 cyclic weights,
unitarity, structural zeros, real orthogonal subfamilies, and operator actions.
Large-time floating-point rank breakdown is documented rather than hidden.

`tests/test_geronimus_spectral_measure.m` validates the analytic infinite
constant-Verblunsky Geronimus e1 decomposition (continuous arc and possible
point mass), Poisson smoothing through the Schur/Carathéodory function,
finite-CMV weighted versus equal-mass trace differences, and smooth
Rogers–Szegő controls. These are not finite continuous spectral densities.

Spherical-cap Slepian tests (`tests/test_spherical_slepian.m`) verify real
orthogonality, analytic degree-one cap integrals, addition-theorem Shannon
trace, eigenvector ordering/signs, cap boundary cases, and the separation
between concentration factors and orthogonal-matrix eigenphases.

Fourier-Laguerre ball tests (`tests/test_fourier_laguerre_ball.m`) verify
independent generalized Laguerre radial function evaluations and physical
r^2dr quadrature normalization, zonal Legendre angular quadrature, square
radial-zonal Kronecker orthogonality, multi-RHS matrix-free forward/adjoint
and full-spherical rectangular sampling isometries, including P=L=1.

The finite Ablowitz–Ladik suite (`tests/test_abl_ladik_cmv.m`) tests exact
prescribed eigenphases, e1 weight evolution, the analytic finite Lax
commutator against independent time differences for harmonics k=1,2,3,
real Schur/dmKdV coefficient ODE reduction, fixed boundary phase, operator
forward/adjoint and multiple RHS. Normalized finite trace moments and
cyclic-vector moments are tested separately.

`tests/test_etna_historical_examples.m` verifies ETNA §7.1 common-draw
core columns, terminal factor, complex and real seeded realizations and
matrix-free adjoints. The modern §7.5 inverse-OPUC analogue is tested
against prescribed equispaced eigenphases, forward construction residual,
real conjugate-paired weights and finite trace/cyclic measures. This
is not a bitwise reproduction of the historical Fortran test setup.

`tests/test_discrete_fractional_fourier.m` independently checks Candan's
parity-ordered Hermite basis against DFT for odd and even orders, including
the double oscillator eigenvalue at n divisible by 4 and the omitted
index n-1. Regression tests also verify F(0)=I, F(1)=DFT, 4-periodicity,
index additivity, adjoint inversion, real symmetric number-operator
Hamiltonian, exact trace moments, and deterministic column-sign rules.

`tests/test_rational_arnoldi_circle.m` independently checks orthonormal
weighted rational samples and exact finite unitary multiplication closures,
the n-by-(n-1) Hessenberg pencil identity U*K=H, infinity/polynomial
reductions, inside/outside finite poles, nearly coincident poles,
conjugate-reciprocal pole pairs, positive e1 weights versus uniform
trace phases, multiple node orders and explicit rational breakdown.


`tests/test_gauss_jacobi_fast_subsets.m` compares FFT-based Gauss-
Chebyshev I (DCT-III) and II (DST-I) forward and adjoint operations
with the dense positive-leading Golub-Welsch transform for small
odd/even n, checks independent trigonometric reference matrices,
multi-RHS and complex-vector round trips, and uses n=16384 to
verify that the operator supports large sizes without a dense n-by-n
setup. General Jacobi and associated-Legendre parameters are explicitly
tested to reject the unsupported fast mode, not silently fall back to
quadratic/cubic operations.




`tests/test_unitary_hessenberg_pencil.m` tests independent
unitarity/upper-Hessenberg form of both strict factors, correct
prescribed finite subdiagonal projective poles, exact zero/infinite
(reducible) cases, near-reducible poles, generalized eigenvalue
agreement with eig(B'*A), seeded reproducibility, real conjugate
eigenvalue pairing, and n=4096 complex multiple-RHS forward/adjoint
O(n) core application. The auxiliary
`tests/test_hessenberg_unitary_generalized.m` tests the broader
Hessenberg–unitary-Hessenberg pencil with fixed poles but demonstrably
nonunitary numerator and quotient; this remains outside the matrix catalog.

`tests/test_real_spherical_harmonic_rotation.m` checks the
Condon–Shortley complex-to-real unitary basis, exact real
SO(3) orthogonality and inverse, arbitrary Euler-angle group
composition, singular Euler limits, spectral multiplicities
L-|m| and trace characters, and blockwise matrix-free
forward/adjoint actions for complex and real harmonics with
multiple RHS. Sparse explicit matrices are independently
compared to the cached O(L^3) block operators.

`tests/test_gauss_large_order_stability.m` exercises
Gauss–Hermite and generalized Gauss–Laguerre transforms with
alpha as close as -0.999, quadrature orders through n=384,
Sturm endpoint orientation, Gram/diagonalization residuals,
forward/adjoint multiple RHS round trips, double underflow in
tail weights versus finite log Christoffel weights, and
independently normalized physical Hermite/Laguerre degree-0
functions. It explicitly rejects an unimplemented `'operator'`
mode and tests that finite Gaussian nodes are not described
as eigenphases of Q.


`tests/test_spin_spherical_harmonic_transform.m` independently
compares the Wigner–Jacobi recurrence with dense
`wigner_rotation` exponentials for integer m,n and ell
through 6, checks strict square fixed-(s,m) Gauss–Jacobi
orthogonality, exact full-sphere rectangular Gram
A'*A=I and coefficient recovery for spins 0,+/-1,+/-2,
complex multi-RHS FFT forward/adjoint comparisons,
Condon–Shortley conjugation across opposite spin,
scalar Y_00/Y_11 normalization and a larger L=24
operator with no explicit sampling matrix.
