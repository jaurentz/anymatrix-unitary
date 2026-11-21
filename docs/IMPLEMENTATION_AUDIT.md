# Matrix-free implementation audit

**Revision:** 21 November 2026. This inventory accompanies the matrix-free
function-handle refactor and subsequent generator additions. It distinguishes
implemented `'operator'` availability from fast **operator actions**, and both from
structured mathematical descriptions that do not themselves yield fast
actions. MATLAB numeric matrices remain the default Anymatrix output.

## Implemented and shared

All tabulated costs use uppercase `N` for the actual matrix order. Let `g`
denote the number of gates, `m` the truncation order for exponential actions,
and `RHS` the number of right-hand sides. Physical input sizes may differ.

| Matrix families | Kernel | Memory without explicit matrix | Action per column |
| --- | --- | --- | --- |
| OPUC / Verblunsky coefficient variants (CMV, Hessenberg, inverse, twisted) | `um_core_apply` | O(N) | O(N) |
| `phase_circulant`, `cycle_propagator`, `kicked_harper`, `kicked_rotor`, `quantum_baker` | `um_fft_apply` | O(N) | O(N log N) |
| `cycle_shift`, `weighted_shift`, `diagonal_phase` | `um_monomial_apply` | O(N) | O(N) |
| `grover_operator` | `um_reflection_apply` | O(N) | O(N) |
| `haar_wavelet` | `um_wavelet_apply` | O(1) coefficients | O(N) |
| `floquet_ising`, `selfdual_kicked_ising` | `um_spin_apply` | O(N) phases, N=2^n | O(N log N) |
| `dual_unitary_circuit`, `xxz_floquet` | `um_circuit_apply` | O(g) gates | O(gN), N=2^n |
| `cayley_jacobi` | `um_cayley_apply` | O(N) sparse factors (typically) | O(N) |
| `kicked_bose_hubbard`, `kicked_top` | `um_kick_tridiag_apply` + `um_tridiag_exp_apply` (`expmv`/fallback) | O(N) coefficient data and O(N·RHS) workspace | O(N·m·RHS), m is Chebyshev truncation order |
| `jaynes_cummings` | `um_jaynes_apply` | O(N) excitation block coefficients | O(N·RHS) |
| `chalker_coddington` | `um_cc_geometry` + `um_scattering_apply` | O(N) scattering ports, phases, and blocks | O(N·RHS) |
| `dirac_automaton` | `um_dirac_apply` | O(N) state, constant local parameters | O(N·RHS) |
| `discrete_fourier`, `orthogonal_dct`, `discrete_sine`, `discrete_hartley`, `walsh_hadamard`, `affine_fourier`, `finite_zak` | Shared FFT, DCT/DST, Walsh, chirp and Zak kernels | O(N) or less operator data, O(N) FFT workspace | O(N log N·RHS) |
| `orthogonal_butterfly` | `um_butterfly_apply` | O(N) rotation angles | O(N log N·RHS) |
| `haar_tensor_power`, `haar_block_tensor_power`, `fourier_laguerre_ball`, `tensor_permutation` | `um_tensor_apply` / `um_tensor_permute_apply` | Small factors/permutation; full state O(N) | Factor-size dependent; O(N) for permutation |
| `daubechies_wavelet` | `um_daubechies_apply` / `um_wavelet_packet_apply` | Fixed taps, levels/tree | O(N·RHS) for standard decimated wavelets; packet tree dependent |
| `gauss_jacobi_transform` at alpha=beta=+/-1/2 **only** | `um_gauss_jacobi_cheb_apply` | O(N) workspace/signs | O(N log N·RHS) |
| `unitary_hessenberg_pencil` (strict unitary case) | `um_hessenberg_pencil_apply` | O(N) cores and diagonal phases | O(N·RHS) |
| `optical_mesh` (full universal mesh) | `um_optical_mesh_apply` | O(N^2) couplers | O(N^2·RHS), **not a fast full mesh** |

The table is a summary, not an exhaustive list. In this document `N`
denotes the actual order of the unitary matrix, while physical dimensions,
qubit counts, ranks, and local block sizes are stated separately.
A few operators have restricted parameter domains: in particular,
`gauss_jacobi_transform` rejects `'operator'` outside the two Chebyshev
cases. `optical_mesh` implements an operator handle, but full universality
retains quadratic work and parameters. The source of truth for **registered**
operator availability is [metadata/catalog.json](../metadata/catalog.json).
No fast claim follows from the flag alone; the corresponding MATLAB source,
forward/adjoint tests, and parameter conditions must be checked.

The final `'operator'` argument returns `@(X,flag) ...`; valid flags are
`'notransp'`, `'transp'` (conjugate transpose), `'ctransp'`,
and `'size'`. Multiple right-hand sides are supported. Materialize
an $N\times N$ matrix using `op(eye(N),'notransp')` only when memory allows.

The Hermitian-tridiagonal exponential defaults to MATLAB's **`expmv`**
(R2023b and later), using a sparse tridiagonal matrix. For R2020b–R2023a,
a Chebyshev–Bessel fallback with numerical tail selection remains available.
No uniformly $O(N)$-time or certified arbitrary-tolerance claim is made.
The *explicit* Bose–Hubbard and kicked-top generators call MATLAB **`expm`**
on a sparse Hermitian-tridiagonal input; the result is generally dense.
Matrix-free mode never constructs that full exponential.
Jaynes–Cummings and Chalker–Coddington retain their sparse numeric defaults.

## Already economical as MATLAB sparse arrays

`orthogonal_blocks`, `split_step_walk`, `dirac_automaton`,
`jaynes_cummings` and `chalker_coddington` have $O(N)$ nonzeros. There is no compelling memory advantage from replacing
their normal sparse return values. Their sparse adjoints already use MATLAB
`A'` correctly. Unlike `orthogonal_blocks` and `split_step_walk`,
`dirac_automaton` also exposes a separate O(N) forward/adjoint
`'operator'` action; its catalog flag must remain true.

CMV is also five-diagonal in its usual representation, but the same OPUC
coefficients may generate structurally full Hessenberg or snake-shaped
matrices. These are served by the same shared core kernel rather than
a CMV-only implementation.

## Dense or restricted cases without a verified fast full operator

Full DPSS and arbitrary graph-Fourier eigenvector matrices are dense and
have no verified generic fast complete-basis operator in this repository.
A fast leading-Slepian-subspace projection or a special graph topology
must not be confused with a generic fast complete eigenvector transform.
The general Gauss--Jacobi, Gauss--Hermite, Gauss--Laguerre and
associated-Legendre evaluation matrices use dense algorithms here;
**only** the two Gauss--Chebyshev specializations noted above have
FFT-fast operator actions. For details see
[STRUCTURED_OPERATORS.md](STRUCTURED_OPERATORS.md) and
[MATRIX_FREE.md](MATRIX_FREE.md).

## Deferred specialized operators

| Family | Obstacle | Proposed route |
| --- | --- | --- |
| `driven_kitaev` | Open-boundary BdG exponentials are dense even though the Hamiltonians are banded | Periodic translation-invariant case: FFT/2-by-2 momentum blocks; general open boundaries: Krylov/Chebyshev exponential action with quantified tolerance |
| `julia_unitary` | Generic contraction defect square roots are dense | Specialize scaled-unitary/permutation/circulant contractions with scalar defect operators; general contraction remains dense or requires structured square-root actions |
| `polar_unitary` | Generic matrix polar factor is dense, and an arbitrary source has no fast action | Preserve dense output; implement specialized normal/circulant sources only with a verifiable structure flag |
| `haar_unitary`, `haar_orthogonal` | Generic Haar samples have no exact short structured representation | Keep dense Gaussian QR construction. Random butterfly substitutes do not have the Haar law |

## General shapes beyond CMV

Watkins, *Some Perspectives on the Eigenvalue Problem*, *SIAM Review* **35**
(1993), 430–471, DOI: https://doi.org/10.1137/1035090,
develops the unit-circle Lanczos/Laurent-polynomial
perspective. The current core engine handles any single-core insertion
ordering given by a `twists` vector (snake-shaped factorizations) and a
terminal boundary phase. Additional block, repeated-core, or multi-snake
layouts require a general factor schedule and explicit mathematical
normalization; they are not claimed as implemented just because they
are based on isometric Arnoldi.

## Validation rules

1. Test `U(X,'notransp') == A*X` and `U(X,'transp') == A'*X` for complex
   matrices and multiple right-hand sides, with numerical tolerances.
2. Independently test `U'*U=I` and the adjoint inner-product identity.
3. Confirm finite boundary phases, parity, and factor order.
4. Never infer an exact finite eigenphase law from an asymptotic limit.
5. Distinguish `n` physical sites from `2^n` Hilbert-space dimension.

See [matrix-free tests](../tests/test_matrix_free.m),
[quantum-operator tests](../tests/test_quantum_matrix_free.m), and
[operator interface](MATRIX_FREE.md).
