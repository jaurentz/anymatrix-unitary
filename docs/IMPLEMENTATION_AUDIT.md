# Matrix-free implementation audit

**Revision:** 13 November 2026. This inventory accompanies the matrix-free
function-handle refactor. It distinguishes fast **operator actions** from
structured mathematical descriptions that do not themselves yield fast
actions. MATLAB numeric matrices remain the default Anymatrix output.

## Implemented and shared

| Matrix families | Kernel | Memory without explicit matrix | Action per column |
| --- | --- | --- | --- |
| All ten OPUC / Verblunsky variants (CMV, Hessenberg, inverse, twisted) | `um_core_apply` | O(n) | O(n) |
| `phase_circulant`, `cycle_propagator`, `kicked_harper`, `kicked_rotor`, `quantum_baker` | `um_fft_apply` | O(n) | O(n log n) |
| `cycle_shift`, `weighted_shift`, `diagonal_phase` | `um_monomial_apply` | O(n) | O(n) |
| `grover_operator` | `um_reflection_apply` | O(n) | O(n) |
| `haar_wavelet` | `um_wavelet_apply` | O(1) coefficients | O(n) |
| `floquet_ising`, `selfdual_kicked_ising` | `um_spin_apply` | O(2^n) phases | O(n 2^n) |
| `dual_unitary_circuit`, `xxz_floquet` | `um_circuit_apply` | O(n) local gates | O(depth n 2^n) |
| `cayley_jacobi` | `um_cayley_apply` | O(n) sparse factors (typically) | O(n) |
| `kicked_bose_hubbard`, `kicked_top` | `um_kick_tridiag_apply` + `um_tridiag_exp_apply` (`expmv`/fallback) | O(n) coefficient data and O(n·RHS) workspace | O(n·m·RHS), m is Chebyshev truncation order |
| `jaynes_cummings` | `um_jaynes_apply` | O(n) excitation block coefficients | O(n·RHS) |
| `chalker_coddington` | `um_cc_geometry` + `um_scattering_apply` | O(n) scattering ports, phases, and blocks | O(n·RHS) |

The final `'operator'` argument returns `@(X,flag) ...`; valid flags are
`'notransp'`, `'transp'` (conjugate transpose), `'ctransp'`,
and `'size'`. Multiple right-hand sides are supported. Materialize
an n-by-n matrix using `op(eye(n),'notransp')` only when memory allows.

The Hermitian-tridiagonal exponential defaults to MATLAB's **`expmv`**
(R2023b and later), using a sparse tridiagonal matrix. For R2020b–R2023a,
a Chebyshev–Bessel fallback with numerical tail selection remains available.
No uniformly O(n)-time or certified arbitrary-tolerance claim is made.
The *explicit* Bose–Hubbard and kicked-top generators call MATLAB **`expm`**
on a sparse Hermitian-tridiagonal input; the result is generally dense.
Matrix-free mode never constructs that full exponential.
Both new sparse families retain their sparse numeric defaults.

## Already economical as MATLAB sparse arrays

`orthogonal_blocks`, `split_step_walk`, `dirac_automaton`,
`jaynes_cummings` and `chalker_coddington` have O(n) nonzeros. There is no compelling memory advantage from replacing
their normal sparse return values. Their sparse adjoints already use MATLAB
`A'` correctly.

CMV is also five-diagonal in its usual representation, but the same OPUC
coefficients may generate structurally full Hessenberg or snake-shaped
matrices. These are served by the same shared core kernel rather than
a CMV-only implementation.

## Deferred specialized operators

| Family | Obstacle | Proposed route |
| --- | --- | --- |
| `driven_kitaev` | Open-boundary BdG exponentials are dense even though the Hamiltonians are banded | Periodic translation-invariant case: FFT/2-by-2 momentum blocks; general open boundaries: Krylov/Chebyshev exponential action with quantified tolerance |
| `julia_unitary` | Generic contraction defect square roots are dense | Specialize scaled-unitary/permutation/circulant contractions with scalar defect operators; general contraction remains dense or requires structured square-root actions |
| `polar_unitary` | Generic matrix polar factor is dense, and an arbitrary source has no fast action | Preserve dense output; implement specialized normal/circulant sources only with a verifiable structure flag |
| `haar_unitary`, `haar_orthogonal` | Generic Haar samples have no exact short structured representation | Keep dense Gaussian QR construction. Random butterfly substitutes do not have the Haar law |

## General shapes beyond CMV

Watkins, *Some Perspectives on the Eigenvalue Problem*, SIAM Review 35
(1993), 430–471, develops the unit-circle Lanczos/Laurent-polynomial
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
