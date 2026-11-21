# Matrix-free operator convention

The Anymatrix `unitary` group retains numeric MATLAB matrices as its default
result. Selected families also accept a final output argument `'operator'`
(or `'implicit'`), returning a plain MATLAB function handle:

```matlab
A = anymatrix('unitary/kicked_rotor',16,5,.2,0);
U = anymatrix('unitary/kicked_rotor',16,5,.2,0,'operator');
x = randn(16,3)+1i*randn(16,3);
y = U(x,'notransp');  % A*x
z = U(x,'transp');    % A'*x, conjugate transpose
z = U(x,'ctransp');   % same
n = U([],'size');     % [16 16]
A = U(eye(16),'notransp');   % explicit, possibly large
```

**Important:** The iterative-solver convention `'transp'` means the
**conjugate transpose** of a complex matrix, not MATLAB's nonconjugating
`A.'`. All actions support multiple right-hand sides and preserve the
factor order under adjoints. No custom MATLAB matrix class is required.

## Shared kernels

The table distinguishes **stored operator data** from per-vector work.
Here `N` is the actual matrix order, `q` is the qubit count when `N=2^q`,
`g` is the number of gates or cores, and `d,t` describe a
tensor-product dimension `N=d^t`. Workspace for an arbitrary vector is
at least `O(N)` when a full state vector must be read. Complexity claims
are per right-hand side, after any separately described setup.

| Kernel(s) | Families / representation | Stored operator data | Action per RHS |
| --- | --- | --- | --- |
| `um_core_apply` | CMV, unitary Hessenberg/inverse-Hessenberg and twisted OPUC | `O(N)` coefficients/cores | `O(N)` |
| `um_monomial_apply` | Diagonal phases, permutations, Weyl–Heisenberg shifts | `O(N)` or less | `O(N)` |
| `um_fft_apply` | Circulants, CAZAC/Fresnel/Hilbert, kicked rotor/Harper, baker | Typically `O(N)` | `O(N log N)` |
| `um_dft_apply`, `um_dct_apply`, `um_dst_apply`, `um_hartley_apply`, `um_affine_apply`, `um_zak_apply` | DFT, cosine/sine, Hartley, affine Fourier and Zak | `O(N)` or less; FFT workspace `O(N)` | `O(N log N)` |
| `um_walsh_apply` | Normalized Walsh–Hadamard | `O(N)` or less | `O(N log N)` |
| `um_wavelet_apply`, `um_daubechies_apply`, `um_wavelet_packet_apply` | Haar, fixed-tap orthonormal wavelets, wavelet packets | Filters plus transform schedule | `O(N)` for standard fixed-tap decimated wavelets; packets depend on tree |
| `um_butterfly_apply` | Orthogonal recursive rotation butterfly | `O(N)` angles | `O(N log N)` |
| `um_gate_apply`, `um_circuit_apply` | Local quantum gates and complete circuits | Gate descriptions `O(g)` | `O(gN)` for full `q`-qubit state |
| `um_spin_apply` | Floquet Ising families | State-length phase data, typically `O(N)` | Depends on spin layout; commonly `O(qN)` |
| `um_reflection_apply`, `um_entangled_reflection_apply` | Rank-one reflections and maximally entangled phase updates | `O(N)` state/vector parameters | `O(N)` |
| `um_tensor_apply` | Haar tensor powers and radial-zonal Fourier–Laguerre product | `O(td^2)` small factors in the uniform case | `O(tdN)`; arbitrary vector remains length `N` |
| `um_tensor_permute_apply` | Tensor-factor permutations | Small permutation descriptor | `O(N)` |
| `um_dirac_apply` | Periodic Dirac quantum automaton | Constant local parameters plus order | `O(N)` |
| `um_cayley_apply` | Tridiagonal Hermitian Cayley transform | `O(N)` factors | `O(N)` after banded factorization |
| `um_tridiag_exp_apply`, `um_kick_tridiag_apply` | Tridiagonal Hamiltonian exponential actions | `O(N)` coefficients and `O(N)` workspace | Accuracy/time-dependent; no uniform `O(N)` bound |
| `um_jaynes_apply` | Jaynes–Cummings excitation blocks | `O(N)` local coefficients | `O(N)` |
| `um_cc_geometry`, `um_scattering_apply` | Chalker–Coddington network | `O(N)` graph and scattering data | `O(N)` |
| `um_gauss_jacobi_cheb_apply` | **Only** Gauss–Chebyshev I/II weighted evaluation matrices | `O(N)` signs / FFT workspace | `O(N log N)`; other Jacobi parameters reject `'operator'` |
| `um_hessenberg_pencil_apply` | Strict independently unitary Hessenberg pairs; quotient `B'*A` | `O(N)` core and phase data | `O(N)` |
| `um_optical_mesh_apply` | Unrestricted Reck/Clements lossless optical mesh | `O(N^2)` coupler parameters | `O(N^2)`; valid operator API, **not** a fast full mesh |

This table summarizes selected verified kernels and is not an exhaustive
generator inventory. The machine-readable source for the registered
`'operator'` interface is [metadata/catalog.json](../metadata/catalog.json);
`matrix_free_supported=true` means an operator output is available for
the documented parameter regime, **not** that its application is
subquadratic or that it is available for all parameter values.
For complete classification of sparse/compact/fast/dense structures,
see [STRUCTURED_OPERATORS.md](STRUCTURED_OPERATORS.md), and for status
by implementation see [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md).

The specialized `um_core_order` encodes the parity-dependent placement of the
finite CMV terminal phase and arbitrary snake-shaped insertion sequences.
Unreduced unitary Hessenberg matrices are **not** treated as sparse matrices
for the purposes of efficient matvecs. The OPUC descriptor supplies the
Verblunsky coefficients and core ordering; the kernel is shared.

`dual_unitary_circuit`, `xxz_floquet`, `floquet_ising`, and
`selfdual_kicked_ising` permit larger spin counts in implicit mode than in
the dense default. The state vector still has exponential length `2^n`.
The dense caps and other original parameter conventions remain in place.

The original eight supported signal-processing generators accept the
shared final `'operator'` argument; additional transforms are covered below.
DCT I–IV use FFT-based actions even for complex multiple RHS. CAZAC,
Fresnel, and Hilbert directly reuse the circulant FFT kernel, while
Weyl–Heisenberg shifts reuse `um_monomial_apply`.

The DFT, DST I–IV, and Walsh–Hadamard families additionally support
fast actions: `um_dft_apply` forwards to the shared Fourier kernel,
`um_dst_apply` uses odd FFT extension and shared DCT routines, and
`um_walsh_apply` uses fast Walsh butterflies. The DPSS/Slepian and
general graph Fourier transforms intentionally **do not offer an
`'operator'` output**: their explicit eigenvector matrices have no
verified generic lower-storage or lower-complexity **full-basis** implicit
action in this repository. This is not a mathematical impossibility claim.

Already economical sparse families (e.g., the cyclic shifts and the split-step
walk) retain their sparse defaults. The Dirac automaton also offers a direct
matrix-free forward/adjoint interface through `um_dirac_apply`, with O(n)
arithmetic and no assembled sparse matrix. General Haar and polar-factor
matrices remain dense. The kicked-top generator now reuses the Hermitian-tridiagonal exponential
kernel developed for the two-mode Bose–Hubbard dimer. In MATLAB R2023b and later,
its implicit mode assembles a sparse Hermitian-tridiagonal generator and
uses MathWorks `expmv` separately for each right-hand side. The repository's
R2020b minimum is preserved by a Chebyshev–Bessel fallback if `expmv`
is unavailable. The fallback is selected using a numerical Bessel-tail
criterion, not a certified tolerance bound. With either method, arithmetic
cost depends on the spectral time range; no uniformly linear-time claim
is made.

The **explicit** Bose–Hubbard and kicked-top paths also construct sparse
tridiagonal generators, then call MATLAB `expm`. Their exponential matrices
are generally **dense** and must be stored as such. Sparse generator storage
does not make `expm(H)` a sparse matrix.

Jaynes–Cummings and Chalker–Coddington retain their economical **sparse
numeric defaults**. Their optional implicit modes use analytic two-state
blocks and directed-node scattering, respectively; neither requires
materializing even the sparse matrix. In both, the action and adjoint
require O(n·RHS) arithmetic and working memory (plus O(n) operator data).
The driven-Kitaev BdG family still needs a validated high-order exponential
action before an implicit interface can be promised.

## Testing

Run `test_run` or `runtests('tests/test_matrix_free.m')`. Tests cover
complex adjoints, multiple right-hand sides, implicit/explicit agreement,
unitarity, and the Saraceno Fourier convention. The existing test suite
provides independent structural and spectral checks for the same matrices.
`tests/test_quantum_matrix_free.m` checks the shared tridiagonal exponential
against independent dense EXPM, as well as multi-RHS forward/adjoint equivalence
for the three added families and the kicked top.

## Haar tensor products and exact ensemble spectral references

The new `haar_tensor_power(d,t,seed,'operator')` and
`haar_block_tensor_power(d,t1,t2,seed,'operator')` preserve the seed and
allow dimensions larger than the explicit order-4096 cap. Their forward and
conjugate-transpose actions apply the small Haar factors along tensor axes.
They require the full *state vector* of length d^t, but no dense d^t-by-d^t
Kronecker matrix. The `tensor_permutation(d,pi,'operator')` mode similarly
reshapes and permutes tensor dimensions, without assembling a sparse matrix.

The terminal `'spectral'` output for Haar SU(n), unimodular COE and CSE
returns the exact *ensemble expectation*, including finite-size deviations
from the uniform phase law. The COE n=2 exceptional case has normalized
density |sin(theta)|/4. Tensor-power `'spectral'` output instead needs a
fixed SEED and returns exact eigenphase products and moments of that
*particular sample*, computed from only its d-by-d base factors. These are
not ensemble means. No nonuniform finite-n ensemble density is asserted
for the tensor-power constructions.

The square radial-zonal `fourier_laguerre_ball` transform uses `um_tensor_apply`
with precomputed dense Gaussian radial (P by P) and zonal Legendre (L by L)
factors. Applying its P*L orthogonal tensor costs O(P*L*(P+L)) per RHS and
stores O(P^2+L^2) coefficients after O(P^3+L^3) setup; there is no claim
of a linear- or FFT-time algorithm. Its separately exposed full-spherical
sampling output is rectangular rather than unitary.

## Gauss-Jacobi: exact FFT subset only

`gauss_jacobi_transform(n,alpha,beta,'operator')` supplies a real
O(n log n) FFT-based, complex-linear forward/adjoint **only** for
alpha=beta=-1/2 (Gauss-Chebyshev I) and alpha=beta=+1/2
(Gauss-Chebyshev II). At ascending Gaussian quadrature nodes,
their positive-leading orthonormal polynomial evaluation matrices
are precisely `DCTIII*diag((-1).^k)` and
`DSTI*diag((-1).^k)`, respectively. The shared
`um_gauss_jacobi_cheb_apply` composes existing FFT kernels with
O(n)-storage alternating signs. The default dense generator remains
unchanged and is used as the independent numerical reference.
For **all other** Jacobi exponents, the operator mode fails explicitly;
a cached dense eigensystem is not a fast matrix-free implementation.
The exact transform is orthogonal in arithmetic; the floating-point FFT
round trip is numerically approximate. See issue #39 for general fast
Jacobi and associated-Legendre transforms.




## Strict core-factored unitary Hessenberg pencils

`unitary_hessenberg_pencil(n,poles,seed,variant,'operator')`
returns a function handle for U=B'*A, with two independent
unitary upper-Hessenberg adjacent-core products A,B and N-1
user-prescribed projective pole ratios. The implementation
`um_hessenberg_pencil_apply` applies the forward quotient
with an A core sweep, O(n) left diagonal phases, and the
B adjoint core sweep; adjoint actions reverse the factor
order and conjugate phases. **Both cost O(n)** arithmetic
and use O(n) stored parameters (not dense A or B).
The real orthogonal specialization is provided by
`variant='real'`. Exactly zero/infinite poles require
a reducible numerator/denominator factor respectively.
The optional `'pencil'` mode explicitly materializes
A,B and the unitary quotient for modest n. The related
`um_hessenberg_unitary_generalized` helper returns a
**nonunitary** generalized problem and MUST NOT be treated
as a unitary matrix-free generator.

## Wigner direct-sum SO(3) rotations: blockwise, not fully fast

For `spherical_harmonic_rotation(L,alpha,beta,gamma,'operator')`
and `real_spherical_harmonic_rotation(...,'operator')`,
`um_wigner_direct_sum_apply` caches one independently computed
unitary/real orthogonal matrix of order 2*ell+1 per degree ell<L.
The full representation has order N=L^2; blockwise forward,
adjoint and complex multi-RHS applications need O(L^3)=O(N^(3/2))
operations and the same order of operator storage, rather than
O(N^2) for a generic dense full matrix. **Setup still costs O(L^4)**
via dense per-block matrix exponentials, and this does NOT implement
a fast O(ell log ell) single-degree Wigner D. The real basis is
unitarily conjugated from complex Condon–Shortley harmonics using
`um_wigner_real_basis`, not `real(D)`. The real family also
provides `'basis'` changes and `'spectral'` exact phase counts.
