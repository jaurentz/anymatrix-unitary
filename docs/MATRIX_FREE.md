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

| Kernel | Families / mathematics | Storage |
| --- | --- | --- |
| `um_core_apply` | CMV, Hessenberg, inverse-Hessenberg, arbitrary twisted OPUC core chains and all Szegő coefficient families | O(n) |
| `um_monomial_apply` | Cyclic/permutation shifts, weighted shifts, diagonal phases, and Weyl–Heisenberg time-frequency shifts | O(n) |
| `um_fft_apply` | Circulants, CAZAC/Fresnel chirps, completed Hilbert transform, tight-binding cycle, kicked Harper/rotor, quantized baker | O(n) |
| `um_gate_apply`, `um_circuit_apply` | Two-site qubit gates, XXZ and dual-unitary Floquet circuits | O(1) per local gate |
| `um_spin_apply` | Floquet Ising chains | O(2^n) phase vector |
| `um_reflection_apply` | Grover rank-one reflection | O(n) temporary memory |
| `um_entangled_reflection_apply` | Maximally entangled rank-one phase update | O(d^2) state |
| `um_tensor_apply` | Repeated or independent-block Haar Kronecker products | O(d^t) state; O(t*d^(t+1)) arithmetic per RHS |
| `um_tensor_permute_apply` | Schur--Weyl tensor factor permutations by reshape/permute | O(d^t) state |
| `um_dirac_apply` | Two-component periodic Dirac automaton, shift and local mass mixing | O(n) for multiple RHS |
| `um_wavelet_apply` | Haar wavelet multilevel transform | O(n) |
| `um_dct_apply` | DCT I–IV using FFT even extensions | O(n) |
| `um_hartley_apply` | Hartley involution from normalized DFT and adjoint | O(n) |
| `um_affine_apply` | Affine Fourier: cached chirp factors and shared Fourier kernel | O(n) |
| `um_zak_apply` | Polyphase reshape and short FFTs | O(n) |
| `um_dft_apply` | Reuses normalized Fourier FFT, with explicit size contract | O(n) |
| `um_dst_apply` | DST-I odd FFT extension and DST-II–IV via DCT | O(n) |
| `um_walsh_apply` | Fast Walsh butterflies with sequency permutation | O(n) |
| `um_cayley_apply` | Cayley transform of a Jacobi recurrence using cached sparse triangular solves | O(n) for tridiagonal factors, modulo fill |
| `um_tridiag_exp_apply`, `um_kick_tridiag_apply` | Bose–Hubbard hopping and kicked-top spin rotations; diagonal phase followed by an Hermitian-tridiagonal exponential | O(n) coefficients; O(n·RHS) workspace |
| `um_jaynes_apply` | Analytic cavity-QED excitation blocks with uncoupled edge phases | O(n) |
| `um_cc_geometry`, `um_scattering_apply` | Disjoint two-port scattering networks, CC directed-link topology | O(n) |

The specialized `um_core_order` encodes the parity-dependent placement of the
finite CMV terminal phase and arbitrary snake-shaped insertion sequences.
Unreduced unitary Hessenberg matrices are **not** treated as sparse matrices
for the purposes of efficient matvecs. The OPUC descriptor supplies the
Verblunsky coefficients and core ordering; the kernel is shared.

`dual_unitary_circuit`, `xxz_floquet`, `floquet_ising`, and
`selfdual_kicked_ising` permit larger spin counts in implicit mode than in
the dense default. The state vector still has exponential length `2^n`.
The dense caps and other original parameter conventions remain in place.

All eight signal-processing generators accept the shared final `'operator'` argument.
DCT I–IV use FFT-based actions even for complex multiple RHS. CAZAC,
Fresnel, and Hilbert directly reuse the circulant FFT kernel, while
Weyl–Heisenberg shifts reuse `um_monomial_apply`.

The DFT, DST I–IV, and Walsh–Hadamard families additionally support
fast actions: `um_dft_apply` forwards to the shared Fourier kernel,
`um_dst_apply` uses odd FFT extension and shared DCT routines, and
`um_walsh_apply` uses fast Walsh butterflies. The DPSS/Slepian and
general graph Fourier transforms intentionally **do not offer an
`'operator'` output**: their explicit eigenvector matrices have no
generic lower-storage or lower-complexity implicit action.

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
