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

Already economical sparse families (e.g., the cyclic shifts, the Dirac
automaton, and the split-step walk) retain their sparse defaults and do not
need redundant custom matvec algorithms. General Haar and polar-factor
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
