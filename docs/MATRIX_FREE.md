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
| `um_fft_apply` | Circulants, tight-binding cycle, kicked Harper/rotor, quantized baker | O(n) |
| `um_gate_apply`, `um_circuit_apply` | Two-site qubit gates, XXZ and dual-unitary Floquet circuits | O(1) per local gate |
| `um_spin_apply` | Floquet Ising chains | O(2^n) phase vector |
| `um_reflection_apply` | Grover rank-one reflection | O(n) temporary memory |
| `um_wavelet_apply` | Haar wavelet multilevel transform | O(n) |
| `um_cayley_apply` | Cayley transform of a Jacobi recurrence using cached sparse triangular solves | O(n) for tridiagonal factors, modulo fill |

The specialized `um_core_order` encodes the parity-dependent placement of the
finite CMV terminal phase and arbitrary snake-shaped insertion sequences.
Unreduced unitary Hessenberg matrices are **not** treated as sparse matrices
for the purposes of efficient matvecs. The OPUC descriptor supplies the
Verblunsky coefficients and core ordering; the kernel is shared.

`dual_unitary_circuit`, `xxz_floquet`, `floquet_ising`, and
`selfdual_kicked_ising` permit larger spin counts in implicit mode than in
the dense default. The state vector still has exponential length `2^n`.
The dense caps and other original parameter conventions remain in place.

Already economical sparse families (e.g., the cyclic shifts, the Dirac
automaton, and the split-step walk) retain their sparse defaults and do not
need redundant custom matvec algorithms. General Haar and polar-factor
matrices remain dense. The BdG driven-Kitaev and kicked-top families still
need dedicated, validated high-order exponential actions before promising
a fast implicit implementation.

## Testing

Run `test_run` or `runtests('tests/test_matrix_free.m')`. Tests cover
complex adjoints, multiple right-hand sides, implicit/explicit agreement,
unitarity, and the Saraceno Fourier convention. The existing test suite
provides independent structural and spectral checks for the same matrices.
