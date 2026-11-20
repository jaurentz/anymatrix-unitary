# Anymatrix Unitary

A MATLAB collection of finite **unitary and orthogonal test matrices**, designed as
an external group for [Anymatrix](https://github.com/north-numerical-computing/anymatrix).
The repository is organized around reproducible mathematical constructions,
known spectral information, and executable regression tests. It contains
89 generators; additional benchmark families are under development.

**Status:** development collection with MATLAB-only tests. Automated MATLAB
regression and Anymatrix integration tests run on GitHub Actions after every
push or pull request. See the [Actions tab](https://github.com/jaurentz/anymatrix-unitary/actions)
for the current branch status.

A structured registry of verified DOI and archive/publisher permalinks is available in [metadata/references.json](metadata/references.json).

## Clone and test locally

Requires **MATLAB R2020b or later**, with no additional MATLAB toolboxes for
these generators or tests.
The tridiagonal matrix-free action uses MATLAB `expmv` on R2023b and later;
on R2020b–R2023a the shared kernel falls back to Chebyshev–Bessel.
For explicit Bose–Hubbard and kicked-top matrices, `expm` receives a sparse
Hermitian-tridiagonal generator but produces a generally dense propagator.

```bash
git clone https://github.com/jaurentz/anymatrix-unitary.git
cd anymatrix-unitary
```

From MATLAB, with the current folder set to the cloned repository:

```matlab
addpath(pwd)
results = test_run;             % full test suite; throws on any failure
```

Or run specific tests:

```matlab
results = runtests(fullfile(pwd,'tests'));
assertSuccess(results)
```

> A *local* commit does not trigger GitHub Actions until pushed to GitHub.
> The workflow runs on every **push** and **pull request**, and can be started
> manually in the Actions tab.

## Matrix generator interfaces

The collection currently has **89 public generators**. Each has a complete
entry with formulas, spectral information, applications, parameter conventions,
and references in the [matrix encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).
The [mathematical notes](docs/MATHEMATICS.md) explain the shared OPUC engine.
The [vibron ESQPT example](examples/vibron_esqpt_sector_gaps.m) compares U(3) angular-momentum and U(2) Lipkin parity gaps.

| Category | Generators |
|---|---|
| Prescribed spectra and circulants | `cycle_shift`, `weighted_shift`, `diagonal_phase`, `phase_circulant`, `cycle_propagator`, `orthogonal_blocks` |
| OPUC / Szegő / core matrices | `finite_cmv`, `opuc_matrix`, `jacobi_szego`, `bernstein_szego`, `circular_beta`, `circular_jacobi`, `rogers_szego`, `geronimus_cmv`, `periodic_cmv`, `bernstein_szego_finite`, `lee_yang_ising`, `scattering_zipper`, `rational_opuc`, `unitary_rank_one_update`, `etna_hessenberg` |
| Random and signal-processing | `haar_unitary`, `haar_orthogonal`, `haar_wavelet`, `daubechies_wavelet`, `haar_special_unitary`, `haar_unimodular_coe`, `haar_unimodular_cse`, `haar_tensor_power`, `haar_block_tensor_power` |
| Signal-processing transforms | `orthogonal_dct`, `discrete_hartley`, `cazac_circulant`, `fresnel_transform`, `affine_fourier`, `finite_zak`, `time_frequency_shift`, `periodic_hilbert`, `discrete_fourier`, `discrete_sine`, `walsh_hadamard`, `dpss_transform`, `graph_fourier`, `blaschke_potapov`, `orthogonal_butterfly`, `complex_hadamard`, `krawtchouk_transform`, `hahn_transform`, `dual_hahn_transform` |
| Riemannian geometry | `sphere_geodesic_transport` |
| Quantum maps, walks, and algorithms | `split_step_walk`, `dirac_automaton`, `grover_operator`, `quantum_baker`, `quantum_cat_map`, `kicked_harper`, `kicked_top`, `extended_lipkin_propagator`, `vibron_2d_esqpt`, `agassi_propagator`, `floquet_ising`, `selfdual_kicked_ising`, `driven_kitaev`, `kicked_rotor`, `dual_unitary_circuit`, `xxz_floquet`, `kicked_bose_hubbard`, `jaynes_cummings`, `chalker_coddington`, `unitary_almost_mathieu`, `szegedy_walk`, `coined_graph_walk`, `metric_quantum_graph` |
| Matrix functions and dilations | `cayley_jacobi`, `polar_unitary`, `julia_unitary`, `entangled_reflection` |
| Quantum group actions | `tensor_permutation` |
| Nonlinear unitary dynamical systems | `lohe_matrix_flow` |
| Strict unitary eigenvalue pencils | `mobius_unitary_pencil`, `unitary_delay_pencil`, `block_cmv_pencil`, `floquet_factor_pencil`, `cayley_hermitian_pencil`, `left_weighted_unitary_pencil` |

Generators return **actual finite square matrices by default**. Selected
generators also accept `'operator'` as their final argument and return a
function handle supporting forward and conjugate-transpose actions, multiple
right-hand sides, and dimension queries. See [matrix-free documentation](docs/MATRIX_FREE.md).
`split_step_walk`, `dirac_automaton`, `driven_kitaev`, and `julia_unitary` return order `2*n`
for an input physical size `n`. `floquet_ising` and `selfdual_kicked_ising` return order `2^n`, capped at 11 spins for dense assembly. `dual_unitary_circuit` and `xxz_floquet` also return order `2^n`, capped at 10 spins. `driven_kitaev` is a 2n-dimensional BdG operator, not the full Fock-space propagator. `kicked_bose_hubbard(n)` has order `n+1` for n bosons; `jaynes_cummings(n)` has order `2*(n+1)` for a photon cutoff n; `chalker_coddington(nx,ny)` has order `4*nx*ny`. `lee_yang_ising(n)` has order `2*n` and represents the zeros of a classical Ising partition function through a twisted Floquet-CMV matrix. `scattering_zipper(n,q)` has order `n*q` (n even) and `unitary_almost_mathieu(n)` has order `2*n`; irrational flux in the latter produces a finite periodically closed window, not a literally periodic irrational lattice. `szegedy_walk(n,P)` has order `n^2` (use `operator` for large n). `coined_graph_walk(n,A)` and `metric_quantum_graph(n,L,k)` have order twice the number of undirected edges. The latter exposes analytic `T(k)=I-U(k)` through `secular` mode; its metric wavenumbers are not eigenphases of one fixed matrix. `agassi_propagator(j,epsilon,g,V,t,particles)` evolves a fermionic Agassi model in full Fock space (order `2^(4*j)`) or a fixed-particle-number sector; `j=1` reproduces the paper's four-orbital simulation. `vibron_2d_esqpt(N,xi,t)` has order `(N+1)*(N+2)/2`, with real-symmetric tridiagonal Hamiltonian sectors indexed by angular momentum `ell=-N:N`; a selected sector is available as `vibron_2d_esqpt(N,xi,t,ell)`. The `hamiltonian` mode scales to larger N and the solvable `xi=0,1` limits have exact finite eigenphases. `extended_lipkin_propagator(n,lambda,alpha,t)` has order `n+1`, with a sparse symmetric pentadiagonal Hamiltonian exposed by terminal `hamiltonian` mode and a generally dense unitary propagator. `quantum_baker` requires even order;
`haar_wavelet` requires power-of-two order. MATLAB R2020b or newer is required.

**2026 spectral and implicit tensor modes:** Try
`S=haar_special_unitary(6,[],'spectral')` for the exact finite-order
ensemble-mean density, `S=haar_unimodular_coe(2,[],'spectral')` for the
exceptional |sin(theta)|/4 law, or
`O=haar_tensor_power(2,13,19,'operator')` to apply an order-8192 unitary
without building its dense matrix. Seeded tensor-product `'spectral'` mode
returns exact eigenvalue products and normalized moments for the chosen
sample. The `tensor_permutation` generator also supports implicit tensor
index reordering. See [matrix-free documentation](docs/MATRIX_FREE.md).

**2026 quantum-information additions:** Haar SU(n) and the unimodular COE/CSE ensembles have finite-order analytic expected eigenphase laws, in contrast to ordinary CUE's uniform mean density. Tensor powers, tensor permutations and the entangled reflection have analytic sample spectra or moments. The existing `time_frequency_shift(n,p,q)` already provides Weyl--Heisenberg 1-design members; it has not been duplicated. See [the encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).

**Spectral reference data:** exact eigenvalues are provided or described for
several families; other families have ensemble laws, limiting densities, or
no universal law. The encyclopedia distinguishes these cases. Do not infer a
limiting eigenphase distribution from unitarity alone.

## Eigenvalue-problem subclass

Six new generators, together with the existing `finite_cmv`, represent strict pencils `A-z*B` for
which `B\\A` is unitary in the ordinary Euclidean inner product.
Each returns the finite numeric unitary by default, with
terminal `'pencil'` mode returning a structure `P` with fields `P.U`, `P.A`,
and `P.B`. The normal two-output `[U,properties]` contract is unchanged.
The structural
criterion is `A*A' == B*B'` (up to floating-point error); the
coefficients need not individually be unitary.

~~~matlab
U = unitary_delay_pencil(3,[2 3 4]');
P = unitary_delay_pencil(3,[2 3 4]','pencil');
assert(norm(full(P.A*P.A'-P.B*P.B'),'fro') < 1e-11);
assert(norm(full(U-P.B\\P.A),'fro') < 1e-11);
S = mobius_unitary_pencil(40,.8,'spectral');
~~~

The separate [eigenvalue-problem registry](metadata/eigenvalue_problems.json)
also records transcendental quantum-graph equations, rational inner
functions, matrix paraorthogonal polynomials and NLEVP palindromic
benchmarks. These are **not** falsely advertised as unitary matrices:
reciprocal or circle-valued eigenvalues alone do not imply a strict
unitary pencil. See the [encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).

## Analytic spectral references

Every public generator also accepts terminal 'spectral' mode and returns a
classified MATLAB structure with exact finite spectra, limiting trace laws,
ensemble expectations, or underlying OPUC weights, where analytically known.
This mode does not construct a matrix or advance the random-number generator.

~~~matlab
S = anymatrix('unitary/dirac_automaton',100,.45,'spectral');
theta = linspace(-pi,pi,500);
plot(theta,S.limit.density(theta));
~~~

See [spectral reference interface](docs/SPECTRAL_REFERENCE.md) for the
standardized data contract, distinctions between trace and vector measures,
and complete positional argument rules.

## Optional matrix-free operators

```matlab
U = anymatrix('unitary/phase_circulant',4096,linspace(-2,2,4096),'operator');
X = randn(4096,2);
Y = U(X,'notransp');  % U*X
Z = U(Y,'transp');    % U'*Y (conjugate transpose)
A = U(eye(4096),'notransp'); % only if the explicit matrix fits memory
```

The group uses one shared core kernel for the Verblunsky-based models, one
FFT kernel for Fourier-factor models, and shared local-gate kernels for
many-body circuits, and the kicked top, Bose-Hubbard dimer, Jaynes-Cummings blocks and Chalker-Coddington network. Default numeric output and legacy argument ordering
remain unchanged. See [operator convention and coverage](docs/MATRIX_FREE.md).

## Measure-first Jacobi--Szegő construction

The authoritative object for this subfamily is a normalized probability
measure, rather than a matrix representation:

```matlab
mu = szego_jacobi_measure(0,0);      % Legendre weight
alpha = szego_verblunsky(mu,15);   % interior coefficients alpha_0,...,alpha_14
[U, coefficients] = szego_matrix(mu,16,'cmv',1);
H = szego_matrix(mu,16,'hessenberg',1);
L = szego_matrix(mu,16,'inverse_hessenberg',1);
T = szego_matrix(mu,16,'twisted',1);
[z,w] = szego_quadrature(U);

assert(norm(full(U'*U-speye(16)),'fro') < 1e-12)
```

The underlying Jacobi weight is proportional to `(1-t)^a (1+t)^b` on
`[-1,1]`, `a,b > -1`. The symmetric Szegő lift under
$x=2t=z+z^{-1}$ has the explicit Verblunsky sequence

```math
\begin{aligned}
\alpha_{2k}&=\frac{b-a}{2k+a+b+2},\\
\alpha_{2k+1}&=-\frac{a+b+1}{2k+a+b+3}.
\end{aligned}
```

For a finite `n`-by-`n` realization, the first `n-1` coefficients are
retained and an independent unimodular `beta` closes the system. The
`szego_quadrature` output `w_j=|v_j(1)|^2` describes the cyclic-vector
spectral measure, **not** the equally weighted trace spectral measure.
See [docs/MATHEMATICS.md](docs/MATHEMATICS.md) for the normalization,
Geronimus relations, and core-ordering conventions.

### Representation conventions

`'hessenberg'`, `'inverse_hessenberg'`, `'cmv'`, and `'twisted'` share
the same **interior Verblunsky coefficients** and finite boundary phase.
`'twisted'` accepts `n-2` binary insertion bits, with `1` meaning prepend
the next elementary core and `0` meaning append it. A default alternating
pattern is used when omitted. All four representations are isospectral and
have matching distinguished-vector moments in the stated convention, but
are not identical matrices. They may have different sparsity patterns.

The shared implementation is in root-level `js_*` MATLAB helpers; mathematical
interface functions and individual historical Anymatrix generators are at
the repository root. This directory layout is compatible with Anymatrix's
remote-group installation, which clones the repository inside its own
`unitary/private` directory.

## Additional real and complex OPUC measures

Five finite-matrix families share the canonical measure-first pipeline:

```matlab
% Real orthogonal CMV from a wrapped Gaussian.
R = rogers_szego(40,0.6,'cmv',-1);

% Generalized circular Jacobi weight with complex coefficients.
U = circular_jacobi(40,0.3,0.7,'cmv',exp(1i*0.2));

% Constant and periodic Verblunsky parameters (possibly arc/gap spectra).
G = geronimus_cmv(40,0.25+0.15i);
P = periodic_cmv(40,[0.2; 0.3i]);

% Rational Bernstein-Szego weight from arbitrary finite Schur data.
B = bernstein_szego_finite(40,[0.25; 0.2i; -0.1]);
```

All five provide corresponding `szego_*_measure` descriptors accepted by
`szego_verblunsky` and `szego_matrix`. Real interior coefficients with
real terminal phase give real orthogonal matrices; otherwise the
finite matrices are complex unitary. See the encyclopedia for support,
weights, and literature provenance. The measure at the distinguished
vector is not generally the equal-weight eigenvalue-counting measure.


## Citation-led rational and butterfly additions (November 2026)

Three new families complement polynomial OPUC, lossless filter-bank and
real-orthogonal benchmarks.

```matlab
% Rational OPUC finite unitary model; exterior poles can include Inf.
R = rational_opuc(24,[2*ones(22,1);Inf], ...
    .2*ones(23,1),exp(1i*.3));
S = rational_opuc(24,[],[],[],'spectral');

% Small rational-inner matrix evaluated at a frequency on the circle.
B = blaschke_potapov(12,[.2+.1i;-.3+.2i],[],.5);
Bop = blaschke_potapov(12,[],[],[],[],'operator');

% N-by-N real orthogonal butterfly with O(N log N) implicit action.
Q = orthogonal_butterfly(64,2*pi*rand(63,1));
Qop = orthogonal_butterfly(8192,[],'operator');
```

`rational_opuc` uses an operator Möbius transform of a finite unitary
Hessenberg recurrence matrix. The exterior poles and recurrence parameters
jointly define the rational construction; its finite trace law is generally
unknown. `blaschke_potapov` returns a fixed-order matrix that is unitary
**pointwise in frequency**, not a finite bank acting across a signal's
frequency samples. `orthogonal_butterfly` uses a separate angle at each
internal node of a recursive binary tree; random angles are not, in
general, Haar random orthogonal matrices. The latter two support
forward/adjoint matrix-free actions. See the encyclopedia for proofs,
conventions and primary references.

## Hadamard, multiplicative perturbation, and ETNA QR benchmarks

```matlab
% Diţă-deformed complex Hadamard; the modulus of every entry is 1/sqrt(N).
H = complex_hadamard(24,4);
Hf = complex_hadamard(8192,128,[],'operator');

% Rank-one unitary perturbation of an upper Hessenberg core model.
U = unitary_rank_one_update(64,pi/3);
Uf = unitary_rank_one_update(8192,pi/3,[],[],[],'operator');

% Types I, II (real), III (four near-decoupled blocks) in ETNA (2015).
E = etna_hessenberg(40,'I',2026);
R = etna_hessenberg(40,'II',2026);
C = etna_hessenberg(40,'III',2026,0.001);
Co = etna_hessenberg(8192,'III',2026,1e-3,'operator');
```

Diţă's family uses a block Fourier tensor product with independent
diagonal unit phases; its generic eigenphase distribution is not fixed.
The rank-one update is multiplicative, **not** a generic additive
unitary-plus-rank-one matrix. Its default choice preserves the upper
Hessenberg shape, and a zero-interior-coefficient case has exact finite
eigenphases. The ETNA Type III generator is a reproducible controlled
four-block near-degeneracy variant of Example 7.4, **not** an exact
recreation of a specific published random sample. All three expose
terminal `spectral` mode and matrix-free `operator` mode. See
[the encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md) for references and
precise constructions.


## Quantum cat, Krawtchouk and periodic D4 wavelet matrices

~~~matlab
C = quantum_cat_map(37,2,2);                % hyperbolic torus map
Cf = quantum_cat_map(4096,2,2,'operator');  % fast chirp-FFT-chirp
K = krawtchouk_transform(33,.3);           % symmetric orthogonal
Sk = krawtchouk_transform(33,.3,'spectral');% exact +/-1 eigenvalues
D = daubechies_wavelet(64);                % periodic D4/db2
Df = daubechies_wavelet(4096,[],'operator');% linear-time filter bank
~~~

The cat-map generator uses an even-diagonal classical hyperbolic
automorphism [a,1;a*d-1,d], with positive even a,d and a*d>=4.
This avoids parity issues with the common classical [2,1;1,1]
matrix. Its global metaplectic phase is specified, and no
generic analytic eigenphase law is claimed.

The Krawtchouk transform is built from ordered, signed eigenvectors
of a binomial Jacobi matrix. The Jacobi eigenvalues (grid nodes
0:N-1) are NOT eigenvalues of the transform, which is a symmetric
orthogonal involution. The construction uses a dense eigensolver
rather than a fast polynomial recurrence.

The Daubechies D4/db2 transform explicitly uses periodic boundary
closure, including the length-2 final stage. It differs from Haar
and provides O(n) forward/adjoint operations.


## Hahn and dual Hahn transforms, and spherical geodesic transport

~~~matlab
H = hahn_transform(20, .4, .9);       % polynomial Hahn transform
D = dual_hahn_transform(20, .4, .9);  % exactly H' (dual quadratic grid)

R = sphere_geodesic_transport(7);   % SO(7) shortest-geodesic lift
Fr = sphere_geodesic_transport(10000,[],[],'operator');
Sr = sphere_geodesic_transport(7,[],[],'spectral'); % exact phase pair
~~~

The Hahn family is defined using the self-adjoint Hahn difference
operator. Its **Jacobi** eigenvalues are quadratic in polynomial
degree, but its **orthogonal transform** eigenvalues are generally
different and not available in closed form. The dual Hahn transform
uses the matching transpose basis: the dual polynomial grid is
quadratic, and matrix orthogonality is preserved exactly in the
normalization conventions. Both generators use a dense symmetric
eigensolver, not a fast polynomial transform.

The spherical family records the ambient SO(n) rotation whose
restriction to tangent spaces is **Levi-Civita parallel transport**
along the unique shortest geodesic from x to y on the unit sphere.
Because source and destination tangent spaces differ, the ambient
orthogonal extension—not an arbitrarily chosen tangent-frame matrix—
is the well-defined finite matrix in the database. The spectrum is
the exactly known pair exp(±i theta), where theta=acos(x'*y),
and n-2 copies of 1. Antipodal points are excluded.

## Install as an Anymatrix external group

Install the upstream Anymatrix package first and add its root folder to the
MATLAB path. Then run:

```matlab
anymatrix('groups','unitary','jaurentz/anymatrix-unitary');
U = anymatrix('unitary/jacobi_szego',16,0,0,'cmv',1);
anymatrix('test','unitary');
```

Anymatrix requires MATLAB R2020b or later. Its GitHub remote update path has
historically run `git pull` from the current working directory instead of
the group worktree. If an update fails, update the repository in
`<anymatrix-root>/unitary/private/` using `git -C <path> pull --ff-only`,
then run `anymatrix('scan')`.

## Example gallery

```matlab
% Exact spectrum, but dense circulant eigenvectors.
U = phase_circulant(24,linspace(-2,2,24));

% Sparse random circular-beta matrix; beta=2 has the CUE eigenvalue law.
C = circular_beta(64,2,123);

% Matrix of order 2*n with a known gapped dispersion relation.
D = dirac_automaton(40,0.45);

% Orthogonal finite transform from real signal processing.
W = haar_wavelet(32);

% Floquet spin, self-dual chaos and topological BdG examples.
T = kicked_top(9,3.2,pi/2);
I = floquet_ising(4,0.65,0.43);
S = selfdual_kicked_ising(4,[0.1,-0.2,0.3,0.4]);
K = driven_kitaev(6,struct('mu1',0,'mu2',2.5));

% Quantum rotor and exactly unitary two-site circuit constructions.
R = kicked_rotor(48,5);
D = dual_unitary_circuit(6,2,.34,[.37,-.61,.29,.83]);
X = xxz_floquet(6,.6,pi/3);

% Kicked bosonic dimer and cavity-QED/quantum-network matrices.
BH = kicked_bose_hubbard(12,1,.4,1,1,.17);
BHop = kicked_bose_hubbard(100,1,.4,1,1,.17,'operator');
JC = jaynes_cummings(12,.8,.35,1,1);
CC = chalker_coddington(4,3,pi/4,2*pi*rand(48,1));

% Nonuniform Szego probability measure and a CMV realization.
mu = szego_bernstein_measure(0.3+0.2i);
B = szego_matrix(mu,24,'cmv',1);
```

### Spectral-band examples

The [OPUC and spectral-gap quickstart](examples/quickstart.m) combines
measure-first Jacobi--Szego construction, a **Dirac two-band model**, and a
matrix-free periodic CMV action. The
[unit-circle spectral-band examples](examples/unit_circle_spectral_bands.m) cover:

- Exactly prescribed finite spectral arcs (via `phase_circulant`).
- The Dirac automaton's two gapped bands, exact finite eigenphases and
  analytic limiting **equal-weight trace density**.
- A real period-two CMV discriminant and its **infinite-operator** band
  endpoints, plus the Geronimus constant-coefficient essential arc.
- An unwrapped tight-binding propagator with a known arcsine phase law
  and a periodic split-step walk with two exact quasienergy branches.

From the repository root, run:

```matlab
addpath(pwd);
run('examples/quickstart.m');
plotSpectralBands = true;    % use false for batch runs
run('examples/unit_circle_spectral_bands.m');
```

For periodic CMV and Geronimus, the theoretical bands describe the
**infinite operator**; finite eigenvalues depend on the terminal phase
and can occur in gaps. These examples never equate OPUC quadrature
weights with the empirical eigenvalue-counting measure. The showcase
plots finite eigenphases separately from theoretical limiting bands
or densities.

## Signal-processing matrix-free examples

The eight signal-processing generators below return numeric matrices by
default or a function handle with a final `'operator'` argument, following the
[matrix-free convention](docs/MATRIX_FREE.md).

```matlab
C = cazac_circulant(13,5);            % unitary Chu circulant
D = orthogonal_dct(16,2);            % real orthogonal DCT-II
H = discrete_hartley(16);           % real orthogonal involution
F = fresnel_transform(16);          % chirp circulant
A = affine_fourier(16,.01,.03);      % chirp-DFT-chirp
Z = finite_zak(16,4);               % Zak polyphase transform
W = time_frequency_shift(16,3,5);   % Weyl-Heisenberg monomial
T = periodic_hilbert(16,1,-1);      % completed Hilbert transform

applyF = fresnel_transform(4096,'operator');
X = randn(4096,3)+1i*randn(4096,3);
Y = applyF(X,'notransp');
X2 = applyF(Y,'transp');            % conjugate-transpose action
```

The conventional discrete Hilbert filter with zero DC/Nyquist
multipliers is not unitary on the entire space. Our finite completed
version assigns unimodular signs. Signal-frequency variables must
not be confused with eigenphases of the transform matrix.
### Fourier, sine, Walsh, Slepian, and graph transforms

The normalized DFT, orthonormal DST I–IV, natural/sequency-ordered
Walsh–Hadamard, DPSS/Slepian, and graph Fourier provide complementary
signal-processing transform families.

```matlab
S = discrete_sine(24,4);           % real orthogonal DST-IV
W = walsh_hadamard(32,'sequency'); % order by zero-crossing count
P = dpss_transform(48,.18);       % Slepian concentration eigenbasis
G = graph_fourier(16);            % path-graph Laplacian eigenvectors
```

DFT, DST, and Walsh support efficient `'operator'` multiplication. **DPSS
and general graph Fourier do not**: they return explicit orthogonal matrices,
because caching dense eigenvectors in a function handle is not a meaningful
matrix-free implementation. Laplacian or concentration eigenvalues are
different from the eigenvalues of the resulting transform matrix.


### Nonlinear unitary dynamics (Lohe matrix oscillators)

The `lohe_matrix_flow` family returns a **unitary matrix on a nonlinear
coupled trajectory**, rather than an eigenbasis or a fixed unitary
propagator. The numerical method evolves all oscillators together using
synchronous second-order exponential midpoint steps; every step preserves
unitarity up to floating-point roundoff. Use a smaller step size for
convergence studies. The terminal `'ensemble'` mode returns all evolved
oscillator matrices in a 3-D array.

```matlab
U = lohe_matrix_flow(4,3,1,1);      % order 4; three oscillators
All = lohe_matrix_flow(4,3,1,1,[],[],1,32,'ensemble');
err = norm(U'*U-eye(4),'fro');       % near machine precision
W = lohe_matrix_flow(4,3,1,1,[],[],2,64); % oscillator 2, 64 steps
S = lohe_matrix_flow(4,[],[],[],[],[],[],[],'spectral');
% S.finite.kind = 'unknown': no generic analytical eigenphase law.
```

Unlike the ordinary matrix exponential with a fixed Hermitian generator,
these systems are coupled and nonlinear; they are not intrinsically a
structured unitary eigenvalue problem. For source models and discretization
see [the matrix encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).

### Riemannian geometry and parallel-transport matrices

The geometry collection reuses the existing sphere-transport generator and adds four distinct constructions. It exposes `sphere_geodesic_transport` (shortest-arc rank-two
rotation), `grassmann_geodesic` (ambient unitary lift with exact angular
spectrum), `geodesic_holonomy` (closed spherical geodesic polygon),
`connection_cycle` (sparse discrete-connection block permutation), and
`lie_geodesic_transport` (real orthogonal parallel transport on u(n)).

```matlab
A = sphere_geodesic_transport(1000);               % rank-two sphere transport rotation
S = grassmann_geodesic(60, [.3;.8], [], 'spectral');
H = geodesic_holonomy(5,eye(5,3));      % exact triangle eigenphase +/-pi/2
F = connection_cycle(4096, [], 'operator');
L = lie_geodesic_transport(10, [], 1); % order 100, sparse, known phases
```

The holonomy of a spherical triangle has an exact area-driven trace law;
more general geodesic polygons are not automatically labeled exact.
The connection cycle is a proposed benchmark based on documented discrete
connection angles, not a built-in mesh geometry solver. See
[the encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md) for definitions, provenance,
and finite-spectrum qualifications.

## Automatic tests

- **MATLAB unit tests** under [`tests/`](tests/) check normalizations,
  Geronimus identities, coefficient formulas, all twisted/core forms,
  spectral and cyclic-vector moment invariants, structural sparsity,
  exact reference spectra, random reproducibility, and input validation.
- The Anymatrix `test_run.m` entry point invokes the **same suite**.
- [GitHub Actions](.github/workflows/matlab-tests.yml) runs MATLAB on each
  push and pull request. It also clones upstream Anymatrix and tests the
  installed local Git checkout as an external group.
- `tests/test_independent_stieltjes.m` independently reconstructs the Jacobi
  recurrence from Gauss--Jacobi nodes and weights using MATLAB only.

The public [unitary-spectral-density](https://github.com/jaurentz/unitary-spectral-density)
repository additionally tests this matrix library with Fourier moments,
weighted Gauss--Szego quadrature, and Poisson-resolvent algorithms. Its CI runs
against both a pinned submodule and the latest Anymatrix main commit every day.
When upstream MATLAB CI passes on `main`, this repository's workflow can
immediately dispatch the tested SHA downstream. This requires the one-time
upstream Actions secret `SPECTRAL_DENSITY_DISPATCH_TOKEN`, a fine-grained token
restricted to the downstream repository with **Contents: Read and write**.
Until that secret is configured, the dispatch job skips with a notice and
the downstream daily CI remains active. See the
[development guide](https://github.com/jaurentz/unitary-spectral-density/blob/main/docs/DEVELOPMENT.md).

**CI status** is displayed in the repository's *Actions* tab. The status of
any unexecuted workflow should not be interpreted as a passing test.

## Scope and mathematical roadmap

The OPUC branch is **measure → Verblunsky coefficients → finite closure →
representation**, leaving room for future isometric/rational Arnoldi and
non-core representations. Historical CMV/Hessenberg interfaces share their
numerical core implementation. The encyclopedia is the authoritative
human-readable record of included families, reference spectra and citations.
Future candidates are labeled as unimplemented, rather than silently treated
as available matrices.

## License

BSD 2-Clause; see [LICENSE](LICENSE). Original references are cited in
individual generators and [docs/MATHEMATICS.md](docs/MATHEMATICS.md).
