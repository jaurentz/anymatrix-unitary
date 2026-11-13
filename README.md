# Anymatrix Unitary

A MATLAB collection of finite **unitary and orthogonal test matrices**, designed as
an external group for [Anymatrix](https://github.com/north-numerical-computing/anymatrix).
The repository is organized around reproducible mathematical constructions,
known spectral information, and executable regression tests. It contains
50 generators; additional benchmark families are under development.

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

The collection currently has **50 public generators**. Each has a complete
entry with formulas, spectral information, applications, parameter conventions,
and references in the [matrix encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).
The [mathematical notes](docs/MATHEMATICS.md) explain the shared OPUC engine.

| Category | Generators |
|---|---|
| Prescribed spectra and circulants | `cycle_shift`, `weighted_shift`, `diagonal_phase`, `phase_circulant`, `cycle_propagator`, `orthogonal_blocks` |
| OPUC / Szegő / core matrices | `finite_cmv`, `opuc_matrix`, `jacobi_szego`, `bernstein_szego`, `circular_beta`, `circular_jacobi`, `rogers_szego`, `geronimus_cmv`, `periodic_cmv`, `bernstein_szego_finite` |
| Random and signal-processing | `haar_unitary`, `haar_orthogonal`, `haar_wavelet` |
| Signal-processing transforms | `orthogonal_dct`, `discrete_hartley`, `cazac_circulant`, `fresnel_transform`, `affine_fourier`, `finite_zak`, `time_frequency_shift`, `periodic_hilbert`, `discrete_fourier`, `discrete_sine`, `walsh_hadamard`, `dpss_transform`, `graph_fourier` |
| Quantum maps, walks, and algorithms | `split_step_walk`, `dirac_automaton`, `grover_operator`, `quantum_baker`, `kicked_harper`, `kicked_top`, `floquet_ising`, `selfdual_kicked_ising`, `driven_kitaev`, `kicked_rotor`, `dual_unitary_circuit`, `xxz_floquet`, `kicked_bose_hubbard`, `jaynes_cummings`, `chalker_coddington` |
| Matrix functions and dilations | `cayley_jacobi`, `polar_unitary`, `julia_unitary` |

Generators return **actual finite square matrices by default**. Selected
generators also accept `'operator'` as their final argument and return a
function handle supporting forward and conjugate-transpose actions, multiple
right-hand sides, and dimension queries. See [matrix-free documentation](docs/MATRIX_FREE.md).
`split_step_walk`, `dirac_automaton`, `driven_kitaev`, and `julia_unitary` return order `2*n`
for an input physical size `n`. `floquet_ising` and `selfdual_kicked_ising` return order `2^n`, capped at 11 spins for dense assembly. `dual_unitary_circuit` and `xxz_floquet` also return order `2^n`, capped at 10 spins. `driven_kitaev` is a 2n-dimensional BdG operator, not the full Fock-space propagator. `kicked_bose_hubbard(n)` has order `n+1` for n bosons; `jaynes_cummings(n)` has order `2*(n+1)` for a photon cutoff n; `chalker_coddington(nx,ny)` has order `4*nx*ny`. `quantum_baker` requires even order;
`haar_wavelet` requires power-of-two order. MATLAB R2020b or newer is required.

**Spectral reference data:** exact eigenvalues are provided or described for
several families; other families have ensemble laws, limiting densities, or
no universal law. The encyclopedia distinguishes these cases. Do not infer a
limiting eigenphase distribution from unitarity alone.

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
