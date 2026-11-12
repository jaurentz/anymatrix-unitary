# Anymatrix Unitary

A MATLAB collection of finite **unitary and orthogonal test matrices**, designed as
an external group for [Anymatrix](https://github.com/north-numerical-computing/anymatrix).
The repository is organized around reproducible mathematical constructions,
known spectral information, and executable regression tests. It starts with
22 generators; additional benchmark families are under development.

**Status:** development collection with MATLAB-only tests. Automated MATLAB
regression and Anymatrix integration tests run on GitHub Actions after every
push or pull request. See the [Actions tab](https://github.com/jaurentz/anymatrix-unitary/actions)
for the current branch status.

## Clone and test locally

Requires **MATLAB R2020b or later**, with no additional MATLAB toolboxes for
these generators or tests.

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

The collection currently has **22 public generators**. Each has a complete
entry with formulas, spectral information, applications, parameter conventions,
and references in the [matrix encyclopedia](docs/MATRIX_ENCYCLOPEDIA.md).
The [mathematical notes](docs/MATHEMATICS.md) explain the shared OPUC engine.

| Category | Generators |
|---|---|
| Prescribed spectra and circulants | `cycle_shift`, `weighted_shift`, `diagonal_phase`, `phase_circulant`, `cycle_propagator`, `orthogonal_blocks` |
| OPUC / Szegő / core matrices | `finite_cmv`, `opuc_matrix`, `jacobi_szego`, `bernstein_szego`, `circular_beta` |
| Random and signal-processing | `haar_unitary`, `haar_orthogonal`, `haar_wavelet` |
| Quantum maps, walks, and algorithms | `split_step_walk`, `dirac_automaton`, `grover_operator`, `quantum_baker`, `kicked_harper` |
| Matrix functions and dilations | `cayley_jacobi`, `polar_unitary`, `julia_unitary` |

Generators return **actual finite square matrices**, not matrix-free objects.
`split_step_walk`, `dirac_automaton`, and `julia_unitary` return order `2*n`
for an input physical size `n`. `quantum_baker` requires even order;
`haar_wavelet` requires power-of-two order. MATLAB R2020b or newer is required.

**Spectral reference data:** exact eigenvalues are provided or described for
several families; other families have ensemble laws, limiting densities, or
no universal law. The encyclopedia distinguishes these cases. Do not infer a
limiting eigenphase distribution from unitarity alone.

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
`x=2*t=z+z^(-1)` has the explicit Verblunsky sequence

```text
alpha_(2k)   =  (b-a)/(2k+a+b+2)
alpha_(2k+1) = -(a+b+1)/(2k+a+b+3).
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

% Nonuniform Szego probability measure and a CMV realization.
mu = szego_bernstein_measure(0.3+0.2i);
B = szego_matrix(mu,24,'cmv',1);
```

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
