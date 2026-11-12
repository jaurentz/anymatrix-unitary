# Anymatrix Unitary

A MATLAB collection of finite **unitary and orthogonal test matrices**, designed as
an external group for [Anymatrix](https://github.com/north-numerical-computing/anymatrix).
The repository is organized around reproducible mathematical constructions,
known spectral information, and executable regression tests. It starts with
eight generators; it is not yet an exhaustive catalog.

**Status:** initial development release. MATLAB is the sole implementation
and testing environment. GitHub Actions executes the MATLAB tests after
pushes and on pull requests; consult the Actions tab for the latest result.

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

| ID in Anymatrix | Generator | Structural/spectral point |
|---|---|---|
| `unitary/cycle_shift` | `cycle_shift(n)` | Sparse permutation; exact roots of unity |
| `unitary/diagonal_phase` | `diagonal_phase(n,theta)` | Sparse prescribed eigenphases |
| `unitary/cycle_propagator` | `cycle_propagator(n,t)` | Circulant propagator; analytically known phases |
| `unitary/finite_cmv` | `finite_cmv(n,alpha)` | Sparse finite CMV matrix |
| `unitary/jacobi_szego` | `jacobi_szego(n,a,b,rep,beta,twists)` | Jacobi/Szegő family; several representations |
| `unitary/opuc_matrix` | `opuc_matrix(n,alpha,beta,rep,twists)` | General finite unitary core product |
| `unitary/split_step_walk` | `split_step_walk(n,theta1,theta2)` | Sparse real orthogonal quantum walk |
| `unitary/haar_unitary` | `haar_unitary(n,seed)` | Haar unitary, reproducible seeded draws |

**Do not interpret `n` uniformly across physical models**: `split_step_walk(n)`
returns a `2*n`-by-`2*n` matrix. The other current generators return `n`-by-`n`.

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

The shared implementation is in root-level `js_*` helper files; mathematical
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

## Automatic tests

- **MATLAB unit tests** under [`tests/`](tests/) check normalizations,
  Geronimus identities, coefficient formulas, all twisted/core forms,
  spectral and cyclic-vector moment invariants, structural sparsity,
  exact reference spectra, random reproducibility, and input validation.
- The Anymatrix `test_run.m` entry point invokes the **same suite**.
- [GitHub Actions](.github/workflows/matlab-tests.yml) runs MATLAB on each
  push and pull request. It also clones upstream Anymatrix and tests the
  installed local Git checkout as an external group.
- `tests/test_independent_stieltjes.m` constructs a Gauss--Jacobi quadrature
  rule using a symmetric tridiagonal Jacobi matrix, independently recovers
  its recurrence by the discrete Stieltjes procedure, and then checks the
  resulting inverse Geronimus coefficients against the explicit formulas.
  No Python, SciPy, or NumPy runtime is required.

**CI status** is displayed in the repository's *Actions* tab. The status of
any unexecuted workflow should not be interpreted as a passing test.

## Scope and mathematical roadmap

The OPUC branch is **measure -> Verblunsky sequence -> finite closure ->
representation**, leaving room for future isometric/rational Arnoldi and
non-core representations. Other families will use shared construction
layers where appropriate while retaining historical public names. The
research survey and bibliography remain in project source documents;
additional families are not treated as implemented until tested.

## License

BSD 2-Clause; see [LICENSE](LICENSE). Original references are cited in
individual generators and [docs/MATHEMATICS.md](docs/MATHEMATICS.md).
