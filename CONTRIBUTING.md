# Contributing

We welcome mathematically precise, reproducible unitary/orthogonal matrix
families. Before proposing a new public generator, check whether it is a
new mathematical family, a different representation of an existing family,
or a benchmark parameter preset. Prefer shared implementations for core
factorizations, recurrence coefficients, and matrix functions.

Each public generator must:

1. Return a square unitary/orthogonal MATLAB array (with explicit physical
   dimension conversions, where relevant).
2. Accept an argument-free call returning `[]` and a list of recognized
   Anymatrix properties; the property list should describe the family
   accurately, including parameter-dependent variants.
3. Document dimension, parameters, finite-boundary conventions, analytical
   spectral information, applications and primary references.
4. Have tests in `tests/` covering unitarity, dimension, edge cases and
   independent reference data wherever possible.
5. Add a mathematically complete entry in
   [`docs/MATRIX_ENCYCLOPEDIA.md`](docs/MATRIX_ENCYCLOPEDIA.md),
   with verifiable primary references and qualifications on spectral laws.
6. Add a machine-readable entry to `metadata/catalog.json` and update
   `Contents.m` and the README index.
7. Not introduce additional toolbox dependencies without documenting them.

Run `test_run` in MATLAB before submitting a pull request. GitHub Actions
must pass before merging. Random generators should support reproducible
seeds without altering the caller's MATLAB RNG state.

Do not conflate vector spectral measures (Szegő/Gauss--Szegő weights) with
the equal-weight trace spectral measure, and do not assume truncating an
infinite unitary operator preserves unitarity.

## Spectral claims, operator contracts, and source conventions

For every claim of a known spectrum, state whether it is an exact finite
multiset, a limiting density, or a random ensemble law. Preserve phase and
boundary conventions, coefficient indexing, and physical dimension. For an
OPUC family, place explicit weights or coefficient formulas in the common
measure/Verblunsky layer and reuse the multiplication representation code.

In complexity statements, let `N` denote the **actual square matrix order**;
specify any different physical site count, qubit count, polynomial degree,
local dimension, rank, bandwidth, depth, and number of right-hand sides.
Use standard NLA big-O bounds such as `O(N)`, `O(N log N)`, and
`O(N^2)`, stating setup costs and parameter assumptions separately.
Describe subquadratic algorithms in words rather than little-o notation.
A bound of `O(N^2)` alone does not establish a fast algorithm.

Distinguish these independent claims:
- The matrix is explicitly sparse, with its actual nonzero count.
- The matrix has a compact exact factorization or mathematical description.
- A fast forward **and adjoint** action is implemented and verified for the
  specified parameters and accuracy.
- The final output accepts the optional `'operator'` mode. This alone does
  **not** imply subquadratic complexity: a full optical mesh, for instance,
  retains `O(N^2)` couplers and work.
- A sparse input Hamiltonian or recurrence does not imply that its matrix
  exponential or eigenvector matrix is sparse or fast.
- An approximation/compression need not preserve exact unitarity.

Use `metadata/catalog.json` for registered generator names and
`matrix_free_supported`, `operator_kernel` and `operator_note` for
current operator coverage. The last two fields require interpretation:
`matrix_free_supported=true` may cover a restricted parameter subset
(e.g., the two Chebyshev cases of `gauss_jacobi_transform`).
The operator API and `'transp'`/adjoint convention are defined by
[the matrix-free documentation](docs/MATRIX_FREE.md); the distinctions among
structural classes and unimplemented fast algorithms are described by
[the structured-operators taxonomy](docs/STRUCTURED_OPERATORS.md).
Do not claim a new kernel is implemented without checking MATLAB source and tests.

Documentation must be self-contained and use GitHub-supported `$...$` for
inline mathematics and fenced `math` blocks for display mathematics; do
not refer to unpublished working notes or private project files.

GitHub's mathematics renderer restricts some otherwise valid LaTeX macros.
Do not use `\operatorname` or `\DeclareMathOperator` in rendered
mathematics. For named operators use, for example,
`\mathop{\mathrm{diag}}` instead of `\operatorname{diag}`.
Keep inline mathematics separate from immediately following letters: write
`$n$-th` instead of `$n$th`, and prefer `circular β` rather than
`circular-$\beta$` in headings and prose. Inspect both the GitHub preview
and the Markdown source.

Example of fenced display mathematics:

````markdown
```math
\alpha_0=0
```
````

Run `test_math_markup` and `test_catalog` as well as the relevant
matrix-free numerical tests before merging.
