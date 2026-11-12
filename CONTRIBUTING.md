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

## Spectral claims and source conventions

For every claim of a known spectrum, state whether it is an exact finite
multiset, a limiting density, or a random ensemble law. Preserve phase and
boundary conventions, coefficient indexing, and physical dimension. For an
OPUC family, place explicit weights or coefficient formulas in the common
measure/Verblunsky layer and reuse the multiplication representation code.

Documentation must be self-contained and use GitHub-supported `$...$` for
inline and `$$...$$` for display mathematics; do not refer to unpublished
working notes or external private project files.
