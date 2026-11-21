# Intrinsic parallel transport on normal-metric real flag manifolds

This note documents `real_flag_parallel_transport.m` (issue #55).
It describes a **real orthogonal map between tangent spaces**, not an
ambient frame rotation and not a complex-unitary flag matrix.

## Flag geometry and normalization

Choose integers `k_1,...,k_r>=1`, with `r>=2` and sum `n`.
The real (unoriented) partial flag manifold is

```math
F=SO(n)/S(O(k_1)\times\cdots\times O(k_r)).
```

Each flag is represented by a complete tuple of mutually orthogonal
symmetric projectors `P_a=G P_a(0) G^T`, where each `P_a(0)` selects
the `a`th coordinate block and `G\in SO(n)`. The quotient normal
metric is induced by

```math
\langle X,Y\rangle=\tfrac12\mathrm{tr}(X^TY)
\quad\text{on }\mathfrak{so}(n).
```

The Lie algebra of the isotropy subgroup consists of within-block
skew matrices. Its orthogonal complement `m` consists of skew
matrices with **zero within-block diagonal blocks**, of dimension

```math
d=\sum_{a<b}k_a k_b=\tfrac12(n^2-\sum_a k_a^2).
```

For every cross-block pair `i<j`, define
`F_ij=e_i e_j^T-e_j e_i^T`. These form an orthonormal basis of
`m` for the metric above. Given any horizontal skew `K\in m`,
the normal-metric geodesic is `G(t)=\exp(tK)`; its flag
projectors are `P_a(t)=G(t)P_a(0)G(t)^T`.

## Explicit finite orthogonal parallel-transport matrix

Write the moving orthonormal tangent basis as
`G(t) F_ij G(t)^T` at the group level, or equivalently the
tuples `[G(t)[F_ij,P_a(0)]G(t)^T]_a` at the projector level.
The Levi–Civita connection in this moving basis is the
horizontal projection of the group bracket:

```math
S_{(i,j),(a,b)}
=-\tfrac12\big([K,F_{ab}]\big)_{ij}.
```

This matrix is real skew-symmetric; therefore

```math
Q(t)=\exp(tS),\qquad Q(t)^TQ(t)=I_d
```

is an *intrinsic* orthogonal parallel-transport map between the
initial and final flag tangent fibers. This follows from the
naturally reductive normal homogeneous connection; no dense polar
projection is involved.

The finite eigenphases of `Q(t)` are **unknown in general**.
Do not substitute the eigenvalues of `K`, the Gram matrix of
the projectors, or arbitrary pointwise principal angles.

## Independent projector connection test

The tuple-of-projectors realization has the induced inner product

```math
\langle X,Y\rangle_P
=\tfrac14\sum_a\mathrm{tr}
\big([X,P_a]^T[Y,P_a]\big)
=\tfrac12\mathrm{tr}(X^TY)
\quad (X,Y\in\mathfrak m).
```

It gives a second, independently implemented connection formula.
For basis elements `F_i,F_j`, with
`E_{j,a}=[F_j,P_a]` and
`E'_{j,a}=[K,E_{j,a}]`:

```math
A_{ij}=\tfrac14\sum_a\mathrm{tr}
\big([F_i,P_a]^T [K,[F_j,P_a]]\big).
```

Direct algebra yields `A=-S`. Thus the independently
computed parallel transport `\exp(-tA)` agrees with
`Q=\exp(tS)` in the same moving coordinates. The tests
verify this connection, the projector tangent metric, normal
metric isometry, reverse time, and parallel geodesic velocity.

## Special cases

- **Two blocks:** This is the real Grassmann symmetric space.
  `[m,m]` lies in the vertical Lie algebra, so `S=0` and
  `Q(t)=I_d` in the moving geodesic frame; the geometric
  transport is not identity in fixed ambient coordinates.
- **Complete flag (all block sizes one):** The isotropy Lie
  algebra is zero, and `Q` reduces to the bi-invariant
  `SO(n)` connection. In the moving frame the transported
  skew matrix is `E(t)=exp(-tK/2)E(0)exp(tK/2)`.
- **Three or more blocks:** The projected bracket is
  generally nonzero; genuine nontrivial orthogonal
  tangent-space matrices arise.

## Implementation scope and complexity

The MATLAB prototype forms a dense `d\times d` projected
connection and evaluates `expm(tS)`. Its assembly is
`O(d²)` for each fixed explicit `K`, its dense exponential
requires `O(d³)`, and its storage/action per RHS is
`O(d²)`. Hard limits are `3<=n<=32,d<=256`.

This is **not** the low-rank fast matrix exponential action of
Nguyen–Sommer, nor a general non-normal invariant flag metric, a
complex flag, a discrete flag random walk, or a new eigenphase law.
It shares the normal homogeneous Lie quotient philosophy with
`stiefel_parallel_transport` but has a different isotropy subgroup,
tangent dimension, and off-block horizontal basis.

## References

- [168] Nguyen and Sommer, *Parallel Transport on Matrix Manifolds
  and Exponential Action*, SIAM J. Matrix Anal. Appl. (2025).
  https://doi.org/10.1137/24M1684189
- For normal homogeneous flag connection formulas,
  *Effective formulas for the geometry of normal homogeneous spaces.
  Application to flag manifolds*, https://arxiv.org/abs/2302.14810
- Issue [#55](https://github.com/jaurentz/anymatrix-unitary/issues/55).
