# Exact geodesic-flow Koopman blocks on the round sphere (issue #61)

**Disposition:** documented exact construction; no new algebraically distinct family recommended.

## Geometry and invariant measure

Let `S²` carry the round metric of radius one and let `T¹S²` be its unit tangent bundle. For `x in S²` and `v in T_x S²` of unit length, define the positively oriented frame

```math
G(x,v)=[\,x\ \ v\ \ x\times v\,]\in SO(3).
```

The unit-speed great-circle geodesic flow is

```math
x(t)=\cos(t)x+\sin(t)v,\qquad
v(t)=-\sin(t)x+\cos(t)v .
```

Consequently, the frame follows **right multiplication**:

```math
G(t)=G(0)R_z(t),\qquad
R_z(t)=\begin{bmatrix}
\cos t&-\sin t&0\\
\sin t&\cos t&0\\
0&0&1
\end{bmatrix}.
```

The Liouville measure on `T¹S²` corresponds, after normalization, to Haar measure on `SO(3)`. The Koopman pullback `(K_t f)(G)=f(GR_z(t))` is therefore unitary on `L²(SO(3))`. A *generic* finite projection of an infinite-dimensional Koopman operator need not be unitary; the subspaces below **are exactly invariant**.

## Finite Peter–Weyl blocks and eigenphases

Let `D^ell(G)` denote the standard integer-degree irreducible Wigner representation, dimension `2ell+1`, with basis indices `-ell,...,ell`. Matrix coefficients `D^ell_{mn}(G)`, with fixed left index `m` and varying right index `n`, span a `(2ell+1)`-dimensional invariant subspace:

```math
(K_t D^\ell_{mn})(G)
=D^\ell_{mn}(GR_z(t))
=\sum_{k=-\ell}^{\ell}D^\ell_{mk}(G)D^\ell_{kn}(R_z(t)).
```

Thus the finite Koopman matrix, with the **right-index coefficient columns**, is `D^ell(R_z(t))` (or its transpose if the basis/operator convention is transposed). Under the convention implemented by `wigner_rotation`,

```math
D^\ell(R_z(t))
=\operatorname{diag}(e^{-imt})_{m=-\ell}^{\ell}
=\texttt{wigner\_rotation(ell,t,0,0)}.
```

The eigenphases are **exactly** `-m t mod 2pi`. Their multiplicity is one in each fixed-left-index block, and `2ell+1` in the entire degree-`ell` Peter–Weyl matrix-coefficient space (dimension `(2ell+1)^2`). The latter representation is `I_(2ell+1) \otimes D^ell(R_z(t))`, up to basis ordering. It can be assembled from the existing generator without a new algorithm or parameter class.

For `ell=1`, the coordinate observables `G_{ab}` give the **real Cartesian representation** `R_z(t)`, unitarily equivalent to the degree-one Wigner matrix. One explicit normalized complex basis uses columns

```math
B=\begin{bmatrix}
1/\sqrt2&0&-1/\sqrt2\\
-i/\sqrt2&0&-i/\sqrt2\\
0&1&0
\end{bmatrix},
\qquad
B^*R_z(t)B
=\operatorname{diag}(e^{it},1,e^{-it}).
```

This directly checks the sign and basis conventions.

The geodesic flow on a sphere of radius `r` changes `t` to `t/r` for unit-speed geodesics.

## Why this does not warrant a new MATLAB matrix generator

- The degree blocks **already are** `wigner_rotation(ell,t,0,0)`; their direct sums and multiple copies use existing `spherical_harmonic_rotation` and standard Kronecker/block constructions. A new named constructor would be an application wrapper rather than a distinct matrix family.
- The flat-torus flow `x(t)=x(0)+tv` likewise gives diagonal phases and duplicates `diagonal_phase` at the matrix level.
- On generic negatively curved compact manifolds, the Koopman spectrum typically has continuous components and the needed finite invariant subspaces should **not** be inferred from these highly symmetric examples.
- The linearized Poincaré map of a closed geodesic is symplectic, not generally Euclidean orthogonal/unitary; it is a different object from the Koopman restriction.

## Reproducible regression

`tests/test_spherical_geodesic_koopman.m` verifies the explicitly defined frame flow, the real `ell=1` basis similarity, exact spectral phase/multiplicity formulas and the existing Wigner generator at various angles and orders. No approximation/polar correction is used.

## References

- Standard Peter–Weyl theorem and Wigner `D` representation on `SO(3)`; the repository's `wigner_rotation.m` and `spherical_harmonic_rotation.m` implement these exact irreducible blocks.
- Existing finite Koopman issue #61: https://github.com/jaurentz/anymatrix-unitary/issues/61
- A general warning on finite invariant Koopman subspaces: https://pmc.ncbi.nlm.nih.gov/articles/PMC4769143/
