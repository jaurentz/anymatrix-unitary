# Finite unitarity audit: Hankel, RBF, Zernike, and spherical Gauss–Laguerre transforms

**Research status:** Issue [#29](https://github.com/jaurentz/anymatrix-unitary/issues/29), 21 November 2026. This document identifies finite matrices that are *provably* square orthogonal/unitary, separates those matrices from fast approximations and rectangular evaluation operators, and specifies what would constitute a nonredundant Anymatrix generator. It does **not** claim implementation of a new fast transform.

## Acceptance principle

A square complex matrix `Q` qualifies if its **finite definition** yields `Q^*Q=I`, not merely if its continuous basis is orthogonal or if a numerical algorithm is approximately inverse to another. For a matrix `A` that is not unitary, `polar(A)` is a *derived* matrix with potentially different fast structure and spectrum. Store operator setup, action, adjoint and spectral reference separately. Eigenvalues of quadrature nodes, differential operators or kernel matrices must not be reported as eigenphases of `Q`.

## 1. Zernike radial blocks: a genuine exact finite orthogonal family

Fix nonnegative integer azimuthal order `m` and radial polynomial degree `p`. Set

```math
x=2r^2-1,\qquad R_{m+2k}^{m}(r)=r^m P_k^{(0,m)}(2r^2-1),\quad k=0,\ldots,p-1.
```

The radial inner product (before any angular normalization) obeys

```math
\int_0^1 r^{2m}P_j^{(0,m)}(2r^2-1)P_k^{(0,m)}(2r^2-1)r\,dr
=2^{-(m+2)}\int_{-1}^{1}(1+x)^m P_j^{(0,m)}(x)P_k^{(0,m)}(x)\,dx.
```

Take the `p`-node Gauss–Jacobi rule `(x_j,w_j)` for weight `(1+x)^m` and the orthonormal polynomials `p_k^{(0,m)}` for that same weight. The **square** matrix

```math
Q_{jk}=\sqrt{w_j}\,p_k^{(0,m)}(x_j),\qquad 1\le j\le p,\quad 0\le k<p
```

satisfies `Q^TQ=I_p` exactly in exact arithmetic: each product `p_kp_l` has degree at most `2p-2`, and Gaussian quadrature is exact through `2p-1`. The factors `r_j^m` and `r_jdr_j` cancel after the corresponding physical radial weighting is incorporated; using bare values `R_n^m(r_j)` is **not** an orthogonal matrix.

**Repository consequence:** The matrix is exactly the already implemented `gauss_jacobi_transform(p,0,m)` mathematical family (subject to checking the generator's argument order). A Zernike wrapper would supply meaningful disk-coordinate metadata but would *not* introduce a new algebraically distinct matrix. Its available fast `operator` mode is not automatically inherited for general integer `m`: issue [#39](https://github.com/jaurentz/anymatrix-unitary/issues/39) remains relevant.

A full, mixed-`m` Zernike transform on an arbitrary disk sample grid is usually rectangular and is **not** made unitary by this fixed-`m` argument. A claim of a square full-disk transform needs an explicitly specified set of angular nodes, radial nodes, sampling weights, truncation/bandlimit, and an exact discrete orthogonality proof.

## 2. Spherical Gauss–Laguerre radial blocks: another exact Gaussian-quadrature specialization

For fixed spherical harmonic degree `ell\ge0`, put `x=r^2`. In the natural radial measure `r^2dr`, radial factors of the form `r^ell L_k^{ell+1/2}(r^2)e^{-r^2/2}` have

```math
r^{2\ell}e^{-r^2}r^2dr=\tfrac12 x^{\ell+1/2}e^{-x}dx.
```

Consequently, if `x_j,w_j` are `p`-point generalized Gauss–Laguerre nodes and weights for `x^{ell+1/2}e^{-x}`, and `l_k^{(ell+1/2)}` are orthonormal polynomials for that weight, then

```math
Q_{jk}=\sqrt{w_j}\,l_k^{(ell+1/2)}(x_j)
```

is a **square exact orthogonal** matrix in exact arithmetic, again by Gaussian quadrature exactness. This is a radial coefficient-to-normalized-sample transform, not automatically the complete three-dimensional spherical Gauss–Laguerre evaluation map. The latter also contains angular harmonics and often uses oversampling/non-square grids.

**Repository consequence:** This radial block is a parameter specialization of `gauss_laguerre_transform`, with `alpha=ell+1/2`. It is not a distinct new matrix class unless the application supplies a reproducible genuinely different finite coupling or a fast transform interface; stable large-order application belongs with [#40](https://github.com/jaurentz/anymatrix-unitary/issues/40).

## 3. RBF and kernel matrices: avoid a false unitary classification

For distinct sites and a strictly positive-definite real Gaussian radial kernel, the finite Gram matrix `K` is real symmetric positive definite. Thus its **orthogonal polar factor is exactly `I`**, not a nontrivial RBF-derived rotation. The common interpolation/evaluation matrix and its numerical inverse are not orthogonal transforms. A mathematically valid alternative is `U(t)=exp(itK)`, which is unitary for real `t`; however, `K` is generally dense, `U` need not inherit a fast kernel matvec, and a truncated polynomial/Krylov approximation need not be exactly unitary. This is an application-specific instance of the Hermitian-exponential construction, not by itself a new fast-generator result. Repeated eigenvalues create a noncanonical eigenbasis, so simply returning an arbitrary eigensolver eigenbasis is not a reproducible parameterized generator without an orientation convention.

## 4. Discrete Hankel and Bessel transforms: research unresolved

The following distinctions must be checked **from the published finite formulas**, including endpoint handling and weights:

1. Bessel-zero (Lemoine) discrete Hankel formulations and quasi-discrete Hankel formulations use orthogonality relations that can be *asymptotic* or approximate at a prescribed finite order. A claimed exact matrix needs an identity for the actual finite weighted sampling matrix, not merely its continuous Bessel-function limit.
2. A fast Bessel/Hankel transform may approximate the **same** finite mathematical matrix, or an entirely different interpolation/expansion operator. Its cost and approximation error must be stated independently from the exact matrix definition.
3. Bessel functions of half-integer order have trigonometric reductions, but a resulting DST/DCT specialization is algebraically redundant unless it offers distinct physical parameter or boundary data.
4. Do not run `polar` silently to force finite unitarity: the dense polar factor generally forfeits the original fast action and gives a different matrix.

**No general exact finite orthogonal Hankel/Bessel generator is certified by this note.** Before implementing one, reproduce its normalization and derive `Q^*Q=I` or quantify `Q^*Q-I` as a finite-order approximation.

## 5. Proposed regression and metadata checklist

- For each fixed Zernike `m` or SGL `ell`, compare a separately weighted physical evaluation of the radial functions with the existing Gaussian-polynomial `Q`; test multiple orders and parameter extremes.
- Validate `norm(Q'*Q-eye(p))`, adjoint round trips, and multiple right-hand sides. The exact theoretical identity is separate from double-precision error; no claim of exact floating-point orthogonality.
- For any full disk/3D transform, verify its actual dimension and exact discrete Gram identity; reject naive square cropping of oversampled grids.
- Label spectrum of the **finite normalized evaluation matrix** as `unknown` unless independently derived. Do not substitute radial differential-operator eigenvalues or quadrature nodes.
- Report explicit construction cost, any operator setup cost, per-RHS apply cost, storage, and errors. An `O(p log p)` interpolation/connection algorithm is not automatically an `O(p log p)` operator for this `Q`.
- Add wrapper aliases only if there is a real application benefit, with references to the existing Jacobi/Laguerre families and to issues #39/#40. Do not duplicate a mathematically identical generator to inflate the catalog.

## Source pointers

- J. Townsend, *A fast algorithm for the discrete Hankel transform*, SIAM J. Sci. Comput., https://doi.org/10.1137/151003106.
- A. E. Siegman, *Quasi fast Hankel transform*, Opt. Lett. 23 (1998), 409–411, https://doi.org/10.1364/OL.23.000409. (Check bibliographic attribution and exact finite formula against the publisher before adding to the reference registry.)
- M. Slevinsky, **FastTransforms** project, https://mikaelslevinsky.github.io/FastTransforms/.
- *Fast Fourier Transforms for Spherical Gauss-Laguerre Basis Functions*, https://arxiv.org/abs/1604.05140.
- NIST Digital Library of Mathematical Functions, Jacobi and Laguerre orthogonality and Gaussian quadrature: https://dlmf.nist.gov/18.
