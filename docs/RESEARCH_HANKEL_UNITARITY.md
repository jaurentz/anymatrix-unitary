# Finite discrete Hankel transforms: an exact-unitarity exclusion test

**Status (21 November 2026).** Research note for [issue #29](https://github.com/jaurentz/anymatrix-unitary/issues/29). This note does **not** register a new unitary generator. In particular, near-orthogonal finite Bessel-zero discretizations must not acquire the Anymatrix `unitary` or `orthogonal` properties without a separate exact finite proof.

## Finite kernel and its actual Gram defect

For integer Bessel order $\nu\ge 0$, let $j_{\nu,k}$ be the *k*-th strictly positive zero of $J_\nu$, and write $c=j_{\nu,n+1}$. Consider the real symmetric $n\times n$ Bessel-zero transform

```math
T_{jk}=\frac{2J_\nu(j_{\nu,j}j_{\nu,k}/c)}
 {c\,|J_{\nu+1}(j_{\nu,j})J_{\nu+1}(j_{\nu,k})|},
\qquad 1\le j,k\le n.
```

This differs from the signed-denominator convention in Baddour--Chouinard only by diagonal sign similarity and therefore has the **same Gram singular values and Frobenius unitarity defect**. It is the symmetric normalized quasi-DHT (not GSL's unnormalized sampling kernel). The scaling of the original forward and inverse *continuous Hankel approximations* must not be conflated with Euclidean orthogonality of the finite sample matrix.

Crucially, in **Baddour, *The Discrete Hankel Transform*, §4.4**, the underlying discrete Bessel orthogonality identity is explicitly said to hold **exactly in the infinite-order limit**; finite orders are approximations, with a finite numerical error. Later statements calling the resulting finite matrix orthogonal are thus not an exact finite-dimensional algebraic proof.

Independent double-precision calculations using SciPy's Bessel zeros and Bessel values give $\|T^T T-I\|_F$:

| $\nu$ | $n=3$ | $n=5$ | $n=8$ |
|:--:|--:|--:|--:|
| 0 | 9.1437849317e-6 | 3.6355997265e-6 | 1.3165558111e-6 |
| 1 | 2.1386000003e-5 | 9.0140614707e-6 | 3.4418923694e-6 |
| 2 | 8.5653186341e-5 | 3.7946946823e-5 | 1.5153930963e-5 |

These discrepancies are **orders of magnitude above machine precision**. The executable MATLAB regression in `tests/test_hankel_nonunitarity.m` independently brackets Bessel zeros with `fzero` and checks symmetry, nonzero Gram defects and near-involution for the specified convention. The numerical table is evidence against strict finite orthogonality, not a general exact formula for the errors.

### Half-integer exception and duplication warning

For $\nu=1/2$, $J_{1/2}(x)=\sqrt{2/(\pi x)}\sin x$ and $j_{1/2,k}=k\pi$. The same absolute-denominator normalization reduces **exactly** to

```math
T_{jk}=\sqrt{\frac{2}{n+1}}\sin\frac{jk\pi}{n+1},
```
the orthonormal DST-I matrix. This is a genuine finite orthogonal transform but **not a new generic Bessel family**, and should be deduplicated against the existing DST / Chebyshev-II coverage. Other special-order collapses need individual verification.

### Orthogonalization is a different family, not a repair to the source

For nonsingular finite $T$, the polar factor $Q=T(T^T T)^{-1/2}$ is mathematically orthogonal. When $T$ is real symmetric, $Q$ is a symmetric involution with eigenvalues $\pm1$, assuming no zero eigenvalue. But this is **a derived polar-corrected transform**, not the original finite Bessel quadrature kernel. A generic dense polar computation destroys any promised fast-asymptotic apply complexity, so the current library should not silently use it as a replacement. The general `polar_unitary` constructor already covers the abstract polar operation; duplicating it here would add no independently justified fast Hankel algorithm.

### Distinguish other fast-transforms questions

- **Townsend (2015):** Fast approximation to an order-zero discrete Hankel transform using asymptotics and FFTs. Its computational complexity is not evidence that the chosen finite Bessel-zero matrix is exactly unitary; this source should not be labeled as a fast application of an exact orthogonal Bessel transform without an explicit identification and proof.
- **Guizar-Sicairos--Gutiérrez-Vega (1998):** Quasi-discrete Hankel transfer used in optical propagation, requiring its own precise kernel and Gram check.
- **RBF:** A real positive-definite Gaussian Gram matrix is *not* orthogonal in general. The separate database `periodic_rbf_propagator` is exactly unitary because it exponentiates an FFT-diagonalized real symmetric *periodic* kernel.
- **Zernike:** `fourier_zernike_disk` is exactly unitary as a **modal-radius hybrid**, not as physical samples on a shared disk grid.
- **Spherical Gauss--Laguerre:** General spherical sampled harmonic transforms are often rectangular. A proposed square angular/radial representation requires a genuine discrete orthogonality proof, not merely a quadrature construction that is an isometry into an oversampled grid.

## Remaining acceptance work in #29

Research a truly distinct, exactly finite orthogonal Bessel construction, if one exists; or mark the generic Bessel-zero and quasi-DHT kernels as exclusions. Reconcile Townsend's fast algorithm with the *identical* target finite matrix and with explicit forward/adjoint errors, and investigate SGL square-unitary cases and any nontrivial additional RBF subclasses. An accepted generator must have source-matched formulas, tests, Anymatrix properties, catalog metadata, eigenphase descriptors and honest setup/multiplication complexity.

## References

1. N. B. Baddour and U. Chouinard, *Theory and operational rules for the discrete Hankel transform*, J. Opt. Soc. Am. A **32** (2015), 611--622. https://doi.org/10.1364/JOSAA.32.000611
2. N. B. Baddour, *The Discrete Hankel Transform* (2019), especially §4.4 (finite/infinite Bessel orthogonality). https://www.intechopen.com/chapters/65719
3. A. Townsend, *A Fast Analysis-Based Discrete Hankel Transform Using Asymptotic Expansions*, SIAM J. Numer. Anal. **53** (2015), 1897--1917. https://doi.org/10.1137/151003106
4. M. Guizar-Sicairos and J. C. Gutiérrez-Vega, *Computation of quasi-discrete Hankel transforms of integer order for propagating optical wave fields*, Opt. Lett. **23** (1998), 409--411. https://doi.org/10.1364/OL.23.000409
5. D. Lemoine, *A note on orthogonal discrete Bessel representations*, J. Chem. Phys. **118** (2003), 6697--6699. https://doi.org/10.1063/1.1559677
