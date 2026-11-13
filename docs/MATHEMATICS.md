# Mathematical and verification notes

## Independent paths from Jacobi measure to Verblunsky coefficients

The prototype deliberately retains **two distinct paths**:

1. Jacobi `(a,b)` to the explicit Jacobi–Szegő formula for `alpha_k`.
2. Jacobi `(a,b)` to its orthonormal three-term Jacobi coefficients on
   `[-2,2]`, then **inverse Geronimus** endpoint Riccati recurrences to `alpha_k`.

It also computes the **forward Geronimus** mapping back to Jacobi recurrence
coefficients. This avoids making a self-consistency test from two copies of the
same formula.

For an orthonormal real-line polynomial basis `p_k` on `[-2,2]`, write

`x*p_k = a_(k+1)*p_(k+1) + b_(k+1)*p_k + a_k*p_(k-1)`.

The direct Geronimus formulas, with `alpha_(-1)=-1` and zero contribution from
any `alpha_(-2)` term in the `k=0` identity, are

```math
 a_{k+1}^2=(1-\alpha_{2k-1})(1-\alpha_{2k}^2)
           (1+\alpha_{2k+1}),
```
```math
 b_{k+1}=(1-\alpha_{2k-1})\alpha_{2k}
  -(1+\alpha_{2k-1})\alpha_{2k-2}.
```

The inverse endpoint Riccati algorithm sets `u_0=2-b_1`, `v_0=2+b_1` and
for `k>=1`

```math
 u_k=2-b_{k+1}-a_k^2/u_{k-1},\quad
 v_k=2+b_{k+1}-a_k^2/v_{k-1}.
```

Then

```math
 \alpha_{2k}=(v_k-u_k)/(v_k+u_k),\qquad
 \alpha_{2k-1}=1-(u_k+v_k)/2 \quad(k\ge1).
```

These use the `x=2*cos(theta)` normalization, so factors of 2 must not be
silently imported from formulas written for Jacobi measures on `[-1,1]`.
The initial positive Jacobi off-diagonal is

```math
 a_1^2=\frac{16(a+1)(b+1)}{(a+b+2)^2(a+b+3)}.
```

The code treats this term separately to avoid a removable singularity when
`a+b=-1`. Both recurrence paths are useful for regression but the inverse
Riccati version may lose relative accuracy for extreme parameters/large order.

## Discrete quadrature versus trace spectral measure

For a finite `n` and closure beta, compute normalized eigenpairs
`U v_j = z_j v_j`, and define `w_j=abs(v_j(1))^2`. The weighted measure
`sum_j w_j delta_(z_j)` matches the original circle measure's moments for
`z^k`, `|k|<=n-1`, but is not identical to the original continuous measure.
In contrast, `n^(-1) sum_j delta_(z_j)` is the empirical trace spectral measure,
which is relevant to spectral density and need not have the same moments.

Theorem 4.8 of Bultheel–Cantero–Cruz-Barroso (2015) shows that finite snake
representations from the same Schur coefficients and terminal parameter share
both eigenvalues and the first component of normalized eigenvectors.

## Source references

- B. Simon, *OPUC on One Foot*, Bull. Amer. Math. Soc. 42 (2005), 431–460,
  esp. Section 11, https://web.ma.utexas.edu/mp_arc/c/05/05-78.pdf
- NIST DLMF, *Polynomials Orthogonal on the Unit Circle*, §18.33,
  https://dlmf.nist.gov/18.33
- A. Bultheel, M. J. Cantero, R. Cruz-Barroso,
  *Matrix methods for quadrature formulas on the unit circle: A survey*,
  J. Comput. Appl. Math. 284 (2015), 78–100, especially Section 4.4
  and Theorem 4.8. DOI: https://doi.org/10.1016/j.cam.2014.11.002
- M. J. Cantero, L. Moral, L. Velázquez,
  *Five-diagonal matrices and zeros of orthogonal polynomials on the unit circle*,
  Linear Algebra Appl. 362 (2003), 29–56.
  DOI: https://doi.org/10.1016/S0024-3795(02)00457-3
- R. Killip, I. Nenciu, *Matrix models for circular ensembles*,
  Int. Math. Res. Not. (2004), 2665–2701, §6 for Jacobi recurrences.
  DOI: https://doi.org/10.1155/S1073792804141597
- N. J. Higham, M. Mikaitis, *Anymatrix: An Extensible MATLAB Matrix Collection*,
  Numer. Algorithms 90 (2022), 1175–1196.
  DOI: https://doi.org/10.1007/s11075-021-01226-2

## Canonical finite multiplication-matrix pipeline

For a measure on the unit circle, the design principle is

```math
\text{probability measure}\;\longrightarrow\;
(\alpha_0,\alpha_1,\ldots)\;\longrightarrow\;
(\alpha_0,\ldots,\alpha_{n-2};\beta)\;\longrightarrow\;U_n.
```

The interior coefficients satisfy $|\alpha_j|<1$; the boundary closure
$\beta\in\mathbb T$ is independent of the infinite measure. With the
specified terminal parameter, different admissible snake-shaped core
orderings give isospectral finite matrices and the same $e_1$ spectral
moments, even though their sparsity patterns are different. The module
`js_core_product.m` owns this representation logic; the measure and
coefficient modules never construct a second copy of the core formulas.

The entry points `szego_jacobi_measure`, `szego_bernstein_measure` and
`szego_verblunsky` define the measure and coefficient layers.
`szego_matrix(mu,n,representation,beta)` closes and assembles the finite
matrix. The historically named matrix generators are thin adapters.

## Normalized Jacobi–Szegő measure

For $a,b>-1$, let the real-line probability measure be

```math
 d\nu(t)=\frac{(1-t)^a(1+t)^b}{2^{a+b+1}B(a+1,b+1)}\,dt,
 \qquad -1<t<1.
```

Under $t=\cos\theta$, its conjugation-symmetric circle lift has density
(with respect to $d\theta$)

```math
 \frac{d\mu}{d\theta}=
 \frac{(1-\cos\theta)^{a+1/2}(1+\cos\theta)^{b+1/2}}
 {2^{a+b+2}B(a+1,b+1)}.
```

The explicit Verblunsky coefficients are

```math
 \alpha_{2k}=\frac{b-a}{2k+a+b+2},\qquad
 \alpha_{2k+1}=-\frac{a+b+1}{2k+a+b+3}.
```

The Jacobi real-line recurrence uses the rescaled variable
$x=2t=z+z^{-1}$. This accounts for the factor-of-two differences between
Jacobi polynomial formulas stated on $[-1,1]$ and the direct Geronimus
relations above.

## One-parameter Bernstein–Szegő measure

For a complex parameter $a\in\mathbb D$, the Poisson/Szegő probability
measure in our convention has density

```math
 d\mu_a(e^{i\theta})=\frac{1-|a|^2}
 {2\pi|e^{i\theta}-\overline{a}|^2}\,d\theta.
```

Its moments and Verblunsky sequence are known explicitly:

```math
 \int z^k\,d\mu_a(z)=\overline{a}^{\,k}\quad(k\geq0),
 \qquad \alpha_0=a,\quad \alpha_j=0\;(j\geq1).
```

This is the first family added beyond Jacobi to the measure-first dispatcher,
`js_verblunsky`. Its finite CMV, Hessenberg, inverse Hessenberg and twisted
realizations are all constructed through the same shared core function.

## Random circular-$\beta$ measures

The Killip–Nenciu distribution generates a *finite random* discrete spectral
measure, rather than starting from one predetermined continuous weight.
For matrix order $n$ and Dyson parameter $\beta>0$, independent interior
coefficients have uniform argument and

```math
 |\alpha_j|^2\sim\mathop{\mathrm{Beta}}\!\left(
 1,\frac{\beta(n-j-1)}2\right),\qquad j=0,\ldots,n-2.
```

The terminal $\beta_{\mathrm{terminal}}$ is uniform on $\mathbb T$ (not
the positive Dyson parameter). A sparse CMV realization then has precisely
the circular-$\beta$ eigenangle law in finite dimension. The CUE eigenvalue
law occurs when $\beta=2$, but the sparse CMV matrix is not itself Haar
unitary in its entries. The implementation samples radial squares by an
inverse CDF using MATLAB's core `rand` function.

## Orthogonal versus complex unitary representation

For real Verblunsky coefficients and closure $\beta\in\{\pm1\}$, the
finite CMV, Hessenberg, inverse-Hessenberg and twisted matrices are real
orthogonal. Complex Verblunsky coefficients or terminal phases generally
produce complex unitary matrices. The generic interfaces therefore declare
`unitary` but not `orthogonal` as an unconditional Anymatrix property.

Additional representations motivated by isometric or rational Arnoldi
require further basis/pole parameters and cannot always be described by
simply permuting the same adjacent cores. They will use the shared measure
and coefficient layers with new representation backends.

## Real and complex OPUC expansions (November 2026)

The classification is by the *interior* Verblunsky coefficients, independently
of whether the underlying positive measure has a closed-form density.
A real coefficient sequence with terminal beta = +/-1 gives a real
orthogonal finite CMV matrix; nonreal interior coefficients or a complex
terminal phase generally give a complex unitary matrix. The terminal phase
is always independent of the infinite measure. The same OPUC data can be
realized as CMV, Hessenberg, inverse Hessenberg, or twisted core products.

### Generalized circular Jacobi

For a > -1/2 and eta real, define the positive measure on (0,2*pi) by

```math
d\mu(\theta)=Z^{-1} e^{(\pi-\theta)\eta}
 (4\sin^2(\theta/2))^a\,d\theta.
```

For b=a+i*eta the complex Verblunsky sequence is

```math
\alpha_k=-\frac{(b)_{k+1}}{(1+\overline{b})_{k+1}},
\qquad k\geq 0.
```

The code evaluates this formula by successive ratios to avoid large
factorials and normalizes the positive weight by numerical quadrature.
At eta=0, alpha_k=-a/(k+a+1) is real (ordinary circular Jacobi).
This is a deterministic circle weight, not a draw from the circular
Jacobi *random matrix ensemble*. See Simon (2005) and Bourgade,
Nikeghbali, Rouault (2009), doi:10.1093/imrn/rnp092.

### Rogers--Szego (wrapped Gaussian)

For 0<q<1, let s=-log(q) and wrap N(0,s) onto the circle. Then

```math
\alpha_k=(-1)^k q^{(k+1)/2},\qquad k\geq 0.
```

The density is evaluated as a wrapped normal sum. Its exact moments
are integral exp(i*k*theta) dmu = q^(k^2/2), for integer k.
Reference: Simon (2005), Example 1.6.5; DLMF 18.33.

### Constant and periodic Verblunsky coefficients

The constant (Geronimus) measure uses alpha_k=a, |a|<1, for all k.
The spectral measure may have a proper absolutely continuous arc and
an additional atom. The periodic family takes a finite coefficient
pattern in the open disk and repeats it indefinitely; the associated
two-sided/infinite CMV theory gives spectral bands and gaps.
These statements concern the infinite measures. The code produces
finite *unitary* n-by-n closures, not principal submatrices of
infinite CMV operators. Reference: Simon (2005), Parts 1--2.

### Finite Bernstein--Szego data

For alpha_0,...,alpha_(m-1) in the disk and alpha_k=0 for k>=m,
write Phi_m^* for the reversed monic Szego polynomial. Then

```math
\frac{d\mu}{d\theta} =
\frac{\prod_{j=0}^{m-1}(1-|\alpha_j|^2)}
 {2\pi|\Phi_m^*(e^{i\theta})|^2}.
```

The code evaluates the reversed polynomial through the Szego recurrence,
so it accepts arbitrary finite complex coefficients and includes the
existing one-coefficient Poisson/Bernstein--Szego measure as a special
case. A zero-length coefficient list denotes the uniform measure.

All five descriptors are implemented by public SZEGO_*_MEASURE functions;
`szego_verblunsky` dispatches their coefficient formulas, and
`szego_matrix` performs finite unitary closure.

**Spectral interpretation.** The finite e1 spectral measure is a positive
Gauss--Szego quadrature measure; its nodes with *equal* weights are not,
in general, distributed according to the original OPUC weight. An
absolutely continuous circle weight need not give a nonuniform limiting
unweighted eigenvalue density.
