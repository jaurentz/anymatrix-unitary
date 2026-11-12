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

\[
 a_{k+1}^2=(1-\alpha_{2k-1})(1-\alpha_{2k}^2)
           (1+\alpha_{2k+1}),
\]
\[
 b_{k+1}=(1-\alpha_{2k-1})\alpha_{2k}
  -(1+\alpha_{2k-1})\alpha_{2k-2}.
\]

The inverse endpoint Riccati algorithm sets `u_0=2-b_1`, `v_0=2+b_1` and
for `k>=1`

\[
 u_k=2-b_{k+1}-a_k^2/u_{k-1},\quad
 v_k=2+b_{k+1}-a_k^2/v_{k-1}.
\]

Then

\[
 \alpha_{2k}=(v_k-u_k)/(v_k+u_k),\qquad
 \alpha_{2k-1}=1-(u_k+v_k)/2 \quad(k\ge1).
\]

These use the `x=2*cos(theta)` normalization, so factors of 2 must not be
silently imported from formulas written for Jacobi measures on `[-1,1]`.
The initial positive Jacobi off-diagonal is

\[
 a_1^2=\frac{16(a+1)(b+1)}{(a+b+2)^2(a+b+3)}.
\]

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
