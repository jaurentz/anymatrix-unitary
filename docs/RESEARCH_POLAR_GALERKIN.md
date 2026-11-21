# Shen–Tran–Leok (2022): source-accurate polar Galerkin stage kernel

**Status:** implemented and tested the stage/endpoint **kinematic**
equations only. The full variational polar-decomposition
**symplectic** method remains open in
[issue #98](https://github.com/jaurentz/anymatrix-unitary/issues/98).
This code is deliberately an internal helper, NOT a new Anymatrix
matrix family or a claimed symplectic integrator.

## Source and exact equations implemented

X. Shen, K. Tran and M. Leok, *High-order symplectic Lie group methods
on SO(n) using the polar decomposition*, Journal of Computational
Dynamics **9** (2022), 529–551,
[DOI 10.3934/jcd.2022003](https://doi.org/10.3934/jcd.2022003).

The paper constructs a constrained Galerkin *discrete action* and
obtains the **full** Hamiltonian VPD equations (24a)–(24f).
Its (24b)–(24c) are the following kinematic subsystem:

~~~math
U_i = \mathbb P\left(g_0+h\sum_{j=1}^{s}a_{ij}U_j\Omega_j\right),
\qquad
g_1 = \mathbb P\left(g_0+h\sum_{i=1}^{s}b_iU_i\Omega_i\right).
~~~

Here $\mathbb P$ is the orthogonal polar factor of a nonsingular
matrix in $GL_+(n)$. The code in
**um_polar_galerkin_stages.m** implements precisely these equations,
and **um_polar_galerkin_coefficients.m** provides the source
Table 1 one-, two-, and three-stage Gauss–Legendre coefficients.
The stage skew matrices $\Omega_i$ are *supplied*, NOT calculated
from the missing Hamiltonian derivative (24f). Thus the
$2s$-order label for the Gauss–Legendre tableaux is **not**
a claim about the full dynamics of this partial implementation.

## Numerical contract

Input $g_0$ must already lie in $SO(n)$ and all $\Omega_i$ must
be skew-symmetric. Stage polar factors are evaluated by SVD
and their arguments must remain well-conditioned and in
$GL_+(n)$; negative-determinant or singular inputs cause errors,
not a silent orientation repair.

The fixed-point iteration updates all stages simultaneously and
checks the **true source residual** after convergence. Default
tolerance is $10^{-12}$; a failed nonlinear solve raises an
error. Endpoint and all stages are returned in $SO(n)$ to
floating-point accuracy. Tests independently verify Eqs. (24b,c),
the Gauss Table 1 identities, frame covariance under left
$SO(n)$ multiplication, and local error on a constant-angular-
velocity reference flow $g(t)=g_0\exp(t\Omega)$.

The method uses dense SVD per stage, not fast matrix-free polar
operations. No exact finite eigenphase law or symplectic map
is inferred from its rotation outputs.

## Crucial missing physics and invariant tests

The complete variational polar-decomposition integrator also
solves for internal momenta $\mu_i$, final cotangent momentum,
the multiplier $\Lambda$, the adjoints of the polar derivative
and associated Lyapunov equations in (24a), (24d), (24e),
(24f), (17) and (19). **No momentum is updated here** and
no phase-space symplecticity, Lie–Poisson preservation or exact
energy conservation is claimed.

The paper's SO(3) dipole-on-a-stick numerical experiments
use the Hamiltonian and initial momentum in Section 5.
Their published order and energy-error figures depend
on solving **the coupled full system**, and cannot be
reproduced by the kinematics-only helper.

## Remaining #98

Implement the coupled VPD nonlinear solve; supply Hamiltonian,
inertia, derivatives and momentum input; test full
symplectic/reduced Lie–Poisson invariants; reproduce the
actual paper's n=3 benchmarks; extend to n>=4 and register
the resulting SO(n) configuration increment as a new
unitary matrix family. Keep track of energy **error** rather
than claiming exact energy preservation.

**References:** [Published source](https://doi.org/10.3934/jcd.2022003),
[arXiv 2201.10768](https://arxiv.org/abs/2201.10768), §§3–5 and Table 1.
