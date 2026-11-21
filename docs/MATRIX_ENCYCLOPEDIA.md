# Encyclopedia of unitary and orthogonal test matrices

**Scope:** all generators currently implemented by this repository. Each entry
states the construction, the dimension convention, spectral information,
structure, an application or useful numerical pathology, and a literature
reference. The existence of a reference does not mean that its authors
introduced the particular default parameters chosen here.

The canonical machine-readable inventory is [`metadata/catalog.json`](../metadata/catalog.json).
The mathematics of the measure-first OPUC construction and its independent
Geronimus checks is in [`MATHEMATICS.md`](MATHEMATICS.md). This encyclopedia is
part of the public repository and is intended to stand alone.

## General conventions

For a unitary $U\in\mathbb C^{N\times N}$, eigenvalues are written
$z_j=e^{i\theta_j}$. Throughout this encyclopedia, **$N$ denotes the actual
matrix order**. A MATLAB generator parameter such as $n$ may instead describe
physical sites, spins, modes or polynomial coefficients; the relationship to
$N$ is specified for each family. The **empirical spectral measure** and its moments are

```math
\mu_U^{\mathrm{trace}}=\frac1N\sum_{j=1}^N\delta_{z_j},\qquad
m_k(U)=\frac1N\mathop{\mathrm{tr}}(U^k).
```

For a unit vector $v$, the different **vector spectral measure** is

```math
\mu_{U,v}=\sum_{j=1}^N |q_j^*v|^2\delta_{z_j},
```

where $q_j$ are normalized eigenvectors (counting multiplicities with
appropriate spectral projections). For finite OPUC quadrature matrices,
$v=e_1$ is the distinguished cyclic vector. **Gauss–Szegő weights are
vector weights, not equal-mass empirical weights.** A finite spectrum is
atomic; a continuous spectral density is a limiting or explicitly smoothed
object. Eigenphase distributions on the circle must account for wrapping
modulo $2\pi$.

**Anymatrix usage:** each matrix entry is available by its filename, such as
`jacobi_szego(32,0,0,'cmv',1)`, or as `anymatrix('unitary/jacobi_szego',...)`
after installing the external group. Optional inputs and defaults are defined
by the header of the linked MATLAB function. Some generators return order
$2n$ for a physical system of size $n$; such cases are explicitly marked.

## A. Prescribed spectra, circulants, and finite dynamics

### 1. Cyclic shift — [`cycle_shift.m`](../cycle_shift.m)

**Order:** $N=n$. **Structure:** real, orthogonal, sparse permutation,
circulant. **Spectrum:** precisely the $n$-th roots of unity.

```math
P_ne_j=e_{j+1}\quad(e_{n+1}=e_1),\qquad P_n^n=I,\qquad
\frac1n\mathop{\mathrm{tr}}(P_n^k)=\mathbf{1}_{n\mid k}.
```

Useful as an exact uniform-phase reference and to detect phase aliasing in
Fourier moment methods. It is a single-cycle permutation rather than a
generic permutation with multiple cycle lengths. **Reference:** [1].

### 2. Phase-weighted cyclic shift — [`weighted_shift.m`](../weighted_shift.m)

**Order:** $N=n$. **Structure:** sparse monomial unitary; becomes a real signed
permutation only for phases $0$ or $\pi$. For real weights $\phi_j$,

```math
Ue_j=e^{i\phi_j}e_{j+1}\quad(e_{n+1}=e_1),\qquad
U^n=e^{i\sum_j\phi_j}I.
```

Its spectrum comprises the roots of $z^n=e^{i\sum_j\phi_j}$. A useful
counterexample to the idea that changing the entries of a sparse unitary
necessarily changes its relative eigenvalue spacings: changing the phases
only rotates the perfectly spaced cycle spectrum. **Reference:** [1].

### 3. Prescribed diagonal phases — [`diagonal_phase.m`](../diagonal_phase.m)

**Order:** $N=n$. **Structure:** sparse diagonal unitary. For real phase inputs
$\theta_j$,

```math
D=\mathop{\mathrm{diag}}(e^{i\theta_1},\ldots,e^{i\theta_n}).
```

The spectrum, algebraic multiplicities, and trace moments are known exactly.
Use clustered phases, open gaps, isolated atoms, or an irregular deterministic
phase mesh to test spectral estimators independently of the eigensolver.
The default gives equally spaced roots of unity. **Reference:** [1].

### 4. Unitary circulant with prescribed phases — [`phase_circulant.m`](../phase_circulant.m)

**Order:** $N=n$. **Structure:** dense circulant with FFT-diagonalization.
For the normalized DFT $F$ and real phases $\theta_j$,

```math
U=F^*\mathop{\mathrm{diag}}(e^{i\theta_j})F.
```

It has **the same eigenvalues as the diagonal control**, but generally
nontrivial eigenvectors and dense entries. This pairs well with entrywise or
matvec algorithms, exposing accidental dependence on eigenvector localization.
Our implementation uses an inverse FFT of the Fourier eigenvalue samples.
**Reference:** [1].

### 5. Periodic tight-binding propagator — [`cycle_propagator.m`](../cycle_propagator.m)

**Order:** $N=n$. **Structure:** dense circulant. Write $H=P_n+P_n^*$; then

```math
U(t)=e^{-itH},\quad
z_k=e^{-2it\cos(2\pi k/n)},\quad k=0,\ldots,n-1.
```

The infinite-chain energy law is arcsine on $[-2,2]$, and its time-scaled
pushforward onto the circle produces a known nonuniform phase law. Large $|t|$
wraps the energy bands and can create overlapping singularities. This is a
strong test for circular spectral densities with analytic reference phases.
The $n=1,2$ cases retain the literal $H=P+P^*$ definition. **Reference:** [2].

### 6. Real orthogonal rotation blocks — [`orthogonal_blocks.m`](../orthogonal_blocks.m)

**Order:** $N=n$. **Structure:** sparse real block diagonal, orthogonal and
normal. Each $2\times2$ block is

```math
R(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix},\qquad
\sigma(R)=\{e^{i\theta},e^{-i\theta}\}.
```

An odd-dimensional matrix also has a chosen $+1$ or $-1$ terminal block.
Prescribe clusters, repeated conjugate pairs, determinant parity and eigenvalue
multiplicities. This is a controlled test for real-normal eigenvalue
algorithms and conjugation-symmetric phase density. **References:** [3, 4].

## B. OPUC, Szegő measures, and core-factorized representations

### Finite OPUC Christoffel measure change — [`opuc_christoffel.m`](../opuc_christoffel.m)

Given an n-atomic positive finite measure
`mu=sum_j w_j delta_(z_j)` on the unit circle and any finite
complex `zeta`, the positive Christoffel transformation is

```math
d\widetilde\mu(z)
=\frac{|z-\zeta|^2}{\sum_jw_j|z_j-\zeta|^2}\,d\mu(z).
```

The routine forms this new discrete measure and delegates inverse
unitary CMV/Hessenberg construction to `inverse_opuc`; it does
**not** implement a CMV Darboux matrix factorization. The two
operations should not be identified. `'operator'` reuses adjacent
Verblunsky core actions and their Hermitian adjoints.

If zeta is unequal to **every** original support point,
all transformed weights remain strictly positive. Therefore
the **finite trace spectrum** (the uniform measure on the
support) is unchanged, even though the distinguished-vector
measure, Verblunsky parameters and actual CMV matrix change.
In contrast, if zeta equals a support point exactly, that
atom is removed and the finite matrix's order decreases.
There is no claim of preserved finite eigenphases or order
under that support loss. If all mass is annihilated, the
transformation is invalid. Close or nearly coincident nodes
may make inverse Szego reconstruction ill conditioned and
are rejected when numerical breakdown is detected.

The additional metadata output records the old and new support
and weights, the positive normalization and removed indices;
the terminal `'spectral'` mode exposes the new trace and
cyclic measures separately. A finite atomic Christoffel
transformation does **not** automatically give a faithful
construction for a continuous infinite OPUC measure.
**References:** Simon, *OPUC* [5]; Cantero and collaborators,
“Darboux transformations for CMV matrices,”
[DOI](https://doi.org/10.1016/j.aim.2016.03.042).

### Unit-circle polynomial / companion-like unitary realization — [`circle_polynomial_unitary.m`](../circle_polynomial_unitary.m)

**Order:** $N=n$. Input is an admissible polynomial of degree n in
descending monomial coefficients, with n **simple** unit-circle roots,
and optionally positive cyclic spectral weights. For a polynomial,
the root ordering is not intrinsic; WEIGHTS are therefore paired
with roots sorted by increasing principal phase angle in [-pi,pi].
Such a polynomial
is necessarily **self-inversive**, but self-inversiveness alone
does *not* ensure circle-rootedness (for example, z²−3z+1).
The generator checks both coefficient reciprocity and numerical
root unimodularity, rejects repeated or poorly resolved roots,
and limits the numerical `roots`-based procedure to moderate n.

A generic monomial companion matrix is **not** generally unitary.
Instead this generator passes the computed circle-root nodes and
positive weights to `inverse_opuc` and obtains a **strictly
unitary** finite Hessenberg (default) or CMV realization.
The characteristic polynomial agrees, up to rounding in the
root computation, with the monic supplied polynomial; the
Hessenberg realization is built from adjacent Verblunsky
unitary cores rather than a nonunitary companion similarity.
This is a constructive example of a paraorthogonal-polynomial
unitary representation. The terminal unimodular closure is
computed from the root product, and all interior Verblunsky
coefficients have modulus less than one.

`'operator'` uses the existing adjacent-core fast actions and
Hermitian adjoint. `'pencil'` reuses the strict block-CMV
unitary generalized pencil, which shares the spectrum but
need not equal the Hessenberg matrix. `'parameters'` returns
the normalized input coefficients, computed nodes, Schur
parameters and boundary phase.

The **trace** eigenphase measure weights every root equally,
whereas the prescribed **cyclic** measure uses the independent
positive weights for e1. `'spectral'` provides these two
distinct finite atomic measures, but its root phases come
from a numerical `roots` calculation, not a symbolic
closed-form polynomial root formula. General monomial,
Jacobi-basis and rational-basis **nonunitary** colleagues or
companions are not automatically included in this class.

**Applications:** structure-preserving unitary eigenvalue
benchmarks and finite paraorthogonal eigenvalue problems.
**References:** Simon, *OPUC* [5], unit-circle quadrature [9],
and structured companion background
[DOI](https://doi.org/10.1016/0024-3795(94)90190-2).

### Inverse spectral CMV / Hessenberg — [`inverse_opuc.m`](../inverse_opuc.m)

**Order:** $N=n$. Input is a collection of n **distinct** unit-circle
points z_j and strictly positive numbers w_j, normalized so
sum(w_j)=1. This defines a finite **cyclic** spectral measure
mu_e1=sum_j w_j delta_(z_j). The inverse finite Szego procedure
computes the n-1 interior Verblunsky coefficients and a
unimodular terminal boundary beta, and the resulting standard
CMV or unitary Hessenberg core product is formed by the existing
`js_core_product`. An arbitrary monomial companion matrix
is **not** asserted to be unitary.

The inverse uses evaluations of the Szego monic polynomial
Phi_k and its reversed polynomial Phi_k^* on the support,
with inner product weighted by w_j:

```math
c_k=\frac{\sum_j w_j\,\overline{\Phi_k^*(z_j)}\,z_j\Phi_k(z_j)}
{\sum_j w_j|\Phi_k^*(z_j)|^2},\quad
\alpha_k=\overline{c_k},\quad
\Phi_{k+1}(z)=z\Phi_k(z)-c_k\Phi_k^*(z).
```

A common normalization at every step improves scaling. Since
the measure has n atoms with nonzero weights, `|alpha_k|<1`
for k<n-1 in exact arithmetic. The terminal coefficient has
exact phase `beta=(-1)^(n+1)*conj(prod_j z_j)`. Finite CMV
factors are unitary, independently of the reconstruction method.
`'parameters'` exposes the computed reflection coefficients.
Near-breakdown or unresolved coincident nodes are **rejected**
rather than silently labeled reliable; arbitrarily clustered
phases are not claimed to be numerically stable.

The finite *trace* measure is the uniform atomic measure
`(1/n)sum_j delta_(z_j)`, while the **distinguished-vector**
measure is `sum_j w_j delta_(z_j)`, matching the squared
first entries of normalized eigenvectors. `'spectral'`
exposes the two independent measures as `S.finite` and
`S.cyclic`, with their distinct Fourier moments.
`'operator'` reuses the existing OPUC adjacent-core
matrix-free actions, forward and Hermitian adjoint.

**Application:** inverse unitary eigenvalue problems, rational
quadrature on the circle and experimental counterexamples
for randomized trace-versus-vector spectral estimators.
References: Simon, *OPUC* [5], and unit-circle matrix
quadrature [9].

### 7. Finite CMV — [`finite_cmv.m`](../finite_cmv.m)

**Order:** $N=n$. **Structure:** sparse five-diagonal unitary. From
$\alpha_0,\ldots,\alpha_{n-2}$ in the open unit disk and a terminal
$\beta\in\mathbb T$, use the local core

```math
\Theta(\alpha)=
\begin{pmatrix}\overline{\alpha}&\sqrt{1-|\alpha|^2}\\
\sqrt{1-|\alpha|^2}&-\alpha\end{pmatrix}
```

and alternating even/odd block-diagonal factors $\mathcal{C}=\mathcal{L}\mathcal{M}$.
The finite matrix closes with the terminal boundary phase. Its eigenvalues
are the zeros of a paraorthogonal polynomial on the circle. Real interior
coefficients and a real closure produce an orthogonal CMV matrix.
**References:** [5–7].

### 8. General Verblunsky multiplication matrix — [`opuc_matrix.m`](../opuc_matrix.m)

**Order:** $N=n$. **Representations:** `hessenberg`, `inverse_hessenberg`,
`cmv`, or `twisted`. The shared function `js_core_product.m` assembles all
forms from the same interior coefficients and boundary phase. A twisted or
snake ordering is specified by a binary sequence of prepend/append decisions
for successive adjacent cores.

For this convention, the finite matrices have the same eigenvalues and the
same distinguished-$e_1$ spectral moments. **This is more than a statement
that arbitrary permutations of matrix factors preserve eigenvalues:**
the construction implements the special snake-shaped equivalence class.
Numerical challenges include clusters, nearly-reducible core sequences and
comparison of sparse and dense representations. **References:** [8–10].

### 9. Jacobi–Szegő / Geronimus lift — [`jacobi_szego.m`](../jacobi_szego.m)

**Order:** $N=n$. **Source measure:** the Jacobi probability measure
$(1-t)^a(1+t)^b\,dt$ on $[-1,1]$, $a,b>-1$, lifted by $t=\cos\theta$
into a conjugation-invariant circle measure. The explicit coefficients are

```math
\alpha_{2k}=\frac{b-a}{2k+a+b+2},\qquad
\alpha_{2k+1}=-\frac{a+b+1}{2k+a+b+3}.
```

Take the first $n-1$ and append a unit-modulus closing parameter. The
same measure can produce Hessenberg, inverse Hessenberg, CMV or twisted
matrices. Classical choices include Chebyshev I ($a=b=-1/2$), Legendre
($a=b=0$), Chebyshev II ($a=b=1/2$), and Gegenbauer
($a=b=\lambda-1/2$). The finite Gauss–Szegő nodes and weights are not
the original continuous measure. **References:** [5, 7, 11, 12].

### 10. Bernstein–Szegő / Poisson measure — [`bernstein_szego.m`](../bernstein_szego.m)

**Order:** $N=n$. **Source:** an explicit absolutely continuous probability
measure on the unit circle, with parameter $|a|<1$ and density (with respect
to $d\theta$)

```math
\frac{d\mu}{d\theta}=
\frac{1-|a|^2}{2\pi|e^{i\theta}-\overline{a}|^2}.
```

The coefficients in the stated OPUC convention are $\alpha_0=a$ and
$\alpha_k=0$ for $k\geq1$, giving the exact moment
$\int z^k\,d\mu(z)=\overline{a}^{\,k}$ for $k\geq0$.
This provides an analytically transparent **nonuniform** weight with
rapidly terminating Verblunsky data, a useful dual check of the measure
and multiplication matrix. All four existing representation choices work.
**References:** [5, 12].

### 11. Circular β random CMV — [`circular_beta.m`](../circular_beta.m)

**Order:** $N=n$. **Structure:** sparse five-diagonal random unitary.
The Killip–Nenciu construction chooses independent interior coefficients
with uniform angle and

```math
|\alpha_{j}|^2\sim\mathop{\mathrm{Beta}}\!\left(
1,\frac{\beta(n-j-1)}{2}\right),\qquad 0\leq j\leq n-2,
```

and an independent uniform terminal phase. The eigenangles have joint
density proportional to $\prod_{j<k}|e^{i\theta_j}-e^{i\theta_k}|^\beta$.
At $\beta=2$ the **eigenvalue law** is CUE, although the matrix is not
Haar-distributed on $U(n)$. This tests level repulsion while preserving
$O(N)$ sparse structure. The RNG-seeded option restores the caller's state.
No Statistics Toolbox is required. **Reference:** [13].

## C. Random and orthogonal signal-processing matrices

### 12. Haar unitary — [`haar_unitary.m`](../haar_unitary.m)

**Order:** $N=n$. **Structure:** dense, random. Phase-adjusted QR of a
complex Gaussian matrix gives Haar distribution on $U(n)$.
For the circular unitary ensemble (CUE), the **expected** one-point phase
density is uniform, but finite samples exhibit correlated phases with
repulsion; individual draws are not equally spaced. This is an ensemble
benchmark rather than a deterministic known-spectrum matrix. **Reference:** [14].

### 13. Haar orthogonal — [`haar_orthogonal.m`](../haar_orthogonal.m)

**Order:** $N=n$. **Structure:** dense, real orthogonal. Sign-adjusted QR of
an $N\times N$ real Gaussian matrix produces Haar distribution on $O(N)$.
Its eigenvalues obey complex-conjugate pairing and determinant/parity
constraints. **Haar on $O(n)$ is not Dyson's COE.** This distinction matters
when testing random matrix universality claims. **Reference:** [14].

### 14. Haar wavelet transform — [`haar_wavelet.m`](../haar_wavelet.m)

**Order:** $N=n=2^p$. **Structure:** sparse, real orthogonal, multiresolution.
The first analysis step forms adjacent sums/differences divided by
$\sqrt2$; subsequent levels act only on the coarse half. The `levels`
parameter ranges from zero (identity) to $p$ (full decomposition).
It is a signal-processing test of real orthogonality, highly nonuniform
sparsity and multiscale eigenvectors. **The eigenphases of a finite analysis
matrix are not the wavelet frequencies.** **Reference:** [15].

## D. Quantum maps, walks, and algorithms

### 15. Split-step quantum walk — [`split_step_walk.m`](../split_step_walk.m)

**Order:** $N=2n$ for $n$ sites. **Structure:** sparse, real orthogonal in the
selected rotation-coin convention. The periodic one-step operator is

```math
U=S_-C(\theta_2)S_+C(\theta_1).
```

It has a two-state internal coin and periodic conditional shifts; the
basis alternates coin states site by site. Parameters change the two
quasienergy branches and potential gaps. The current periodic uniform-coin
version has no physical edge, so **edge states require a separate open or
interface construction**. **Reference:** [16].

### 16. Periodic Dirac quantum cellular automaton — [`dirac_automaton.m`](../dirac_automaton.m)

**Order:** $N=2n$ for $n$ sites. **Structure:** sparse complex unitary.
With mass $m\in[0,1]$, $\nu=\sqrt{1-m^2}$ and periodic cyclic shift $P_n$,

```math
U=\begin{pmatrix}\nu P_n&-imI\\-imI&\nu P_n^*\end{pmatrix},\qquad
\cos\omega_k=\nu\cos(2\pi k/n).
```

The exact eigenvalues are $e^{\pm i\omega_k}$. For $0<m<1$, the large-$n$
phase density has gapped support and integrable band-edge singularities:

```math
\rho(\theta)=\begin{cases}
\frac{|\sin\theta|}{2\pi\sqrt{\nu^2-\cos^2\theta}},
  &|\cos\theta|<\nu,\\
0,&\text{otherwise}.
\end{cases}
```

Thus this is a particularly useful nonuniform, exact-density benchmark.
At $m=0$ it reduces to two oppositely directed shifts; at $m=1$ its
spectrum is concentrated at $\pm i$. **References:** [17, 18].

### 17. Grover search iterate — [`grover_operator.m`](../grover_operator.m)

**Order:** $N=n\geq2$. **Structure:** dense real orthogonal product of two
reflections. Mark the first $M$ basis states, $1\leq M<n$, and let
$|s\rangle=n^{-1/2}(1,\ldots,1)^T$. Then

```math
G=(2|s\rangle\langle s|-I)(I-2\Pi_M),\qquad
\sin^2\alpha=M/n.
```

The exact spectrum is $e^{\pm2i\alpha}$ (one each), $+1$ with multiplicity
$M-1$, and $-1$ with multiplicity $n-M-1$. Two nontrivial eigenphases
may have only $2/n$ empirical spectral mass, so a density estimator can
miss them. **Reference:** [19].

### 18. Quantum baker map — [`quantum_baker.m`](../quantum_baker.m)

**Order:** $N=n$ (with $n$ even). **Structure:** dense, product of Fourier transforms.
The Balazs–Voros quantization is

```math
B_n=F_n^*\mathop{\mathrm{diag}}(F_{n/2},F_{n/2}).
```

The optional Saraceno variant uses half-integer/antiperiodic Fourier phases
and a different parity-symmetry convention. Eigenphase statistics are
parameter/size/convention dependent; this matrix **does not inherit the
four-atom eigenvalue spectrum of $F_n$**. **References:** [20, 21].

### 19. Kicked Harper torus map — [`kicked_harper.m`](../kicked_harper.m)

**Order:** $N=n$. **Structure:** dense explicit matrix, FFT-factorized
mathematical action. With $D_x(K)$ and $D_p(L)$ defined by diagonal
phases $e^{-iK\cos(2\pi j/n)}$ and $e^{-iL\cos(2\pi j/n)}$,

```math
U=D_p(L)F_n^*D_x(K)F_n.
```

This fixed finite-torus convention is exactly unitary. Quasienergy
statistics, clustering and band structure depend on $K,L,n$ and the torus
quantization; statements about multifractality in irrational limits must
not be transferred indiscriminately to every finite test matrix.
**Reference:** [22].

## E. Matrix functions, polar factors, and dilations

### 20. Cayley image of a Jacobi matrix — [`cayley_jacobi.m`](../cayley_jacobi.m)

**Order:** $N=n$. **Structure:** dense complex unitary, derived from a real
symmetric tridiagonal Jacobi matrix $J_n$ with Jacobi weight parameters
$a,b>-1$. For scale $t>0$,

```math
U=(J_n-itI)(J_n+itI)^{-1},\qquad
\lambda_U=\frac{x-it}{x+it}\quad(x\in\sigma(J_n)).
```

This is an explicit Möbius pushforward of a classical real-line spectrum,
providing controlled phase clustering near $+1$ and $-1$ as the scale
changes. **It is not the same finite matrix as the Jacobi–Szegő lift.**
**Reference:** [23].

### 21. Nonsingular polar factor — [`polar_unitary.m`](../polar_unitary.m)

**Order:** $N=n$. **Structure:** generally dense. If $A=L\Sigma R^*$ is a
full SVD of a nonsingular matrix, the uniquely determined unitary polar
factor is $U=LR^*$. For real $A$, $U$ is orthogonal. It is also a
nearest unitary matrix to $A$ in every unitarily invariant norm.

We restrict the generator to numerically nonsingular inputs: for singular
$A$, an SVD-based square unitary extension is not unique and depends on
choices in nullspaces. Useful for polar algorithms, structure restoration,
Procrustes problems and perturbation analysis. **Reference:** [24].

### 22. Julia/Halmos unitary extension — [`julia_unitary.m`](../julia_unitary.m)

**Order and assumptions:** A square contraction produces a block unitary of
twice its order:

```math
T\in\mathbb{C}^{n\times n},\qquad
\lVert T\rVert_2\leq 1,\qquad
\mathcal{U}(T)\in\mathbb{C}^{2n\times 2n}.
```

**Structure:** generally dense, block unitary. Its defect operators are

```math
D_T=(I-T^*T)^{1/2},\qquad D_{T^*}=(I-TT^*)^{1/2}.
```

The Julia/Halmos extension is

```math
\mathcal{U}(T)=
\begin{pmatrix}
T & D_{T^*} \\
D_T & -T^*
\end{pmatrix}.
```

The input contraction occupies the leading diagonal block of the extension.
**This is not automatically a unitary power dilation:** taking the leading
block of higher powers of the extension does not, in general, reproduce the
corresponding powers of the input contraction. With the projection onto the
first block, this distinction can be written as

```math
P=\begin{pmatrix}I_n&0\end{pmatrix},\qquad
P\mathcal{U}(T)^kP^*\neq T^k
\quad\text{in general for }k\geq 2.
```

Applications include conservative system matrices, passive control, and
dilation theory. **Reference:** [25].

## F. Further real and complex OPUC matrix generators

### 23. Generalized circular Jacobi — [`circular_jacobi.m`](../circular_jacobi.m)

**Order:** $N=n$. **Structure:** unitary, sparse five-diagonal by default.
For real $a>-1/2$ and $\eta\in\mathbb R$, the positive circle weight is
proportional to $e^{(\pi-\theta)\eta}|1-e^{i\theta}|^{2a}$ for
$0<\theta<2\pi$. With $b=a+i\eta$, the Verblunsky sequence is

```math
\alpha_k=-\frac{(b)_{k+1}}{(1+\overline{b})_{k+1}}.
```

At $\eta=0$ all coefficients are real; otherwise they are generally
complex. Supports an explicitly prescribed nonuniform *vector* measure,
but does not promise the same limiting unweighted phase density.
The finite boundary $\beta$ is independent of the measure.
**References:** [5, 12, 27].

### 24. Rogers--Szego wrapped Gaussian — [`rogers_szego.m`](../rogers_szego.m)

**Order:** $N=n$. **Structure:** real orthogonal if the terminal boundary is
real. For $0<q<1$, the circle measure is a wrapped Gaussian with variance
$-\log q$, and its Verblunsky coefficients satisfy

```math
\alpha_k=(-1)^k q^{(k+1)/2}.
```

Its positive circular density is smooth and the $k$-th Fourier moment is
$q^{k^2/2}$. This is a nontrivial parameter-dependent real-coefficient
OPUC benchmark with a known cyclic-vector measure. **References:** [5, 12].

### 25. Constant Geronimus CMV — [`geronimus_cmv.m`](../geronimus_cmv.m)

**Order:** $N=n$. **Structure:** unitary five-diagonal by default;
orthogonal for real $a$ and real $\beta$. Use $\alpha_k=a$ for
all $k\ge0$, $|a|<1$. In the infinite model the measure can have a
proper absolutely continuous arc and a point mass. This generator
tests finite approximations to arc gaps with complex coefficients.
**References:** [5, 7].

### 26. Periodic Verblunsky CMV — [`periodic_cmv.m`](../periodic_cmv.m)

**Order:** $N=n$. **Structure:** sparse unitary. Repeat a nonempty
complex or real coefficient pattern $a_0,\ldots,a_{p-1}$ with all
$|a_j|<1$; close the first $n-1$ entries using a terminal unit
phase. The infinite periodic theory has bands and gaps governed by
a transfer-matrix discriminant. Boundary effects and finite eigenvalues
must be kept distinct from the infinite spectral measure.
**References:** [5, 7].

### 27. Finite generalized Bernstein--Szego — [`bernstein_szego_finite.m`](../bernstein_szego_finite.m)

**Order:** $N=n$. **Structure:** unitary, sparse five-diagonal by default.
Specify a finite complex coefficient list $a_0,\ldots,a_{m-1}$
strictly inside the disk and set $\alpha_k=0$ for $k\ge m$.
With the associated reversed polynomial $\Phi_m^*$, the positive
weight is

```math
\frac{d\mu}{d\theta}=
\frac{\prod_{j=0}^{m-1}(1-|a_j|^2)}
{2\pi|\Phi_m^*(e^{i\theta})|^2}.
```

The previous one-coefficient Poisson family is a special case.
These rational weights give arbitrarily chosen real or complex
Schur parameters with an exact and evaluable cyclic-vector density.
**References:** [5, 12].

### Lee-Yang Ising zeros via twisted Floquet-CMV — [`lee_yang_ising.m`](../lee_yang_ising.m)

**Order:** $N=2n$ for `n` Ising spins. **Application:** zeros of the partition
function of a one-dimensional ferromagnetic Ising model with periodic spin
boundary conditions. Given positive dimensionless magnetic couplings
`p_j=J_j/(k_B*T)`, set the `2*n` Verblunsky parameters

```math
\alpha_{2j}=\exp(-2p_{j+1}),\qquad \alpha_{2j+1}=0,
\quad 0\leq j<n.
```

They define an **extended periodic** CMV operator, whose finite Floquet
restriction at the *twisted* boundary phase `theta=pi/2` is a sparse
unitary `U=F_{2n}(pi/2)=L*M`. This is **not** the standard finite-CMV
paraorthogonal truncation: the cyclic boundary core uses the nonterminal
coefficient `alpha_{2n-1}=0` and a nontrivial complex twist. The
strict pencil `M-z*L'` is returned by terminal `'pencil'` mode.

The exact finite spectrum is the set of zeros, counted with multiplicity,
of the Laurent polynomial

```math
Z_n(z)=\sum_{s\in\{-1,+1\}^n}
\exp\!\left(\sum_{j=1}^n p_j s_j s_{j+1}\right)
z^{\sum_{j=1}^n s_j},\qquad s_{n+1}=s_1.
```

Thus `z^n Z_n(z)` has degree `2n` and contains only even powers.
Unlike the Floquet spin-chain propagators elsewhere in this collection,
this is a **classical statistical-mechanics partition-function** problem
whose complex-field zeros are realized as eigenvalues of a unitary matrix.

**Presets:** Fibonacci couplings `a->ab, b->a`, with
`(p_a,p_b)=(2/3,1/100)` as in Example 5.2 of [70]; uniform,
alternating, and explicitly supplied positive coupling vectors. The
Fibonacci limit can exhibit Cantor-type gaps and gap labels. The
`'spectral'` metadata does not claim a universal analytic trace density.
The construction uses the periodic Floquet boundary formulas and
propositions 2.5/3.2 of [70], with their connection to Schur parameter
pencils discussed in [69]. **References:** [69, 70].

### General multichannel scattering zippers — [`scattering_zipper.m`](../scattering_zipper.m)

**Order:** $N=nq$, with even `n>=2` sites and `q>=1` channels per
site. Marin and Schulz-Baldes [71] specify arbitrary local
`2q-by-2q` unitary scattering matrices `S_2,...,S_n` whose
upper-right `q-by-q` blocks are nonsingular, together with
independent `q-by-q` endpoint unitaries. The assembled operator is

```math
Z_n=V_n W_n,\quad
V_n=\mathrm{diag}(S_2,S_4,\ldots,S_n),\quad
W_n=\mathrm{diag}(U_{\mathrm L},S_3,S_5,\ldots,S_{n-1},U_{\mathrm R}).
```

Every allowed local scattering matrix has the representation
(Proposition 1 of [71])

```math
S(a,X,Y)=\begin{bmatrix}
a&(I-aa^*)^{1/2}X\\
Y(I-a^*a)^{1/2}&-Ya^*X
\end{bmatrix},
\quad \|a\|_2<1,\quad X,Y\in U(q).
```

The independent mixing factors `X,Y` and boundary scatterers make this
a strictly more general **parameterized finite model** than the
previously implemented canonical block-CMV factors. Its `q=1`
specializations can be gauge-equivalent to CMV models and are not
counted separately. The default generator has nontrivial deterministic
channel mixing; callers may pass arbitrary admissible local matrices.

`'pencil'` mode returns `A=W_n`, `B=V_n^*`, so that
`B^{-1}A=Z_n`. The construction has a sparse `O(Nq)` description.
Finite eigenvalues lie on the circle; no universal analytic trace
density is asserted. The source develops transfer-matrix methods,
matrix-valued spectral measures, and finite/periodic Prüfer-phase
eigenvalue counting. **Reference:** [71].

## G. Floquet spin chains and driven fermionic systems

### 28. Quantum kicked top — [`kicked_top.m`](../kicked_top.m)

**Order:** $N=n=2j+1$ (integer or half-integer spin $j>0$). **Structure:**
dense spin rotation times diagonal phase twist. In the $J_z$ basis
$m=j,j-1,\ldots,-j$, define $J_y=(J_+-J_+^*)/(2i)$ and

```math
U_F=\exp\left(-\frac{i\kappa}{2j}J_z^2\right)\exp(-ipJ_y).
```

The nonzero superdiagonal of $J_+$ is $\sqrt{k(n-k)}$ for $1\le k<n$.
Both factors are unitary. At $p=0$ the eigenvalues are
$\exp(-i\kappa m^2/(2j))$; at $\kappa=0$ they are $\exp(-ipm)$.
Otherwise the quasienergies vary with torsion and rotation.
Applications: quantum chaos and dense small-to-medium spectral tests.
The spin rotation can also be applied without forming a dense matrix
using final output `'operator'`. This uses the same sparse Hermitian-
tridiagonal generator and `expmv` action as the Bose–Hubbard dimer
(MATLAB R2023b+), with a Chebyshev–Bessel fallback on older releases.
Explicit mode uses `expm` of a sparse tridiagonal generator, whose
exponential is generally dense. No uniformly $O(N)$-time action is claimed.
**References:** [28, 29].

### 29. Integrable kicked Ising chain — [`floquet_ising.m`](../floquet_ising.m)

**Order:** $N=2^n$ for $n$ spins. **Structure:** a diagonal Ising layer
followed by a Kronecker product of single-spin rotations:

```math
U_F=\exp\left(-iJ\sum_rZ_rZ_{r+1}\right)
    \exp\left(-ih\sum_{r=1}^nX_r\right).
```

The bond sum includes $(n,1)$ for periodic boundaries and excludes it
for open boundaries. The second layer is
$(\cos h\,I-i\sin h\,X)^{\otimes n}$.
In the computational basis (spin $r$ corresponds to bit $r-1$),
the first layer has exact diagonal entries
$\exp(-iJ\sum_rz_rz_{r+1})$, $z_r=\pm1$.
With $n=2$, the periodic convention counts both directed bonds.
The **full many-body eigenphase measure** is not the single-particle
quasiparticle measure. Dense assembly costs $O(4^n)$ memory and is
capped at $n=11$ sites. **References:** [30, 31].

### 30. Self-dual kicked Ising chain — [`selfdual_kicked_ising.m`](../selfdual_kicked_ising.m)

**Order:** $N=2^n$ for $n\ge3$ periodic spins, capped at $n=11$
for explicit dense storage. With real longitudinal fields $b_r$,

```math
U_F=\exp\left(-\frac{i\pi}{4}\sum_rX_r\right)
    \exp\left[-\frac{i\pi}{4}
      \left(\sum_rZ_rZ_{r+1}+\sum_rb_rZ_r\right)\right].
```

The first factor is a tensor product of $2\times2$ unitary kicks;
the second is diagonal. For zero longitudinal fields, the two layers
are reversed relative to the $J=h=\pi/4$ integrable Ising operator,
so the eigenvalue multiset is the same ($AB$ and $BA$ are similar
for invertible square matrices). Inhomogeneous longitudinal fields
provide the interacting model used in quantum-chaos research.
Random-matrix spectral-form-factor results require explicit
disorder and limiting assumptions and **do not** imply a universal
finite-size phase density. Default $b_r=0.3\sin r$ is deterministic.
**Reference:** [32].

### 31. Driven Kitaev chain, BdG representation — [`driven_kitaev.m`](../driven_kitaev.m)

**Order:** $N=2n$ for $n$ sites in Nambu basis
$(c_1,\ldots,c_n,c_1^\dagger,\ldots,c_n^\dagger)$;
**not** the $2^n$-dimensional Fock-space propagator.
Let

```math
H_{\mathrm{BdG}}(\mu)=
\begin{pmatrix}h(\mu)&\Delta\\ \Delta^T&-h(\mu)^T\end{pmatrix},
\qquad
U_F=e^{-i\tau_2H_{\mathrm{BdG}}(\mu_2)}
    e^{-i\tau_1H_{\mathrm{BdG}}(\mu_1)}.
```

Here $h_{rr}=-\mu$, nearest-neighbor $h_{r,r+1}=-t$,
$\Delta_{r,r+1}=\delta$ and $\Delta_{r+1,r}=-\delta$;
optional periodic closure preserves the oriented pairing bonds.
Parameters are supplied as a scalar MATLAB options struct.
For real parameters and
$C=\begin{psmallmatrix}0&I\\I&0\end{psmallmatrix}$,
$C\overline{U_F}C=U_F$, so phases occur in conjugate pairs.
Certain open-boundary regimes feature Majorana edge modes at
quasienergies $0$ or $\pi$; their presence is not universal.
Periodic $n=2$ is disallowed to avoid double-bond ambiguity.
**References:** [33, 34].

## H. Further quantum Floquet maps and Yang–Baxter circuits

### 32. Quantum kicked rotor — [`kicked_rotor.m`](../kicked_rotor.m)

**Order:** $N=n$, with normalized negative-exponent Fourier matrix $F_n$.
**Structure:** dense but explicitly factored into unitary Fourier and diagonal
stages. On the finite Fourier grid, define

```math
U_F=D_pF_n^*D_xF_n,\qquad
(D_p)_{kk}=\exp(-i\hbar k^2/2),\qquad
(D_x)_{jj}=\exp[-i(K/\hbar)\cos(2\pi(j+\eta)/n)],
```

for $j,k=0,\ldots,n-1$. The momentum labels are **uncentered** and
$\eta$ is a position-grid offset, not a general Bloch-twist parameter.
Default $\hbar=2\pi/n$, $K=5$, and $\eta=0$.
All factors are unitary, independently of finite $n$. The zero-kick
spectrum is known exactly:

```math
\sigma(U_F|_{K=0})=\{\exp(-i\hbar k^2/2):0\le k<n\}.
```

For nonzero kick strength, quasienergy statistics, clustering, and
localization-like behavior depend on dimension and the quantization
convention; no universal finite-size density is asserted.
This is a **specified finite-grid model**, not an arbitrary principal
truncation of the infinite kicked-rotor operator, and should not be
identified with every possible twisted torus quantization.
Applications: quantum chaos, phase aliasing and DFT-based matrix tests.
**Reference:** [35].

### 33. Dual-unitary brickwork Floquet circuits — [`dual_unitary_circuit.m`](../dual_unitary_circuit.m)

**Order:** $N=2^n$ for $n$ qubits; dense matrices are capped at $n=10$.
Each elementary two-qubit gate is

```math
G(\chi,a,b,c,d)=
(R_y(a)\otimes R_z(b))\,S\,
e^{-i\chi Z\otimes Z}\,
(R_z(c)\otimes R_y(d)),
```

where $S$ swaps the two qubits and $R_{y,z}$ denote single-qubit
rotations. This choice is **dual-unitary**: $G^*G=I_4$ and its
space–time reshuffling is also unitary. In explicit indices, define
$\widetilde G_{(a,c),(b,d)}=G_{(a,b),(c,d)}$, where the indices range over
$\{0,1\}$. Dual unitarity of the *local gate* does **not** mean
that the full many-body circuit is itself a dual-unitary $4\times4$
matrix or has a fixed one-point spectral law.

A Floquet period applies disjoint odd bonds first and even bonds
second, then repeats this period a nonnegative integer number of times.
Open boundaries are available for all $n\ge2$; periodic wrap bonds
require even $n\ge4$. The default is periodic for even $n\ge4$ and
open otherwise. For $n=2$ and open boundaries, the matrix is a single
two-qubit gate when depth is one; at zero depth it is the identity.
Spatially homogeneous fixed-angle gates are deterministic and do not
automatically inherit random-matrix claims derived for disordered
dual-unitary circuit ensembles.

Applications: structured many-body spectral problems, local
space–time-duality checks, products of low-order gates, and spectral-form-
factor experiments with appropriately chosen randomness.
**References:** [36, 37].

### 34. Integrable XXZ Floquet circuits — [`xxz_floquet.m`](../xxz_floquet.m)

**Order:** $N=2^n$ for $n$ spin-$1/2$ sites, capped at $n=10$ for
dense assembly. For real rapidity $u$, anisotropy $0<\gamma<\pi$,
and $\Delta=\cos\gamma$, define

```math
\rho=\sqrt{\sinh^2 u+\sin^2\gamma},\qquad
\check R(u,\gamma)=\frac{1}{\rho}
\begin{pmatrix}
\sinh(u+i\gamma)&0&0&0\\
0&i\sin\gamma&\sinh u&0\\
0&\sinh u&i\sin\gamma&0\\
0&0&0&\sinh(u+i\gamma)
\end{pmatrix}.
```

This is a normalized, **unitary** checked six-vertex
$R$-matrix. It also satisfies the braid-form Yang–Baxter equation,

```math
\check R_{12}(u)\check R_{23}(u+v)\check R_{12}(v)
=\check R_{23}(v)\check R_{12}(u+v)\check R_{23}(u).
```

The finite Floquet step is $U_F=U_{\rm even}U_{\rm odd}$, where the
odd and even bond layers comprise copies of this gate. For periodic
boundaries, $n$ must be even and at least four, with oriented wrap bond
$(n,1)$. Open boundaries are an explicitly specified finite-circuit
extension; integrable reflection-boundary conditions are **not**
automatically claimed for that case.

The gate and hence the global evolution preserve total $Z$ magnetization,
so sector-resolved spectra can differ substantially from the full
eigenphase measure. At $u=0$ the local gate equals $iI_4$ and the entire
operator is a known scalar phase: $i^nI$ (even periodic $n$) or
$i^{n-1}I$ (open $n$). This gives exact spectral and degeneracy checks.
Default $u=0.6$, $\gamma=\pi/3$; default boundaries are periodic
for even $n\ge4$ and open otherwise.

Applications: Yang–Baxter-integrable interacting Floquet operators,
conservation-law sectors, and test cases contrasting integrable and
chaotic many-body spectral statistics. **References:** [38, 39].

### Magnetic unitary almost-Mathieu walk — [`unitary_almost_mathieu.m`](../unitary_almost_mathieu.m)

**Order:** $N=2n$, for `n` cyclic positions with two spin states per
position. From [72], equations (2.3)–(2.6), with `lambda1,lambda2`
in `[0,1]`, magnetic flux frequency `phi` and phase `theta`
(in cycles), define `W=S_lambda1 Q` where

```math
Q_j=\begin{bmatrix}
\lambda_2\cos(t_j)+i\sqrt{1-\lambda_2^2}&-\lambda_2\sin(t_j)\\
\lambda_2\sin(t_j)&\lambda_2\cos(t_j)-i\sqrt{1-\lambda_2^2}
\end{bmatrix},
\qquad t_j=2\pi(j\phi+\theta).
```

The shift/reflection `S_lambda1` maps a `+` spin at site `j`
to `lambda1` times the forward `+` spin plus
`sqrt(1-lambda1^2)` times the on-site `-` spin, and maps a
`-` spin to `lambda1` times the backward `-` spin minus
`sqrt(1-lambda1^2)` times the on-site `+` spin.
Both factors are exactly unitary.

**Infinite-volume regimes:** subcritical `lambda1>lambda2`,
critical `lambda1=lambda2>0`, and supercritical
`lambda1<lambda2`, with the absolutely continuous / singular
continuous / pure-point spectral classifications subject to
the frequency/phase hypotheses in [72]. Critical infinite
irrational-frequency spectra can have zero-measure Cantor structure.
A finite matrix has only discrete eigenvalues, so these claims
must not be applied literally to finite eigenphase histograms.

The implemented finite matrix closes a window of `n` sites
periodically, preserving exact unitarity. For rational `phi=p/q`
with `q` dividing `n`, the result is a periodic-cell model.
For irrational flux it is a **periodically closed finite window**,
not a translationally periodic infinite operator.
`'pencil'` returns factors `A=Q,B=S^*`. Compared with the
existing real constant-coin `split_step_walk`, this benchmark
has complex magnetic quasiperiodic coins and two independent
couplings. **Reference:** [72].

### Agassi nuclear pairing and particle-hole propagator — [`agassi_propagator.m`](../agassi_propagator.m)

**Order:** $N=2^{4j}$ in the full Fock space of `4j` fermionic orbitals,
or `nchoosek(4j,N)` in a selected invariant `N`-particle sector.
The original Agassi model combines two levels (upper/lower), each with
`2j` degenerate magnetic substates, in a genuine **fermionic** model.
The Hamiltonian of Pérez-Fernández et al. [75], Eqs. (1)-(6), is

```math
H=\varepsilon J_0-g\sum_{\sigma,\sigma'=\pm1}
A_\sigma^*A_{\sigma'}
-\frac{V}{2}(J_+^2+J_-^2),
```

with `J_+=sum_m c_(+,m)^* c_(-,m)`,
`J_-=J_+^*`, `J_0=(N_+-N_-)/2` and
`A_sigma^*=sum_(m>0)c_(sigma,m)^*c_(sigma,-m)^*`.
All fermionic creation/annihilation operators use the
Jordan--Wigner sign convention in an occupation-number basis.
The Hamiltonian is sparse Hermitian and commutes with total
fermion number `N_++N_-`. Its propagator

```math
U(t)=\exp(-itH)
```

is unitary and generally dense. Parameter-dependent pairing `g`
and particle-hole interactions `V` generate competing superconducting
and deformed nuclear phases. These are structurally distinct from
the bosonic extended-Lipkin model and the quantum spin-chain
Floquet propagators already implemented.

**Finite benchmarks:** `j=1` gives `4` orbitals and `16` full Fock
states, matching the four-qubit simulation in [75]. That paper also
gives an explicit Pauli representation (Eqs. (14)-(17)).
The implemented fermionic and Pauli Hamiltonians are isospectral
up to an unimportant constant energy shift in the respective
occupation conventions, checked independently in the regression tests.
Specifying `particles=2` for `j=1` restricts to a `6`-state invariant
subspace. `j=2` gives up to `256` full Fock states or `70` states
at half filling.

When `g=V=0`, the eigenvalues are known exactly from diagonal
occupations: `exp(-it*epsilon*(Nplus-Nminus)/2)`. For `t=0`,
all eigenvalues equal one. Generic finite spectral data are otherwise
marked unknown. The `'hamiltonian'` mode returns sparse `H`,
and `'spectral'` provides the exact solvable limits without calling
`eig` on the full propagator. The initial implementation restricts
explicit unitary materialization to `j<=2` to avoid dense memory
costs. **Reference:** [75].


### Extended Lipkin model: collective two-level nuclear dynamics — [`extended_lipkin_propagator.m`](../extended_lipkin_propagator.m)

**Order:** $N=n+1$ for a fixed population of `n>=1` identical two-level bosons,
with basis indexed by `k=0,...,n` occupied `t` bosons.
The extended Lipkin Hamiltonian of Baid et al. [74], Eqs. (4)-(5), is

```math
H_{n,\lambda,\alpha}=(1-\lambda)n_t-\frac{\lambda}{n}Q_\alpha^2,
\qquad
Q_\alpha=s^*t+t^*s+\alpha\,t^*t .
```

The operators are the two bosonic modes `s` and `t` restricted to the
fixed-total-occupation sector. In the stated basis,

```math
(n_t)_{kk}=k,\qquad
(Q_\alpha)_{kk}=\alpha k,\qquad
(Q_\alpha)_{k,k+1}=(Q_\alpha)_{k+1,k}
=\sqrt{(k+1)(n-k)}.
```

Consequently `H` is **real symmetric pentadiagonal**, with only
`|i-j|<=2` nonzero. The generator returns the **unitary**
`U(t)=exp(-itH)`, which is generally dense. The optional terminal
`'hamiltonian'` mode returns the sparse symmetric matrix `H` without
forming the exponential.

The dimensionless parameters `0<=lambda<=1` and `alpha` determine
the nuclear shape-transition phase diagram. The large-`n` mean-field
critical line of [74] is `lambda_c=1/(5+alpha^2)`.
At `alpha=0` the symmetric-to-deformed transition is second-order;
for nonzero `alpha` it is generically first-order. These transitions
do **not** imply a genuine nonanalytic phase change in any fixed finite
matrix, though finite-size spectral precursors are numerically relevant.

Two exactly solvable finite benchmarks are explicitly available:

- `lambda=0`: `H=n_t`, so the phases are
  `exp(-it*k)`, `k=0,...,n`.
- `lambda=1, alpha=0`: `H=-4*S_x^2/n`, so
  `exp(4it*m^2/n)`, `m=-n/2,...,n/2`.
- `t=0`: `U=I` with all eigenvalues one.

These eigenphase formulas are returned lazily by `'spectral'` mode;
the generic trace spectrum is deliberately marked unknown.
The model is **not** a duplicate of `kicked_top` (kicked spin rotations)
or `kicked_bose_hubbard` (two-step Floquet dimer): here the physical
dynamics is governed by a time-independent nuclear shape-transition
Hamiltonian with first- and second-neighbor couplings. **Reference:** [74].


### Two-dimensional U(3) vibron ESQPT and symmetry-sector clustering — [`vibron_2d_esqpt.m`](../vibron_2d_esqpt.m)

**Order:** $N=(n+1)(n+2)/2$ in the fixed-total-boson-number `n`
Hilbert space with one scalar boson `sigma` and two circular bosons
`tau_+`, `tau_-`. Alternatively, specify conserved angular momentum
`ell=n_+-n_-` for an individual `floor((n-|ell|)/2)+1` dimensional
sector. The physical model describes bending vibrations and an
excited-state quantum phase transition (ESQPT). It is a **U(3) model**,
not another U(2) Lipkin parameterization.

The model in [76], Eqs. (5)-(7), is

```math
H_N(\xi)=(1-\xi)n_\tau+\frac{\xi}{n}
\left[n(n+1)I-W^2\right],\qquad
n_\tau=n_++n_-,\quad 0\leq\xi\leq 1,
```

where `W^2=(D_+D_-+D_-D_+)/2+ell^2`,
`D_+=sqrt(2)*(tau_+^* sigma-sigma^* tau_-)` and
`D_-=D_+^*`. In a fixed `ell` sector use
`t=n_tau=|ell|,|ell|+2,...,n` and
`s=n-t`, `p=(t+ell)/2`, `m=(t-ell)/2`:

```math
(W^2)_{tt}=2st+2s+t+\ell^2,\qquad
(W^2)_{t,t+2}=-2\sqrt{s(s-1)(p+1)(m+1)}.
```

Thus the sector Hamiltonian `H_(n,ell)` is **real symmetric
tridiagonal**; the full Hamiltonian is block diagonal in
`ell=-n,...,n`. We supply its unitary evolution
`U(t)=exp(-itH)` in the same sector order. Its full matrix
is sparse block diagonal, although each nontrivial exponential block
is typically dense. `'hamiltonian'` exposes the sparse structured
matrix without computing the exponential.

**Exact finite limits:** At `xi=0`, the energies are `t=n_tau`.
At `xi=1`, the SO(3) dynamical-symmetry label
`omega=n,n-2,...,|ell|` gives

```math
E_{n,\ell,\omega}=\frac{n(n+1)-\omega(\omega+1)}{n},
\qquad z_{n,\ell,\omega}(t)=e^{-itE_{n,\ell,\omega}}.
```

These formulas provide exact eigenphases and trace moments in
`'spectral'` mode; `t=0` gives the identity for any `xi`.
The generic finite eigenphase law is not claimed analytically known.

**Clustered eigenphases:** For `xi>1/5` the broken-symmetry regime
has an excited-state separatrix. Pairs of low-lying energies from
`ell=0` and `ell=1` sectors approach one another **algebraically**
as `n` increases. The comparable even/odd parity splittings of the
U(2) Lipkin model can decrease **exponentially** (ultimately below
floating-point resolution); see Fig. 2 of [76]. The `ell=+q`
and `ell=-q` sector matrices themselves are exactly identical in
this simplified Hamiltonian, a distinct symmetry degeneracy from the
near-degeneracy of the `ell=0,1` sectors. Those distinctions matter
for structured eigensolvers and spectral-density estimators.

The example [`examples/vibron_esqpt_sector_gaps.m`](../examples/vibron_esqpt_sector_gaps.m)
compares the sector spectra, unit-circle eigenphases, and finite-size
cross-sector gaps without presupposing an asymptotic density.
Independent tests construct `D_+` in the three-boson occupation
basis and verify the tridiagonal blocks and solvable SO(3) limits.
**Reference:** [76].

## I. Bosonic, cavity-QED, and quantum Hall network operators

### 35. Kicked Bose–Hubbard dimer — [`kicked_bose_hubbard.m`](../kicked_bose_hubbard.m)

**Order:** $N=n+1$ for $n$ identical bosons on **two sites**, in the
fixed-occupation Fock basis $|k,n-k\rangle$ ($0\le k\le n$).
The two Hermitian generators are

```math
H_{\rm hop}=-J(a_1^\dagger a_2+a_2^\dagger a_1),\qquad
H_{\rm kick}=\frac{g}{2}\sum_{r=1}^{2}n_r(n_r-1)
                 +\frac{b}{2}(n_1-n_2).
```

The Floquet operator in the stated order is

```math
U_F=\exp(-i\tau_{\rm kick}H_{\rm kick})
     \exp(-i\tau_{\rm hop}H_{\rm hop}).
```

The kick is diagonal with entries

```math
E_k=\frac{g}{2}\bigl[k(k-1)+(n-k)(n-k-1)\bigr]
      +\frac{b}{2}(2k-n),
```

and the hopping Hamiltonian is real symmetric tridiagonal, with
adjacent off-diagonal entries $-J\sqrt{(k+1)(n-k)}$.
Thus both factors are exactly unitary mathematically.
In explicit mode the Hermitian tridiagonal hopping generator is assembled
as a MATLAB sparse matrix and passed to `expm`; the exponential itself
is generally dense. In matrix-free mode MATLAB `expmv` (R2023b+) acts on
that sparse generator without forming the exponential, with a
Chebyshev–Bessel fallback for older MATLAB releases.
This is a number-conserving **many-boson** operator of order $n+1$,
rather than a $2^n$ spin-chain matrix.

When $g=b=0$ the **complete exact spectrum** is

```math
\sigma(U_F)=
\bigl\{\exp[-i\tau_{\rm hop}J(2k-n)]:k=0,\ldots,n\bigr\}.
```

When $J=0$ the eigenvalues are exactly
$\exp(-i\tau_{\rm kick}E_k)$. For $n=0$ the vacuum sector
is the scalar unitary $[1]$.
Applications: Bose–Einstein condensates in double wells, quantum
chaos, Floquet tunneling, self-trapping and avoided quasienergy
crossings. Related to spin-$n/2$ kicked tops via Schwinger bosons,
but the interaction and phase conventions are specified independently.
The general interacting eigenphase law is parameter-dependent.
**References:** [40, 41]; the alternate periodically biased dimer
model of [42] is related but uses a different drive.

### 36. Jaynes–Cummings cavity-QED propagator — [`jaynes_cummings.m`](../jaynes_cummings.m)

**Order:** $N=2(n+1)$ for a photon cutoff at occupation $n$.
The **interleaved atom–photon basis** is
$|g,0\rangle,|e,0\rangle,\ldots,|g,n\rangle,|e,n\rangle$.
In units $\hbar=1$, write

```math
H=\omega_c a^\dagger a+\frac{\omega_a}{2}\sigma_z+
g(a\sigma_++a^\dagger\sigma_-),\qquad
\Delta=\omega_a-\omega_c.
```

The model compresses the **Hermitian Hamiltonian** to levels
$0,\ldots,n$ before taking $U(t)=\exp(-itH_n)$.
This avoids the nonunitarity of principal submatrices of an infinite
unitary. The excitation blocks act on $|e,k\rangle,|g,k+1\rangle$
for $k=0,\ldots,n-1$:

```math
H_k=\omega_c(k+\tfrac12)I_2+
\begin{pmatrix}
\Delta/2 & g\sqrt{k+1}\\
g\sqrt{k+1} & -\Delta/2
\end{pmatrix}.
```

Set $r_k=\sqrt{(\Delta/2)^2+g^2(k+1)}$.
The implemented analytic block exponential is

```math
e^{-itH_k}=e^{-it\omega_c(k+1/2)}
\left[\cos(tr_k)I_2-
i\frac{\sin(tr_k)}{r_k}
\begin{pmatrix}\Delta/2&g\sqrt{k+1}\\
g\sqrt{k+1}&-\Delta/2\end{pmatrix}\right],
```

using the continuous $r_k=0$ limit.
The **complete finite spectrum** consists of

```math
e^{-it[\omega_c(k+1/2)\pm r_k]},\quad 0\le k<n,
\qquad e^{it\omega_a/2},\quad
e^{-it(n\omega_c+\omega_a/2)}.
```

The last two are the isolated boundary states $|g,0\rangle$ and
$|e,n\rangle$, respectively. The upper state is an artifact of
**finite photon truncation** and not isolated in the infinite model.
The propagator is sparse with two-by-two blocks, making it useful
for exact trace moments, eigenphase clusters and repeated degeneracies.
**Reference:** [43].

### 37. Chalker–Coddington directed network — [`chalker_coddington.m`](../chalker_coddington.m)

**Order:** $N=4n_xn_y$, with directed links $(x,y)$ on the periodic
$2n_x\times2n_y$ grid, ordered by the MATLAB index
$1+(x\bmod 2n_x)+2n_x(y\bmod 2n_y)$.
The implementation follows the **one-step oriented network**
in Eq. (3.3) of [46], not an arbitrary coined-walk substitute.
The two types of $2\times2$ scattering nodes are

```math
S_e=\begin{pmatrix}
\cos\theta_e&\sin\theta_e\\
-\sin\theta_e&\cos\theta_e
\end{pmatrix},\qquad
S_o=\begin{pmatrix}
\sin\theta_o&\cos\theta_o\\
-\cos\theta_o&\sin\theta_o
\end{pmatrix}.
```

An even node maps incoming $(2j,2k)$, $(2j+1,2k-1)$
to outgoing $(2j+1,2k)$, $(2j,2k-1)$.
An odd node maps incoming $(2j+2,2k+1)$, $(2j+1,2k)$
to outgoing $(2j+1,2k+1)$, $(2j+2,2k)$.
All coordinates are taken modulo the even lattice dimensions.
Every input link and output link belongs to exactly one node;
hence the global scattering matrix $S$ is unitary. Assigning
independent input-link phases gives

```math
U=S\,\mathrm{diag}(e^{i\phi_1},\ldots,e^{i\phi_{4n_xn_y}}).
```

The result has at most two nonzero entries per row and column.
A scalar angle applies to both node types; a two-vector gives
independent $\theta_e,\theta_o$. By default phases are zero;
for Anderson-localization experiments, supply random phases explicitly.
At angles in $\{0,\pi/2\}$ the matrix is a phase-weighted permutation.
At $\theta_e=\theta_o=\pi/4$, both node types are balanced.

**Exact spectrum constraint:** with
$C_{(x,y),(x,y)}=(-1)^{x+y}$, the directed network alternates
between the two link parities. Thus

```math
CUC=-U,\qquad \sigma(U)=-\sigma(U),\qquad
\mathrm{tr}(U^{2r+1})=0\quad(r\ge0).
```

This $\pi$-paired eigenphase symmetry holds with arbitrary link phases
and node angles, and gives a particularly useful exact benchmark for
unit-circle spectral-density estimators.
Applications: integer quantum Hall plateau transitions, localization,
critical scattering and sparse unitary operators. No universal
finite-size one-point eigenphase density is asserted.
**References:** [44–46].

### Szegedy walks for finite Markov chains — [`szegedy_walk.m`](../szegedy_walk.m)

**Order:** $N=n^2$ for an `n-by-n` row-stochastic matrix `P`.
In the ordered basis `|i,j>`, introduce the isometries

```math
A e_i=\sum_j\sqrt{P_{ij}}\,|i,j\rangle,\qquad
B e_j=\sum_i\sqrt{P_{ji}}\,|i,j\rangle.
```

The Szegedy quantum evolution [26] is the product of reflections

```math
W=(2BB^*-I)(2AA^*-I),\qquad
D=A^*B,\quad D_{ij}=\sqrt{P_{ij}P_{ji}}.
```

Let `sigma_j` be the singular values of the real-symmetric
discriminant `D`. Every `0<sigma_j<1` produces
`exp(+/-2i*acos(sigma_j))`; every zero singular value produces two
`-1` eigenvalues. The remaining eigenvalues are `+1`. In particular,
the entire spectrum can be obtained without diagonalizing the generally
dense `n^2-by-n^2` unitary matrix. A lazy exact finite spectral
reference evaluates the singular values of `D` only when needed.
The `'operator'` mode applies the two reflections directly in
`O(n^2)` time and memory for dense data. This is not another
finite-CMV realization or a Grover search reflection on marked states.
**Reference:** [26].

### Coined quantum walks on arbitrary finite graphs — [`coined_graph_walk.m`](../coined_graph_walk.m)

**Order:** $N=2m$, the number of directed arcs of an undirected graph with
`m` edges. Vertex `v` has `d_v` outgoing arcs and a local unitary
coin `C_v`. Use the canonical arc ordering sorted by source and
destination. The shift `S` reverses each arc, so `S^2=I`. Thus

```math
U=SC,\qquad C=\bigoplus_{v\in V}C_v,\qquad S|v,w\rangle=|w,v\rangle.
```

The default weighted Grover coin is

```math
C_v=2s_vs_v^*-I,\qquad
(s_v)_w=\sqrt{A_{vw}\big/\sum_u A_{vu}},
```

where `A` is a symmetric nonnegative edge-weight matrix with zero
diagonal and no isolated vertices. `'fourier'`, `'identity'`, and
custom vertex-local unitary coins are supported. Uniform Grover coins
on an unweighted connected graph admit an explicit spectral
relationship to the normalized adjacency/random-walk matrix; general
local coins need not. `'operator'` applies the block coins and
arc reversal without assembling the product; `'pencil'` exposes
the strict unitary pair `A=C`, `B=S^*`. This is genuinely
different from fixed one-dimensional walk topologies. **Reference:** [73].

### Parameterized metric quantum-graph scattering maps — [`metric_quantum_graph.m`](../metric_quantum_graph.m)

**Order:** $N=2m$ directed bonds on a finite undirected metric graph.
`L` is a symmetric zero-diagonal matrix whose positive off-diagonal
entries are physical edge lengths; zeros indicate missing edges. Using
outgoing directed bonds, the real-wavenumber bond-scattering map is

```math
U(k)=D(k)CR,\qquad
D(k)=\mathrm{diag}(e^{ik\ell_b}),\quad
R|v,w\rangle=|w,v\rangle .
```

`R` reverses directed bonds; `C` is the direct sum of the
vertex-local unitary scattering matrices. For Neumann--Kirchhoff
boundary conditions, `C_v=2J_{d_v}/d_v-I`. For any real `k`,
`D(k)`, `C`, and `R` are unitary, hence so is `U(k)`.
Arbitrary unitary vertex scatterers may also be supplied.

**Two distinct spectral questions** must be separated:

- At a fixed real `k`, `eig(U(k))` contains unit-circle eigenvalues
  of the finite bond-scattering map, giving a valid unitary matrix test.
- The physical metric-graph wavenumbers are the *real zeros* of
  `det(I-U(k))=0`, in which `k` itself is the spectral variable.
  In general this is a transcendental nonlinear eigenvalue problem.

`'operator'` supplies matrix-free unitary action at fixed `k`.
`'secular'` supplies analytic function handles `U(z)`,
`T(z)=I-U(z)`, and `det(T(z))` for complex `z`. The analytic
continuation `U(z)` is **not** generally unitary when `z` is
nonreal. Unlike `unitary_delay_pencil`, incommensurate real
lengths have no general finite-dimensional polynomial lift. This
generator reuses directed-arc and local-coin components, but its
metric dependence and nonlinear secular problem are new.
**References:** [58, 68].


## J. Non-quantum signal-processing transforms

This section gives independent finite matrix constructions. They are distinct
from the FFT, wavelet and paraunitary filter-bank families discussed earlier;
where a generator specializes Fourier conjugation, the spectral
redundancy is explicit. Frequencies of a **signal** are not to be mistaken for
eigenphases of the **transform matrix**.

### Orthonormal cosine transforms — [`orthogonal_dct.m`](../orthogonal_dct.m)

**Order:** $N=n$. The parameter `kind` selects DCT-I through DCT-IV.
For type II, with zero-based indices $j,k$, the matrix is

```math
U_{kj}=\sqrt{\frac{2}{n}}\,c_k
\cos\left(\frac{\pi k(j+1/2)}{n}\right),
\qquad c_0=2^{-1/2},\quad c_k=1\ (k>0).
```

Type III is its transpose; type IV uses
$\sqrt{2/n}\cos(\pi(k+1/2)(j+1/2)/n)$. Type I uses
$\sqrt{2/(n-1)}c_kc_j\cos(\pi kj/(n-1))$ for $n\geq2$,
with $c_0=c_{n-1}=1/\sqrt2$.

All variants are **real orthogonal** with FFT-based implicit actions; I and IV are also symmetric
involutions. Eigenphases of DCT-II/III generally depend on order and
are unrelated to the transform-frequency sample points. Uses:
compression, image coding, cosine-basis spectral methods. [52]

### Discrete Hartley — [`discrete_hartley.m`](../discrete_hartley.m)

**Order:** $N=n$. The orthonormal kernel is

```math
H_{jk}=\frac{\cos(2\pi jk/n)+\sin(2\pi jk/n)}{\sqrt n}.
```

This is real, symmetric, orthogonal, and $H^2=I$.
The entire spectrum is in $\{-1,+1\}$; exact multiplicities may be
determined from $n$ and the trace. Useful as a two-atom spectrum in
real signal processing, independent of the Hadamard construction. [47]

### Chu/CAZAC unitary circulants — [`cazac_circulant.m`](../cazac_circulant.m)

**Order:** $N=n$. Choose an integer root $r$ coprime to $n$.
For $k=0,\ldots,n-1$, define the constant-amplitude sequence as follows.

For **even** $n$,

```math
c_k=\exp\left(-\frac{\pi i r k^2}{n}\right).
```

For **odd** $n$,

```math
c_k=\exp\left(-\frac{\pi i r k(k+1)}{n}\right).
```

Let $d$ be the residue of $j-k$ modulo $n$. The unitary circulant is

```math
U_{jk}=\frac{c_d}{\sqrt{n}}.
```

Chu's perfect periodic autocorrelation yields $U^*U=I$,
with **constant-modulus entries** $|U_{jk}|=n^{-1/2}$.
The exact eigenvalues, including multiplicities, are
$z_\ell=n^{-1/2}\mathop{\mathrm{fft}}(c)_\ell$.
This is a restricted, application-derived subfamily of
unitary circulants (not a distinct general matrix class).
Applications: channel synchronization, radar/sonar pulse coding. [48]
### Discrete Fresnel transform — [`fresnel_transform.m`](../fresnel_transform.m)

**Order:** $N=n$. For $d=(j-k)\bmod n$ and
$\delta=0$ (even $n$) or $\delta=1/2$ (odd $n$), set

```math
U_{jk}=\frac{e^{-i\pi/4}}{\sqrt n}
\exp\left(\frac{i\pi(d+\delta)^2}{n}\right).
```

Both parity choices are circularly periodic and exactly unitary.
It is circulant; its eigenvalues are the FFT of its first column.
OCDM uses this discrete Fresnel transform for chirp multiplexing.
The matrix is unitary independently of the physical channel. [49]

### Affine Fourier / AFDM modulation — [`affine_fourier.m`](../affine_fourier.m)

**Order:** $N=n$. For real parameters $c_1,c_2$, define

```math
A=D(c_2)F_nD(c_1),\qquad
D(c)_{kk}=e^{-2\pi i c k^2},\quad k=0,\ldots,n-1.
```

where $F_n$ is the normalized negative-exponent DFT.
This is an exactly unitary chirp-DFT-chirp multiplication.
At $c_1=c_2=0$ it reduces to the DFT, whose spectrum has
at most four atoms; away from that point the phase law can differ.
Application: affine frequency-division multiplexing in
delay-Doppler channels. Do not conflate modulation $A$ with the
generally **nonunitary** wireless channel. [50]

### Finite Zak transforms — [`finite_zak.m`](../finite_zak.m)

**Order:** $N=n=pq$, $p\mid n$. Use the indexing $j+pk$
for $0\leq j<p$, $0\leq k<q$. Apply a normalized $q$-point
DFT to the second coordinate for each fixed first coordinate:

```math
Z_{n,p}=F_q\otimes I_p.
```

This is unitary; $p=1$ gives the DFT, and $p=n$ gives $I_n$.
It is a particular finite normalized Zak convention, used in
Gabor analysis and polyphase signal processing. Its spectrum
is exactly that of $F_q$, each eigenvalue repeated $p$ times. [51]

### Finite Weyl–Heisenberg shifts — [`time_frequency_shift.m`](../time_frequency_shift.m)

**Order:** $N=n$. Define $T e_j=e_{j+1\bmod n}$ and
$M e_j=e^{2\pi i j/n}e_j$, with zero-based indices.
The generator returns $U=M^pT^q$ for integer $p,q$.
It is a sparse monomial unitary, and its exact eigenvalues
are determined by the cycles of translation $q$ and the
product of modulation weights around each cycle. Uses:
Gabor time-frequency shifts and Weyl–Heisenberg systems. [51,53]

### Periodic Hilbert transform completion — [`periodic_hilbert.m`](../periodic_hilbert.m)

**Order:** $N=n$. Start from the DFT multiplier $-i$ at
positive frequencies and $+i$ at negative frequencies.
Replace the conventional *zero* DC multiplier (and, for
even $n$, zero Nyquist multiplier) by signs $s_0,s_\pi\in\{-1,+1\}$.
The resulting multiplier $d$ satisfies
$|d_j|=1$ and the conjugation symmetry needed for a real matrix.

```math
U=F_n^*\mathop{\mathrm{diag}}(d_0,\ldots,d_{n-1})F_n.
```

Hence $U$ is **real orthogonal and circulant**, with exact spectral
atoms among $\{-i,+i,-1,+1\}$. The zero-ended conventional periodic
Hilbert filter is **not** unitary on the whole finite-dimensional
space. The endpoint completion is a benchmark choice, not a claim
that all conventional Hilbert transforms are orthogonal. [52]

### Normalized discrete Fourier — [`discrete_fourier.m`](../discrete_fourier.m)

**Order:** $N=n$. Define the negative-exponent DFT by

```math
F_{jk}=\frac{e^{-2\pi i jk/n}}{\sqrt n},\quad 0\leq j,k<n.
```

Then $F^*=F^{-1}$ and $F^2$ is the reversal permutation
$j\mapsto-j\bmod n$. Therefore $F^4=I$ and all eigenvalues
belong to $\{1,-1,i,-i\}$. The default materialization is
dense; operator mode reuses the existing normalized FFT kernel.
This is a distinct generator interface but the same DFT
already appearing implicitly in several earlier matrix families.
Applications: Fourier spectral analysis, convolution, OFDM. [1,52]

### Orthonormal sine transforms — [`discrete_sine.m`](../discrete_sine.m)

**Order:** $N=n$. Type II uses the kernel

```math
S_{kj}=\sqrt{\frac2n}\,c_k
\sin\left(\frac{\pi(k+1)(j+1/2)}{n}\right),
\quad c_{n-1}=\frac1{\sqrt2},\quad c_k=1\ (k<n-1).
```

Type III is the transpose of type II. Type I has
$\sqrt{2/(n+1)}\sin(\pi(k+1)(j+1)/(n+1))$ and
type IV has $\sqrt{2/n}\sin(\pi(k+1/2)(j+1/2)/n)$.
All four are real orthogonal; types I and IV are symmetric
involutions. Their operator mode uses an odd FFT extension
or the existing fast DCT kernel with signs/reversals, avoiding
dense materialization. Applications: odd-boundary spectral
methods, transform coding, and PDE discretizations. [54]

### Walsh–Hadamard transforms — [`walsh_hadamard.m`](../walsh_hadamard.m)

**Order:** $N=n=2^m$. Natural Sylvester order satisfies

```math
H_1=[1],\qquad
H_{2n}=\frac1{\sqrt2}
\begin{pmatrix}H_n&H_n\\H_n&-H_n\end{pmatrix}.
```

Hence $H_n$ is real, orthogonal, symmetric, involutory,
and every entry has magnitude $1/\sqrt n$. Sequency order
permutes the rows into increasing zero-crossing count using
bit-reversed Gray indices; the permutation generally changes
the eigenvalues of the *transform matrix*, though its singular
values remain 1. The matrix-free kernel uses a fast Walsh
butterfly and row permutations, with $O(N)$ storage.
Applications: spread-spectrum signals, coding and image
processing. [55]

### Slepian/DPSS orthogonal basis — [`dpss_transform.m`](../dpss_transform.m)

**Order:** $N=n$, with half-bandwidth $0<W<1/2$.
Let the time-band concentration matrix $B$ have entries

```math
B_{jk}=\begin{cases}
2W,&j=k,\\
\dfrac{\sin(2\pi W(j-k))}{\pi(j-k)},&j\neq k.
\end{cases}
```

A commuting symmetric tridiagonal matrix has
$T_{jj}=((n-1)/2-j)^2\cos(2\pi W)$ and
$T_{j,j+1}=(j+1)(n-j-1)/2$, using zero-based $j$.
The generator forms its orthonormal eigenvectors in descending
concentration order and returns their transpose. This defines
a real orthogonal transform; its eigenphases are **not**
the concentration eigenvalues of $B$. The matrix is **explicit
only**: the generic DPSS eigenbasis does not admit a validated
fast $O(N)$-storage transform merely by wrapping its dense
eigenvectors in a function handle. Applications: multitaper
estimation, finite-band concentration and leakage suppression. [56]

### Graph Fourier transform — [`graph_fourier.m`](../graph_fourier.m)

**Order:** $N=n$. For a symmetric, nonnegative zero-diagonal
adjacency $A$, let $L=\mathop{\mathrm{diag}}(A\mathbf{1})-A$.
If $L=V\Lambda V^T$, sorted by nondecreasing eigenvalues,
the transform is $U=V^T$. It is real orthogonal. The generator
uses a path graph by default and fixes repeated eigenspaces
by vertex-order projection/orthonormalization. **Graph
Laplacian frequencies are not the eigenvalues of $U$**.
This is **explicit only**: a generic weighted graph has no
automatic fast graph-Fourier transform. Applications: graph
filtering and irregular-domain data analysis. [57]


## K. Unitary eigenvalue problems and strictly unitary pencils

An eigenvalue problem may possess a circle spectrum without being a unitary
matrix eigenproblem. We use the **strict Euclidean convention** for the
regular pencil `A-zB`, with nonsingular `B`:

```math
U=B^{-1}A\in U(N)
\quad\Longleftrightarrow\quad AA^*=BB^*.
```

This is distinct from `A*B=B*A`, `A'*A=B'*B`, or reciprocal eigenvalue
pairing. All seven formulations below are **Anymatrix-compatible numeric unitary
matrices by default**, including the existing `finite_cmv`. Terminal `'pencil'` mode returns
`P.A`, `P.B`, and `P.U`, without changing the standard two-output
Anymatrix generator interface. The extra
[eigenvalue-problem registry](../metadata/eigenvalue_problems.json) separately
classifies nonlinear functions and external nonunitary examples.

### Strict unitary Möbius pencil — [`mobius_unitary_pencil.m`](../mobius_unitary_pencil.m)

**Order:** $N=n$. Let `P=cycle_shift(n)`, `|a|<1`, and set

```math
A=P-aI,\quad B=I-\overline a P,\quad U=B^{-1}A.
```

The spectrum is explicitly `(omega-a)/(1-conj(a)*omega)` for the `n`th
roots of unity `omega`. The limiting **equal-weight trace density** is

```math
\rho_a(\theta)=\frac{1-|a|^2}
 {2\pi|1+\overline a e^{i\theta}|^2},
\qquad m_k=(-a)^k\ \ (k\geq0).
```

A synthetic benchmark derived from disk automorphisms, with a known
nonuniform limit and very concentrated phases as `|a|` approaches 1.
This is not the weighted OPUC Poisson measure of a different finite CMV
matrix. **Reference:** [5] (background); formula and pencil specialization
are derived from the displayed construction.

### Lossless integer-delay network pencil — [`unitary_delay_pencil.m`](../unitary_delay_pencil.m)

**Order:** $N=\sum_j\ell_j$, not the number `n` of channels.
For positive integer delays `ell_j` and unitary scattering `S`,
the polynomial nonlinear secular equation is

```math
T(z)=I-\mathop{\mathrm{diag}}(z^{\ell_j})S,
\qquad \det T(z)=0.
```

Form a state-space lift by cyclically shifting each delay line and
mixing its terminal coordinates through `S^*`. If `P` is the combined
delay permutation and `C` the terminal scattering operator, the finite
unitary lift is `U=PC`; `A=C, B=P^*` form a strict sparse pencil.
The unitary lift has exactly the polynomial secular zeros including
multiplicity. Distinguish this **commensurate integer-delay** model from
transcendental quantum graphs with incommensurate metric lengths.
**Applications:** quantum-graph scattering; lossless feedback delay
networks. **References:** [58, 60].

### Scalar CMV factor pencil — [`finite_cmv.m`](../finite_cmv.m)

**Order:** $N=n$; interior Verblunsky coefficients in the open disk, final
coefficient unimodular. Take the existing finite-CMV factorization
`C=LM`, with sparse alternating unitary `2-by-2` core blocks and
terminal closure. The existing `finite_cmv` generator exposes the generalized coefficients
`A=M`, `B=L^*` via terminal `'pencil'` mode, satisfying `B\\A=C`.
Its circle spectrum consists of scalar paraorthogonal polynomial zeros;
no duplicate generator is registered.
**References:** [6, 7].

### Finite block-CMV matrix, matrix spectral measure and strict pencil — [`block_cmv_pencil.m`](../block_cmv_pencil.m)

**Order:** $N=nq$. An independent q-by-q contractive matrix Verblunsky
coefficient alpha_j is supplied for each of n-1 core interfaces.
The q-by-q terminal beta is unitary. The alternating block factors
L and M are products of two-fiber matrices

```math
\Theta(\alpha_j)=
\begin{bmatrix}
\alpha_j^*&\rho_j^L\\
\rho_j^R&-\alpha_j
\end{bmatrix},\quad
\rho_j^L=(I-\alpha_j^*\alpha_j)^{1/2},\quad
\rho_j^R=(I-\alpha_j\alpha_j^*)^{1/2}.
```

The two defect matrices are generally **not equal**, and
noncommuting alpha_j cannot be simplified to scalar coefficients.
The numeric default is already the **standalone finite unitary
matrix U=L M**; this is not merely a generalized pencil.
The strict pencil option gives (A,B)=(M,L^*) with B^{-1}A=U.
No duplicate matrix generator is registered. q=1 agrees
with `finite_cmv` using the identical boundary phase convention.

The new `'operator'` output applies even then odd two-fiber cores
and their boundary blocks in O(N q) work per right-hand side;
the adjoint reverses both factor layers and uses the
Hermitian adjoints of the noncommuting cores. It supports
multiple complex right-hand sides without constructing U.

The `'quadrature'` mode (moderate N only) computes an eigendecomposition
U V=V diag(z_j) and returns the **matrix-valued first-block
cyclic spectral weights**

```math
W_j=V_{1:q,j}V_{1:q,j}^*,\qquad
\sum_jW_j=I_q,\qquad
\sum_jz_j^kW_j=E^*U^kE,\quad E=\begin{bmatrix}I_q\\0\end{bmatrix}.
```

These weights are positive semidefinite q-by-q matrices,
not the scalar trace weights 1/(nq). The trace moments are
(1/(nq))sum_j z_j^k. In a repeated-eigenvalue cluster,
only the sum of the matrix weights across a full eigenspace
is basis invariant. This optional **numerical quadrature**
does not constitute an analytic closed-form trace phase law;
the generic terminal `'spectral'` descriptor remains unknown.
**References:** Damanik–Pushnitski–Simon [59], and [67].

### Floquet factor pencil — [`floquet_factor_pencil.m`](../floquet_factor_pencil.m)

**Order:** $N=2n$. With periodic two-state partial shifts `Splus,Sminus`
and real orthogonal on-site coin rotations `C1,C2`,

```math
U=S_-C_2S_+C_1,\quad
A=S_+C_1,\quad B=(S_-C_2)^T.
```

The generalized problem is exactly the previously implemented
`split_step_walk`, without forming an unstructured linearization.
The pencil is real orthogonal with sparse factors. **Reference:** [16].

### Hermitian Cayley pencil — [`cayley_hermitian_pencil.m`](../cayley_hermitian_pencil.m)

**Order:** $N=n$. From the same real symmetric Jacobi recurrence `H`
used by `cayley_jacobi`, take `A=H-itI`, `B=H+itI`, `t>0`.
Then `AA^*=BB^*` and the spectrum is the Cayley image of the known
real Jacobi eigenvalues. The pencil coefficients are sparse tridiagonal
although the quotient is generally dense. **Reference:** [23].

### Ill-scaled strict unitary pencil — [`left_weighted_unitary_pencil.m`](../left_weighted_unitary_pencil.m)

**Order:** $N=n$. Given the cyclic permutation `P`, choose a positive
diagonal `B=diag(d_j)` with endpoint ratio `kappa>=1` and set
`A=BP`. Then `B^{-1}A=P` is real orthogonal even when `A` and
`B` are far from orthogonal. The eigenvalues are exactly roots of
unity; coefficient conditioning is controlled independently.
This is a **synthetic stress test**, not an application family
attributed to the background reference [4].

### Nonlinear and reciprocal-only benchmarks (not Anymatrix unitary generators)

The [machine-readable eigenvalue-problem registry](../metadata/eigenvalue_problems.json)
also tracks seven additional cases:

- **Transcendental quantum graphs:** `I-diag(exp(ik*lengths))*S`
  is unitary-derived for real `k`; the physical spectral parameters
  are roots of a secular determinant, not the spectrum of one fixed
  matrix. For general incommensurate lengths no finite polynomial
  linearization follows. [58, 68]
- **Blaschke–Potapov rational inner functions:** unitary on `|z|=1`,
  but the nonlinear secular function itself is not a fixed unitary
  matrix. Poles and normalization require a rational NEP interface. [61, 62]
- **Scalar paraorthogonal companions:** use the explicit finite CMV
  unitary realization; do not call a generic Frobenius companion
  unitary merely because its zeros lie on the circle. [5, 66]
- **NLEVP `gen_tpal2` and `gen_tantipal2`:** prescribed circle
  roots of real palindromic/antipalindromic quadratic polynomials.
  Public reference code is in `ftisseur/nlevp`; these are **not**
  automatically strict unitary pencils. [63]
- **Surface acoustic waves:** palindromic polynomial formulations
  include propagating and evanescent modes; reciprocal pairing is
  weaker than circle confinement. [65]
- **Rail vibration `railtrack2`:** sparse, low-rank, large
  palindromic polynomial from the NLEVP collection, with potentially
  zero/infinite eigenstructure and roots away from the circle. [63, 64]

The library deliberately does **not** return these external nonlinear
functions as numeric unitary matrices. Their representation and reference
links remain accessible through the registry.

## L. Unimodular random matrices and quantum-information test families (November 2026)

The seven functions in this section are explicit *square unitary* test generators.
A Choi matrix of a quantum channel is positive semidefinite, not usually
unitary; a Stinespring isometry is rectangular unless the environment is
trivial. Such objects provide context but are not added to this matrix catalog.
The cited application papers motivate or explicitly contain the constructions;
do not attribute the proposed benchmark parameterizations to them without
checking their exact definitions.

### Haar special-unitary ensemble — [`haar_special_unitary.m`](../haar_special_unitary.m)

**Construction.** Draw a Haar matrix `Q` from U(n) by phase-corrected complex
Gaussian QR; set `U=diag(conj(det(Q)),1,...,1)*Q`. Haar measure pushes forward to
SU(n). The matrix is dense, of order n, and has determinant one. A seed is
optional; its use does not alter the caller's RNG state. For n=1, U=1.

**Exact statistical reference.** For n>=2, the normalized *ensemble-mean*
one-point eigenphase density is
`rho_n(theta)=[1-(-1)^n*(2/n)*cos(n*theta)]/(2*pi)`. This is not the empirical
density of one draw. For n=2, it is `sin(theta)^2/pi`. Compare a large sample
average with this function; a single matrix has atomic eigenphases. The
terminal `'spectral'` mode returns this normalized mean density and its
exact Fourier moments, without drawing a random matrix. [84,85]

### Unimodular circular orthogonal ensemble — [`haar_unimodular_coe.m`](../haar_unimodular_coe.m)

**Construction.** With `Q` Haar in SU(n), return `U=Q*Q.'` (*ordinary
transpose*, not conjugate transpose). It is symmetric, complex unitary,
determinant one, and distributed on the compact symmetric space SU(n)/SO(n).
This is **not** a real-orthogonal Haar sample.

**Reference density.** For n>=3, write
`a_n=sqrt(pi)*factorial(n-1)/(2^(n-1)*gamma(n/2+3/2)*gamma(n/2+1))`.
The normalized ensemble-mean density is
`rho(theta)=[1-(-1)^n*a_n*cos(n*theta)]/(2*pi)`.
Compute `a_n` via log-gamma to avoid overflow for large n. **Important exception:**
for n=2, the normalized density is exactly `abs(sin(theta))/4`, *not* the
preceding cosine expression (Nishigaki, footnote 1). For n=1 the single
eigenphase is the atom at zero. These cases are handled explicitly by the
`'spectral'` interface. [84]

### Unimodular circular symplectic ensemble — [`haar_unimodular_cse.m`](../haar_unimodular_cse.m)

**Construction.** Draw Haar `Q` in SU(2n), set
`J=[0 I;-I 0]`, and return `U=Q*J*Q.'*J'`.
The order is **2n**, the matrix is unitary, determinant one, and self-dual under
`A^D=J*A.'*J'`. Each eigenphase occurs twice (Kramers degeneracy);
this is the quotient SU(2n)/Sp(2n), not Haar sampling of Sp(n).

**Reference density.** For n>=2 the normalized ensemble-mean density of its
n distinct doubly repeated phases is
`[1-(-1)^n*b_n*cos(n*theta)+c_n*cos(2*n*theta)]/(2*pi)`, where
`b_n=(2n)!!/((2n-1)!!*n)`, `c_n=2/((2n-1)*n)`.
For n=1 this ensemble is the deterministic identity of order 2; use its point
mass at phase zero rather than the generic continuous formula. All three
unimodular generators provide these laws in the `'spectral'` output; the
`'ensemble'` field is never confused with one random matrix's finite atoms. [84]

### Maximally entangled rank-one phase reflection — [`entangled_reflection.m`](../entangled_reflection.m)

**Construction.** For physical local dimension d, let
`v=vec(I_d)/sqrt(d)` and
`U=I_(d^2)+(exp(i*phi)-1)*v*v'`. At the default phase `phi=pi`,
`U=I-2vv'` is the explicit local-certification counterexample of [77],
Eq. (15). Order d^2; sparse identity plus an O(d^2)-nonzero rank-one term.
The `operator` mode applies the rank-one update in O(d^2) per RHS.
The `'spectral'` output includes exact finite eigenvalues, moments, and
atomic weights (including the degenerate phase-zero case).

**Exact spectrum.** The exceptional eigenvalue is `exp(i*phi)` once, and 1
occurs d^2-1 times. At phi=pi and d>=2 the ordinary numerical range contains
zero, whereas the product-state numerical range is the real interval
`[1-2/d,1]`. This separation is informative for quantum certification but
does not change the empirical trace spectral measure.

### Lossless optical interferometer meshes — [`optical_mesh.m`](../optical_mesh.m)

**Order:** $N=n$ optical modes. The Reck triangular geometry and the
Clements rectangular geometry use n(n-1)/2 consecutive two-port,
nearest-neighbor lossless couplers. Their adjacency schedules are
exposed by selecting `layout='reck'` or `layout='clements'`.
For a selected pair (j,j+1), the core is

```math
G_j(\theta,\phi)=
\begin{bmatrix}
e^{i\phi}\cos\theta&-\sin\theta\\
e^{i\phi}\sin\theta&\cos\theta
\end{bmatrix},
\qquad 0\leq\theta\leq\pi/2.
```

For parameters in application order, the complete transfer matrix
is `U=D G_M ... G_1`, where `D=diag(exp(i*external))`
is the output phase-shifter bank. The Reck implementation applies
sequential descending adjacent pairs in triangular sweeps; the
Clements implementation uses n alternating disjoint odd/even
nearest-neighbor columns. Every factor is strictly unitary.
These topologies describe *ideal lossless* finite-mode models; no
claims about scattering loss, experimental error or gain are made.

The `'operator'` mode retains the O(N^2) coupler parameters and
multiplies vectors or multiple right-hand sides in O(N^2) arithmetic
per right-hand side, with the exact reversed core order on the adjoint.
At `theta=pi/4` each coupler is 50:50, and `theta=pi/2` gives
phase-signed routing/permutation matrices. Arbitrary input parameters
give exactly unitary products but an **unknown general analytical
trace eigenphase law**. With all couplings turned off, U is diagonal
and the `'spectral'` mode supplies its exact phase multiset and moments.
The spectrum of a local beamsplitter must not be mistaken for that
of the full network, nor should lossy networks be called unitary.

This provides photonic-core structure absent from a generic
`haar_unitary` or `complex_hadamard` matrix; no stable
decomposition of an arbitrary target unitary is claimed by the
parameterized forward model. References:
Reck et al. (1994) [118], Clements et al. (2016) [119].

### Bosonic symmetric-power lift — [`bosonic_symmetric_power.m`](../bosonic_symmetric_power.m)

For a one-particle V in U(d) and p indistinguishable bosons in d
modes, Sym^p(V) is a unitary matrix of order N=binomial(d+p-1,p)
on the normalized occupation-number basis (n_1,...,n_d), with
nonnegative occupancies summing to p. The mode-occupation ordering
is the lexicographic order of the underlying stars-and-bars separators.
For output occupation m and input n, repeat the rows of V according
to m and its columns according to n. The transition amplitude is
the permanent of this p-by-p repeated-index matrix, divided by
sqrt(prod(m_j!)*prod(n_j!)). These factorial normalizations follow
from normalized bosonic creation operators and are essential for
strict unitarity. P=0 gives [1], p=1 gives V, and d=1 gives V^p.

The exact N trace-eigenphases are all products of p one-particle
eigenvalues with **repetition allowed**, counted by occupation
multiplicity. For integer k the normalized trace moment is the
complete homogeneous symmetric polynomial
h_p(z_1^k,...,z_d^k)/N, evaluated by an O(d*p) recurrence.
Thus a large many-body spectral reference is available without
forming the N-by-N matrix. The explicit Ryser-permanent method
is intentionally restricted to N<=100 and p<=7, and there is
**no claimed generic fast matrix-free operator**.

This is a proper symmetric subrepresentation of V^(tensor p);
it differs from the fermionic exterior power, in which occupancy
is limited to 0 or 1 and determinants replace permanents. This
family provides exact phase-multiplicity tests for linear optics
and bosonic many-body problems.

**Reference [117]:** S. Scheel, *Permanents in linear optical
networks* (2004), [arXiv](https://arxiv.org/abs/quant-ph/0406127).

### Fermionic exterior-power lift — [`fermionic_exterior_power.m`](../fermionic_exterior_power.m)

Given V in U(d) and a number p of conserved fermions, 0<=p<=d,
the exterior-power representation U=wedge^p V acts on the antisymmetric
p-particle Hilbert space of dimension N=binomial(d,p). In the canonical
increasing occupied-subset basis (I,J), the matrix element is the
determinant of the p-by-p minor V(I,J). Antisymmetrization of the full
tensor action proves U is strictly unitary, and Cauchy--Binet gives
wedge^p(VW)=(wedge^p V)(wedge^p W). P=0 gives [1], p=1 gives V,
and p=d gives [det(V)].

Its exact finite trace-eigenphase multiset consists of products of
p **distinct** one-particle eigenvalues, with multiplicities counted
by p-element subsets. For integer k, the normalized trace moment is
the elementary symmetric polynomial e_p(z_1^k,...,z_d^k)/N, which
is computable in O(d*p) per k without constructing the large matrix.
The `spectral` mode offers this exact reference, while the explicit
minor construction is deliberately capped at N<=256. There is
**no claimed fast matrix-free operator** for arbitrary V. This is a
proper invariant subrepresentation, not `haar_tensor_power` on
dimension d^p, and it excludes number-nonconserving Bogoliubov
transformations. Application: fermionic optics and Slater determinants.

S. Bravyi, "Lagrangian representation for fermionic linear optics,"
*Quantum Information and Computation* **5** (2005), 216–238 [116].

### Repeated Haar tensor power — [`haar_tensor_power.m`](../haar_tensor_power.m)

Given a Haar V in U(d), form `U=kron(V,...,V)` with t identical factors on
the full tensor-product space of order d^t. The exact eigenvalues are products
of t eigenvalues of V (with repetitions), and
`tr(U^k)/d^t=(tr(V^k)/d)^t`. The repeated base unitary induces special
algebraic degeneracies. This is **not** `Sym^t(V)` or an exterior power.
The default function explicitly materializes matrices only through order 4096.
The optional `'operator'` output uses O(d^t) state-vector memory and repeated
small-factor actions, allowing larger dimensions without a dense Kronecker
matrix. With a fixed random seed, `'spectral'` computes eigenvalues of the
small d-by-d base matrix and returns exact full-space phase products/moments;
it does not claim a nonuniform ensemble-mean density. The repeated-unitary
setting is present in [81].

### Independent-block Haar tensor power — [`haar_block_tensor_power.m`](../haar_block_tensor_power.m)

Use independent Haar V1,V2 in U(d) to form
`V1^(tensor t1) kron V2^(tensor t2)`. Its order is d^(t1+t2), with exact
normalized moments `[tr(V1^k)/d]^t1*[tr(V2^k)/d]^t2`. Same-seed output is
reproducible with independent sequential draws. Contrast the *same* U on
both blocks in `haar_tensor_power`, as in the comparison of Haar-averaged
measurement models [81]. Explicit materialization is capped at order 4096. The
`'operator'` mode applies only small factors, while `'spectral'` (requiring
a seed) computes full-space spectral data from the two small matrices.

### Tensor-factor permutation — [`tensor_permutation.m`](../tensor_permutation.m)

For a permutation pi of t factors, `P_pi` maps
`e_(j1) kron ... kron e_(jt)` to
`e_(j_pi(1)) kron ... kron e_(j_pi(t))` on order d^t.
It is a sparse **real orthogonal** matrix. Write c(pi^k) for the number of
cycles of pi^k. Then `tr(P_pi^k)=d^c(pi^k)` exactly. At t=2 with pi=[2 1],
the swap operator has +1 multiplicity d(d+1)/2 and -1 multiplicity
d(d-1)/2. This is a concrete representation-theoretic specialization of
generic permutations relevant to twirling and Schur--Weyl duality. The
`'operator'` mode reshapes and permutes each RHS with no sparse matrix;
`'spectral'` provides exact normalized moments `d^(c(pi^k)-t)` at large
orders, as well as an optional finite eigenvalue list derived from basis
orbits rather than a numerical eigensolver. [79,81]

### Finite Weyl--Heisenberg 1-design (already implemented)

No duplicate generator is introduced. Existing
[`time_frequency_shift.m`](../time_frequency_shift.m) yields
`W_(a,b)=Z^a X^b`; enumerate a,b=0,...,d-1 to obtain an exact unitary
1-design (global phases are immaterial for the twirl):
`(1/d^2)*sum_(a,b) W_(a,b)*A*W_(a,b)' = tr(A)/d*I`.
This is a useful exact averaging test, not a new family [79,80].

## M. Nonlinear unitary matrix flows

### Lohe coupled matrix oscillators — [`lohe_matrix_flow.m`](../lohe_matrix_flow.m)

**Definition and construction.** The returned unitary has order $N=d$. For `NOSC=p` coupled oscillators of matrix
order `d`, evolve `U_j(t)` on `U(d)` according to

```math
\dot U_j=\left[-iH_j+\frac{\kappa}{2}
\left(\bar U U_j^*-U_j\bar U^*\right)\right]U_j,
\qquad \bar U=\frac1p\sum_{\ell=1}^p U_\ell,
\qquad H_j=H_j^*.
```

The bracketed generator is skew-Hermitian, so each exact trajectory remains
unitary. The sample matrices are generally dense, with time-dependent
eigenphases. Coupling makes the evolution **nonlinear**; there is no single
fixed unitary propagator that simultaneously represents all initial conditions.
The family is a matrix-valued, non-Abelian generalization of the Kuramoto
synchronization model [86,87].

**Numerical implementation.** We use a *synchronous exponential midpoint*
method: from all current `U_j` compute all skew-Hermitian generators,
predict each half-step using `expm(h*A_j/2)*U_j`, then recompute the common
mean and apply `expm(h*A_j_half)*U_j` from the original states. This is
second-order accurate for smooth trajectories and keeps each step on the
unitary group up to floating-point roundoff. It is **not** an exact
closed-form solution. Distinct geometric integrators for Lohe systems have
also been studied [87]. Complexity is O(`STEPS*NOSC*d^3`), and the
implementation caps `d<=128`, `NOSC<=64`. Default Hermitian frequencies
are deterministic; default initial unitary matrices are independently
seeded Haar draws and do not consume caller random state.

**Parameters and outputs.** `lohe_matrix_flow(d,NOSC,kappa,t,H,U0,which,steps)`
returns oscillator `which` at time `t`. Custom `H` and `U0` are
`d`-by-`d`-by-`NOSC` arrays of Hermitian and unitary matrices,
respectively. Terminal `'ensemble'` returns all evolved oscillator states
as a single 3-D array instead of a single selected matrix.
Terminal `'spectral'` returns an **unknown** general finite counting
measure, not a fabricated analytic spectral density.

**Special cases and useful invariants.**
- At `kappa=0`, `U_j(t)=exp(-i*t*H_j)*U_j(0)` exactly.
- With identical `H_j` and synchronized initial data, all states remain
  synchronized and evolve under the common Hermitian propagator.
- At `d=1`, the model reduces to scalar Kuramoto phases; the symmetric
  two-oscillator zero-frequency case admits an exact solution for benchmarking
  global time-stepping error.
- If `H_j=0` and initial `U_j(0)` are real orthogonal, the flow remains
  real orthogonal. Right multiplication by any common fixed unitary `R`
  commutes with the full dynamics: `U_j(t;U0*R)=U_j(t;U0)*R`.
- **Do not assume `SU(d)` invariance:** the Lohe coupling need not be
  traceless, and individual determinants can evolve. A unitary matrix
  family does not automatically lie in a special-unitary subgroup.

**Spectral benchmarking.** Compare the time-dependent empirical measures
`(1/d)*sum_k delta_(lambda_k(U_j(t)))`, and check unitarity, spectral
continuity, and numerical sensitivity to step size. No generic exact
eigenphase law, bulk density or OPUC measure is claimed.


## N. Citation-led rational and butterfly constructions (November 2026)

### Pole-parameterized rational OPUC finite unitary model — [`rational_opuc.m`](../rational_opuc.m)

**Construction.** Supply exterior poles `poles` of length `n-1`
(`abs(poles)>1` or `Inf`), strict-disk recurrence coefficients
`schur` of length `n-1`, and a unimodular closing parameter `beta`.
Set the diagonal normal strict contraction
`A=diag([0;1./conj(poles)])` with `1/Inf=0`, and
`D=diag(sqrt(1-abs(diag(A)).^2))`. If `H` is the finite unitary
Hessenberg core matrix from the recurrence, define

```math
U=A+D H (I+A^*H)^{-1}D.
```

This is an operator Möbius transform (equivalent to the normalized
linear-fractional representation for a normal strict contraction).
It is exactly unitary when `H` is unitary and `norm(A)<1`. All-infinite
poles recover the polynomial Hessenberg model. The poles are
**recurrence parameters**, not eigenvalues or singularities of `U` as a
fixed numeric matrix. The finite boundary choice determines a rational
paraorthogonal approximation; generic equal-weight eigenphase laws
are unknown. With poles near the unit circle the linear solve may
become ill-conditioned. This is not claimed to be a universal finite
rational-Krylov pencil or to preserve an independently fixed
polynomial OPUC measure.

**Applications:** orthogonal rational functions and rational quadrature,
pole-sensitive eigenphase tests, structured unitary Möbius models.
The explicit matrix uses an N-by-N solve; no matrix-free option is
claimed. **References:** [88], [5], and the 2017 SIGMA ORF survey.

### Multichannel paraunitary FIR and finite periodic operators — [`paraunitary_fir.m`](../paraunitary_fir.m)

**Size and conventions.** L temporal polyphase blocks with d channels
per block yield a finite dL-by-dL unitary. The M unit-norm
projection directions v_j and constant d-by-d unitary U0 define
a degree-M matrix Laurent polynomial

```math
E(z)=U_0 S_M(z)\cdots S_1(z),\qquad
S_j(z)=I_d+(z^{-1}-1)v_j v_j^*,\qquad |z|=1.
```

Every S_j(z) acts as z^(-1) on span(v_j) and as identity
on its orthogonal complement. Consequently E(z) is
**paraunitary** for every point of the unit circle, and
its determinant has winding number -M. `'frequency'`
returns one small d-by-d E(exp(i*omega)), and
`'polyphase'` returns the exact Laurent coefficient
array `E(z)=sum_r E_r z^(-r)` and winding -M.
The tap autocorrelations satisfy
`sum_r E_r^* E_(r+k) = delta_(k,0) I_d`.

The actual finite periodic critically sampled polyphase model
uses circular delay `v_j^*x_(t-1)` of the projected channel
amplitude in each stage, then applies U0 at every time block.
The action is exactly unitary under **circular boundary closure**
for all positive L, even L=1 (where the delay is the identity).
The `'operator'` mode performs O(dLM) arithmetic per right-hand
side using rank-one local updates and applies the Hermitian
adjoint with reversed factors and inverse circular shifts.

The finite operator is block circulant in time: the
normalized temporal DFT reduces it to the direct sum
of the L **small frequency responses**

```math
E(e^{2\pi i k/L}),\quad k=0,\ldots,L-1.
```

Hence its full finite trace-eigenphase measure is exactly
the union of the eigenvalues of these L small d-by-d blocks,
with multiplicities. `'spectral'` derives eigenvalues and
moments from that reduction without diagonalizing the dL matrix.
The eigenphases of a **single response** do *not* give the
finite operator spectrum, and the continuous response's
frequency curves are not an independently defined trace law.
The sampled blocks' small eigenvalues are computed numerically,
not asserted to have general closed forms.

This is distinct from the rational pointwise
`blaschke_potapov` and the predefined two-channel
`daubechies_wavelet` analysis. The factored construction
gives exact alias cancellation and perfect reconstruction
for ideal, lossless **critically sampled** polyphase filters.

**Applications:** multirate signal processing, filter banks,
lossless lattice sections, block-circulant spectral density
benchmarks. See Daubechies [15] and Vaidyanathan [120].

### Rank-one Blaschke–Potapov products — [`blaschke_potapov.m`](../blaschke_potapov.m)

**Construction.** For strict-disk poles `a_j`, normalized vectors `v_j`,
frequency `z=exp(1i*omega)`, and an initial unitary `U0`, define

```math
B(z)=U_0\prod_{j=1}^m
 \left(I+\left(\frac{z-a_j}{1-\overline{a_j}z}-1\right)v_jv_j^*\right).
```

Each factor has eigenvalues 1 on `v_j^perp` and the unimodular
Blaschke factor on `span(v_j)`, so `B(z)` is exactly unitary
for `abs(z)=1`. Projectors need not commute. The determinant
satisfies `det(B(z))=det(U0)*prod(b_j(z))`. The generic
matrix eigenphases need not be the individual factor phases.

**Applications:** lossless multichannel scattering, rational all-pass
filter banks, phase winding, noncommuting rank-one updates.
The `operator` mode applies factors in O(n*m) plus the cost of `U0`;
the default identity `U0` permits a genuinely low-storage action.
These are **pointwise** unitaries of order n, not automatically the
finite analysis matrix of a sampled filter bank. **Reference:** [61].

### Recursive orthogonal butterflies — [`orthogonal_butterfly.m`](../orthogonal_butterfly.m)

**Construction.** For `n=2^k`, prescribe `n-1` real angles indexed
in preorder across a binary tree. With recursively generated
`B_left, B_right` of order n/2,

```math
B_n=
\begin{bmatrix}
 cI&sI\\
 -sI&cI
\end{bmatrix}
\begin{bmatrix}
 B_{\rm left}&0\\
 0&B_{\rm right}
\end{bmatrix},
\quad c=\cos\theta,\ s=\sin\theta,\quad B_1=1.
```

Every node is orthogonal with positive determinant, hence `B_n` is
real orthogonal with determinant 1, a factored O(N log N) action,
and conjugation-symmetric eigenphases. A particular distribution of
angles (e.g., independent uniform) defines its random ensemble;
it is not generally Haar on all of `SO(n)`. The matrix-free
operator supports forward and adjoint multiple-RHS actions.

**Applications:** large-scale orthogonal eigensolver benchmarks,
spectral-density estimation with exact real symmetry, randomized
linear algebra and structured fast transformations. **Reference:** [89].

## O. Complex Hadamard and structured Hessenberg benchmark families (November 2026)

### Butson Hadamard phase tables — [`butson_hadamard.m`](../butson_hadamard.m)

A Butson matrix BH(n,q), with $N=n$, is an unnormalized $N\times N$ matrix
with entries among the q-th roots of unity, whose Gram matrix
is exactly nI. The normalized matrix U=H/sqrt(n) is therefore
strictly unitary. The generator constructs reproducible **integer
exponent tables** E with H(j,k)=exp(2*pi*i*E(j,k)/q), and exposes
them with `'exponents'` output, separately from its ordinary
normalized unitary matrix.

**Tao's isolated BH(6,3).** The supported 6-by-6 construction
has q=3 and an explicitly tabulated symmetric integer phase
table (entries 0,1,2). Every pair of distinct rows has exactly
two phase differences in each residue class modulo 3, so
the complex Hadamard Gram identity follows from
1+omega+omega^2=0 for omega=exp(2*pi*i/3).
This is the isolated S6 matrix cataloged by Tadej and
Zyczkowski [92], and **not** another Diţă deformation or Fourier
block matrix.

**Paley type-I real BH(p+1,2).** For prime p congruent to
3 (mod 4), set Q(i,j)=chi(i-j) for i~=j (quadratic
character modulo p) and Q(i,i)=0. The bordered sign matrix

```math
H=\begin{bmatrix}1&\boldsymbol{1}^T\\
\boldsymbol{1}&Q-I_p\end{bmatrix}
```

satisfies H H^T=(p+1)I, since Q is skew-symmetric and
Q Q^T=pI-J. The construction gives real orthogonal
normalized matrices of orders 4, 8, 12, 20, 24 and
other admissible prime-plus-one sizes, including orders
not covered by tensor powers of the Sylvester transform.
The sign/exponent table is exact, not a floating-point
optimization output. Source: Paley [121].

**Spectra.** The root-of-unity restriction on entries
does **not** imply the unitary eigenvalues are q-th
roots of unity. Generic finite trace eigenphase laws
are accordingly recorded as unknown. Left/right
Hadamard equivalence operations (row/column phases
and permutations) do not generally give a unitary
similarity and need not preserve eigenvalues.
No general fast implicit action is asserted.

### Diţă-deformed Fourier-block complex Hadamard — [`complex_hadamard.m`](../complex_hadamard.m)

**Construction:** For n=m r with normalized Fourier matrices F_m and
F_r, and r-by-m real phase angles t_{a,j}, set D_j=diag(exp(i*t(:,j))).
Define

```math
U=(F_m\otimes I_r)\,\mathrm{diag}(D_1F_r,\ldots,D_mF_r).
```

The product is unitary and *every entry* of U has modulus 1/sqrt(n).
Thus it is a normalized complex Hadamard matrix. Distinct angle choices
yield Diţă-type deformations, generally not simply Fourier matrices.
The default choice is deterministic and nontrivial when both factors
have order at least 2. For a prime n, the default m=1 recovers F_n.
The `operator` mode applies two FFTs and diagonal phases with
O(N log N) complexity, including multiple RHS and adjoints. Only
the undeformed Fourier corner provides an asserted analytic finite
eigenphase law. A complex Hadamard matrix need not be an involution.

**Application:** constant-modulus unitary transforms, structured
randomized linear algebra, signal coding, mutually unbiased bases,
and dense normal eigenvalue tests with highly non-generic entries.

**References:** [91], [92]. Equimodularity is a property of matrix
entries, not an assertion about eigenphase uniformity.

### Multiplicative rank-one unitary perturbations — [`unitary_rank_one_update.m`](../unitary_rank_one_update.m)

Let H be a unitary upper Hessenberg matrix encoded by Schur/Verblunsky
coefficients. For a unit vector v and real phase phi, set

```math
U=H\left(I+(e^{i\phi}-1)vv^*\right).
```

This preserves unitarity and ensures rank(U-H)<=1 and
det(U)=exp(i*phi)*det(H). By default v=e_n, so only the last column
of H changes and U is still upper Hessenberg. Arbitrary v gives
another valid unitary, but usually loses the Hessenberg constraint.
This is multiplicative *unitary-to-unitary* rank-one perturbation,
not the additive nonnormal companion-matrix model H+uv^*.
The latter belongs to a different matrix class.

When all interior coefficients vanish and v=e_n, the exact finite
eigenvalues are n-th roots of exp(i*phi)*conj(beta), for boundary beta.
For generic coefficients and vectors, finite trace eigenphases remain
unknown. The matrix-free kernel applies the rank-one projector and
unitary core product in O(N).

**Application:** eigenphase interlacing and motion, rank-sensitive
unitary eigensolver benchmarks, circular-ensemble perturbation theory,
and exact-cycle regression problems. **References:** [90], [5].

### ETNA core-factor unitary Hessenberg tests — [`etna_hessenberg.m`](../etna_hessenberg.m)

This generator exposes the Type I, II and III constructions discussed
in Section 7 of Aurentz–Mach–Vandebril–Watkins [8]:

- **Type I:** Each core has first column
  `[(-1)^(k-1)*m_k*p_k; sqrt(1-m_k^2)]`, with
  `m_k` sampled independently in (0,1) and `p_k` on the unit circle.
  The unitary second column is fixed by the repository core convention.
- **Type II:** Real core first columns `[cos(t_k);sin(t_k)]`,
  `t_k` sampled from (0,pi); the chosen orthogonal completion follows
  the same core convention. The result is real orthogonal Hessenberg.
- **Type III:** Four identical Type-I *interior* core sequences,
  coupled at three block boundaries by a real core with small
  `rho=coupling` (default 0.001), and terminal phase -1. This
  provides an exactly fourfold block spectrum as coupling tends to
  zero, with clusters split by small coupling at finite n. The
  coupling phase is fixed by construction instead of being
  sampled as in the original paper.

The published Example 7.4 reported a difficult real instance
with sixteen eigenvalues clustered near -1. **The implementation
does not claim to reproduce that same sample or that specific
cluster.** Its Type III gives a reproducible family with the
near-reducible structure of the benchmark; a user-controlled seed
and coupling allow systematic parameter studies. Coupling values too
small to resolve in double-precision Schur parameters are rejected.

All three variants are specified by n-1 core parameters and a
unimodular closure, have O(N) implicit matvecs, and generate finite
unitary upper Hessenberg matrices. The seeded generator temporarily
changes the random stream and restores it. No general sample
eigenphase density is asserted.

**Application:** fast unitary QR/QZ, double shifts on real orthogonal
Hessenberg matrices, difficult clustered eigenvalues, and
structure-preserving spectral approximation. **Reference:** [8].


## P. Quantum cat, Krawtchouk and periodic Daubechies (November 2026)

### Hyperbolic quantum cat map — [`quantum_cat_map.m`](../quantum_cat_map.m)

The classical map is the integral hyperbolic torus automorphism

```math
S=\begin{bmatrix}a&1\\ad-1&d\end{bmatrix},
\qquad a,d\in 2\mathbb Z_{>0},\quad ad\ge 4,\quad \det S=1.
```

The specific parity-compatible quantization, at Hilbert-space
dimension n, is the chirp-Fourier-chirp propagator

```math
U_{jk}=\frac{1}{\sqrt n}\exp\left(
\frac{i\pi}{n}(d j^2-2jk+a k^2)\right),\qquad 0\le j,k<n.
```

The default classical map is [2,1;3,2], rather than the
standard textbook [2,1;1,1] with a different parity convention.
For even diagonal a and d, the kernel is periodic under j or k
shifts by n. The matrix is unitary: its phases are diagonal
unimodular chirps surrounding the normalized DFT. Its overall
metaplectic phase is fixed to 1, an explicit **projective phase
convention**. Other global phases rotate the finite eigenphases
even though they have identical adjoint dynamics.

The implementation does not quantize arbitrary SL(2,Z) maps,
nor does it assert a uniform finite eigenphase trace measure.
In operator mode, both forward and adjoint matvecs use one FFT
and O(N) chirp multiplications, hence O(N log N) complexity.

**Application:** quantized hyperbolic dynamics, quantum chaos,
semiclassical spectral statistics and unitary FFT-chirp benchmarks.
**Reference:** [93], Hannay and Berry (1980).

### Real orthogonal Krawtchouk polynomial transform — [`krawtchouk_transform.m`](../krawtchouk_transform.m)

For matrix order $N=n=m+1$ and $0<p<1$, define a real symmetric Jacobi matrix J,
with x=0,...,m,

```math
J_{xx}=p(m-x)+(1-p)x,\qquad
J_{x,x+1}=J_{x+1,x}=
-\sqrt{p(1-p)(x+1)(m-x)}.
```

The eigenvalues of J are exactly 0,...,m. Order its
real orthonormal eigenvectors by these eigenvalues and fix the
signs by requiring a positive first-row component. Their matrix
Q is the symmetric orthogonal Krawtchouk transform; its first
column contains the square roots of the binomial probability
weights. With the indicated sign/ordering convention,

```math
Q^T=Q,\quad Q^2=I,\quad
\sigma(Q)=\{1,-1\},
```

with eigenvalue multiplicities ceil(n/2) for +1 and
floor(n/2) for -1. Thus the **Jacobi nodes 0,...,m are
NOT the eigenvalues of Q**. The finite normalized trace
spectral measure of Q consists of two atoms at 0 and pi,
not a smooth distribution on the circle.

The first version forms Q by dense symmetric eigendecomposition
at O(N^3) work. It does not claim a fast Krawtchouk algorithm.
The first-component sign convention can become numerically
delicate at high polynomial degrees or extreme binomial
probability parameters; validation rejects unresolved zero
first components.

**Applications:** exact-moment unitary test problems, binomial
signal representations and discrete orthogonal polynomial transforms.
**Reference:** [94], Feinsilver and Schott.

### Periodic Daubechies wavelets and wavelet-packet trees — [`daubechies_wavelet.m`](../daubechies_wavelet.m)

**Order:** $N=n=2^m$, m>=0. The legacy four-argument function constructs
periodized orthogonal multilevel approximation-only wavelets.
ORDER=2,3,4 selects db2/db3/db4 (filter lengths 4,6,8), with db2
the unchanged default. The analysis low-pass filter h uses the documented
Daubechies reconstruction-tap convention; the high-pass filter is
g(k)=(-1)^k*h(L-1-k). For each block of even length m,

```math
a_j=\sum_{k=0}^{L-1}h_k x_{(2j+k)\bmod m},\qquad
d_j=\sum_{k=0}^{L-1}g_k x_{(2j+k)\bmod m}.
```

Optional fifth argument PACKET may be 'dwt' (the old transform), 'full'
(split all bands), 'none' (identity), or a vector of binary-heap labels
of internal split nodes. Label 1 is the root and children of node j
are 2j, 2j+1. Every split node must also have its ancestors split.
The coefficient ordering is a deterministic left-to-right dyadic
leaf ordering, and arbitrary incomplete trees are allowed.

Alternatively ORDER may be a finite real even-length low-pass QMF tap
vector with unit norm, sum sqrt(2) and vanishing correlation with its
positive even shifts. This condition supplies an orthogonal periodic
two-channel analysis even when the filter wraps repeatedly on m=2
or other short dyadic blocks. Each split stage is orthogonal, so all
packet-tree products are orthogonal. The operator and its adjoint
support complex, multi-RHS input with O(n L levels) arithmetic.
Reflection/symmetric boundary extensions are not this family.

The finite trace-eigenphase law of the resulting square matrix
is *not* its filter frequency response and is generally unknown;
the empty tree is I, with exact eigenphases all +1. These transforms
serve as contrasting unitary spectral-density benchmarks for Haar,
db2 and db4 packet trees, with exact finite orthogonality.

**Application and reference:** finite multiresolution and
wavelet-packet filter banks, Daubechies, *Ten Lectures on Wavelets*
[15], https://doi.org/10.1137/1.9781611970104.


## Q. Hahn/dual Hahn transforms and spherical parallel transport (November 2026)

### Normalized Hahn polynomial transform — [`hahn_transform.m`](../hahn_transform.m)

For matrix order $N=n$, let $m=n-1$ denote the maximal polynomial degree.
For real parameters alpha,beta>-1, define
Hahn polynomials, normalized by Q_k(0)=1, as

```math
Q_k(x;\alpha,\beta,m)
={}_3F_2\!\left(
\begin{matrix}-k,\; k+\alpha+\beta+1,\;-x\\
\alpha+1,\;-m\end{matrix};1\right),
\quad 0\le k,x\le m.
```

The explicit weight is positive on the finite integer grid:

```math
w(x)={\alpha+x\choose x}{\beta+m-x\choose m-x}.
```

A numerically reproducible approach builds the symmetric
self-adjoint difference matrix J from the coefficients

```math
B(x)=(x+\alpha+1)(x-m),\qquad
D(x)=x(x-\beta-m-1),
```

with

```math
J_{xx}=-B(x)-D(x),\quad
J_{x,x+1}=J_{x+1,x}=-\sqrt{B(x)D(x+1)}.
```

Here B(x),D(x) are nonpositive on their respective interior
indices, so J is real symmetric. Its spectrum is the
**quadratic degree grid**
`lambda_k=k(k+alpha+beta+1)`. The eigenvectors, arranged
with increasing lambda and a positive first component,
form the weighted Hahn evaluation matrix V, where rows
index x and columns index k. The generator returns
`U=V'`, so U has degree k in its rows and grid x in its
columns. Consequently

```math
U^T U=I,\qquad J U^T=U^T
\,\mathrm{diag}\!\bigl(k(k+\alpha+\beta+1)\bigr).
```

**Spectral distinction:** the quadratic values are eigenvalues
of the self-adjoint **difference operator J**. They are not,
in general, eigenvalues of the orthogonal matrix U. Therefore
this database does not assign the polynomial grid as
eigenphases of U. The small, deterministic `eig(J)`
construction is O(N^3); no fast Hahn transform is claimed.
Eigenvector first-component sign setting can become
ill-conditioned for extreme orders or parameters.

**Applications:** finite orthogonal polynomial analysis,
special-function quadrature, numerical tests with
nonuniform discrete polynomial weights. **Reference:** [95],
Koekoek, Lesky and Swarttouw; DLMF 18.22.

### Dual Hahn transform on quadratic spectral grid — [`dual_hahn_transform.m`](../dual_hahn_transform.m)

The dual Hahn polynomial
`R_x(lambda(k);alpha,beta,m)` is defined by

```math
R_x(\lambda(k);\alpha,\beta,m)
={}_3F_2\!\left(
\begin{matrix}-x,\;-k,\;k+\alpha+\beta+1\\
\alpha+1,\;-m\end{matrix};1\right),
\quad \lambda(k)=k(k+\alpha+\beta+1).
```

The Hahn/dual Hahn duality is exactly

```math
Q_k(x;\alpha,\beta,m)
=R_x(\lambda(k);\alpha,\beta,m).
```

Under the matching positive weight and orthonormal
normalization conventions, their **analysis transforms**
are adjoints:
`dual_hahn_transform(n,alpha,beta) =
hahn_transform(n,alpha,beta)'`.
This guarantees exact algebraic transpose duality without
duplicating an independent, potentially inconsistent weight
computation. The dual Hahn transform is real orthogonal,
and therefore its eigenvalues are complex conjugates of
those of Hahn's transform; its quadratic **polynomial grid**
is not itself a unit-circle eigenphase law.

The current implementation reuses the shared Hahn
Jacobi eigensolver (O(N^3)). **Applications:** quadratic-grid
discrete spectral expansions, exactly normalized dual bases,
polynomial transform pair benchmarks. **Reference:** [95];
DLMF 18.21, "Duality of Hahn and Dual Hahn."

### Spherical shortest-geodesic parallel transport — [`sphere_geodesic_transport.m`](../sphere_geodesic_transport.m)

Given two unit vectors x,y on the unit sphere S^(n-1) with
non-antipodal endpoints, let

```math
c=x^Ty,\quad
K=yx^T-xy^T,\quad
R=I+K+\frac{K^2}{1+c}.
```

R is the **unique minimal ambient rotation** in SO(n) that
maps x to y and is the identity on the orthogonal complement
of `span(x,y)`. Its restriction to
`T_x S^{n-1}\to T_y S^{n-1}` is exactly the Levi-Civita
parallel transport along the unique shortest geodesic,
which on tangent vectors v satisfying x^Tv=0 has the
closed form

```math
P_{x\to y}v
=v-\frac{y^Tv}{1+x^Ty}(x+y).
```

The **tangent-space map is not an endomorphism** until the
source and destination tangent frames are identified.
This generator stores R, its canonical square ambient
orthogonal lift, not an arbitrary frame-coordinate
identification. The antipodal case is excluded because
no unique shortest great circle joins x to -x, and near
antipodal points are rejected to avoid division by
an ill-conditioned `1+x^Ty`.

Writing `theta=acos(x^Ty)`, the finite exact eigenvalues are

```math
e^{i\theta},\quad e^{-i\theta},\quad
\underbrace{1,\ldots,1}_{n-2},\qquad
m_k=\frac{n-2+2\cos(k\theta)}n.
```

These give an **exact discrete trace spectral measure**
with at most three distinct atoms. At x=y, R=I and the
measure is one atom at zero. The ambient rotation has
rank(R-I)<=2, so the implicit forward/adjoint kernel
applies in O(N) per right-hand side. This is an
application-derived real orthogonal matrix family with
an exact finite spectral reference, not a random rotation
or an arbitrary product of Givens matrices.

**Applications:** Riemannian optimization, spherical
statistics, geometric numerical integration, geodesic
transport of vectors and exact finite moment benchmarks.
**Reference:** [96], Absil, Mahony and Sepulchre.
The ambient formula is verified by the tangent-space
transport identity in the regression tests.

## R. Families reserved for future versions (not implemented)

- **Additional ETNA benchmark cases:** the exact Example 7.1 common-parameter
  and Example 7.5 inverse-eigenvalue presets are not yet implemented.
  Types I, II and a controlled Type III are already available in
  `etna_hessenberg`. [8]
- **Rational and isometric Arnoldi:** pole sequences, orthonormal rational
  basis conventions, unitary closures and appropriate generalized matrices/
  pencils. Not every finite compression is unitary. [9, 10]
- **Other OPUC measures:** positive Pastro and generalized $q$-weight
  specializations, inverse spectral finite atomic measures, Aleksandrov
  families, skew-shift parameters, and matrix-valued/block CMV. [5, 12]
- **Additional signal-processing matrices:** other Daubechies orders,
  wavelet packets, paraunitary polyphase banks and lossless scattering.
  Boundary and frequency conventions must be specified before admitting a
  finite matrix as unitary. [15]
- **Larger quantum models:** Szegedy graph walks, quantum-graph scattering,
  Floquet Ising models, kicked tops, cat maps and interferometer meshes.
  Different symmetry sectors and truncations cannot be conflated. [26]

### Graph connection-Laplacian unitary propagators — [`graph_connection_propagator.m`](../graph_connection_propagator.m)

**Order:** $N=nq$, with n graph vertices and q complex fiber
coordinates per vertex. For every oriented representative edge
(i,j), the supplied q-by-q unitary connection link R_ij
maps the **fiber at j to the fiber at i**. A nonnegative
undirected weight w_ij determines the Hermitian graph
connection Laplacian:

```math
(L_\nabla x)_i=\sum_{j\sim i}w_{ij}(x_i-R_{ij}x_j),\quad
R_{ji}=R_{ij}^*,\qquad
x^*L_\nabla x=\sum_{\{i,j\}}w_{ij}\|x_i-R_{ij}x_j\|^2\ge0.
```

The finite **unitary** matrix in this family is
`U(t)=expm(-1i*t*L_nabla)`. The sparse positive-semidefinite
Laplacian itself is **not** unitary. `'laplacian'` returns
its complete sparse matrix, oriented links, weights and
gauge convention. For moderate N, `'matrix'` computes the
dense true matrix exponential, which is unitary to numerical
rounding and has eigenphases `exp(-it*lambda_j(L_nabla))`.
For a generic supplied graph no closed-form analytical
trace phase distribution is asserted.

The default N-cycle has unit edge weights, identity interior
links, and a diagonal q-by-q closing holonomy with phases
`phi_a=0.33*a` for a=1,...,q. Each channel's exact
Laplacian eigenvalues and propagator eigenphases are

```math
\lambda_{a,k}=2-2\cos\left(\frac{2\pi k+\phi_a}{n}\right),
\qquad z_{a,k}=\exp(-it\lambda_{a,k}),\quad k=0,\ldots,n-1.
```

This default therefore has **exact finite trace spectral
data** with equal weights 1/(nq) over all channels and sites.
General non-Abelian graph links need not commute and
their finite eigenvalue laws remain unknown analytically.

Under a local unitary frame change G_i, the links transform as
`R'_{ij}=G_i R_{ij}G_j^*` and the full operator is conjugated
by D=diag(G_i): `L'=DLD^*`, `U'=DUD^*`. This is genuine
gauge covariance; no assumption is made that the user-supplied
links arise from Levi–Civita parallel transport on a mesh.

For large graphs, `'cayley_operator'` supplies a separately
**labeled approximation**, the s-fold implicit midpoint
rational map

```math
Q_s=\left[(I+itL_\nabla/(2s))^{-1}
(I-itL_\nabla/(2s))\right]^s.
```

This approximation is exactly unitary in exact arithmetic,
admits sparse cached-LU forward/adjoint actions and converges
to `U(t)` at second order in 1/s for fixed t.
It is **not identical to exp(-itL_nabla)**, and is not registered
as an exact matrix-free action for the primary family.
Generic polynomial/Krylov exponential actions are not falsely
certified as exactly unitary.

**Applications:** discrete vector diffusion, gauge-aware geometric
graphs, connection-aware spectral density estimation, and finite
unitary approximation benchmarks. Distinct from the
single-cycle block permutation `connection_cycle` and the
general block-monomial `bundle_transport` families.

References: Singer and Wu, *Vector diffusion maps and the
connection Laplacian* [123], and Crane, Desbrun and Schröder [98].

### Non-Abelian Berry/Wilczek–Zee holonomy from moving eigenspaces — [`berry_wilson_loop.m`](../berry_wilson_loop.m)

**Order:** $N=r$, the rank of the occupied subspace. The input is a
closed sequence of m>=2 orthonormal d-by-r frames Q_j, with the
final overlap closing Q_m to Q_1. Unlike a collection of arbitrary
gauge links, these frames can represent actual parameter-dependent
eigenspaces of a band Hamiltonian. The generator does **not**
claim to infer a Hamiltonian or its occupied eigenstates; the
user is responsible for supplying physically meaningful frames.

From the adjacent frame overlaps, form

```math
M_j=Q_j^*Q_{j+1},\quad
V_j=\mathrm{polar}(M_j)=U_jZ_j^*,\quad
M_j=U_j\Sigma_jZ_j^*,\quad Q_{m+1}=Q_1,\quad
W=V_1V_2\cdots V_m.
```

Although each M_j is generally **contractive** and its raw
product is not unitary, every V_j is unitary if M_j is
nonsingular, and therefore W is genuinely in U(r).
Singular or poorly conditioned overlaps are rejected instead of
using an arbitrarily chosen polar factor across a rank loss.

Gauge changes Q_j to Q_jG_j imply
M_j to G_j^*M_jG_(j+1) and V_j to G_j^*V_jG_(j+1), so

```math
W\mapsto G_1^*WG_1.
```

Thus the eigenphases, traces and determinant of a closed
holonomy are gauge invariant, while its individual entries
depend on the base-frame convention. The `'operator'` mode
reuses the ordered unitary-link kernel from
`lattice_wilson_loop`, applying the factors in reverse
order to vectors and their adjoints in forward order. This
avoids materializing the full W product and handles complex
multiple RHS.

The rank-one default is a spinor latitude at polar angle
theta=1.1. Its sampled polar-overlap Berry phase converges
quadratically under uniform refinement to
`exp(i*pi*(1-cos(theta)))` with the stated
`M_j=Q_j^*Q_(j+1)` orientation.
For r>=2 the deterministic default has genuinely
noncommuting link matrices; the `'spectral'` descriptor
does not invent a generic analytic non-Abelian phase law.
`'data'` explicitly distinguishes raw overlap products,
singular values, unitary polar links and sampled W.

Wilson eigenphases alone are not automatically bulk Chern
numbers or topological invariants; such conclusions require
a specified gapped band family and the relevant physical
symmetries or global parameter-space assumptions.
This **complex occupied-band** fiber is distinct from
sphere tangent-bundle parallel transport and from
user-supplied lattice-gauge links.

**Reference:** Wilczek and Zee [124].

## Geometry: geodesics, discrete connections, and holonomy

The shortest-arc spherical transport operator is already provided by [`sphere_geodesic_transport.m`](../sphere_geodesic_transport.m); its geometric interpretation is also discussed by Fiori [100]. The construction is not duplicated here.

### Grassmannian geodesic ambient lift — [`grassmann_geodesic.m`](../grassmann_geodesic.m)

**Order:** $N=n$; real orthogonal or complex unitary according to the supplied
orthonormal frame B=[Y,U] of width 2r. With tangent singular angles a_j,
the ambient lift is identity outside span(B) and on span(B) rotates by

```math
C=\begin{bmatrix}\cos(D)&-\sin(D)\\
                  \sin(D)&\cos(D)\end{bmatrix},
\qquad D=\mathrm{diag}(a_1,\ldots,a_r).
```

Its exact eigenvalues are e^{+/-i a_j} and 1 with multiplicity n-2r.
The subspace path Y(t)=Y cos(tD)+U sin(tD) is a Grassmann geodesic; this
**ambient lift is not intrinsic parallel transport of the tangent bundle**.
The low-rank action costs O(nr). See Edelman–Arias–Smith [97].

### Spherical geodesic-loop holonomy — [`geodesic_holonomy.m`](../geodesic_holonomy.m)

**Order:** $N=n\ge3$. Columns of POINTS specify vertices of a closed spherical
geodesic polygon. Multiply shortest-geodesic ambient rotations in path order.
The result fixes the base point and is orthogonal; its restriction to the
base tangent space is the Levi-Civita holonomy.

For a triangle with Gram matrix G of its three unit vertices, the unsigned
spherical excess in the induced oriented 3-plane is

```math
a=2\mathrm{atan2}\left(\sqrt{\det G},
1+G_{12}+G_{23}+G_{31}\right),
```

and the ambient eigenvalues are e^{+/-ia} and n-2 copies of 1.
For polygons with more than three vertices, this generator does not claim
an automatic analytic spectrum; the implicit action takes O(nm) for m legs.
See sphere parallel transport [100] and the classical Gauss–Bonnet theorem.

### Block-permutation bundle transport — [`bundle_transport.m`](../bundle_transport.m)

**Order:** $N=md$ for $m$ sites, a permutation pi of sites, and unitary
d-by-d link matrices R_j. Define U(e_j tensor v) = e_{pi(j)} tensor R_j v.
This is a sparse **unitary** block-permutation (real **orthogonal** for
orthogonal links). For each cycle (j_1,...,j_ell), its holonomy is the
ordered product H = R_{j_ell} ... R_{j_1}. The characteristic polynomial
on that cycle is det(z^ell I_d - H). Thus each eigenvalue of H produces
all ell-th roots, with exact multiplicities. For integer k, its trace
contribution is ell*trace(H^(k/ell)) if ell divides k, and zero
otherwise. Dividing by m*d gives the equal-weight trace moments.

The operator/adjoint apply link matrices and reverse edges in O(m*d^2)
work per right-hand side. A change of site frames G_j replaces a link
with G_{pi(j)} R_j G_j^* and conjugates U by the global block-diagonal
frame matrix. The old `connection_cycle` is the single-cycle d=2
planar-rotation specialization and is retained unchanged.

This is a proposed exact finite benchmark of **supplied** connection
links, not a numerical construction of a Levi-Civita connection from
a surface mesh. Connection background: Singer and Wu, "Vector diffusion
maps and the connection Laplacian," *Comm. Pure Appl. Math.* 65 (2012),
[DOI](https://doi.org/10.1002/cpa.21395).


### Cyclic discrete-connection transport — [`connection_cycle.m`](../connection_cycle.m)

**Order:** $N=2n$; exactly sparse real orthogonal block-permutation matrix.
Each leg carries a prescribed SO(2) connection rotation R(phi_j).
The holonomy is R(sum(phi_j)), and the exact eigenvalues are

```math
\exp\left(i\frac{2\pi k\pm\sum_j\phi_j}{n}\right),
\qquad k=0,\ldots,n-1.
```

This is an **original proposed benchmark** assembled from documented local
discrete-connection rotations; the generator accepts connection angles but
does not compute them from a triangulated mesh. Spectral density of this
particular cyclic model is a pair of interlaced root grids, not an arbitrary
mesh spectrum. Reference on discrete connection angles: Crane–Desbrun–
Schröder [98].

### Lie-group geodesic parallel transport — [`lie_geodesic_transport.m`](../lie_geodesic_transport.m)

**Order:** $N=n^2$ real orthogonal matrix representing transport on the real Lie
algebra u(n) of skew-Hermitian matrices. Along exp(tK), K=i diag(kappa),
the left-trivialized transport is exp(-t ad_K/2). In the real orthonormal
basis consisting of n diagonal and n(n-1)/2 off-diagonal generator pairs,
the result is a sparse direct sum of 2x2 rotations and n singleton 1's.
Exact eigenvalues are e^{+/-it(kappa_j-kappa_k)/2} for j<k, plus n ones.

```math
\frac{1}{n^2}\mathrm{tr}(Q^\ell)
=\frac{1}{n^2}\left|\sum_{j=1}^n
e^{i\ell t\kappa_j/2}\right|^2.
```

This is a canonical-coordinate transport matrix, not multiplication by
exp(tK) on the ambient n-dimensional complex space. The eigenphase-difference
formula connects geometric transport to unitary conjugation but represents
a different real action. See Lie-group parallel transport [99].

### Time-dependent noncommuting SO(n) flows — [`time_skew_magnus.m`](../time_skew_magnus.m)

**Order:** $N=n\ge3$. This benchmark treats the time-dependent nonquantum
flow Y'(t)=A(t)Y(t), with Y(0)=I and **real skew-symmetric**

```math
A(t)=A_0+\gamma\sin(\nu t)A_1,\qquad
[A_0,A_1]\ne0.
```

The default sparse A0 and A1 use first- and second-neighbor skew
links, respectively. For every n>=3 their commutator is nonzero.
Consequently there is no reduction to a single time-independent
matrix exponential for general gamma. This distinguishes the model
from `lie_geodesic_transport`, `polar_unitary`, and
the nonlinear matrix-valued `lohe_matrix_flow`.

For h=T/m and stage midpoint times t_j=(j-1/2)h, the **second-order
commutator-free exponential midpoint** scheme is

```math
Q_j=\exp(hA(t_j)),\qquad U_m=Q_m\cdots Q_1.
```

The second supported scheme is **implicit/Cayley midpoint**,

```math
Q_j=(I-hA(t_j)/2)^{-1}(I+hA(t_j)/2),\qquad
U_m=Q_m\cdots Q_1.
```

Both are second-order methods (not an asserted fourth-order modified
Cayley–Magnus scheme). In exact arithmetic every stage belongs to
SO(n), so both products are exactly orthogonal with determinant +1.
The Cayley operator caches sparse LU factorizations of I±hA/2,
applies each factor and reverse-order adjoint with sparse solves;
complex multi-RHS vectors are supported. The exponential variant
caches dense `expm` factors and does **not** claim exact sparse
exponential actions. Both methods are symmetric in step
reversal, but this does not imply that a nonautonomous flow
is autonomous or that generic SO(2k) stages preserve canonical J.
Indeed, the default n=4 example is **not** symplectic.

An independent high-accuracy ODE integration checks second-order
global convergence for genuinely noncommuting stages. The generic
finite trace **eigenphase law is unknown analytically**.
When gamma=0, stage matrices commute, so the `'spectral'`
mode obtains exact algebraic exponent/Cayley eigenphase maps of
the numerical eigenvalues of A0 (independent of forming U).
This is not a closed-form eigenphase law for the nonautonomous case.

**Applications:** geometric integration, nonautonomous rotation
systems, orthogonality-preserving time-ordered numerical benchmarks.
References: Celledoni, Marthinsen and Owren (2003) [122],
and Hairer, Lubich and Wanner [101].

## Geometric integration: positive-metric unit-circle pencils

### Symplectic-Euler oscillator steps within `geometric_oscillator.m`

The periodic oscillator generator supports two physical first-order
symplectic-Euler schemes: `euler_kick_drift` and
`euler_drift_kick`, with canonical state x=(q,p), symmetric
positive K and h^2*lambda_max(K)<4.

For kick--drift, the physical step and its quadratic invariant are

```math
M_{KD}=\begin{bmatrix}I-h^2K&hI\\-hK&I\end{bmatrix},
\qquad G_{KD}=\begin{bmatrix}K&-hK/2\\-hK/2&I\end{bmatrix}.
```

For drift--kick,

```math
M_{DK}=\begin{bmatrix}I&hI\\-hK&I-h^2K\end{bmatrix},
\qquad G_{DK}=\begin{bmatrix}K&hK/2\\hK/2&I\end{bmatrix}.
```

Both physical maps satisfy M^T J M=J and M^T G M=G.
The `'pencil'` mode returns A=GM, B=G, M, G, the
symmetric energy square root T=G^(1/2), and U=T*M*T^(-1).
The generalized eigenvalues lie exactly on the unit circle.

The Euler orderings differ from velocity Verlet in their physical
steps and cross-term metrics, although their exact phase pairs
coincide: exp(±2i asin(h*omega_j/2)). The default normalized
representative is real orthogonal. `'operator'` reuses the
O(N log N) FFT diagonalization of the periodic stiffness.
For a generic positive metric normalization, do **not** claim
standard-J symplecticity without additional structure.
See Hairer, Lubich and Wanner [101,102].

### Periodic coupled oscillators: Gauss–Legendre and Störmer–Verlet — [`geometric_oscillator.m`](../geometric_oscillator.m)

**Order:** $N=2n$ real orthogonal-symplectic, for $n$ mechanical
displacements. **Application:** vibrational lattice dynamics, periodic
semidiscretizations of wave equations and structure-preserving integration.
The physical positive-stiffness oscillator is

```math
q''+Kq=0,\qquad
K=\omega_0^2 I+\kappa(2I-P-P^T),\qquad
\omega_j^2=\omega_0^2+4\kappa\sin^2(\pi j/n),
```

with cyclic permutation $P$, $\omega_0>0$, $\kappa\geq0$ and
$j=0,\ldots,n-1$. The $n=1,2$ periodic bond multiplicities follow
the same algebraic definition. A time step $h>0$ produces a **physical**
$2n$-by-$2n$ phase-space map $M$ in $x=(q,p)$, which generally
is not Euclidean orthogonal. Nevertheless, it preserves a definite
quadratic invariant,

```math
M^TJM=J,\quad M^TGM=G,\quad G=G^T\succ0,\qquad
J=\begin{pmatrix}0&I\\-I&0\end{pmatrix}.
```

Define $T$ so that $G=T^TT$. The default Anymatrix output
$U=TMT^{-1}$ is **real orthogonal and symplectic**, and its spectrum
equals that of $M$. The chosen $T$ need not itself be symplectic.
The terminal `'pencil'` mode returns $A=GM$, $B=G$, together
with $M,G,T,U$. It obeys

```math
A^TB^{-1}A=B,\qquad
Ax=\lambda Bx\ \Longrightarrow\ |\lambda|=1.
```

**Crucial distinction:** $B^{-1}A=M$ is generally *not* Euclidean
orthogonal, and $AA^T=BB^T$ generally fails. Hence this is a
`circle_spectrum_only` positive-metric pencil, not one of the
repository's `strict_unitary_pencil` examples. No generalized pencil
is passed to Anymatrix by default: the returned matrix $U$ is
genuinely unitary.

**Midpoint and Gauss methods.** The methods `'midpoint'`,
`'gauss2'` and `'gauss3'` give the diagonal Padé stability
functions $R_s(z)=P_s(z)/P_s(-z)$ of orders 2, 4, 6, respectively.
For all $h>0$, they preserve the *original* quadratic energy,
so $G=\mathop{\mathrm{diag}}(K,I)$ and

```math
\sigma(U)=\left\{
R_s(-ih\omega_j),\ R_s(ih\omega_j):0\leq j<n
\right\}.
```

**Stable Störmer–Verlet.** With $h\max_j\omega_j<2$, velocity
Verlet preserves a *modified* quadratic energy

```math
G=\mathop{\mathrm{diag}}\left(K(I-h^2K/4),I\right)\succ0.
```

Its eigenvalues are $e^{\pm i\theta_j}$, with
$\theta_j=2\arcsin(h\omega_j/2)$. An independent terminal
`'polynomial'` mode provides the **real symmetric palindromic**
quadratic eigenvalue problem from the position recurrence,

```math
P(\lambda)=\lambda^2I+\lambda(h^2K-2I)+I,\qquad
P(\lambda)^T=\lambda^2P(\lambda^{-1}).
```

At the stability boundary the modified form loses strict positivity;
above it the Verlet multipliers can leave the circle, so the generator
rejects those step sizes. In particular, *symplecticity alone* never
guarantees unimodular eigenvalues.

**Spectral reference:** the finite multiset of $2n$ conjugate unit
phases is exact for every permitted method, returned without `eig`
by terminal `'spectral'` mode. The periodic chain diagonalizes in
the Fourier basis. The explicit real representative is generally
dense. The terminal `'operator'` mode already provides $O(N\log N)$
FFT-based forward/adjoint actions (including complex multiple
right-hand sides), using only $O(n)$ stored frequency samples. Its symmetric phase measure
must be distinguished from the one-sided $n$-point complex
representation of the positive frequencies. The terminal `'spectral'` descriptor conforms to the shared trace-measure
interface: `finite.eigenvalues()` and `finite.moments(k)` are callable,
while the limiting density is marked unknown.
The formulas and
normalization are specified benchmark specializations of established
geometric-integration theory, not named distributions originating
in [93, 94]. **References:** [93, 94].

### One-dimensional finite-element wave equation with consistent mass — [`geometric_fem_wave.m`](../geometric_fem_wave.m)

**Order:** $N=2n$ real orthogonal-symplectic, from a spatially
discretized elastic/acoustic wave model with $n$ interior nodes.
For the scalar model $q_{tt}-c^2q_{xx}+\omega_0^2q=0$ on $(0,L)$,
use continuous piecewise-linear finite elements on a uniform grid
with homogeneous Dirichlet endpoints, mesh width $d=L/(n+1)$,
wave speed $c>0$ and $\omega_0\geq0$. The consistent finite-element
mass and stiffness matrices are

```math
\mathrm{Mass}=\frac{d}{6}\mathop{\mathrm{tridiag}}(1,4,1),\qquad
K=\frac{c^2}{d}\mathop{\mathrm{tridiag}}(-1,2,-1)+\omega_0^2\mathrm{Mass},
```

with boundary truncation of the off-diagonals. These matrices are
symmetric positive definite and simultaneously diagonalized by the
orthonormal DST-I basis. The *exact finite generalized frequencies*
satisfy

```math
K v_j=\omega_j^2\mathrm{Mass}\,v_j,\quad
\theta_j=\frac{j\pi}{n+1},\quad
\omega_j^2=\omega_0^2+
\frac{6c^2(1-\cos\theta_j)}{d^2(2+\cos\theta_j)},
\quad j=1,\ldots,n.
```

The canonical physical Hamiltonian uses momentum $p=\mathrm{Mass}\,q_t$,

```math
H(q,p)=\frac12 q^TKq+\frac12 p^T\mathrm{Mass}^{-1}p.
```

Its midpoint, Gauss-4 and Gauss-6 steps preserve this original
quadratic energy and are symplectic. Stable velocity Verlet
($h\max_j\omega_j<2$) preserves the positive *modified* form

```math
G=\mathop{\mathrm{diag}}\left(
K-\frac{h^2}{4}K\mathrm{Mass}^{-1}K,\,
\mathrm{Mass}^{-1}\right).
```

The terminal `'pencil'` mode returns the *physical* phase-space
one-step matrix $M$, $G$, $T$, $U=TMT^{-1}$, and
$A=GM,B=G$. Its generalized eigenvalues lie on the unit circle,
but it is not generally a strict Euclidean-unitary pencil.
The default output is the genuinely orthogonal-symplectic $U$.

Most importantly, the **Verlet position recurrence retains the mass
matrix**. Terminal `'polynomial'` mode returns the real symmetric
palindromic quadratic

```math
P(z)=z^2\mathrm{Mass}
+z(h^2K-2\mathrm{Mass})+\mathrm{Mass}.
```

The spectrum consists of $e^{\pm i\theta_j^{(h)}}$, where
$\theta_j^{(h)}=2\arcsin(h\omega_j/2)$ for Verlet;
for the Gauss methods, the phases are the diagonal-Padé images
$R_s(\pm ih\omega_j)$. All eigenvalues are exactly characterized
at finite dimension. `'spectral'` returns the frequencies and
phase multiset without an eigensolve; `'operator'` supplies
$O(N\log N)$ forward and adjoint DST-I actions with $O(n)$
stored frequency samples, including multiple complex right-hand
sides. The explicit $U$ is generally dense. The standard terminal
`'spectral'` descriptor exposes `finite.eigenvalues()` and
`finite.moments(k)` and leaves the limiting density unasserted.

Unlike the periodic oscillator in `geometric_oscillator`, this family
has a **nonidentity consistent finite-element mass matrix**,
Dirichlet boundary conditions, sine rather than Fourier modes,
and a genuinely mass-weighted palindromic eigenvalue polynomial.
The particular uniform Dirichlet test family is a mathematical
specialization of standard finite-element and geometric-integrator
constructions, not a claim of a novel scheme in [93, 95].
**References:** [93, 95].

### Daubechies db3/db4 parameter extensions (same canonical generator)

The existing `daubechies_wavelet.m` also offers the 6-tap db3 and 8-tap db4
orthonormal extremal-phase filters as an optional fourth argument,
`daubechies_wavelet(n,levels,output,order)` with `order=2,3,4`.
The legacy three-argument call remains db2/D4. Periodic extension,
orthonormal QMF coefficients, and the multilevel forward/adjoint kernels
are shared rather than creating separate matrix generators.
The filter coefficients and original phase convention are documented by
Daubechies [15] and the [PyWavelets coefficient table](https://github.com/PyWavelets/pywt/blob/main/pywt/_extensions/c/wavelets_coeffs.template.h).
The Hannay–Berry cat-map and Krawtchouk references from the same research
batch remain [93] and [94], attached to their existing canonical families.

## Q. Quadrature-orthogonal transforms and spherical harmonic rotations (November 2026)

**Scope and caveat.** The following are six mathematically unitary or real
orthogonal finite matrices, not simply nonsingular polynomial coefficient
conversions. The three Gauss transform generators are constructed through a
dense Golub--Welsch symmetric eigensolver; this is *not* a claim that the
current implementation has the fast matvec of Townsend/Slevinsky algorithms.

### General Gaussian transform from Jacobi recurrences — [`gauss_recurrence_transform.m`](../gauss_recurrence_transform.m)

**Order:** $N=n$. Input is the real diagonal b_0,...,b_(n-1) and strictly
positive subdiagonal a_1,...,a_(n-1) of a finite Jacobi tridiagonal
J=diag(b)+diag(a,1)+diag(a,-1). This gives simple real quadrature nodes,
positive Gaussian weights and the orthogonal evaluation transform
Q=V', where columns of V are ascending, first-component-positive
orthonormal eigenvectors of J. Thus Q'*Q=I and Q*J*Q' is the
diagonal matrix of quadrature nodes.

The finite eigenphases of Q are **not** the Jacobi nodes; a generic
analytical trace eigenphase law is unknown. This general generator
reuses the project's `um_gauss_basis` Golub--Welsch routine, extending
Legendre, Jacobi, Hermite and Laguerre from built-in recurrence
coefficients to user-defined classical/nonclassical measures.
It is a dense O(N^3) setup and O(N^2) multiply, not a claimed fast
Slevinsky/Townsend transform. At n=1, Q=[1]. Nonpositive recurrence
subdiagonals are rejected. Numerical conditioning may limit very
large or tightly clustered recurrences.

**Applications:** arbitrary positive spectral measures and their
orthogonal polynomial quadrature transforms; reproducible geometric
and spectral testing. Golub and Welsch [104],
[DOI](https://doi.org/10.1090/S0025-5718-69-99647-1).

### Gauss--Jacobi transform — [`gauss_jacobi_transform.m`](../gauss_jacobi_transform.m)

**Order:** $N=n$. **Parameters:** alpha,beta > -1. The weight on [-1,1]
is `(1-x)^alpha (1+x)^beta`, normalized to probability measure.
The Gauss evaluation matrix is exactly real orthogonal, constructed
through the Jacobi recurrence and Golub--Welsch. The special cases
alpha=beta=0, -1/2, +1/2 give Legendre, Chebyshev I and Chebyshev II;
alpha=beta=lambda-1/2 gives Gegenbauer. The first Jacobi offdiagonal
uses the exact beta-distribution variance formula, avoiding removable
0/0 for alpha+beta=-1. The eigenvalues of the real Jacobi tridiagonal
are quadrature *nodes*, NOT eigenvalues of the orthogonal evaluation
matrix. Generic construction O(N^3); dense multiplication O(N^2).
For the two *exact* Chebyshev specializations, a new **true FFT
operator** returns the very same mathematically defined Gauss evaluation
matrix, **without** constructing or caching it:
`Q_I=DCTIII*diag((-1).^k)` at alpha=beta=-1/2, and
`Q_II=DSTI*diag((-1).^k)` at alpha=beta=+1/2. The alternating signs
follow ascending-node ordering (x_j=-cos(j*pi/(n+1)) for second kind)
and positive-leading polynomial normalization. These forward/adjoint
operators cost O(N log N) per vector and O(n) workspace, via the existing
FFT DCT/DST kernels; roundoff means their explicit numerical products
need only agree to tolerance. For all other Jacobi parameters, including
Legendre alpha=beta=0, `'operator'` is **rejected**, not implemented
using a hidden dense matrix or mislabeled as fast. Genuine Hale--Townsend
or Slevinsky fast connection algorithms remain issue #39. Matrix
eigenphases are unknown, unrelated to Gauss nodes. [104,105,106]

### Gauss--Hermite transform — [`gauss_hermite_transform.m`](../gauss_hermite_transform.m)

**Order:** $N=n$. **Construction:** For the normalized Gaussian weight
`exp(-x^2)/sqrt(pi)` and orthonormal physicists' Hermite polynomials `p_k`,
the weighted Gauss evaluation matrix is

```math
Q_{jk}=\sqrt{w_j}\,p_{k-1}(x_j),\qquad Q^TQ=I.
```

Here the quadrature weights `w_j` sum to one. Compute `Q` as the
transpose of the first-component-positive eigenvector matrix of the
symmetric Jacobi tridiagonal with zero diagonal and offdiagonal
`sqrt(k/2)` for k=1,...,n-1. **Spectral properties:** conjugate-paired
unit-circle eigenvalues; no general closed form for their angles.
Gauss--Hermite nodes are eigenvalues of the *Jacobi* matrix, not of `Q`.
**Computational cost:** O(N^3) setup and O(N^2) explicit multiply.
**Applications:** Hermite spectral approximations and Gaussian quadrature. [104,109]

### Generalized Gauss--Laguerre transform — [`gauss_laguerre_transform.m`](../gauss_laguerre_transform.m)

**Order:** $N=n$. **Parameter:** alpha > -1. The normalized positive measure
has density proportional to `x^alpha exp(-x)` on `[0,infinity)`.
The Jacobi tridiagonal for the positive-leading orthonormal Laguerre
polynomials has diagonal `2*k+alpha+1`, k=0,...,n-1, and offdiagonal
`sqrt(k*(k+alpha))`, k=1,...,n-1. Its ordered eigenvectors yield
the weighted Gauss evaluation matrix `Q` with `Q'*Q=I`.
**Spectral properties:** conjugate eigenphase pairs; quadrature nodes
are NOT eigenvalues of `Q`. **Applications:** radial numerical spectral
methods. **Action:** explicit dense, not fast. [104,106]

### Associated-Legendre specialization — [`gauss_jacobi_transform.m`](../gauss_jacobi_transform.m)

**Not a separate generator:** call `gauss_jacobi_transform(n,m,m)` for nonnegative integer m. The associated Legendre
functions with degrees ell=m,...,m+n-1 are the Jacobi(m,m) orthonormal
polynomials multiplied by `(1-x^2)^(m/2)`, up to an overall phase/sign.
This is the existing Jacobi(m,m) quadrature, not a second mathematical construction. The corresponding n-point Gauss--Jacobi quadrature gives a square real
orthogonal matrix. Its Jacobi tridiagonal has zero diagonal and

```math
a_k=\sqrt{\frac{k(k+2m)}{(2k+2m)^2-1}},\quad k=1,\ldots,n-1.
```

**Boundary:** this fixed-m block is a square orthogonal transform.
The usual all-degrees spherical harmonic transform on a
Gauss--Legendre-longitude grid is generally rectangular. Fast
butterfly/connection transforms are a related *algorithmic* topic;
we do not claim a fast action here. [104,106,107]

### Wigner degree-ell rotation — [`wigner_rotation.m`](../wigner_rotation.m)

**Order:** $N=2\ell+1$, ell=0,1,... . **Parameters:** z-y-z Euler angles
(alpha,beta,gamma) in radians. With m=-ell,...,ell, the generator
builds the complex spherical harmonic representation

```math
D^{(\ell)}(R)=e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z},
\qquad J_z=\mathrm{diag}(-\ell,\ldots,\ell).
```

The middle factor is calculated by exponentiating the real
skew-symmetric matrix representing `-i*Jy`. It is unitary without
evaluating numerically unstable factorial Wigner sums.
For rotation angle delta satisfying
`cos(delta/2)=cos(beta/2)*cos((alpha+gamma)/2)`, the exact eigenvalues
are `exp(-i*m*delta)`, m=-ell,...,ell. **Action:** explicit dense matrix
constructed using `expm`; no fast rotation algorithm claimed. [108]

### Degree-truncated spherical rotation — [`spherical_harmonic_rotation.m`](../spherical_harmonic_rotation.m)

**Order:** $N=L^2$, with degrees ell=0,...,L-1. **Structure:** sparse block
diagonal with `D^(ell)(R)` blocks. The phase `exp(-i*m*delta)`
has multiplicity `L-abs(m)` for m=-(L-1),...,L-1, before merging
coincident phases. An exact Fourier-moment identity is

```math
\mathrm{tr}(D_L(R)^k)=
\left|\sum_{j=0}^{L-1}e^{-ikj\delta}\right|^2.
```

If delta=c/L for fixed 0<c<pi, the normalized eigenphase density
converges to the triangular law `(1-|theta|/c)/c` for |theta|<=c,
and zero elsewhere. This is a derived spectral benchmark, NOT
a claim that the cited rotation paper states the limiting law.
**Applications:** controlled gaps, clusters, exact moments, SO(3)
representation structure. [108]

### Periodic sinc/Dirichlet fractional shift — [`periodic_sinc_shift.m`](../periodic_sinc_shift.m)

**Order:** $N=n$ for any positive integer $n$, with an explicit distinction between odd-
and even-order constructions. **Parameter:** real shift tau in grid samples.
The default call is unchanged: it requires ODD n and gives the real
orthogonal translation with the centered Fourier modes
`k=-(n-1)/2,...,(n-1)/2`. Let F be the normalized negative-exponent DFT:

```math
U_n(\tau)=F^*\,\mathrm{diag}(e^{-2\pi i k\tau/n})F.
```

For odd n, each nonzero frequency has its negative-frequency conjugate,
so the matrix is real, orthogonal and circulant. It is a periodic
Dirichlet-cardinal interpolation matrix, and `U(1)` is the forward cyclic
permutation P. Its eigenvalues are the displayed finite Fourier multipliers.

**Even-order obstruction.** The cyclic permutation is one cycle of length
even n, and therefore `det(P)=(-1)^(n-1)=-1`. Every continuous path
`Q(t)` in the real orthogonal group with `Q(0)=I` has determinant +1
throughout (the determinant maps continuously into the discrete set
`{+1,-1}`). Consequently no continuous real-orthogonal fractional shift
can join I to P for even n, whether or not a group composition law is
required. In the DFT basis, the single real Nyquist mode must have a
continuous real eigenvalue, hence must remain +1 when starting at I,
whereas P acts as -1 on it.

The **explicit complex-unitary variant** uses
`periodic_sinc_shift(n,tau,'matrix','complex')` for EVEN n. Frequencies
are in FFT order with positive Nyquist `k=n/2`: `0,1,...,n/2,
-n/2+1,...,-1`. This defines a continuous one-parameter unitary group,
including the full `U(1)=P`, but matrices are generally complex rather
than real. It is also the corresponding choice of phase parameters in
the existing `phase_circulant` family, rather than an unrelated class.

The **Nyquist-frozen variant** uses
`periodic_sinc_shift(n,tau,'matrix','nyquist-frozen')`. All paired
frequencies retain their multipliers, but the Nyquist multiplier is fixed
at +1. This is a continuous real orthogonal one-parameter group in
SO(n), **not** an interpolation to P. At tau=1 it equals
`P+2vv^T`, where `v_j=(-1)^j/sqrt(n)` for `j=0,...,n-1`.
In particular, the even-order `n=2` frozen family is the identity.

For every admitted variant, `U(tau)U(sigma)=U(tau+sigma)`.
The optional `'operator'` mode applies the matrix and its conjugate
adjoint to multiple columns in O(N log N) using FFTs; `'spectral'`
returns exactly the finite trace-measure eigenphases and moments.
As elsewhere in the library, operator flag `'transp'` means the
conjugate adjoint. No universal continuous limiting phase density
is asserted without specifying how tau scales with n.

**Applications:** periodic interpolation and fractional translation,
Dirichlet kernels, exact eigenphase benchmarks, determinant topology
of O(N), and verification of unitary versus orthogonal numerical
implementations. Existing code uses the same odd-order default and is
unaffected. [1,52]

### Not-yet-admitted related candidates

**Hankel/Bessel** Bessel-zero discretizations often give only
approximately orthogonal finite transform matrices; a fast DHT
does not by itself ensure unitarity. **RBF** positive-definite Gram
matrices are nonunitary and their square polar factors are usually
the identity. **Zernike, disk, triangle, spherical Gauss--Laguerre**
transforms require square and quadrature-exact sampling conventions:
standard grids may instead define rectangular isometries.
See Slevinsky's [FastTransforms](https://mikaelslevinsky.github.io/FastTransforms/)
for related fast *algorithms*, not a blanket assertion of unitarity.

### Stiefel geodesic ambient lift — [`stiefel_geodesic.m`](../stiefel_geodesic.m)

**Order:** $N=n\ge2p$. Orthonormal frames B=[Y,Q] (n-by-2p), skew-Hermitian
p-by-p A, and arbitrary p-by-p R specify tangent H=YA+QR. With

```math
K=\begin{bmatrix}A&-R^*\\R&0\end{bmatrix},
\qquad G(t)=I+B(e^{tK}-I)B^*,
```

the frame geodesic Y(t)=G(t)Y is canonical-metric Stiefel motion.
This extends the Grassmann lift with a nonzero vertical rotation A,
but does **not** model intrinsic Stiefel tangent-bundle parallel transport.
The large unitary spectrum is known from a small 2p-by-2p skew-Hermitian
eigenproblem and n-2p fixed eigenvalues. The generator computes this small
eigenproblem for its spectral reference, not a full $N\times N$ eigensolve (here $N=n$).
Operator mode stores a 2p-by-2p exponential and applies the low-rank update
in O(np) arithmetic. Edelman–Arias–Smith [97].

### Intrinsic canonical Stiefel tangent transport — [`stiefel_parallel_transport.m`](../stiefel_parallel_transport.m)

**Order:** `d=np-p(p+1)/2`, the *real tangent-space dimension*
of `St(p,n)`, for `1<=p<=n`. This is a **genuinely different
square matrix** from the ambient `n`-by-`n` Stiefel geodesic
lift above. The source and target are different tangent spaces;
the matrix is represented in canonical **orthonormal moving frames**.

The round canonical Stiefel metric and geodesic generator are

```math
g_Y(V,W)=\operatorname{tr}(V^T(I-YY^T/2)W),\qquad
K=\begin{bmatrix}A&-B^T\\B&0\end{bmatrix},
\quad Y(t)=e^{tK}\begin{bmatrix}I_p\\0\end{bmatrix},
```

where `A` is real skew-symmetric `p`-by-`p` and `B` is
`(n-p)`-by-`p`. The quotient `SO(n)/SO(n-p)` carries
the normal-homogeneous metric induced by half the Frobenius
form on `so(n)`. Its horizontal space comprises skew
matrices whose lower-right `(n-p)` block is zero.

Use the orthonormal horizontal basis `F_ij=e_i e_j^T-e_j e_i^T`,
for `1<=i<=p, i<j<=n`. At `Y(t)`, the tangent basis vectors
are `e^{tK}F_ij Y(0)`; their canonical metric Gram matrix is
exactly the identity. The quotient Levi–Civita connection, in
these moving coordinates, has constant real skew generator

```math
S_{(i,j),(a,b)}=-\tfrac12\big([K,F_{ab}]\big)_{ij},
\qquad Q(t)=e^{tS},\quad Q(t)^TQ(t)=I_d.
```

The sign follows by projecting the bi-invariant group
connection `E'(t)+[K,E(t)]/2` onto the horizontal space.
Thus `Q(t)` **is intrinsic parallel transport**, not a
tangent retraction, ambient rotation, or arbitrary
orthogonal polar factor. Arbitrary initial orthogonal
frames are obtained by left `SO(n)` equivariance.

**Special cases:** for `p=1`, `[m,m]` is vertical, so
`S=0`: spherical transport is the identity in the moving
geodesic frame. For `p=n`, the quotient has no vertical
subspace and `S=-ad_K/2`, recovering the bi-invariant
orthogonal-group parallel transport. For general `p`,
the geodesic velocity is a fixed vector in these transported
coordinates because `[K,K]=0`.

**Algorithms and scope.** The implementation forms all
`d`-by-`d` structure constants and evaluates a **dense
matrix exponential**, `O(d²)` generator setup,
`O(d³)` exponential, `O(d²)` storage and multiplication
per right-hand side. It is limited to `n<=64,d<=256`.
The improved `O(np²)` small-time and `O(tp³)`
large-time exponential-action method of Nguyen–Sommer [168],
general flag-manifold transport and other Stiefel metrics
are **not** implemented by this prototype.

The tests compare its intrinsic tangent vectors against
an **independently integrated ambient canonical
Levi–Civita ODE**, check tangent constraints,
isometry, time reversal, and parallel geodesic velocity.
Generic eigenphases of the finite **transport** matrix are
not assigned a closed law; the `spectral` mode returns
`unknown`. [97,168]

### Non-Abelian Grassmannian loop holonomy — [`grassmann_holonomy.m`](../grassmann_holonomy.m)

**Order:** $N=r$, from n-by-r orthonormal frames Y_j of a closed subspace
polygon in complex Gr(r,n). On each unique shortest geodesic arc the
canonical tautological-bundle transport, represented in frame coordinates,
is the unitary polar factor of Y_(j+1)^*Y_j (provided this overlap is
invertible). Multiplication in traversal order gives H in U(r).
Changing each frame by an arbitrary gauge Q_j changes the closed
holonomy by Q_1^* H Q_1: eigenphases remain gauge invariant.
This is not Grassmannian Levi-Civita tangent-bundle holonomy.

For r=1 the scalar phase is the normalized Bargmann-invariant product
of pairwise overlaps and is provided exactly without eig. For general
r>=2, the phases can be genuinely non-Abelian; no generic closed-form
phase law is claimed. Source on geometric phases [110] and Grassmann
geometry [111]. The code's default frame polygon is a reproducible
complex test case, not extracted from an application dataset.

### Wigner rotations: canonical spin extension — [`wigner_rotation.m`](../wigner_rotation.m)

The existing Wigner generator also covers half-integer SU(2) spins and the nontrivial 2pi double-cover sign. Its integer-spin version factors through SO(3); both use the same angular-momentum construction, now cross-referenced to Kostelec–Rockmore [112]. No `wigner_geodesic` generator is added.


## Hyperbolic and gauge holonomy (November 2026)

The spherical polygon holonomy is already implemented by [`geodesic_holonomy.m`](../geodesic_holonomy.m), now additionally citing Needham [113]. No duplicate spherical generator is added.

### Hyperbolic geodesic-triangle holonomy — [`hyperbolic_triangle_holonomy.m`](../hyperbolic_triangle_holonomy.m)

Represent the hyperbolic plane H² as the Poincare unit disk
with complex coordinates z_j of three geodesic vertices.
The Riemannian metric has Gaussian curvature K=-1:

```math
ds^2=\frac{4|dz|^2}{(1-|z|^2)^2},\qquad |z|<1.
```

Define the **oriented** hyperbolic triangle area

```math
A=-2\,\arg\!\left[
(1-\overline z_1 z_2)
(1-\overline z_2 z_3)
(1-\overline z_3 z_1)\right].
```

The value lies between -pi and pi for a nondegenerate
geodesic triangle. Independently, the unsigned area
satisfies Gauss-Bonnet's angle defect

```math
|A|=\pi-(\alpha+\beta+\gamma),
```

where the angles are measured in the hyperbolic metric.
The Levi-Civita tangent-space holonomy is the real
2x2 rotation by **theta=-A**, with an oriented basis
fixed by the Poincare disk's complex coordinate:

```math
R(\theta)=
\begin{bmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{bmatrix},
\qquad \theta=-A.
```

Its exact eigenvalues are e^(+/-i theta), so the normalized
finite trace moments equal cos(k theta).
For database scalability, the function accepts **n>=2**,
returning R(theta) direct-summed with I_(n-2). Only the
original **2x2 block is intrinsic** to hyperbolic tangent
transport; higher-dimensional identity coordinates are
spectator benchmark padding, not higher-dimensional
hyperbolic holonomy. The embedded spectrum is the
angle pair and n-2 copies of 1. Both matrix and O(N)
implicit forward/transpose actions are supported.

**Applications:** hyperbolic curvature and area diagnostics,
geodesic holonomy, exact-spectrum normal-matrix tests,
and orientation-sensitive real orthogonal benchmarks.
**References:** [113] (Needham) and [114] (Voight).

### Non-Abelian lattice-gauge Wilson-loop matrices — [`lattice_wilson_loop.m`](../lattice_wilson_loop.m)

Assign a unitary link L_j in U(n) to each directed edge
of a closed lattice contour. Traversing a link in reverse
uses its Hermitian adjoint. For ordered links

```math
W=L_1L_2\cdots L_m.
```

This **matrix-valued Wilson loop** is unitary, but its
eigenphases generally depend on the gauge field and the loop.
The gauge-invariant Wilson observable is typically
tr(W)/n; the trace alone does not give an eigenphase density.

For site-dependent local changes of basis G_j in U(n),
the oriented edge transforms by

```math
L_j\longmapsto G_jL_jG_{j+1}^*,\qquad G_{m+1}=G_1.
```

The intermediate gauges cancel **exactly**, giving

```math
W\longmapsto G_1WG_1^*.
```

Thus its trace, determinant, eigenvalues and normalized
spectral measure are invariant under local gauge changes,
although individual matrix entries are not. This is a
natural gauge-covariant unitary matrix benchmark rather
than an arbitrarily conjugated random unitary matrix.

The default case uses four deterministic noncommuting
special-unitary link matrices, each formed as exp(i H_j)
for a traceless Hermitian H_j. A user may provide arbitrary
N-by-N-by-M **unitary** link arrays, including general U(n)
links with non-unit determinant. The optional "links" mode
returns the link array. The matrix-free mode applies the
ordered factors in O(m n²) per right-hand side; default
link construction itself uses dense matrix exponentials.
A generic exact finite spectral law is not asserted.

**Applications:** lattice Yang-Mills, gauge-invariant
observables, non-Abelian parallel transport, SU(n) spectral
statistics and symmetry-preserving numerical methods.
**Reference:** [115] (Ipp and Müller).


### Heterogeneous gyroscopic/Coriolis oscillators — [`gyroscopic_oscillator.m`](../gyroscopic_oscillator.m)

**Order:** $N=2n$ real orthogonal ($n\geq2$ mechanical displacements). **Application:** conservative gyroscopic vibrations and rotating-structure models. This is a reproducible *benchmark specialization*, not experimental apparatus data.

For $j=1,\ldots,n$, take diagonal positive $M$ with entries $m_j=1+0.2j/n+0.05\sin j$. Let $D$ be the forward-difference matrix with row $j$ equal to $e_{j+1}^T-e_j^T$, and let $E_{j+1,j}=1$ with all other entries zero. Define

```math
K=\omega_0^2M+\kappa D^TD+
0.15\,\mathop{\mathrm{diag}}(1+0.1\cos j),\qquad
C=\gamma(E-E^T).
```

Parameters are $\omega_0>0$, $\kappa\geq0$, $\gamma\geq0$, with $h>0$. For $\gamma>0$ the gyroscopic matrix is nonzero and skew-symmetric, both $M$ and $K$ are sparse SPD, and the mass and stiffness generally do not commute. The physical equation and its **continuous-time gyroscopic QEP** are

```math
M\ddot q+C\dot q+Kq=0,\qquad Q(\lambda)=\lambda^2M+\lambda C+K.
```

The continuous eigenvalues are imaginary and semisimple under the stated positive-energy assumptions. In velocity coordinates $x=(q,v)$,

```math
A_0=\begin{bmatrix}0&I\\-M^{-1}K&-M^{-1}C\end{bmatrix},
\qquad G=\mathop{\mathrm{diag}}(K,M),
\qquad A_0^TG+GA_0=0.
```

Implicit midpoint and two- and three-stage Gauss–Legendre steps use $S=R_s(hA_0)$, with $R_s(z)=P_s(z)/P_s(-z)$ and

```math
P_1(z)=1+z/2,\quad
P_2(z)=1+z/2+z^2/12,\quad
P_3(z)=1+z/2+z^2/10+z^3/120.
```

Consequently $S^TGS=G$, and with the upper Cholesky factor $T^TT=G$, the **returned matrix** $U=TST^{-1}$ is orthogonal. The **operator** mode applies $U$ and its transpose through the sparse rational step matrices, with sparse solves and Cholesky factors; fill-in and conditioning govern performance.

The **pencil** mode exposes $(GS,G)$, the physical $S$, the metric $G$, normalized $U$, and $M,K,C$. Its generalized eigenvalues lie on the unit circle, but it is not in general a strict Euclidean-unitary coefficient pencil. The **quadratic** mode returns coefficients $(K,C,M)$; its continuous eigenvalues $\lambda$ must not be confused with the Gauss multipliers $R_s(h\lambda)$.

The Lagrangian $L=\frac12v^TMv+\frac12v^TCq-\frac12q^TKq$ gives canonical momentum $p=Mv+Cq/2$. The change of coordinates

```math
B=\begin{bmatrix}I&0\\C/2&M\end{bmatrix},\qquad
(BSB^{-1})^TJ(BSB^{-1})=J,\qquad
J=\begin{bmatrix}0&I\\-I&0\end{bmatrix}
```

shows the correct canonical symplecticity. The energy-normalized matrix $U$ need **not** preserve the same standard $J$. Both metric and canonical symplecticity are independently tested.

The generic finite eigenphase law is **unknown analytically**. Numerical QEP roots can provide independent checks; the $\gamma=0$ specialization is a non-gyroscopic limit overlapping the stiffness-only families.

**References:** Qian–Lin [125] (gyroscopic QEP), Hairer–Lubich–Wanner [101] and Ohsawa [126] (geometric integration and quadratic invariants). Kuzmanović Ivičić–Miodragović [127] treat a different negative-stiffness gyroscopic-stabilization regime, not a justification for this positive-energy construction.


## Clifford/stabilizer circuits — [`clifford_circuit.m`](../clifford_circuit.m)

**Matrix order:** $2^q$, where $1\le q\le20$ is a qubit count;
dense assembly is limited to $q\le10$. The computational basis
is $|b_1\cdots b_q\rangle$, with the first qubit most significant.
This finite family consists of ordered words in standard
**Clifford gates** H, S, S-dagger, Pauli X/Y/Z, CNOT, CZ and SWAP.
A MATLAB N-by-2 cell list gives a gate name and its ordered site
indices, with chronological convention $U=G_m\cdots G_1$.
It reuses the existing `um_circuit_apply` and
`um_gate_apply` local-statevector algorithms. It is
not equivalent to the existing generic `dual_unitary_circuit`:
ordinary dual-unitary two-qubit gates need not be Clifford.

The global-phase convention is fixed rather than quotienting the
actual matrices by arbitrary scalars:

```math
H=\frac1{\sqrt2}\begin{bmatrix}1&1\\1&-1\end{bmatrix},
\qquad S=\mathop{\mathrm{diag}}(1,i),\qquad
Y=\begin{bmatrix}0&-i\\i&0\end{bmatrix}.
```

The first listed CNOT site is the **control** qubit; its local basis
is $|00\rangle,|01\rangle,|10\rangle,|11\rangle$.
No arbitrary continuous global scalar phase is inserted.
Clifford gate words normalize the Pauli group: each
$UPU^*$, for a Hermitian Pauli tensor product $P$, is another
Hermitian Pauli tensor product with sign $+1$ or $-1$.
The projective Clifford group is finite. With the specified
finite-phase gate generators, all eigenvalues are roots of unity;
however, their exact **finite multiset** for a general word is
not inferred from the group order, and no CUE eigenphase density
or level-spacing law is asserted.

**Exactly known trace spectra for elementary presets**
with $N=2^q$, allowing idle tensor factors:

| Circuit preset | Exact finite eigenvalues and multiplicities |
| --- | --- |
| identity | $+1$ with multiplicity $N$ |
| hadamard, pauli_x, pauli_y, pauli_z | $+1$ and $-1$, each multiplicity $N/2$ |
| phase | $1$ and $i$, each multiplicity $N/2$ |
| phase_dagger | $1$ and $-i$, each multiplicity $N/2$ |
| cnot, cz, swap ($q\ge2$) | $+1$ multiplicity $3N/4$, $-1$ multiplicity $N/4$ |

The first one or two qubits carry each named gate; the rest are idle.
Any custom one-gate word has the corresponding exact trace spectrum,
regardless of the chosen valid sites. The default bell circuit
performs H on site 1 followed by CNOT(1,2) and creates a Bell state
from $|00\rangle$. The generic multi-gate spectrum is conservatively
marked **unknown analytically** even though the explicit matrix
could be diagonalized for small q. Exact presets expose lazy
finite eigenvalues, probability-weighted atoms and Fourier
moments $m_k=N^{-1}\mathop{\mathrm{tr}}(U^k)$.

Operator mode supports forward and adjoint local-gate actions
with $O(L2^q)$ operations for L elementary gates, using
$O(2^q)$ state storage per right-hand side, without a
dense $4^q$-entry matrix. Its regression tests independently
assemble small test circuits, check exact phase multiplicities,
and verify the Clifford Pauli normalizer property.

**Applications:** stabilizer-state simulation, Clifford twirling,
fault-tolerant quantum computing and highly atomic non-CUE
spectral-density benchmarks.

**References:** Dehaene and De Moor [128] describe the binary
quadratic and stabilizer structure of the Clifford group.
Gross, Nezami and Walter [129] develop representation-theoretic
properties and related applications. The simple finite phase
multiplicities in the table follow directly from the canonical
gate matrices.



## Fibonacci-anyon fixed-charge braid representations — [`fibonacci_braid.m`](../fibonacci_braid.m)

**Application:** topological quantum computation, non-Abelian braid
representations and structured unit-circle test matrices. This generator
does **not** build an ordinary tensor product of qubits: it acts on
the *fusion space* of n identical Fibonacci anyons, all of charge
$\tau$, with a fixed total charge $\mathbf{1}$ (vacuum) or $\tau$.
The fusion rule is $\tau\times\tau=\mathbf{1}+\tau$.
In the left-associated fusion-tree basis, label the charges obtained
after fusing the first j anyons by $a_j\in\{\mathbf1,\tau\}$,
$a_1=\tau$, $a_n$ the selected total charge. Only admissible
fusion paths are retained, ordered lexicographically with
$\mathbf1=0,\tau=1$.

For $n\ge3$ the **exact matrix dimensions** are $F_{n-1}$ in
the vacuum sector and $F_n$ in the tau sector, where
$F_1=F_2=1$. This count comes from the admissible fusion paths,
not the ambient Hilbert space of n qubits. The generator
supports $n\le20$; dense assembly is limited to fusion dimension
at most 512.

In the standard real F-gauge, with $\varphi=(1+\sqrt5)/2$,
the nontrivial recoupling block and counterclockwise exchange
phases are

```math
F=\begin{bmatrix}
\varphi^{-1}&\varphi^{-1/2}\\
\varphi^{-1/2}&-\varphi^{-1}
\end{bmatrix},\quad
R_{\mathbf1}^{\tau\tau}=e^{-4\pi i/5},\qquad
R_{\tau}^{\tau\tau}=e^{3\pi i/5}.
```

The generator $\sigma_1$ is diagonal, its phase determined by
the first two anyons' intermediate fusion charge $a_2$.
For $\sigma_i$ with $i\ge2$, inspect the neighboring path
labels $a_{i-1},a_{i+1}$: the only **two-dimensional** mixing
case is $a_{i-1}=a_{i+1}=\tau$, where the block in
$a_i\in\{\mathbf1,\tau\}$ is
$F\mathop{\mathrm{diag}}(R_{\mathbf1},R_\tau)F$.
All other admissible cases are scalar R-channel phases.
The resulting sparse matrices satisfy the Artin braid relations:

```math
\sigma_i\sigma_{i+1}\sigma_i=
\sigma_{i+1}\sigma_i\sigma_{i+1},\qquad
\sigma_i\sigma_j=\sigma_j\sigma_i\quad(|i-j|\ge2).
```

A signed integer vector is a braid word: positive index i
applies $\sigma_i$, negative index $-i$ applies
$\sigma_i^*$, and the first entry acts first. The empty
word is identity. Both ordinary and opposite chirality are
supported; opposite chirality complex conjugates the R-phases
and hence the braid matrices, **not** a mere basis gauge
change. Any unitary gauge change of the fusion-tree basis
conjugates the represented matrices and preserves each word's
eigenvalues. The overall phase and handedness conventions
must be retained when comparing sources.

**Finite spectral annotation.** The identity has only eigenvalue
$+1$ with multiplicity the sector dimension. Each single positive
$\sigma_i$ has exactly the phases $R_{\mathbf1},R_\tau$
with multiplicities equal to the counts of fusion paths whose
$a_2$ label is respectively vacuum or tau. Single negative
generators have the conjugate phases. A zero multiplicity
channel is omitted. These facts give exact normalized trace
moments without constructing a dense matrix. Arbitrary longer
braid words are **not** assigned a generic closed-form finite
phase law; Fibonacci braid images can be noncommuting and
spectrally much richer than the finite Clifford circuit model.

Each sparse generator has at most two nonzeros per column.
The `operator` output applies braid factors chronologically
and reverse-order adjoint factors using sparse matrix-vector
products; it does not form the full dense product.
The `generators` mode also exposes the complete fusion-path
basis, F/R symbols, chirality and sparse $\sigma_i$ factors.
Tests verify dimension recurrences, unitarity, nontrivial
Yang–Baxter and distant-commutation relations, analytic
three- and four-anyon blocks, exact single-generator
spectral multiplicities, inverse words and the two chiralities.

**References:** Trebst, Troyer, Wang and Ludwig [130] give
the fusion, F and R constructions and braid matrices.
Carnahan, Zeuch and Bonesteel [131] discuss multi-anyon
braid words in topological quantum computation.



## Flux-normalized resonant multichannel scattering — [`resonant_scattering.m`](../resonant_scattering.m)

**Application and distinction.** This is a genuine **open**
finite-dimensional resonant-cavity scattering benchmark, distinct
from the closed directed-bond Floquet operator
`metric_quantum_graph` and the staggered unitary
`scattering_zipper`. At a specified **real** energy
$E$, incoming/outgoing amplitudes lie in an orthonormal,
unit-flux basis of n asymptotic open channels. Every open
channel is retained; there is no absorption, channel threshold
or hidden inelastic channel. The internal cavity has $m\le n$
states described by Hermitian $H\in\mathbb C^{m\times m}$,
and a full-row-rank coupling
$W\in\mathbb C^{m\times n}$.

Set $\Gamma=WW^*>0$ and adopt a wide-band coupling convention
in which $\Gamma$ is energy-independent (the factor of $\pi$
in some Mahaux–Weidenmüller conventions is absorbed into W).
The on-shell scattering matrix is

```math
A(z)=zI_m-H+i\Gamma,\qquad
S(E)=I_n-2iW^*A(E)^{-1}W,\quad E\in\mathbb R.
```

This choice fixes the free/background scattering matrix to $I_n$.
It is a dimensionally explicit finite version of the
Mahaux–Weidenmüller resolvent formula used in transport theory
[132,134]. Because $A(E)-A(E)^*=2iWW^*$,
the resolvent identity proves
$S(E)^*S(E)=I_n$ for every real E.
The full-row-rank assumption gives $\Gamma>0$ and excludes
singular real-energy bound/dark internal states in this model.
For **real symmetric** H and **real** W,
$S(E)^{\mathsf T}=S(E)$ (reciprocity/time-reversal symmetry),
even though $S(E)$ is generally not real.

**Resonances, derivative and Wigner–Smith delay.**
The effective internal operator $H_{\mathrm{eff}}=H-iWW^*$
is dissipative and non-Hermitian, with poles of the analytic
continuation $S(z)$ at its eigenvalues in the *strict lower
half-plane* (under the full-row-rank assumption). Its complex
poles are not the real eigenvalues of H and are **not**
unit-circle eigenvalues of a real-energy S(E). The `model`
mode returns both the pole list and an evaluable meromorphic
continuation. It must not be evaluated *at* a pole. Direct
differentiation at a real energy gives

```math
\frac{dS}{dE}=2iW^*A(E)^{-2}W,\qquad
Q(E)=-iS(E)^*\frac{dS}{dE}
=2W^*A(E)^{-*}A(E)^{-1}W\succeq0.
```

The Hermitian **Wigner–Smith time-delay matrix** Q has
proper-delay eigenvalues in energy-inverse units with
$\hbar=1$ [133]. If $m<n$, Q has rank m and therefore
$n-m$ zero delay eigenvalues. The `time_delay`
mode returns Q, **not** a unitary. The `derivative`
mode returns the analytic energy derivative of S, also **not**
a unitary. For general coupled H, there is no claimed
closed-form finite scattering-phase law. Nor is the
deterministic benchmark automatically a Haar/CUE or
Poisson-kernel ensemble: those are separate probabilistic
transport models with coupling assumptions [132].

**Exactly known finite trace spectra in decoupled cases.**
When H and $\Gamma$ are both diagonal in the supplied internal
basis (equivalently, the coupling rows are pairwise orthogonal),
write $h_j=H_{jj}$ and $\gamma_j=\Gamma_{jj}>0$.
The n on-shell scattering eigenvalues are exactly

```math
\lambda_j(E)=
\frac{E-h_j-i\gamma_j}{E-h_j+i\gamma_j},
\quad j=1,\dots,m;\qquad
\lambda_{m+1}=\cdots=\lambda_n=1.
```

Thus the `spectral` mode returns the **exact finite
trace measure** (with equal eigenvalue weights $1/n$)
and exact Fourier moments without calling eig(S).
For example, with one resonant open channel, its eigenphase
derivative is the Lorentzian delay
$2\gamma/((E-h)^2+\gamma^2)$. In the default coupled
real-tridiagonal cavity with diagonal W, the finite
trace phase law is conservatively marked *unknown analytically*;
it is not replaced by a numerical eigendecomposition.

**Observed channel blocks and loss of flux.** Any submatrix
$B(E)$ obtained by deleting the same incoming and outgoing
channels from the full unitary satisfies $B^*B\preceq I$,
but is generally a strict **contraction**. Such a submatrix
must not be registered as another unitary family. The tests
check exact flux conservation for the complete matrix and
strict flux loss for a chosen incomplete observation, along
with reciprocity, real-energy delay positivity, numerical
differentiation, resonance locations and unwrapped
single-channel eigenphase curves.

**Defaults and dimensions.** The default m=n cavity uses
reproducible real-symmetric nearest-neighbor couplings and
positive diagonal coupling strengths. User-supplied complex
Hermitian H and rectangular full-row-rank W are permitted.
Dense on-shell S has order n (n<=512); no unsupported
matrix-free acceleration is claimed. This is a constant
wide-band finite model, **not** a microscopic waveguide
with energy-dependent channel thresholds.

**References:** Beenakker [132] reviews multichannel
quantum-transport scattering models and explains CUE versus
Poisson-kernel ensemble assumptions. Brouwer, Frahm and
Beenakker [133] give the Wigner–Smith proper-delay framework.
Erdős, Krüger and Nemish [134] state a precise finite
Mahaux–Weidenmüller scattering formula and its quantum-dot
transport setting. The exact diagonal-channel phase identity
above is elementary algebra applied to that fixed convention.



## Symmetry-protected cluster-spin Floquet drive — [`cluster_floquet.m`](../cluster_floquet.m)

**A reproducible ideal spin-chain representative, not a generic
Floquet topology classification.** The matrix order is $2^n$ for
$n=4,\ldots,20$ spins. Site 1 is the most-significant computational
qubit, $Z|0\rangle=|0\rangle$. The drive is an ordered product
of two local unitary layers with a fixed global-phase convention:

```math
\begin{aligned}
K_j&=Z_{j-1}X_jZ_{j+1},\\
U_F(\delta)&=P(\delta)
 \exp\left(-i\sum_{j\in J}\theta_j K_j\right),\\
P(\delta)&=\prod_{j=1}^n
 X_j\exp(-i\delta X_j)
 =i^n\prod_{j=1}^n\exp[-i(\pi/2+\delta)X_j].
\end{aligned}
```

Here $J=\{2,\ldots,n-1\}$ for **open** boundaries, so both
boundary cluster stabilizers are deliberately absent.
For **periodic** boundaries use every $j=1,\ldots,n$ with
indices taken cyclically; n must be even for the
bipartite symmetry. The real parameters $\theta_j$ may be
inhomogeneous/disordered. The real pulse detuning $\delta$
is zero in the exactly solvable limit. The scalar
$i^n$ fixes the global quasienergy convention: at $\delta=0$,
$P(0)=\prod_jX_j$. One must not drop this phase when checking
individual eigenphase multiplicities.

The two commuting parity symmetries are

```math
G_{\mathrm{o}}=\prod_{j\ \mathrm{odd}}X_j,\qquad
G_{\mathrm{e}}=\prod_{j\ \mathrm{even}}X_j,\qquad
[U_F(\delta),G_{\mathrm{o}}]=[U_F(\delta),G_{\mathrm{e}}]=0.
```

The ideal cluster stabilizers $K_j$ commute mutually,
and $G_{\mathrm{o}},G_{\mathrm{e}}$ commute with both drive
layers for every $\delta$. For **open** chains at
$\delta=0$, the two strictly boundary-supported observables
$Z_1$ and $Z_n$ commute with the bulk cluster layer but
anticommute with $P(0)$. Consequently

```math
U_F(0)Z_1=-Z_1 U_F(0),\qquad
U_F(0)Z_n=-Z_n U_F(0).
```

Thus a quasienergy eigenphase $\lambda$ is paired exactly
with $-\lambda$, separated by $\pi$ modulo $2\pi$.
This is an **operator-level edge diagnostic**, not an assumption
about a limiting random-matrix law. Other edge operators
$X_1Z_2$ and $Z_{n-1}X_n$ also anticommute with the ideal
drive. Sector-resolved pairing is visible because $Z_1$
reverses the odd symmetry eigenvalue while leaving the
even one fixed; for even n, $Z_n$ similarly reverses
the even symmetry eigenvalue. With periodic boundaries,
local $Z_1$ no longer commutes with the bulk K terms
and these exact edge identities do not apply.

**Exact finite trace spectra in the commuting ideal limit.**
For the open chain, choosing every interior stabilizer
sign $k_j\in\{-1,1\}$ independently gives

```math
\lambda_{k,\pm}=
\pm\exp\left(-i\sum_{j=2}^{n-1}\theta_jk_j\right),
\qquad\text{each repeated twice}.
```

There are $2^{n-2}$ interior sign configurations and
four edge states per configuration; hence the multiset
contains precisely $2^n$ eigenvalues, counting multiplicities.
The normalized **global trace** Fourier moments are

```math
m_\ell^{\mathrm{open}}=
\begin{cases}
0,&\ell\text{ odd},\\
\prod_{j=2}^{n-1}\cos(\ell\theta_j),&\ell\text{ even}.
\end{cases}
```

For even periodic chains, use the commuting circuit
$C=\prod_j CZ_{j,j+1}$, which conjugates
$X_j$ into $K_j$ and satisfies
$\prod_j K_j=\prod_jX_j=P(0)$. Consequently the exact
finite eigenvalues and trace moments are

```math
\begin{aligned}
\lambda_k^{\mathrm{per}}&=
\left(\prod_{j=1}^n k_j\right)
\exp\left(-i\sum_{j=1}^n\theta_jk_j\right),\\
m_\ell^{\mathrm{per}}&=
\begin{cases}
\prod_j\cos(\ell\theta_j),&\ell\text{ even},\\
\prod_j[-i\sin(\ell\theta_j)],&\ell\text{ odd}.
\end{cases}
\end{aligned}
```

These are **finite sample trace laws** derived from the
commuting stabilizer algebra. They are not densities of an
infinite-dimensional system and not a universal CUE law.
Although the ideal open-chain trace spectrum has exact
$\pi$ symmetry, this global fact **alone** does not identify
an edge-localized mode; the anticommuting local operators
and the symmetry sectors supply the missing localization
information. A feature supported on a fixed number of
states can have negligible normalized trace mass as the
many-body dimension grows, so the global trace measure and
edge-resolved probes should not be conflated.

**Symmetry-preserving detuning.** A nonzero $\delta$
rotates each spin slightly away from the ideal $\pi/2$
pulse while preserving both $\mathbb Z_2$ symmetries.
It generally destroys the *exact* boundary anticommutation
and can produce finite-chain deviations from exact
$\pi$ pairing. This implementation does not infer
thermodynamic topological stability, localization or a
phase transition solely from such finite numerical data.
Bulk disorder in the commuting $\theta_j K_j$ terms,
by contrast, leaves ideal $\delta=0$ pairing exact.

**Efficient construction.** Let C be the product of
nearest-neighbor CZ gates (including the closing bond
for periodic chains). Then
$\exp(-i\sum_{j\in J}\theta_jK_j)=C[\prod_{j\in J}\exp(-i\theta_jX_j)]C$.
The generator implements the chronological CZ, X-rotation,
CZ and pulse gate sequence using existing
`um_circuit_apply`/`um_gate_apply`.
It supports dense assembly for $n\le10$ and local
matrix-free forward/adjoint actions for $n\le20$
with $O(N log N)$ state-vector arithmetic per RHS.
Its `spectral` mode uses the exact finite
multisets and moments **only** at $\delta=0$; otherwise
it reports an unknown finite analytic law without
substituting a numerical eigendecomposition.

This is a **distinct** spin-chain Floquet family from
the existing Ising, self-dual Ising, XXZ and BdG
generators. Its $\mathbb Z_2\times\mathbb Z_2$
cluster symmetry and exact open-edge operator algebra
are central to the construction.

**References and scope:** Potter, Morimoto and Vishwanath
[135] and Else and Nayak [136] establish the broader
Floquet symmetry-protected topology framework. Zhang,
Jiang, Deng and collaborators [137] demonstrate a
related driven spin-chain Floquet-SPT phenomenon
with $\mathbb Z_2\times\mathbb Z_2$ symmetry.
Our particular commuting-cluster/pulse specialization,
its fixed phase convention and finite trace moment
identities are **derived here**; they should not be
misidentified as exact verbatim equations or disorder
protocols taken from those papers.



## Moser–Veselov rigid-body increments in SO(n) — [`moser_veselov.m`](../moser_veselov.m)

**Application and mathematical distinction.** This family represents
the implicit **discrete Euler–Arnold equations** for a generalized
free rigid body, not merely a prescribed orthogonal path, an ordinary
polar factor or a Cayley transform of a chosen skew matrix. Its output
is a *dynamics-determined* real special-orthogonal body increment of
order $n\ge3$ [138,139]. An independent 2022 high-order
polar-decomposition variational integrator [140] is a distinct
construction and is **not implemented** by this generator.

**Exact convention for attitude, body increment and momentum.**
Let $J=J^T\succ0$ be the positive-definite generalized inertia
matrix, $M_k^T=-M_k$ the body angular momentum, $Q_k\in SO(n)$
the orientation matrix, and $h>0$ the physical step parameter.
In the body-increment convention of Moser–Veselov [138]
and its numerical study by Cardoso–Miraldo [139], solve

```math
\begin{aligned}
hM_k &= F_k^T J-JF_k, \qquad F_k\in SO(n),\\
M_{k+1}&=F_kM_kF_k^T,\\
Q_{k+1}&=Q_kF_k^T.
\end{aligned}
```

The literature commonly absorbs the factor $h$ into the
discrete momentum variable. The equation as printed in [139]
is $M=F^T J-JF$ **without an explicit h**. Here the user-supplied
M is a step-scaled physical-momentum convention and is multiplied
by h so that the increment approaches identity as $h\to0$.
This rescaling is a declared interface convention, not a
change to the published discrete matrix equation.

The orientation update **must** be compatible with the chosen
transpose: if one denotes $X=F^T$ instead, the implicit equation
is $hM=XJ-JX^T$, the momentum update is $M_+=X^TMX$,
and the attitude update is $Q_+=QX$. Using one formula from
each convention without the corresponding transposes is
incorrect.

**Unitary matrix structure and a discrete Noether invariant.**
Since $F_k\in SO(n)$, the returned numerical matrix
satisfies $F_k^TF_k=I$ and $\det F_k=1$. For real F,
complex eigenvalues lie in conjugate pairs on the unit circle,
and for odd n there is at least one eigenvalue $+1$.
In particular $SO(3)$ has a mandatory $+1$ eigenvalue,
so nontrivial higher-dimensional test matrices
($SO(4),SO(5),\ldots$) are included.

The angular momentum is advanced by orthogonal conjugation.
Thus its characteristic polynomial and all trace Casimirs
are invariant; in particular $\|M_{k+1}\|_F=\|M_k\|_F$.
The spatial momentum follows from the consistent pair of updates:

```math
Q_{k+1}M_{k+1}Q_{k+1}^T
=Q_kF_k^TF_kM_kF_k^TF_kQ_k^T
=Q_kM_kQ_k^T.
```

This is a strong independent test of the dynamics.
Crucially, **orthogonality of an orientation increment does not
mean the full map on the cotangent bundle is represented by
an orthogonal Euclidean matrix**. Variational
symplecticity of the full discrete mechanical flow [138]
and Euclidean orthogonality of its configuration factor are
different assertions. The generator does not compute or
claim a full phase-space Jacobian.

**Full Moser–Veselov/Manakov Lax invariants.**
Beyond the conjugacy invariants of M, the original paper [138]
derives an isospectral deformation of a matrix polynomial.
For an arbitrary (real or complex) scalar $z$, define

```math
L_k(z)=hM_k+zJ^2,\qquad A_k(z)=F_k-zJ.
```

The implicit equation and the momentum conjugation law imply
the **exact discrete Lax intertwining relation**

```math
L_{k+1}(z)\,A_k(z)=A_k(z)\,L_k(z).
```

Where $A_k(z)$ is nonsingular, the two matrices
$L_k(z)$ and $L_{k+1}(z)$ are similar; as an algebraic
identity, the entire characteristic polynomial in the
spectral variable is invariant from step to step.
This is a particularly useful **structured matrix pencil**
associated with the rigid-body map, but it is **not**
a generally unitary matrix or a unit-circle eigenvalue
problem. The `model` mode exposes
`lax(z)`, `next_lax(z)`
and `lax_intertwiner(z)`. The tests verify the
Lax residual and the first four trace moments of $L_k(z)$
at real and complex z, independently of the SO(n)
orientation invariants. The scaled variable hM is used
consistently with our physical-step convention.


**Numerical solution and branch policy.**
Solving $F^TJ-JF=hM$ is a genuine nonlinear matrix problem.
It need not admit an admissible solution for arbitrary inertia
and momentum, and solutions need not be globally unique
[138,139]. The implementation is intentionally restricted to
$n=3,\ldots,16$ for reproducible dense Newton arithmetic
and uses a **near-identity Cayley coordinate chart**:
write $X=F^T$ and $K^T=-K$, with

```math
X(K)=(I-K)^{-1}(I+K),\quad
R(K)=X(K)J-JX(K)^T-hM,
```

```math
dX=(I-K)^{-1}(dK)(I+X),\qquad
dR=(dX)J-J(dX)^T.
```

The first-order approximation
$KJ+JK\approx hM/2$ yields a Sylvester-equation
initialization, solved in an orthonormal eigenbasis of J.
Newton solves for independent strict-upper entries of K
using the displayed **analytic Jacobian** and backtracking.
Every trial Cayley iterate is in $SO(n)$, subject to roundoff,
and an unacceptably large residual, numerical singularity,
or inability to decrease the residual causes a
**NoSolution** exception. This is **our local analytic-Newton
solver**, not a claim of implementing the Bregman or Riemannian
gradient algorithms compared in [139].

A useful necessary solvability bound is
$\|hM\|_2\le2\|J\|_2$. The solver rejects inputs violating
this bound before attempting Newton; satisfying the bound
is not sufficient to guarantee a solution. At $M=0$, it
returns the **identity** branch, although other remote
solutions to the nonlinear equation may exist.

**Independent benchmarks.** For a principal-axis $SO(3)$
rotation $F=\mathop{\mathrm{diag}}(R(\phi),1)$ and diagonal
$J=\mathop{\mathrm{diag}}(j_1,j_2,j_3)$, the displayed
published implicit equation reduces algebraically to

```math
hM_{12}=(j_1+j_2)\sin\phi,\qquad
M_{13}=M_{23}=0.
```

The regression test takes $j_1=0.8,j_2=1.3,j_3=2$,
$h=0.24,\phi=0.18$ and recovers the independently
computed rotation with small residual. These
**numbers are a manufactured specialization** of the
published $SO(3)$ Moser–Veselov equation, **not**
a numerical table transcribed from [138] or [139].
Additional tests generate admissible $SO(4)$–$SO(6)$
momenta from independently constructed exact rotations
and recover those rotations without passing them to
the numerical generator. Twelve-step full trajectories
test $Q_kM_kQ_k^T$ and the Casimir $\|M_k\|_F$.
The standard default inertia is a deterministic
positive tridiagonal matrix, and the default body
momentum is a deterministic skew banded matrix;
neither is a claimed measured physical apparatus.

**Modes and trace spectral data.** The default
`matrix` output is the first $N\times N$ increment (here $N=n$),
`model` exposes the verified first step and
its momenta and residuals, and `trajectory`
returns the sequence of increments, body momenta,
orientations and equation residuals.
For generic parameters there is **no asserted exact finite
trace eigenphase multiset**, and `spectral`
reports it as analytically unknown, recording only
general $SO(n)$ conjugate pairing and the odd-n $+1$
eigenvalue. The zero-body-momentum near-identity branch
has $F=I$, so its finite normalized trace measure
is exactly one atom at $1$ (multiplicity n) and every
trace Fourier moment is 1.

Dense nonlinear solves do not admit a genuinely
matrix-free speedup merely by wrapping the computed
F in a matrix-vector handle; this family deliberately
reports `matrix_free_supported=false`.
The spectra of these configuration increments should
not be conflated with characteristic multipliers of
the **full** rigid-body phase-space flow, or with the
eigenvalues of its (generally nonorthogonal) Jacobian.

**References:** [138] is the original Moser–Veselov
discrete Euler–Arnold system, [139] analyzes solvability
and numerical solution of the nonlinear orthogonality-
constrained equation, and [140] presents a
**distinct** high-order polar/Galerkin variational
SO(n) method, retained for future independent work.



## Tensor-product 2D/3D finite-element wave propagators — [`tensor_fem_wave.m`](../tensor_fem_wave.m)

**Nonquantum application and independence from the 1D benchmark.**
This family advances the constant-coefficient anisotropic linear
wave equation on a rectangular two- or three-dimensional domain:

```math
u_{tt}-\sum_{a=1}^{d}c_a^2\,\partial_{x_a}^2u
+\omega_0^2u=0,\qquad
u|_{\partial\Omega}=0,\quad d=2,3.
```

The finite-element basis consists of **tensor products of
continuous piecewise-linear 1D functions**, giving conforming
bilinear elements in 2D and trilinear elements in 3D. The
one-dimensional meshes are uniform with $n_a$ unknown interior
Dirichlet nodal values and spacing $\Delta x_a=L_a/(n_a+1)$.
The total number of physical displacement nodes is
$m=\prod_{a=1}^{d}n_a$. The returned real orthogonal step has
matrix order $N=2m$, because the phase-space state includes
both displacement and its canonical momentum.

**Consistent mass, true Galerkin stiffness, and tensor order.**
For coordinate $a$, the assembled one-dimensional **consistent**
mass $M_a$ and Laplace stiffness $K_a$ are

```math
M_a=\frac{\Delta x_a}{6}
\begin{bmatrix}
4&1\\ 1&4&1\\ &\ddots&\ddots&\ddots\\ &&1&4
\end{bmatrix},
\qquad
K_a=\frac{1}{\Delta x_a}
\begin{bmatrix}
2&-1\\ -1&2&-1\\ &\ddots&\ddots&\ddots\\ &&-1&2
\end{bmatrix}.
```

For example in 2D, using the MATLAB first-index-fastest
flattening $i_x+n_x(i_y-1)$, the physical FE system has

```math
\begin{aligned}
M &= M_y\otimes M_x,\\
K &= c_x^2(M_y\otimes K_x)
    +c_y^2(K_y\otimes M_x)+\omega_0^2 M,\\
M\ddot q+Kq&=0,\qquad
\dot q=M^{-1}p,\quad\dot p=-Kq.
\end{aligned}
```

The 3D construction adds the $z$ coordinate as Kronecker
products in the same order; it is **not** a generic nonuniform
or noncommuting variable-coefficient finite-element matrix
(see issue #44). The tensorized assembly is the usual
Galerkin factorization, not a mass-lumped approximation
[103,141].

**Exactly known finite multiband frequencies.** The
orthonormal Dirichlet sine vectors jointly diagonalize each
1D $M_a,K_a$, with generalized eigen-symbols

```math
\begin{aligned}
\theta_{a,j}&=\pi j/(n_a+1),\\
m_{a,j}&=\frac{\Delta x_a}{3}(2+\cos\theta_{a,j}),\\
\kappa_{a,j}&=\frac{2}{\Delta x_a}(1-\cos\theta_{a,j}),\\
\omega_{\boldsymbol j}^{2}
 &=\omega_0^2+
   \sum_{a=1}^{d}c_a^2\frac{\kappa_{a,j_a}}{m_{a,j_a}}.
\end{aligned}
```

The normal-mode multiplicities are retained exactly by
enumerating every multi-index $\boldsymbol j$; e.g., for a
square 2D grid with identical coordinate lengths and speeds,
the modes $(j,k)$ and $(k,j)$ have the same frequency.
Rectangular anisotropy removes many such coincidences.
An orthonormal Kronecker DST basis $V$ diagonalizes the
physical mass/stiffness simultaneously. Under the
positive energy change of variables the physical wave
step becomes a direct sum of planar rotations.

**Gauss–Legendre and Verlet unitary eigenphases.** Write
$x_{\boldsymbol j}=h\omega_{\boldsymbol j}$ and
$R_s(z)=P_s(z)/P_s(-z)$, the diagonal $[s/s]$ Padé
approximation to the exponential, with $s=1,2,3$
(midpoint, Gauss order four, Gauss order six respectively).
The exact **finite unit eigenphase multiset** is

```math
\big\{
z_{\boldsymbol j},\overline{z_{\boldsymbol j}}:
\boldsymbol j\in\prod_a\{1,\ldots,n_a\}
\big\},\qquad
z_{\boldsymbol j}=R_s(-ih\omega_{\boldsymbol j}).
```

For velocity Verlet the corresponding phase branch is

```math
z_{\boldsymbol j}=
1-\frac{x_{\boldsymbol j}^2}{2}
-i x_{\boldsymbol j}\sqrt{1-\frac{x_{\boldsymbol j}^2}{4}},
\qquad h\max_{\boldsymbol j}\omega_{\boldsymbol j}<2.
```

The strict stability inequality makes the **modified energy
positive definite** and ensures distinct unit-circle
conjugate pairs. It is actively enforced; outside this regime
the generator raises a stability error rather than returning
a matrix incorrectly labeled unitary. Midpoint and all
Gauss choices preserve the **original** positive quadratic
energy, while Verlet preserves a method-dependent
**modified** quadratic form. Gauss phase wrapping is
allowed; no artificial principal-argument restriction is
placed on $R_s$.

The normalized trace spectral measure assigns equal weight
$1/(2N)$ to each eigenvalue (including repeated eigenvalues).
Its exact finite Fourier moments for integer $\ell$ are

```math
m_\ell=\frac1{2N}
\sum_{\boldsymbol j}
\left(z_{\boldsymbol j}^{\ell}
+\overline{z_{\boldsymbol j}}^{\ell}\right).
```

These identities produce the `spectral` mode by reusing
`um_geometric_spectral`, without calling eig(U).
**No limiting density or van Hove singularity is asserted.**
Neither mesh refinement with fixed domain nor a thermodynamic
large-domain limit is specified as a measure-theoretic
convergence statement here.

**Real orthogonal representative versus physical symplectic
step.** With $M,K$ the original FE matrices, let
$S_h$ denote the physical time map on $[q;p]$. In the
common orthonormal DST eigenbasis, the positive quadratic
metric $G=T^TT$ is

```math
G=\begin{cases}
\mathop{\mathrm{diag}}(K,M^{-1}),
&\text{midpoint/Gauss},\\[2pt]
\mathop{\mathrm{diag}}\!\left(
K-\tfrac{h^2}{4}KM^{-1}K,\ M^{-1}\right),
&\text{stable Verlet}.
\end{cases}
```

Then $S_h^T G S_h=G$, and
$U=T S_h T^{-1}$ is real orthogonal. Both $U$ and
the physical $S_h$ preserve the canonical symplectic
form with its appropriate coordinates; the positive
normalization T is **not generally symplectic**.
The `pencil` mode exposes $A=G S_h$ and
$B=G$, so $B^{-1}A=S_h$ is circle-spectral, though not
normally Euclidean orthogonal. It also returns both
assembled Kronecker FE matrices, metric G and T.

For Verlet, elimination of momentum from the physical
position recurrence gives the **symmetric palindromic
matrix quadratic**

```math
P(z)=z^2 M+z(h^2K-2M)+M.
```

The `polynomial` mode returns three sparse
physical coefficient matrices, not the unitary U.
The structure is a genuine mass-weighted QEP and is
not replaced by an unjustified identity-mass companion.
The unit-circle claim requires the strict stable
Verlet parameter range.

**Matrix-free performance and validation.** The
`operator` mode stores $O(m)$ complex
phase samples and applies forward/adjoint
multidimensional **orthonormal DST-I** transforms
one coordinate at a time using the existing FFT
kernel `um_dst1_apply`. For fixed dimension
$d\le3$, arithmetic complexity is $O(N\log N)$ (since $N=2m$)
per right-hand side; no $N$-by-$N$ dense
matrix is formed. Dense matrix/pencil output is
restricted to $m\le128$; explicit sparse Verlet
quadratic output to $m\le10^4$. The tests compare
small-grid mass/stiffness with an **independent
entry-by-entry FE assembly**, check the generalized
eigenvalues and all modes, energy and symplectic
metrics, a separate physical Padé step, the
palindromic QEP companion roots, degeneracies
and anisotropic splitting, phase wrapping and
large 2D/3D matrix-free adjoint reversibility.

**References:** Hairer–Lubich–Wanner [101] for
symplectic Gauss–Legendre/Verlet time integrators,
Cockburn–Du–Sánchez [103] for finite-element wave
Hamiltonian semidiscretization, and
Loveland–Valseth–Lukac–Dawson [141] for
tensor-product finite-element assembly.
The particular exact DST symbol and eigenphase
formulas are elementary deductions from this
special constant-coefficient rectangular model.



## Nonuniform finite-element wave equations with noncommuting mass and stiffness — [`nonuniform_fem_wave.m`](../nonuniform_fem_wave.m)

**Why this family differs from existing 1D/DST and tensor-product
wave generators.** This is a genuine 1D heterogeneous
piecewise-linear finite-element model with a **nonuniform mesh**
and/or **spatially variable wave speed**, so its positive-definite
consistent physical mass $M$ and stiffness $K$ generically
satisfy $MK\ne KM$. A uniform consistent-mass Dirichlet chain has
commuting Toeplitz mass and stiffness and exact DST-I frequencies,
but a nonuniform mesh or variable material coefficients destroys
that simultaneous Euclidean diagonalization. Applying the
homogeneous DST formulas to the heterogeneous case would be wrong.

The physical PDE model is
$u_{tt}-\partial_x(c(x)^2u_x)+\omega_0^2u=0$ on $(0,L)$ with
homogeneous Dirichlet boundaries. There are $n\ge2$ interior
continuous piecewise-linear nodal degrees of freedom. Let the
$n+1$ element lengths be $\ell_e>0$, $\sum_e\ell_e=L$, and
piecewise-constant positive wave speeds $c_e$. The **exact
Galerkin element matrices** are

```math
M_e=\frac{\ell_e}{6}
\begin{bmatrix}2&1\\1&2\end{bmatrix},\qquad
K_e=\frac{c_e^2}{\ell_e}
\begin{bmatrix}1&-1\\-1&1\end{bmatrix}.
```

After standard assembly and elimination of Dirichlet boundary
unknowns, include the onsite term $\omega_0^2 M$ in $K$. Both
physical matrices are sparse, real symmetric, and positive
definite. The generator supplies reproducible nonuniform element
lengths and spatially varying speeds by default, or accepts
user-specified positive element vectors through `opts`.
The `coefficients` mode returns the assembled matrices,
geometry, speed array and measured Frobenius commutator norm.
The independently assembled test confirms the mesh physics
and $\|MK-KM\|_F>0$ for the standard heterogeneous cases.

**Canonical Hamiltonian form and the two separate structures.**
For physical displacement $q$ and canonical momentum
$p=M\dot q$, use $x=[q;p]$. Define

```math
\begin{aligned}
H(x)&=\tfrac12 q^T Kq+\tfrac12 p^TM^{-1}p,\\
A_0&=\begin{bmatrix}0&M^{-1}\\-K&0\end{bmatrix},
\quad
G=\begin{bmatrix}K&0\\0&M^{-1}\end{bmatrix}\succ0,\\
J_c&=\begin{bmatrix}0&I\\-I&0\end{bmatrix}.
\end{aligned}
```

For this **constant** Hamiltonian matrix, two separate
skew-adjointness identities hold:

```math
A_0^TG+GA_0=0,\qquad A_0^TJ_c+J_cA_0=0.
```

Midpoint and the 2-/3-stage Gauss–Legendre methods are
symplectic Runge–Kutta integrators; their rational stability
maps $S_h=R_s(hA_0)$ preserve all quadratic invariants of
this linear conservative problem [101,103]. In particular,

```math
S_h^TGS_h=G,\qquad S_h^TJ_cS_h=J_c.
```

These statements concern the **physical canonical phase-space
step** $S_h$. They do not say that $S_h$ is Euclidean
orthogonal, nor that the mass and stiffness commute.

**Positive-metric orthogonal normalization.** Take Cholesky
factors $R_K^TR_K=K$ and $R_M^TR_M=M$ and define

```math
T=\mathop{\mathrm{diag}}(R_K,R_M^{-T}),\qquad
G=T^TT,\qquad U=TS_hT^{-1}.
```

The returned default real $2n$-by-$2n$ U satisfies
$U^TU=I$, so it is also unitary and has all its
eigenvalues on the unit circle. Its **correct transported
symplectic form** is

```math
\widetilde J=T^{-T}J_cT^{-1},\qquad
U^T\widetilde J U=\widetilde J.
```

There is **no generic claim** that $U^TJ_cU=J_c$:
the chosen energy whitening T is not generally
symplectic. The tests include an explicit numerical
noncommuting counterexample exhibiting the discrepancy,
as well as both true invariants of the physical S.
The `pencil` mode exposes $A=GS_h$, $B=G$,
so $B^{-1}A=S_h$ is a **positive-metric circle-spectral
pencil**, satisfying $A^TB^{-1}A=B$. It is not
necessarily a strict Euclidean-unitary pencil.

**Sparse implicit midpoint and Gauss higher-order actions.**
The canonical generator $A_0$ contains an inverse of M,
but forming $M^{-1}$ explicitly would destroy sparsity.
Instead introduce the **sparse generalized first-order
form**

```math
E=\mathop{\mathrm{diag}}(M,I),\quad
L_0=\begin{bmatrix}0&I\\-K&0\end{bmatrix},\quad
A_0=E^{-1}L_0.
```

For $s=1,2,3$ the Gauss stability function has a
diagonal Padé representation, $R_s(z)=P_s(z)/P_s(-z)$,
with $P_1(z)=1+z/2$,
$P_2(z)=1+z/2+z^2/12$, and
$P_3(z)=1+z/2+z^2/10+z^3/120$.
Factor $P_s(z)=\prod_{j=1}^s(1+\alpha_j z)$
using its exact-degree scalar polynomial roots. One
Gauss step is the commuting product of sparse
generalized Cayley factors

```math
S_h=\prod_{j=1}^{s}
\left(E-\alpha_j hL_0\right)^{-1}
\left(E+\alpha_j hL_0\right).
```

The `operator` mode caches sparse LU
factorizations of the needed generalized factors,
and Cholesky triangular factors for the *energy*
normalization; the adjoint action reverses the
factors and uses their inverses. Complex-conjugate
factors combine to a real map, with no false
interpretation of their intermediate complex matrices
as physical time steps. No dense $M^{-1}$ is
constructed in operator mode; no FFT/DST is
assumed for the heterogeneous mesh. Sparse LU
can still fill in, and neither preprocessing
nor arithmetic is promised to be linear in n.
Dense `matrix`/`pencil` modes are limited
to $n\le128$, while operator and sparse coefficient
modes extend to $n\le5000$ with a numerical SPD
conditioning guard. This is a sparse-solve
benchmark, not a matrix-free FFT diagonalization.

**Finite trace spectra, numerical conditioning, and uniform
consistency check.** In a general heterogeneous problem,
a generalized $K v=\omega^2 M v$ calculation can
approximate the frequencies numerically, and the
Padé map then gives numerical unit-circle phases.
This is **not** a known analytical closed form:
`spectral` returns a standardized
`finite.kind='unknown'` record. By contrast,
if all element lengths equal $L/(n+1)$ and all
wave speeds equal the nominal c, the original
commuting 1D finite-element benchmark is recovered.
Its generalized sine-mode frequencies are

```math
\omega_j^2=\omega_0^2+
\frac{(2c^2/\Delta x)(1-\cos\theta_j)}
     {(\Delta x/3)(2+\cos\theta_j)},\qquad
\theta_j=\frac{\pi j}{n+1}.
```

The uniform `spectral` branch reuses the
standard exact finite trace eigenvalue/moment
descriptor `um_geometric_spectral`, rather
than assuming an analytical law for the heterogeneous
case. All modes obey real conjugate eigenphase
pairing. The tests verify that both zero
inhomogeneity parameters recover the uniform
mass/stiffness and exact finite phases.

As $G$ becomes severely ill conditioned, the
normalization $U=TS_hT^{-1}$ is numerically
sensitive even when the underlying physical
dynamics is well defined. A conservative
`condest` guard for M and K (threshold
$10^{12}$) rejects extreme grading instead of
silently reporting a machine-precision orthogonal
matrix. Tests include an explicit near-degenerate
element mesh and a controlled exception.

**Sources:** Hairer, Lubich and Wanner [101]
establish symplectic Gauss collocation, rational
linear-flow maps, and quadratic-invariant
preservation; Cockburn, Du and Sánchez [103]
discuss finite-element Hamiltonian
semidiscretizations and symplectic time marching.
The particular sinusoidal nonuniform nodal
geometry and elementwise speed pattern here are
deterministic reproducible *specializations*, not
published experimental data. Equations and
identities above follow directly from assembled
Galerkin matrices and the cited Gauss methodology.



## True hyperbolic $H^n$ tangent-bundle holonomy — [`hyperbolic_tangent_holonomy.m`](../hyperbolic_tangent_holonomy.m)

**Scope.** Unlike `hyperbolic_triangle_holonomy`, which
returns the tangent rotation of an $H^2$ triangle embedded by
an explicitly artificial identity direct sum, this family
computes the **intrinsic** Levi–Civita holonomy of a closed
geodesic polygon in genuine hyperbolic $H^n$ for $n=2,\ldots,40$.
The output is an $N$-by-$N$ real $SO(N)$ matrix in an
orthonormal tangent basis at the starting vertex.
For $n\ge3$, non-coplanar polygons can produce
noncommuting, non-Abelian holonomy.

**Ambient indefinite geometry.** The hyperboloid is

```math
H^n=\{x\in\mathbb R^{n+1}:x^T J_Lx=-1,\ x_0>0\},
\qquad J_L=\mathop{\mathrm{diag}}(-1,I_n).
```

Here index 0 denotes the Lorentzian time coordinate.
The tangent metric is positive definite:
$T_xH^n=\{v:x^TJ_L v=0\}$ with
$\langle v,w\rangle_x=v^TJ_Lw$.
A useful globally smooth **canonical tangent frame**
is generated by the Lorentz boost carrying $(1,0,\ldots,0)$
to $x=(x_0,u)$:

```math
E(x)=\begin{bmatrix}
u^T\\I_n+\frac{uu^T}{1+x_0}
\end{bmatrix},\qquad
E(x)^TJ_LE(x)=I_n,\quad E(x)^TJ_Lx=0.
```

The coordinate point $x$ and its frame
$[x,E(x)]$ form a proper orthochronous Lorentz
matrix, which is **not** Euclidean orthogonal.
The generator accepts a real $(n+1)$-by-$m$
array of future-pointing Lorentz-unit vertices
($3\le m\le24$), and automatically closes
the final edge back to its first vertex.

**An exact geometric leg, not a retraction.**
For nonidentical vertices $x,y$, the geodesic
distance is
$d(x,y)=\mathop{\mathrm{arcosh}}(-x^TJ_Ly)$.
Let $c=-x^TJ_Ly=\cosh d>1$.
The exact Levi–Civita parallel transport of any
tangent vector $v\in T_xH^n$ is

```math
\mathcal P_{x\to y}v=
v+\frac{y^T J_Lv}{1+c}(x+y).
```

It maps $T_xH^n$ to $T_yH^n$ and preserves
the positive tangent inner product.
It transports the geodesic velocity to the
endpoint geodesic velocity. The formula appears
in equivalent logarithm-map form in hyperboloid
parallel-transport literature such as Nagano
and collaborators [142].

The square **ambient Lorentz extension** of
this tangent isometry is uniquely determined
by sending $x\mapsto y$ and
$E(x)\mapsto\mathcal P_{x\to y}E(x)$:

```math
B_{xy}=
[\,y,\mathcal P_{x\to y}E(x)\,]\,
J_L[\,x,E(x)\,]^TJ_L.
```

Then $B_{xy}^TJ_L B_{xy}=J_L$ and
$\det B_{xy}=1$. It is in $SO^+(1,n)$,
**not** $O(n+1)$. No Lorentz matrix is
incorrectly classified as a unitary matrix.

**Closed-loop tangent holonomy.** For vertices
$x_1,\ldots,x_m$, let $x_{m+1}=x_1$.
Compose the proper Lorentz boosts in path order:

```math
B_{\gamma}
=B_{x_m x_1}\cdots B_{x_2x_3}B_{x_1x_2},
\qquad
U_\gamma=E(x_1)^TJ_L B_{\gamma}E(x_1).
```

The closed product fixes $x_1$ and its restriction
to the tangent fiber is in **$SO(n)$**.
Changing the orthonormal gauge at the starting
vertex from $E$ to $EW$ conjugates the returned
matrix $U_\gamma\mapsto W^TU_\gamma W$;
its eigenphase multiset is gauge-independent.
Reversing loop orientation yields $U_\gamma^{-1}$.

**Exact finite triangle spectra.** Any three
hyperboloid points span a totally geodesic copy
of $H^2$, even inside high-dimensional $H^n$.
For a nondegenerate geodesic triangle with
hyperbolic side lengths $a,b,c$ and corner
angles $A,B,C$, the positive angle defect is

```math
\begin{aligned}
\cos A&=\frac{\cosh b\cosh c-\cosh a}
             {\sinh b\sinh c},\\
\mathcal A&=\pi-(A+B+C)>0.
\end{aligned}
```

The hyperbolic curvature is $-1$, so the
Levi–Civita holonomy on the geodesic triangle's
2-plane is a signed rotation of magnitude
$\mathcal A$. The remaining $n-2$ normal
tangent directions are fixed. Therefore the
exact **finite equal-weight trace** eigenvalues
and integer Fourier moments are

```math
\begin{aligned}
\{\lambda_j\}_{j=1}^n
&=\{e^{i\mathcal A},e^{-i\mathcal A},
  \underbrace{1,\ldots,1}_{n-2}\},\\
m_\ell&=\frac{n-2+2\cos(\ell\mathcal A)}{n}.
\end{aligned}
```

The `spectral` mode computes the
defect directly from hyperbolic distances,
**without** a numerical eigenvalue decomposition.
For general multi-plane polygons of four or
more vertices, no generic closed finite
eigenphase formula is asserted; the analytical
trace reference is marked `unknown`.
Although an $SO(3)$ matrix always has one
fixed tangent direction, an arbitrary spatial
curve need not lie in any single $H^2$.

**Non-Abelian $H^3$ and reproducibility checks.**
The default $H^3$ polygon has four vertices
and extends in all three spatial directions;
it is not a planar $H^2$ triangle padded
with identity directions. Tests construct
two independent triangles in intersecting
totally geodesic 2-planes of $H^3$; their
tangent holonomy matrices do not commute.
They also verify reverse-path inversion,
$SO(n)$ determinant and orthogonality,
proper Lorentz invariants for *every leg*,
ambient basepoint fixing and tangent-frame
gauge covariance.

For an **independent differential validation**
the tests solve the Levi–Civita connection ODE
along the exact hyperboloid geodesic
$\gamma(t)$ via high-order Runge–Kutta:

```math
\begin{aligned}
\gamma(t)&=
\frac{\sinh[(1-t)d]}{\sinh d}x+
\frac{\sinh(td)}{\sinh d}y,\\
\dot V(t)&=
\langle\dot\gamma(t),V(t)\rangle_L\,\gamma(t).
\end{aligned}
```

The result agrees with the closed leg formula
and transports the geodesic velocity correctly.
For $n=2$, mapping hyperboloid vertices to
Poincaré disk coordinates $z=(x_1+i x_2)/(x_0+1)$
recovers the existing 2D signed-area
`hyperbolic_triangle_holonomy`
*as a genuinely intrinsic two-dimensional result*,
not a padded higher-dimensional surrogate.
The positive triangle area is independently
checked using hyperbolic L'Huilier's formula.

The `model` mode returns the intrinsic
$SO(n)$ holonomy, per-leg positive-definite
tangent rotations, each Lorentz boost, all
canonical tangent frames and their product.
A dense matrix suffices at these orders,
so no unsupported matrix-free speedup
is claimed. The construction is useful as
a benchmark for **non-Abelian geometric phases,
Lorentzian parallel transport and connection
holonomy**, not for generic Lorentz-group
Euclidean-unitary eigenvalue testing.

**References:** [142] Nagano et al.,
ICML 2019 (hyperboloid tangent transport);
[113] Needham, *Visual Differential Geometry
and Forms* (negative curvature, Gauss–Bonnet
and geodesic holonomy). The formulas above
explicitly fix the curvature, sign, metric
and fiber-basis conventions.



## Gauss–Hermite finite harmonic-oscillator spectral flow — [`hermite_oscillator.m`](../hermite_oscillator.m)

**What is and is not a fractional Fourier transform.** The generator
`hermite_oscillator(n,theta,offset,mode)`
builds a finite, **exactly unitary** spectral evolution from
the orthogonal Gauss–Hermite quadrature-evaluation basis
already provided by `gauss_hermite_transform(n)`.
It is an independent matrix family, not a second copy of
the quadrature basis Q, because its input and output are
**weighted values on the quadrature grid**, and it
varies continuously as a unitary one-parameter group.
It is not a claim to implement the Candan–Kutay–Özaktas
**DFT-compatible discrete fractional Fourier transform** [144]
(issue #38), nor Campos–Rico-Melgoza–Chávez's different
asymptotic FFT-based quadrature algorithm [143].

The $n$ real quadrature nodes $x_j$ are the increasing
eigenvalues of the real symmetric Jacobi matrix with
offdiagonal entries $\sqrt{k/2}$, $k=1,\ldots,n-1$,
corresponding to the normalized Gaussian measure
$e^{-x^2}dx/\sqrt{\pi}$. Write $w_j>0$ for quadrature
weights normalized so that $\sum_jw_j=1$, and
$p_k$ for the orthonormal physicists' Hermite
**polynomials**. The existing
Golub–Welsch construction [104] gives

```math
Q_{j,k+1}=\sqrt{w_j}\,p_k(x_j),\qquad
Q^TQ=I_n,\quad k=0,\ldots,n-1.
```

Finite Gaussian quadrature exactly integrates products
of these $n$ polynomials, so Q is an orthogonal matrix,
not a loose numerical truncation of a continuous integral
kernel. The **quadrature-weighted nodal** coordinates
are the chosen physical representation, rather than
samples on an equispaced FFT grid.

**Strict finite unitary spectral family.** Let
$\theta\in\mathbb R$ and $\delta\in\{0,\tfrac12\}$,
where $\delta$ is the optional zero-point energy.
The complex finite matrix is

```math
\begin{aligned}
E_\delta&=\mathop{\mathrm{diag}}(\delta,1+\delta,
                   \ldots,n-1+\delta),\\
H_\delta&=QE_\delta Q^T,\\
U_{\theta,\delta}&=Q
 \mathop{\mathrm{diag}}_{k=0}^{n-1}
 \left(e^{-i\theta(k+\delta)}\right)Q^T
                 =e^{-i\theta H_\delta}.
\end{aligned}
```

The **energy matrix** $H_\delta$ is real
symmetric, and the sample-space evolution is
unitary. The output has order n, not twice
the number of quadrature nodes. It obeys the
exact algebraic identities

```math
\begin{aligned}
U_{0,\delta}&=I,\qquad
U_{\theta,\delta}^*=U_{-\theta,\delta},\\
U_{\theta,\delta}U_{\varphi,\delta}
 &=U_{\theta+\varphi,\delta},\\
\frac{d}{d\theta}U_{\theta,\delta}
 &=-iH_\delta U_{\theta,\delta}.
\end{aligned}
```

With $\delta=0$ (the default), $\theta$ is
$2\pi$-periodic. Including the continuous
harmonic oscillator's ground-state contribution
$\delta=1/2$ makes the unitary
$2\pi$-**antiperiodic** and $4\pi$-periodic.
The half-integer shift is a scalar global phase;
it does not alter normalized power-spectrum
multiplicities but does shift every eigenphase.

**Exactly finite trace phase measure.** The
finite trace spectrum is

```math
\left\{e^{-i\theta(k+\delta)}
          :k=0,\ldots,n-1\right\},\qquad
m_\ell=\frac1n\sum_{k=0}^{n-1}
               e^{-i\ell\theta(k+\delta)},
\qquad \ell\in\mathbb Z.
```

These equal-weight trace moments are computed
from the **specified Hermite energy levels**,
not by calling a numerical eigensolver on U.
Coincident/aliased phases retain their finite
multiplicities (for example $\theta=\pi/2$).
No thermodynamic or continuous spectral
density is claimed for this finite family.
The `spectral` mode records exact
eigenvalue multiset and all integer Fourier
moments in the common reference format.

**Do not confuse the Fourier endpoint.** At
$\theta=\pi/2$, one obtains the finite diagonal
oscillator-phase map in a **Gauss–Hermite
quadrature-evaluation basis**. The normalized
DFT instead has a distinct periodic equispaced
sampling space and DFT-compatible Hermite-like
eigenvectors [144]. Thus in general

```math
U_{\pi/2,0}\ne F_n,\qquad
F_n(j,k)=n^{-1/2}e^{-2\pi i\,jk/n}.
```

The tests verify this *non-equality* for several
nontrivial even and odd matrix orders. In particular
the present model is neither a chosen branch
of $F_n^\alpha$ nor the exactly DFT-compatible
construction requested by issue #38.

**Computation and diagnostic modes.** The
`basis` mode returns the real Q,
increasing Gauss–Hermite nodes, normalized
positive weights, and energy levels.
The `hamiltonian` mode returns
the symmetric $H_\delta$;
`matrix` returns the complex U.
The existing dense Golub–Welsch eigensolver
costs $O(N^3)$ to initialize Q.
Subsequent dense matrix assembly and the
first spectral basis construction likewise
have nontrivial cubic cost; no unwarranted
FFT or matrix-free action is advertised.
The generator is limited to $1\le n\le256$
for reproducible moderate-order benchmarks.
Tests cover quadrature normalization and the
Hermite three-term recurrence, strict unitarity,
the group law, both periodicity conventions,
the Hamiltonian matrix exponential and
derivative, and exact trace moments.

**References and precise scope:** [104]
Golub–Welsch for Gauss-quadrature recurrence
eigenvectors; [143] Campos–Rico-Melgoza–Chávez
for a distinct fast FrFT quadrature; [144]
Candan–Kutay–Özaktas for the DFT-compatible
discrete FrFT. The **particular finite
sample-space spectral construction and
global-phase convention above are declared
here** and do not purport to reproduce
the published fast algorithms.



## Spectral-measure computation: Colbrook and collaborators

References [145–148] supply spectral-measure methods and data-driven Koopman benchmarks, rather than four automatically new finite unitary generators.

- **Resolvent-based density methods [145].** The self-adjoint spectral-measure calculations, including graphene/tight-binding and Dirac examples, provide methodological context and potential *Hermitian-to-unitary* test constructions (e.g., explicitly specified Cayley transforms or Hermitian propagators); their raw operators are **not** unitary matrices.
- **Data-driven Koopman spectral measures [146].** For invariant measures, composition Koopman operators are isometries; they are unitary on the full $L^2$ space when the dynamics is invertible modulo measure-zero sets. ResDMD/EDMD outputs are **not generally unitary finite matrices**, even if the underlying Koopman operator is unitary. Examples in [146] include circle rotations, tent and Gauss maps, pendulum systems and Lorenz dynamics; admissible finite unitary test matrices require an independently justified exact invariant subspace, a bijective finite-state discretization, or a genuinely structure-preserving approximation.
- **Spectral-type benchmarks in OPUC [147].** Section 6.3 treats infinite CMV matrices and their **vector** spectral measures, including Rogers–Szegő (wrapped Gaussian) and constant-Verblunsky Geronimus cases with band support and, in some parameter regimes, an additional point mass. The repository already provides `rogers_szego` and `geronimus_cmv` finite generators. A finite unitary CMV truncation has an atomic spectral measure, so it must **not** be labelled as exactly possessing the continuous or mixed spectral type of its infinite counterpart. Those examples instead motivate carefully verified finite-size convergence and smoothing benchmarks.
- **Measure-preserving finite Koopman approximants [148].** mpEDMD uses a data-dependent **orthogonal Procrustes** construction to enforce measure preservation in a finite-dimensional approximation. In an orthonormalized/Gram-whitened dictionary basis the resulting matrix can be unitary; a raw nonorthogonal-coordinate representation is generally unitary only in its associated positive metric. This suggests a distinct data-derived finite matrix class, provided matrix assembly, conditioning, and the metric are exposed and tested, rather than presenting generic EDMD or ResDMD matrices as unitary.

Keep **vector spectral measures** $\langle v,E_U(\cdot)v\rangle$ distinct from the **normalized trace eigenphase measure** $n^{-1}\sum_j\delta_{\theta_j}$; in general the former is not the latter. Also distinguish certified infinite-operator spectral limits from the finite sample spectrum and from smoothing-induced densities.

## Intrinsic Grassmann tangent Levi–Civita holonomy — [`grassmann_tangent_holonomy.m`](../grassmann_tangent_holonomy.m)

**Different bundles and different matrix sizes.** Three
important unitary constructions on a Grassmann manifold
must not be confused: `grassmann_geodesic` returns
an **ambient $N\times N$** unitary ($N=n$) implementing a prescribed
subspace geodesic; `grassmann_holonomy` computes
the **tautological rank-r U(r)** Berry/Wilczek–Zee
transport around geodesic polygons; this new
`grassmann_tangent_holonomy` returns the
**Levi–Civita transport of arbitrary tangent vectors**.
Its square matrix has order

```math
d=r(n-r),
\qquad
\dim_{\mathbb R}\mathrm{Gr}_{\mathbb R}(r,n)=d,
\qquad
\dim_{\mathbb C}\mathrm{Gr}_{\mathbb C}(r,n)=d.
```

For **real** Grassmannians the tangent representation
is a real $d$-by-$d$ orthogonal matrix. For
**complex** Grassmannians with the canonical Kähler
metric, Levi–Civita transport commutes with
the complex structure and gives a complex
$d$-by-$d$ unitary matrix in complex linear
tangent coordinates. The ambient dimension
and tautological fiber rank are **not**
substitutes for this tangent order.

**Horizontal representatives and orthonormal frames.**
A point is an equivalence class of real or complex
$n$-by-$r$ orthonormal frames $Y$ under
$Y\mapsto YR$, $R\in O(r)$ or $U(r)$.
Let $Q_Y=\mathop{\mathrm{null}}(Y^*)$ be an
orthonormal complement and write any horizontal
tangent as

```math
\Xi=Q_Y B,\qquad
Y^*\Xi=0,\quad
B\in\mathbb F^{(n-r)\times r},\qquad
\langle\Xi,\Eta\rangle=\mathrm{tr}(B^*C).
```

Here $\mathbb F=\mathbb R$ or $\mathbb C$.
The ordered coordinates `vec(B)` use MATLAB
column-major vectorization and form an orthonormal
(real or complex) tangent basis. Two vertex frames
may be joined along their **unique minimizing
geodesic** when all principal angles are strictly
below $\pi/2$; a singular adjacent overlap
is rejected as the Grassmann cut locus.

**Exactly source-matched Grassmann geodesics.**
For adjacent frames $Y,Z$ compute the polar alignment

```math
R=\mathop{\mathrm{polar}}(Z^*Y),\quad
Z_a=ZR,\quad C=Y^*Z_a=C^*>0.
```

Let $C=V\mathop{\mathrm{diag}}(\cos\sigma_j)V^*$.
Then the horizontal initial velocity and its
constant skew(-Hermitian) ambient generator are

```math
\begin{aligned}
\Delta&=(Z_a-YC)\,
 V\mathop{\mathrm{diag}}\!\left(
    \frac{\sigma_j}{\sin\sigma_j}\right)V^*,\\
\Omega&=\Delta Y^*-Y\Delta^*,\qquad
A(t)=e^{t\Omega},\\
Y(t)&=A(t)Y,\qquad Y(1)=Z_a.
\end{aligned}
```

The removable value $\sigma/\sin\sigma=1$
is used at zero, and its small-angle series
is used numerically. The formula is the
Grassmann horizontal principal-angle geodesic
and logarithm of Edelman–Arias–Smith [97] and
Bendokat–Zimmermann–Absil [111]. The ambient
lift $A(1)$ is unitary, but it is **not** the
tangent holonomy recorded in this family.

**Exact Levi–Civita parallel transport.**
For a horizontal tangent $\Xi(0)$ at Y,
define the vector field $\Xi(t)=A(t)\Xi(0)$.
Since $\Omega$ has only off-diagonal blocks
between $\mathrm{span}Y$ and its
orthogonal complement, its derivative
$\dot\Xi(t)=\Omega\Xi(t)$ lies in
$\mathrm{span}Y(t)$. Hence its
horizontal covariant derivative is zero:

```math
\left(I-Y(t)Y(t)^*\right)\dot\Xi(t)=0.
```

This directly proves **Levi–Civita parallel
transport**, not merely projection or
first-order vector transport. It preserves
the tangent Frobenius metric and carries
the geodesic velocity to the endpoint
geodesic velocity.

To express the transported tangent in
the **target's own fixed orthonormal
frames** $(Z,Q_Z)$, use

```math
\begin{aligned}
S&=Q_Z^*A(1)Q_Y\in U(n-r),\\
B_{\mathrm{end}}&=S\,B_{\mathrm{start}}\,R^*,\\
\mathrm{vec}(B_{\mathrm{end}})
 &=\left(\overline R\otimes S\right)
   \mathrm{vec}(B_{\mathrm{start}}).
\end{aligned}
```

For real data replace $U$ by $O$.
Thus each **intrinsic tangent link**
$T=\overline R\otimes S$ is a square
$d$-by-$d$ unitary (real orthogonal when
appropriate). The closed polygon's
Levi–Civita holonomy is the ordered
product of these tangent links, **not**
a product of ambient A matrices or
the product of the rank-r alignment R's.

**Gauge covariance, inverse paths and validation.**
Changing local orthonormal Y frames at the
vertices transforms tangent-coordinate
matrices covariantly. A closed-loop change
at its base is a unitary similarity,
preserving the eigenphases. Reversing a
polygon gives the inverse holonomy. Tests
verify all horizontal constraints and
isometries, inverse-path composition,
base-fiber gauge covariance, and the
parallelism of geodesic velocity.
An **independent** connection ODE on each
geodesic,

```math
\dot\Xi(t)=
-Y(t)\big(\dot Y(t)^*\Xi(t)\big),
```

is integrated with fourth-order Runge–Kutta
and compared to $A(t)\Xi(0)$.
This numerical differential check is
separate from the algebraic unitary
product that defines the generator.

**Complex projective geometry and scalar versus
matrix-valued holonomy.** For
$\mathbb{CP}^{m}=\mathrm{Gr}_{\mathbb C}(1,m+1)$,
the tautological bundle is a **complex
line** with a scalar Berry phase, whereas
the **tangent bundle has complex rank m**.
For $m=1$ ($\mathbb{CP}^1$), the
Fubini–Study metric in the standard
normalization has Gaussian curvature $4$.
With the current sign convention, for a
geodesic polygon of unit line frames
$y_j$ the Berry holonomy and intrinsic
tangent holonomy satisfy the exact scalar
relationship

```math
\beta=\prod_j
\frac{y_{j+1}^*y_j}{|y_{j+1}^*y_j|},
\qquad
U_{\mathrm{tangent}}=\overline{\beta}^{\,2}.
```

The **square-inverse** is essential: using
$\beta$ itself would incorrectly classify
tautological U(1) line holonomy as
Levi–Civita tangent holonomy. In terms
of the Fubini–Study signed area, the
tangent phase has curvature-4 scale,
twice the Berry-line phase magnitude.
The `spectral` mode reports
this **exact finite single-eigenvalue
trace measure** for complex CP¹ without
calling eig(U).

For $m=2$, a closed CP² loop acts on a
**genuine two-dimensional complex
tangent fiber**. Tests construct two
loops occupying different, overlapping
complex directions and verify the
resulting $2$-by-$2$ matrices do **not
commute**, a simple concrete manifestation
of non-Abelian $U(2)$ tangent holonomy.
This is different from two scalar
Berry phases associated with the
tautological line.

**Spectrum, limits and cost.** Except for
the explicitly solvable CP¹ scalar
relationship, a generic Grassmann
geodesic polygon has no claimed
analytical finite trace eigenphase
law. The `spectral` mode marks it
`unknown`, even though the
finite holonomy matrix itself can be
diagonalized numerically for exploratory
studies. No CUE/COE law or limiting
density is asserted. The model's
`model` mode exposes the
per-leg geodesic velocities, ambient
unitary lifts, polar alignments, orthogonal
complement changes, and intrinsic
tangent links, allowing independent
metric and covariance checks.

The implementation uses dense SVD
and matrix exponentials of the
$n$-by-$n$ ambient geodesic generators,
with $2\le n\le16$. The intrinsic
$d$-by-$d$ links are Kronecker products;
**there is no claimed matrix-free
speedup** from this dense construction.
Future fast exponential-action techniques
may change that assessment.

**References:** Edelman–Arias–Smith
[97] for Grassmann geodesics and
horizontal geometry; Bendokat–Zimmermann–
Absil [111] for canonical Grassmann
Levi–Civita transport and the cut locus;
Mukunda–Simon [110] for the contrasting
tautological geometric-phase convention.




## Exact finite classical Koopman dynamics on a modular torus — [`finite_koopman_torus.m`](../finite_koopman_torus.m)

**Classical dynamics, not quantum cat-map quantization.** The finite
phase space is the abelian group $X_m=(\mathbb Z/m\mathbb Z)^2$,
equipped with uniform normalized counting measure, $m\ge2$.
For any integer $2$-by-$2$ matrix $M$ whose determinant is invertible
modulo $m$ and any integer translation $b$, define the affine map

```math
T(x)=Mx+b\pmod m,\qquad
\gcd(\det M,m)=1,\qquad N=|X_m|=m^2.
```

The default $M=\begin{bmatrix}2&1\\1&1\end{bmatrix}$
is a familiar hyperbolic cat automorphism, now **reduced
exactly** on the finite grid. The shear
$\begin{bmatrix}1&1\\0&1\end{bmatrix}$,
quarter-turn $\begin{bmatrix}0&-1\\1&0\end{bmatrix}$,
and a nonzero translation with $M=I$ are alternative
reproducible exact models. The modular invertibility
check rejects non-bijective maps (e.g. determinant two
with even grid modulus).

With the specified coordinates $x=(x_1,x_2)$ and
MATLAB flattened state index
$j=1+x_1+m x_2$ ($0\le x_1,x_2<m$), define the
**finite-state Koopman pullback**

```math
(U_T f)(x)=f(Tx),\qquad
(U_T)_{j,k}=\begin{cases}
1,&k=1+(Tx)_1+m(Tx)_2,\\
0,&\text{otherwise}.
\end{cases}
```

The map is a genuine bijection, so each row and
column contains exactly one nonzero equal to one.
Consequently $U_T$ is a **real orthogonal permutation
matrix**, and hence a strict finite unitary:
$U_T^TU_T=I$. Its order is $N=m^2$, unlike
`quantum_cat_map`, whose order n is the
quantized Hilbert-space dimension and whose
matrix is a complex metaplectic chirp–Fourier–
chirp operator [93].

**Exact finite trace eigenphases by orbit structure.**
The permutation decomposes uniquely into disjoint
cycles $\mathcal C_a$ of lengths $\ell_a$,
$\sum_a\ell_a=N$. Restricting Koopman pullback
to a cycle gives an ordinary cyclic shift, so

```math
\mathrm{spec}(U_T)=\biguplus_a
\big\{e^{2\pi i q/\ell_a}:q=0,\ldots,\ell_a-1\big\}.
```

The **finite equal-weight trace** spectral measure
is the uniform empirical distribution on this
multiset, including multiplicity. It is not
the spectral measure of any selected cyclic
vector or invariant observable. In particular,
for any integer $k$,

```math
m_k=\frac1N\mathrm{tr}(U_T^k)
=\frac{\#\{x\in X_m:T^k(x)=x\}}{N}
=\frac1N\sum_{\ell_a\mid k}\ell_a.
```

The fixed-point formula remains valid for
negative powers and $k=0$, when all states
are fixed. The `spectral` mode
stores the cycle lengths, exact finite
root-of-unity eigenvalues and a function
for moments computed using integer
cycle divisibility, not numerical
diagonalization. An optional `map`
mode returns the forward/inverse
permutations, full cycle decomposition,
and the finite permutation order
(where representable without integer overflow).

**Fourier-character conjugacy, including affine phases.**
The normalized Fourier characters are

```math
\chi_k(x)=\frac1{\sqrt N}
e^{2\pi i\,k^Tx/m},\quad k\in X_m.
```

They form an orthonormal basis of
$\mathbb C^{N}$, and the exact pullback
action on a character is

```math
U_T\chi_k
=e^{2\pi i k^Tb/m}\chi_{M^Tk}.
```

Thus the matrix in Fourier-character
coordinates is again a monomial unitary,
now with **dual frequency permutation**
$k\mapsto M^Tk$ and affine phase. The
`map` diagnostics return both the
frequency map and phases, using the
same positive-exponent character
convention as these formulas.
Independent tests explicitly assemble
a dense 2D discrete Fourier basis for
small moduli and verify the identity
$F^*U_TF=U_{\mathrm{dual}}$.
This is an exact finite equivalence,
not an approximation or a claim that the
matrix has the eigenvector structure
of a quantum metaplectic cat map.

A translation $b=(1,0)$ with $M=I$
consists of **m independent cycles of
length m**. Each cycle is one of the
existing `weighted_shift`
families (with zero phases), after
reordering states. General linear
cat/shear/rotation permutations have
multiple cycle lengths and are not
one globally weighted cyclic shift.

**Operator representation and limitations.**
The default `matrix` mode builds a
sparse $N$-by-$N$ permutation, capped
at $N\le4096$ to avoid excessive
explicit allocations. The
`operator` mode computes both
forward and adjoint operations by
indexing the forward or inverse
permutation. These operations cost
$O(N)$ per right-hand side and preserve
complex-valued inputs exactly up to
copying; they apply to moduli
$2\le m\le1000$ ($N\le10^6$).
The complete cycle decomposition is
$O(N)$ in work and memory.

This finite phase-space permutation
should **not** be confused with the
infinite-dimensional Koopman operator
on $L^2(\mathbb T^2)$ associated
with the continuous hyperbolic toral
automorphism. Fixed finite-modulus
spectra are finite roots of unity and
do not establish a continuous
spectrum, convergence of its
eigenphase laws, or numerical rates
for the infinite problem.
Colbrook–Townsend [146] discuss
rigorous data-driven infinite Koopman
spectra; the present elementary
finite permutation is **not
attributed** to their particular
experiments. That distinction is
also relevant to open issues #61
(geodesic-flow Koopman) and #107
(mpEDMD data-derived metric
unitaries).

**Tests.** The implementation is checked
against independent modular-state
enumeration for cat maps, shears,
rotations and translations, including
a noninvertible rejected example;
exact cycle roots versus numerical
matrix eigenvalues and integer trace
powers; fixed-point moments (including
negative powers); Fourier-character
dual conjugacy on small tori and
complex phases for translations;
multi-RHS forward/adjoint reversibility;
and the decomposition of a translation
into a direct sum of existing cyclic
permutation blocks.

**References:** Hannay–Berry [93]
for the distinct quantization of
classical toral maps; Colbrook–
Townsend [146] for the rigorous
Koopman spectral problem on
infinite-dimensional spaces.


## Data-derived metric-unitary mpEDMD Koopman models — [`koopman_mpedmd.m`](../koopman_mpedmd.m)

**Physical origin.** A bijective measure-preserving transformation
$T$ induces a unitary Koopman composition operator on $L^2(\mu)$,
$\mathcal Kf=f\circ T$. A noninvertible measure-preserving map
instead produces an isometry that need not be onto. Neither
ordinary EDMD nor an arbitrary Galerkin compression is generally
unitary. Colbrook's mpEDMD [148] uses weighted snapshot
data and a Gram-aware orthogonal Procrustes solve
to enforce a finite-dimensional metric isometry.
It is distinct from `polar_unitary`, which has no
trajectory, dictionary or measured inner-product structure.

**Dictionary and snapshot inputs.** The generator accepts
`koopman_mpedmd(n,PX,PY,weights,mode)` where
$m\ge n$ rows of PX and PY are dictionary evaluations
at snapshot pairs $(x_j,y_j=T(x_j))$. All observation
weights $w_j$ are positive; their normalization is fixed by

```math
W=\mathrm{diag}(w_j/\sum_\ell w_\ell),\qquad
G=P_X^*W P_X,\qquad A=P_X^*W P_Y.
```

The measured G must be Hermitian SPD with
moderate condition number; duplicate or nearly
linearly dependent dictionary vectors are rejected.
For $H=G^{1/2}$, form

```math
B=H^{-1}A H^{-1}=U_B\Sigma V_B^*,
\qquad Q=U_BV_B^*,\qquad Q^*Q=I.
```

This Euclidean-unitary Q is the **matrix admitted
to the database**. The associated coefficient-space
Koopman matrix and ordinary EDMD matrix are

```math
K_{\mathrm{mp}}=H^{-1}QH,\qquad
K_{\mathrm{mp}}^*G K_{\mathrm{mp}}=G,\qquad
K_{\mathrm{EDMD}}=G^{-1}A.
```

The first is **G-unitary**, generally not Euclidean
unitary; the latter is not in general G-unitary.
The `metric` mode exposes G, A, B,
H, inverse H, Q, K, EDMD and cross-Gram
singular values. The `data` mode returns
these along with the actual sampled dictionaries,
weights and default snapshot coordinates.
`edmd` returns the ordinary EDMD
comparison, **not** a certified unitary.

**Published algorithm orientation.** Colbrook's public
MATLAB `mpEDMD(G,A)` implementation
(https://github.com/MColbrook/Measure-preserving-Extended-Dynamic-Mode-Decomposition)
uses the SVD of $H^{-1}A^*H^{-1}=B^*$
and computes $V_B U_B^*$ in that adjoint-SVD
convention. This is equivalent to $\mathrm{polar}(B)$
above. Incorrectly reversing the adjoint
would conjugate all phases. A regression test
independently reproduces the author's
SVD orientation and coefficient matrix.

A singular **G** is rejected; if **B** is
rank-deficient, a full SVD still supplies
a unitary Procrustes completion, but its
action in the null space is **nonunique**.
The model exposes a numerical rank and
`rank_deficient_completion` flag
rather than claiming that missing-data
directions have determined phases.

**Exactly soluble default.** The default dictionary
comprises $\psi_k(x)=e^{ikx}$, $0\le k<n$,
sampled at $m=\max(64,2n+7)$ equispaced
points on the circle, with the rotation
$y_j=x_j+\alpha\pmod{2\pi}$, $\alpha=0.37$.
This invariant dictionary has $G=I$ in
exact quadrature arithmetic and

```math
Q=\mathrm{diag}(1,e^{i\alpha},\ldots,
e^{i(n-1)\alpha}),\quad
\lambda_k=e^{ik\alpha},\quad
m_\ell^{\mathrm{trace}}=\frac1n\sum_{k=0}^{n-1}\lambda_k^\ell.
```

The `spectral` mode gives this **exact finite
trace eigenphase/moment law** only for the
default invariant Fourier rotation. With
arbitrary noisy or noninvariant snapshot
dictionaries the *analytic* finite phase
reference is `unknown`. Numerical
eigenvalues of Q do not establish an
infinite Koopman spectral limit [146,148].

**Observable-weighted versus trace phases.**
For any coefficient vector $c\ne0$ in
the original (possibly nonorthogonal)
dictionary, let $v=Hc$. If
$Q=V\Lambda V^*$ is a unitary Schur
diagonalization up to roundoff, then

```math
\mu_c=\sum_j
\frac{|(V^*Hc)_j|^2}{c^*Gc}
\,\delta_{\arg\Lambda_{jj}},\qquad
m_\ell(c)=
\frac{(Hc)^*Q^\ell (Hc)}{c^*Gc}.
```

The callback `spectral.observable_measure(c)`
computes **numerical finite** eigenphases and
positive normalized cyclic-vector weights.
These differ from the equal $1/n$ trace weights.
They are not certified analytic
infinite-volume eigenphases.

**Validation and cost.** The independent tests
reproduce the original mpEDMD SVD recipe,
compare $Q^*Q$ and $K^*GK$ with identities,
exhibit standard EDMD non-isometry for
noisy/noninvariant dictionaries, verify
the default circle phases and moments,
weight rescaling invariance, singular
B nonuniqueness, Gram-rank rejection,
real-data orthogonal outputs, and cyclic
spectral-measure moments. Dense dictionary
dimension is capped at $n\le64$ with at
most 20,000 snapshot pairs. No
unwarranted fast matrix-free
action is claimed.

**References:** [148] M. J. Colbrook,
*The mpEDMD Algorithm for Data-Driven
Computations of Measure-Preserving
Dynamical Systems*, SINUM 2023,
https://doi.org/10.1137/22M1521407;
[146] Colbrook–Townsend, rigorous
Koopman spectral computations, CPAM 2024,
https://doi.org/10.1002/cpa.22125.
Unlike the *exact* finite-torus permutation
Koopman family in issue #108, mpEDMD
approximates dynamics from measured
snapshot dictionaries.



## Electric quantum walks and finite skew-shift CMV-factor models — [`electric_quantum_walk.m`](../electric_quantum_walk.m)

**Physical motivation and honest finite closure.**
The electric quantum walk of Cedzich–Werner [149]
and Yang [150] is a two-component, discrete-time
shift-and-coin dynamics on the **infinite** lattice
$\mathbb Z$, with a phase increasing linearly
with lattice position. Its relation to
skew-shift CMV operators is central
to spectral localization and interval-spectrum
results. The database contains a **finite
cyclic ring** of $n$ sites, so the matrix is
$2n$-by-$2n$ and strictly unitary for all
admissible real field parameters. An
irrational field $\omega$ used on this
finite ring is only a finite window
closed by a unitary boundary link:
it is **not** an infinite irrational
periodic operator. Infinite-volume
Anderson localization or full-circle
spectrum theorems [149,150] cannot
be inferred from a single finite
ring's eigenvalue plot.

**Exact basis, coin, shift, and phases.**
Use sites $j=0,\ldots,n-1$ and
the site-interleaved basis
$|j,-\rangle,|j,+\rangle$, flattened
with indices $2j+1,2j+2$.
For $0\le\gamma\le\pi/2$,
let $a=\cos\gamma$ and
$b=\sin\gamma$. Our real
$SU(2)$ coin, in this ordered
minus/plus convention, is

```math
C_\gamma=
\begin{bmatrix}a&b\\-b&a\end{bmatrix}.
```

The conditional shift S acts as

```math
S|j,-\rangle=|j-1,-\rangle,\qquad
S|j,+\rangle=|j+1,+\rangle,
```

with site positions interpreted modulo n
and a boundary twist $\beta$ **in radians**:
minus crossing the left wrap picks up
$e^{-i\beta}$, plus crossing the
right wrap $e^{i\beta}$.
For electric strength $\omega$ and
global offset $\vartheta$, both
measured in cycles, the electric
potential is the **output-site**
diagonal operator

```math
D|j,\pm\rangle
=e^{2\pi i(\vartheta+j\omega)}|j,\pm\rangle,
\qquad
W=D(I_n\otimes C_\gamma)S.
```

The finite matrix W is unitary because
all three factors are unitary. Explicit
sparse output is supported for
$3\le n\le2048$, and sparse
forward/adjoint/transpose actions
for $n\le10^5$. Complex transpose
and conjugate transpose are
**distinct** when the electric phase
or boundary twist is nontrivial;
the operator implementation handles
both correctly for multiple
right-hand sides.

**Exact cyclic CMV-factor identity.**
Introduce an on-site flip
$F=\bigoplus_{j=0}^{n-1}\begin{bmatrix}0&1\\1&0\end{bmatrix}$.
The conditional shift factors as
$S=FM$, where M swaps the
neighboring odd/even intersite pairs
$|j,+\rangle\leftrightarrow |j+1,-\rangle$ around the ring;
the boundary pair carries
the twist $e^{\pm i\beta}$.
The local identity is

```math
C_\gamma
\begin{bmatrix}0&1\\1&0\end{bmatrix}
=\begin{bmatrix}b&a\\a&-b\end{bmatrix}
=\Theta(b),\qquad
\Theta(\alpha)=
\begin{bmatrix}\overline{\alpha}&\rho\\
\rho&-\alpha\end{bmatrix}.
```

Consequently, with
$L=\bigoplus_{j=0}^{n-1}\Theta(b)$,

```math
W=DLM,\qquad L^*L=M^*M=D^*D=I.
```

This is a **finite, exact
phase-decorated cyclic
even/odd CMV-factorization**,
not an approximate CMV surrogate.
It is independently checked
against $D(I_n\otimes C_\gamma)S$
and its eigenvalues and normalized
trace moments in the tests.
When $\omega=\vartheta=\beta=0$,
$D=I$ and the factors form an
exact standard **cyclic** CMV product
with even Verblunsky coefficients
$\alpha_{2j}=b$ and odd coefficients
$\alpha_{2j+1}=0$. This is not
the usual **half-line finite CMV**
with a terminal unit-modulus
Verblunsky coefficient; no
unjustified identification
with `finite_cmv` is made.

For the infinite-lattice walk,
Yang [150, Lemma 1.1] gives,
after a diagonal unitary gauge,
the **formal skew-shift
Verblunsky parameters**

```math
\alpha_{2j}
=b\,e^{2\pi i(j^2\omega+2j\vartheta)},
\qquad
\alpha_{2j+1}=0,\qquad
\rho_{2j}=a,\quad\rho_{2j+1}=1.
```

The `factors` mode supplies
these coefficients for the
selected site window, plus
the actual D, L, and M
matrices and a flag identifying
the zero-field standard cyclic
CMV case. **The whole-line
diagonal gauge can carry
nontrivial monodromy around
a finite ring**, so these
formal coefficients are not
claimed, for arbitrary
$\omega,\vartheta,\beta$, to
define a boundary-matched
standard finite CMV matrix
with the identical spectrum.
The factorization $W=DLM$
*does* remain exact.

**Two solvable finite trace spectra.**
The `spectral` interface
follows the database's
equal-weight **matrix trace**
convention, not the spectral
measure of a preferred
coin/site state.

For $\gamma=0$, the local coin
is the identity, and W
splits into **two weighted
n-cycles** (minus and plus)
with opposite boundary twist.
Writing

```math
\Phi=2\pi\left(n\vartheta+
\frac{\omega n(n-1)}2\right),
```

the exact finite eigenphase
multiset is

```math
\left\{
\exp\!\left(
\frac{i(\Phi-\beta+2\pi k)}n\right),
\exp\!\left(
\frac{i(\Phi+\beta+2\pi k)}n\right):
k=0,\ldots,n-1
\right\}.
```

At the other extreme
$\gamma=\pi/2$, the coin
is perfectly reflecting.
Every intersite pair
$\{|j,+\rangle,|j+1,-\rangle\}$
becomes an invariant
two-dimensional bond,
with eigenphase pair

```math
\lambda_{j,\pm}
=\pm i\exp\!\left(
\frac{i(\phi_j+\phi_{j+1})}{2}\right),
\qquad
\phi_j=2\pi(\vartheta+j\omega)
\pmod{2\pi}.
```

The boundary twist cancels
within the bond eigenvalues.
Each simple bond eigenvector
is supported on exactly
two basis states and has
inverse participation ratio
$\sum|v_i|^4/(\sum|v_i|^2)^2=1/2$.
This is an elementary **fully
decoupled finite localization
limit**, not evidence for
nontrivial infinite-lattice
Anderson localization.

Both limits have **exact finite
eigenvalue/moment descriptors**
without calling eig(W).
For generic $0<\gamma<\pi/2$
the analytical finite
eigenphase law is explicitly
`unknown`. The
matrix itself can of course
be diagonalized numerically
for finite-size comparisons;
no limiting density or
Cantor/interval theorem
is inferred for these rings.

**Reproducibility and tests.**
The MATLAB tests independently
assemble the physical
shift-and-coin matrix,
including asymmetric wrap
phases; independently assemble
the cyclic even/odd CMV factors;
verify W unitary and
$W=DLM=DCS$; check the
standard-CMV zero-field
reduction, formal skew-shift
coefficients and finite
boundary distinction;
check exact transparent-cycle
and reflecting-bond eigenphases
against independent eigensolvers
and trace powers; test
strong-coin bond localization,
rational and irrational
finite windows, global
phase covariance, multi-RHS
forward/transpose/adjoint
actions, and input rejection.
All examples distinguish
the finite **trace** phase
measure from any infinite
spectral type or spatially
weighted wavepacket measure.

**Sources:** [149] Cedzich–Werner,
*Anderson Localization for Electric
Quantum Walks and Skew-Shift CMV
Matrices* (CMP 2021); [150]
Yang, *Interval Spectrum for
Electric Quantum Walk and Related
Skew-Shift CMV Matrices* (LMP
2025). Both infinite-lattice
theorems are cited for the
model and its skew-shift CMV
connection, **not** as theorems
about arbitrary finite cyclic
truncations.




## Finite integrable isospectral trajectories (November 2026)

### Hessenberg QR/Schur flow — [`schur_hessenberg_flow.m`](../schur_hessenberg_flow.m)

Start with an unreduced unitary upper-Hessenberg matrix with prescribed distinct
circle eigenvalues `z_j` and strictly positive e1 spectral weights `w_j(0)`.
The positive-diagonal QR factorization of `exp(t H0)` defines
`H(t)=Q(t)'*H0*Q(t)`. The exact mathematical flow is unitary,
upper-Hessenberg, and isospectral. Its cyclic weights evolve as
`w_j(t) ∝ w_j(0)*exp(2*t*real(z_j))`. The trace weights remain `1/n`.
Direct QR, inverse OPUC Hessenberg reconstruction, and cross-representation
finite CMV output are exposed separately; the CMV model is *not* claimed to
be the direct Hessenberg QR similarity. Real conjugate-paired input
nodes/equal conjugate-pair weights give an orthogonal subclass. The direct
matrix path uses a dense exponential and QR with potential large-time
loss of accuracy, whereas the operator path uses existing core actions.
This is a finite isospectral benchmark, not a claim of a new spectral
distribution. Sources [5, 151, 152].

### Positive discrete CMV spectral-weight steps — [`discrete_schur_cmv.m`](../discrete_schur_cmv.m)

For distinct fixed support `z_j` on the circle, choose `0<h<1/2` and
integer `k>=0` and multiply every strictly positive weight by
`(1+2*h*real(z_j))^k`, followed by normalization. A finite inverse
CMV construction yields a genuine five-diagonal unitary with *unchanged*
eigenvalues at every step. Evolving e1 cyclic weights must not be
confused with equal-mass trace weights. With `h→0`, `k*h→t`,
the weights approach the continuous Schur-measure deformation
`w_j(t) ∝ w_j(0)*exp(2*t*real(z_j))`. Unlike the documented
Mukaihira–Nakamura integrable discretization, this particular
positive-measure step is **only a separately justified benchmark**;
its integrable Lax/auxiliary-variable formulation is not claimed.
All finite measures remain atomic, including when limiting OPUC
measures have continuous densities. Sources [5, 154, 155].

### Toda Jacobi QR flow and symmetric Cayley-unitary matrices — [`toda_cayley.m`](../toda_cayley.m)

For a real symmetric irreducible tridiagonal Jacobi matrix `J0`,
factor `exp(tau*J0)=Q(tau)*R(tau)` with positive diagonal R and set
`J(tau)=Q(tau)'*J0*Q(tau)`. The exact finite nonperiodic Toda
solution remains Jacobi and isospectral. For an *independent*
Cayley scale `gamma>0` set
`U(tau,gamma)=(J(tau)-i*gamma*I)/(J(tau)+i*gamma*I)`.
The unitary is complex symmetric, with fixed eigenvalues
`(x_j-i*gamma)/(x_j+i*gamma)` for `x_j=eig(J0)`.
The normalized trace measure is fixed but the distinguished e1
weights are `w_j(tau) ∝ w_j(0)*exp(2*tau*x_j)`.
The constructor accepts classical Jacobi recurrence parameters or
an arbitrary positive-subdiagonal tridiagonal recurrence. Dense
QR setup is distinct from cached sparse Cayley solves in operator
mode; numerical large-time QR rank breakdown is detected. This
is a nonlinear isospectral deformation, **not** `exp(i*tau*J0)`.
Sources [23, 153].


## Geronimus CMV spectral-measure benchmark — [`geronimus_cmv.m`](../geronimus_cmv.m) and [`szego_geronimus_measure.m`](../szego_geronimus_measure.m)

The existing finite CMV generator admits **complex** constant Verblunsky
coefficients `a` with `|a|<1`, while its last coefficient is a unit
boundary phase. The infinite e1 OPUC measure has an explicitly known
Lebesgue decomposition, *not* a finite CMV trace-density formula.
Set `theta_a=2*asin(|a|)` and `b=2*arg(1+conj(a))`.
Its absolutely continuous density relative to `dtheta` on
`|theta|>theta_a` (principal circle angle) is

```math
\rho_a(\theta)=\frac{\sqrt{\sin^2(\theta/2)-|a|^2}}
{2\pi|1+a|\,|\sin((\theta-b)/2)|}.
```

If `real(a)+|a|^2>0`, it also has an atom at angle `b` of mass

```math
m_a=\frac{2(\mathrm{Re}\,a+|a|^2)}{|1+a|^2};
```

otherwise `m_a=0`. In the real case `a=0.8` this is an atom of
mass `8/9` at phase zero plus continuous mass `1/9`;
for `a=-0.4` there is no atom. The `a=0` special case has
uniform density `1/(2*pi)`. The measure is implemented through
`js_geronimus_density` and the exact Poisson smoothing
`js_geronimus_poisson`, computed from the constant-coefficient Schur
fixed point and its Carathéodory function, not by sampling a finite CMV.

The `'spectral'` output of `geronimus_cmv` exposes this formula only
in `S.underlying` (infinite distinguished-vector measure), and gives
the *essential* finite-size limiting support in `S.limit`.
No absolutely continuous spectral mass exists at any finite order.
For finite-order comparisons, compute eigenpairs of
`geronimus_cmv(n,a,'cmv',beta)` and distinguish the `e1` weights
`|v_j(1)|^2` from equal-mass trace weights `1/n`. The terminal
phase `beta` changes the discrete nodes, and the fixed-radius
Poisson smoothing of the cyclic measure converges to the infinite
underlying Poisson-smoothed density. The Rogers–Szegő smooth
wrapped-Gaussian measure is a complementary control.
See Colbrook, CMP 384 (2021), Section 6.3, reference [147]; Simon,
OPUC, reference [5].

## Spherical polar-cap Slepian orthogonal eigenbases — [`spherical_slepian.m`](../spherical_slepian.m)

The real bandlimited spherical harmonic space (degrees 0,...,L-1) has
dimension L^2. For a polar cap of colatitude Theta form the real symmetric
concentration matrix K_ab = integral_cap Y_a Y_b dOmega. Its eigenvectors
form an orthogonal synthesis matrix Q; the generator returns analysis Q'.
The separately available `concentration` mode returns K, which is NOT unitary.

Polar-cap symmetry separates the real cosine/sine blocks of each azimuthal
order m. For fixed m, products P_ell^m(x) P_k^m(x) are polynomials of degree
at most 2L-2 in x=cos(theta), so L-point Gauss-Legendre quadrature evaluates
the finite cap integrals exactly in mathematics. The eigenvectors are
sorted by decreasing concentration with deterministic sign convention
(largest-magnitude coefficient positive), preserving separate cos/sin
block bases under repeated eigenvalues.

Concentration factors are eigenvalues of K in [0,1], NOT eigenvalues of
the orthogonal matrix Q'. The latter has a generally unknown phase law.
The exact Shannon number is trace(K)=L^2*(1-cos(Theta))/2. Empty/full
caps return the canonical Q=I; `parameters` includes K, eigenvalues and
real-harmonic basis labels. No general polygon/mask or fast matrix-free
mode is claimed. The family is distinct from 1D DPSS and Wigner rotations.
Source [156].
## Square axisymmetric Fourier-Laguerre transforms on the ball — [`fourier_laguerre_ball.m`](../fourier_laguerre_ball.m)

Leistedt–McEwen's Fourier-Laguerre approach uses radial functions
K_p(r) = sqrt(p!/(p+2)!)*tau^(-3/2)*exp(-r/(2*tau))*L_p^(2)(r/tau),
orthonormal in L2([0,infinity),r^2 dr), and spherical harmonics in angle.
The P generalized Gauss-Laguerre nodes x_j of parameter 2 give radii
r_j=tau*x_j. If w_j integrates x^2 exp(-x), the physical radial
quadrature weights are tau^3 exp(x_j) w_j; weighted radial samples of
K_p give the square orthogonal Qr. The (-1)^p sign converts the existing
positive-leading Gauss-Laguerre transform to conventional L_p^(2).

In the exactly axisymmetric subspace m=0, the L Gauss-Legendre polar
nodes and analytically integrated azimuth give another square orthogonal
Qa with entries sqrt(2*pi*wtheta_j)*Y_(ell,0)(theta_j). The new square
family U=kron(Qa,Qr) has order P*L and exact U'*U=I in mathematical
arithmetic, with radial index fastest. The 'operator' interface applies
the tensor factors without materializing U (cost O(P*L*(P+L)) per RHS
after dense P/L factor setup); it is tensor matrix-free, not FFT-fast.

The 'sampling' mode is a DISTINCT full-sphere tensor grid using L
Gauss-Legendre colatitudes and 2L-1 azimuths, returning a rectangular
weighted matrix A of order [P*L*(2L-1),P*L^2]. Quadrature gives A'*A=I;
except L=1 it is NOT square unitary. The sampling mode is an ancillary
isometric output, not a second public square generator. The full exact
Fourier–Laguerre sampling theorem has no automatic square interpretation.

The radial physical scale tau rescales radii and weights but not the
normalized square matrix. The quadrature nodes are not its eigenphases;
its unit-circle trace spectral law is unknown. This family is not the
Gaussian-weight spherical Gauss–Laguerre (SGL) transform. Sources [157,158].
## Finite defocusing Ablowitz–Ladik / CMV Lax hierarchy — [`abl_ladik_cmv.m`](../abl_ladik_cmv.m)

Nenciu's finite-dimensional CMV Lax pair is an exact unitary,
structure-preserving evolution rather than an arbitrary discretized ODE.
Choose n even, distinct unimodular nodes z_j constrained to give terminal
Verblunsky coefficient -1, and strictly positive e1 weights w_j(0).
For the k-th real Hamiltonian H_k=2 Im Tr(C^k)/k, the exact spectral
solution is w_j(t) proportional to w_j(0)*exp(-2*t*Re(z_j^k)).
The eigenphases z_j and equal-mass trace moments stay fixed, whereas
cyclic spectral weights and interior Verblunsky coefficients evolve.

With X_+ equal to strictly upper triangular part plus half the diagonal,
the finite CMV matrix satisfies dC/dt = [C,B_k] for
B_k=(C^k)_+ - ((C^k)_+)^*. The skew-Hermitian generator ensures
isospectrality through a unitary similarity, and inverse finite OPUC
reconstruction ensures five-diagonal CMV structure at every sampled time.
The 'lax' mode exposes the instantaneous matrix, B, and commutator,
which are independently tested against symmetric time differences.
The k=1 conjugate-symmetric real subclass obeys the scalar Schur
(discrete modified KdV) flow: alpha_j'=(1-|alpha_j|^2)
*(alpha_(j-1)-alpha_(j+1)), with endpoint alphas -1.

Output modes matrix, operator, parameters, spectral and lax all respect
the same prescribed finite spectral data. The operator uses existing
CMV core actions; the Lax diagnostic explicitly materializes dense C^k.
This is not an identification of the real Schur flow with the
complex NLS-type AL Hamiltonian 2 Re Tr(C), nor does it claim
stability of the inverse construction for arbitrarily clustered nodes.
Sources [5, 159, 160]; especially Nenciu (2005), Theorem 3 and
Gekhtman–Nenciu (2008), Eq. (5.7)-(5.8).
## ETNA 2015 Section 7 reproducible examples — [`etna_hessenberg.m`](../etna_hessenberg.m)

Alongside existing Type I/II/controlled Type III, `kind='7.1'`
implements ETNA Example 7.1 with a **single common m and a single
common phase p** sampled per matrix. The core k first column is
 [(-1)^(k-1) p sqrt(1-m^(2/B^k)); sqrt(m^(2/B^k))] 
for k=1,...,n-1, with terminal factor (-1)^(n-1).
Here B (passed as `coupling` for this mode only) is the positive
correlation exponent base, B>=1 (B=2 by default), m uniform
in (0,1), p uniform on the circle. The optional `7.1-real`
preset fixes p=1 and gives a real orthogonal analogue. The
existing independent `Type I` coefficients are intentionally distinct.

`kind='7.5'` and `kind='7.5-real'` give a **modern inverse-OPUC
analogue** of the known-spectrum Example 7.5: equally spaced unit
circle eigenphases, deterministic seeded strictly positive cyclic
weights, and an inverse finite Hessenberg core realization. The
real model pairs conjugate weights. The `parameters` output includes
finite forward eigenphase error for n<=64. This uses our `inverse_opuc`
method, NOT the original Mach--Van Barel--Vandebril software cited
in ETNA Example 7.5. Consequently one must NOT attribute the
historical complex-generator failure for n>100, or its exact
random setup, to this modern realization. The `spectral` mode
reports target nodes, not an independently validated computed error.

These additions cover paper-defined correlated-core distribution
and an application-appropriate modern known-spectrum test, but
faithfully recovering the original historical Example 7.5 code
and reproducing the paper's exact historical Type III sample
remain research tasks under issue #33. Source: ETNA 44 (2015),
pp. 335–336, §7, reference [8].
## Exact Candan discrete fractional Fourier transform — [`discrete_fractional_fourier.m`](../discrete_fractional_fourier.m)

The normalized DFT `F(j,k)=exp(-2*pi*i*j*k/n)/sqrt(n)` is an exact
special value of the Candan–Kutay–Özaktas finite fractional family,
not merely a continuous Fourier-integral sampling approximation.
The real symmetric periodic discrete oscillator

```math
S_{j,j}=2\cos(2\pi j/n)-4,\qquad
S_{j,j+1}=S_{j,j-1}=1\quad\text{(cyclic indices; accumulated if coincident)}
```

commutes with F. Its degeneracy when 4 divides n is resolved
by eigenanalysis on distinct **even** and **odd** parity subspaces.
For each parity block, sort the oscillator eigenvectors in descending
eigenvalue order: the j-th even/odd eigenvector is assigned respectively
index 2j or 2j+1. The resulting real orthonormal V is uniquely fixed
up to column signs, which are set by largest-absolute-entry positivity.

The precise discrete Hermite index set is `0,...,n-1` for odd n,
but **`0,...,n-2,n` for even n**. Thus

```math
F_{\alpha}=V\,\mathrm{diag}(e^{-i\pi\alpha k/2})\,V^*,
\qquad F_0=I,\quad F_1=F,\quad
F_{\alpha}F_{\beta}=F_{\alpha+\beta},\quad
F_{\alpha+4}=F_\alpha.
```

All identities hold in exact arithmetic; finite-precision tests
also verify n divisible by 4, the skipped even-n index, and DFT
eigenvectors explicitly. `basis` exposes V, k, parity and S;
`oscillator` returns symmetric S (NOT a unitary matrix), and
`hamiltonian` returns V diag(k) V* (also NOT a unitary matrix).
`spectral` supplies the exact finite trace eigenphase law and moments
for the unitary F_alpha, without inferring a limiting distribution.
Only explicit dense construction is offered; its O(N^3) eigensolver
setup and O(N^2) storage do NOT constitute a fast matrix-free action.
This is distinct from the existing Gauss-Hermite quadrature oscillator
and fractional Fourier integrals computed by approximate chirp FFTs.
Reference [161], Section III, Equations (9), (15), (22) and Table IV.

## Finite positive-measure rational Arnoldi on the circle — [`rational_arnoldi_circle.m`](../rational_arnoldi_circle.m)

For distinct prescribed nodes z_j on the unit circle with strictly
positive normalized weights w_j, let Z=diag(z_j) and q_1=sqrt(w).
Unlike `rational_opuc`, whose model is an operator Möbius transform
of a unitary Hessenberg Schur recurrence, this family constructs a
discrete-measure orthonormal rational Krylov basis Q explicitly.

At step j, with arbitrary pole p_j outside the circle, inside the
circle, or at infinity, use the stable normalized rational step

```math
v_j=\begin{cases}
Zq_j, & p_j=\infty,\\
Z(I-Z/p_j)^{-1}q_j, & |p_j|>1,\\
(Z-p_jI)^{-1}q_j, & |p_j|<1.
\end{cases}
```

Two-pass modified Gram–Schmidt constructs the complete full N-vector
orthonormal system Q. Every finite step is a rational function of Z;
poles on the unit circle are excluded. If a mixed sequence returns to
an already generated subspace, an explicit rational Arnoldi breakdown
is reported rather than silently completing by arbitrary vectors.

The exact full-space **unitary matrix** is

```math
U=Q^*ZQ, \qquad U^*U=I.
```

Its eigenphases are precisely the prescribed nodes z_j, with
**uniform trace weights** 1/N. The e1-cyclic spectral weights are
the distinct initial w_j: Q(:,1)=sqrt(w) exactly. The true
weighted-sample rational functions are Q(j,k)/sqrt(w_j).

The auxiliary rational Arnoldi Hessenberg pencil H,K has dimensions
N-by-(N-1) and satisfies U K = H. With h_j defined by
v_j=Q*h_j, the recurrence column is respectively

```math
\begin{array}{c|cc}
\mathrm{pole} & H_j & K_j\\\hline
\infty & h_j & e_j\\
|p_j|>1 & h_j & e_j+h_j/p_j\\
|p_j|<1 & e_j+p_jh_j & h_j
\end{array}
```

Thus H,K are **rectangular** recurrence matrices, not individually
unitary, and the conventional reduced $N-1$ order rational compression
does not have to be unitary. The full Q'*Z*Q closure is unitary
because the entire finite discrete L2 measure space is represented.
An infinite-measure orthogonal-rational-function recurrence or
arbitrary structured square-pencil completion is not claimed.

Modes 'matrix', 'basis'/'parameters', 'pencil' and 'spectral' expose
the true distinctions and exact finite trace/cyclic measures.
No 'operator' mode is claimed since forming Q requires $O(N^3)$ work
and caching it requires $O(N^2)$ storage. Sources [9,162].

## Core-factored strict unitary–unitary Hessenberg pencils — [`unitary_hessenberg_pencil.m`](../unitary_hessenberg_pencil.m)

For N>=2, construct TWO independent core-factored, unitary
upper-Hessenberg matrices A and B, then use U=B'*A as the
strictly unitary default matrix. A=diag(d)*H(alphaA,betaA),
B=H(alphaB,betaB), where H is a product of N-1 local
2-by-2 Verblunsky cores plus a unimodular terminal closure,
and d is a diagonal unit-phase vector. The core parameterization
and quotient have O(N) storage and O(N) forward/adjoint actions.
Each core has sine rho=sqrt(1-|alpha|^2); the subdiagonal
entries of H are exactly rho and satisfy Hessenberg sparsity
below the first subdiagonal. The genuine generalized pencil
A-zB is a strict unitary-unitary pair, and

```math
\mathrm{eig}(A,B)=\mathrm{eig}(B^*A),
\qquad (B^*A)^*(B^*A)=I.
```

Prescribe N-1 **projective poles** xi_j=A_(j+1,j)/B_(j+1,j),
not transform eigenvalues. For 0<|xi|<infinity,
rhoA=t min(1,|xi|) and rhoB=t min(1,1/|xi|), with t=.7.
The left phase d_(j+1)=xi/|xi| supplies each pole argument.
Independent seeded Schur phases and terminal phases make real
(orthogonal) and complex benchmark families. For xi=0,
rhoA=0 while rhoB>0; for xi=infinity, rhoB=0 while rhoA>0.
These exactly realizable zero/infinite poles imply a REDUCIBLE
Hessenberg factor, hence are not proper unreduced pairs.
In exact arithmetic the strict finite poles can have any nonzero
modulus; in double precision the generator limits them to
1e-4<=|xi|<=1e4 for numerically resolved subdiagonals and
reports achieved ratios separately from targets.

The `pencil` output returns A,B,U and all pole/proper flags;
`parameters` returns only the compact phase/core data;
`operator` applies B'*A and its adjoint by core sequences
without storing dense A or B. The `spectral` output correctly
marks generic finite trace eigenphases as UNKNOWN: projective
pole ratios are not generalized eigenvalues.

**Separately classified nonunitary problem:** the auxiliary
[`um_hessenberg_unitary_generalized.m`](../um_hessenberg_unitary_generalized.m)
returns A0+gamma*diag((1:N)/(N+1)), B0 for a strict reference
pair A0,B0. The deformation preserves its offdiagonal pole ratios
and Hessenberg structure but makes A and B'*A GENERALLY
NONUNITARY, with generalized eigenvalues not constrained to
the circle. This belongs to the broader Hessenberg–unitary-
Hessenberg RQR/rational-QZ class; it is deliberately excluded
from the square unitary matrix catalog.

The construction is a new **derived benchmark specialization**
motivated by the pole-swapping literature [163,164]. It does
NOT implement the Camps–Mach–Vandebril–Watkins swapping
iterations themselves or claim all RQR problems have unitary
quotients. The separate nonunitary helper is registered as a
generalized problem in `metadata/eigenvalue_problems.json`.

## Real spherical-harmonic SO(3) rotations — [`real_spherical_harmonic_rotation.m`](../real_spherical_harmonic_rotation.m)

**Exactly real orthogonal (in arithmetic), not entrywise real(Wigner D).**
The repository's existing [`wigner_rotation.m`](../wigner_rotation.m)
uses the complex harmonics Y_l^m, ordered m=-l,...,l with
Condon–Shortley phases Y_l^(-m)=(-1)^m conj(Y_l^m).
For m>0, the associated REAL harmonics have the convention

```math
R_l^{m}=(Y_l^{-m}+(-1)^mY_l^m)/\sqrt2,
\quad R_l^{-m}=(Y_l^{-m}-(-1)^mY_l^m)/(i\sqrt2).
```

Define a sparse unitary basis change C_l, one dense but simple
two-column unitary 2-by-2 block for each pair +/-m plus the
unchanged m=0 harmonic. The real degree-l rotation is

```math
D_l^{\mathrm{real}}=C_l^*D_l^{\mathrm{complex}}C_l.
```

This matrix is REAL and orthogonal in exact arithmetic, and
its conjugate spectral pairs and trace characters coincide
with D_l^complex. The new real generator contains all degrees
l=0,...,L-1 in the direct sum (order N=L^2), with exact phases
exp(-i*m*delta) and multiplicities L-|m|, independent of axis.
The real basis convention and complex representation composition
were independently tested through arbitrary 3D Euler products.

For **both** complex `spherical_harmonic_rotation` and the new
real counterpart, `output='operator'` precomputes the per-degree
blocks and applies them without assembling an L^2-by-L^2 matrix.
Per RHS application/storage is O(sum_l (2l+1)^2)=O(L^3)
=O(N^(3/2)), compared with O(N^2)=O(L^4) generic dense action.
However, the present individual degree-l blocks are **still
formed by dense matrix exponentials**. The total precomputation
is O(sum_l (2l+1)^3)=O(L^4), and one irrep's dense multiply
still costs O(l^2); we do NOT claim to implement the
Gimbutas–Greengard pseudospectral fast rotation or S2FFT
Wigner algorithms. Issue #41 remains open for that work.
Source [108,112].

## Large-order Hermite/Laguerre weighted Gauss transforms — stability diagnostics

The existing [`gauss_hermite_transform.m`](../gauss_hermite_transform.m)
and [`gauss_laguerre_transform.m`](../gauss_laguerre_transform.m)
remain **dense** real orthogonal Gauss-evaluation matrices Q in
mathematical arithmetic. Their rows are Gaussian quadrature nodes;
their columns are positive-leading orthonormal polynomials weighted
by the square roots of normalized Gaussian quadrature weights.
These real nodes are NOT eigenphases of the finite orthogonal Q.

At large order, tiny first components of eigenvectors of the real
Jacobi tridiagonal can underflow, so the old rule
`sign(V(1,j))` can spuriously flip polynomial signs. For ascending
simple Jacobi eigenvalues and strictly positive offdiagonals,
Sturm oscillation gives the **exact arithmetic** relation
sign(V(n,j))=(-1)^(n-j)*sign(V(1,j)). The shared
`um_gauss_basis` now orients each column by whichever endpoint
has the larger numerical magnitude, while retaining the original
dense Golub–Welsch eigensystem and positive-leading convention.

The new `output='quadrature'` diagnostic independently computes
Gaussian normalized Christoffel weights via

```math
w_j/\mu_0=1/\sum_{k=0}^{n-1}p_k(x_j)^2,
```

where the three-term polynomial recurrence is dynamically rescaled
and squared sums accumulated with log-sum-exp, O(N^2). Its
`log_normalized_weights` can remain finite when the corresponding
representable double `normalized_weights` underflows to zero.
These are independently **evaluated approximate weights**;
large-order logweights are not given certified relative error.

**Polynomials versus physical functions.** For Hermite polynomials
orthonormal in exp(-x^2)/sqrt(pi), the L2(dx)-normalized Hermite
functions are h_k(x)=pi^(-1/4)*exp(-x^2/2)*p_k(x). Thus
h_k(x_j)=Q(j,k+1)/d_j with

```math
\log d_j=\tfrac14\log\pi+\tfrac12\log(w_j/\sqrt\pi)
+\tfrac12 x_j^2.
```

For generalized Laguerre weight x^alpha exp(-x)/Gamma(alpha+1),
the L2(dx)-normalized associated functions are
h_k(x)=x^(alpha/2)*exp(-x/2)*p_k(x)/sqrt(Gamma(alpha+1)).
Again h_k(x_j)=Q(j,k+1)/d_j with

```math
\log d_j=\tfrac12[\log\Gamma(\alpha+1)
+\log(w_j/\Gamma(\alpha+1))+x_j-\alpha\log x_j].
```

The new `log_function_scaling` reports these positive diagonal
row factors logarithmically; the unweighted physical function
sample matrix is generally NOT orthogonal by the Euclidean norm.
This diagnostic does not claim to implement Webb's full stabilized
physical Hermite transform [109] or the fast Hermite/Laguerre
asymptotic/butterfly methods. O(N^3) setup and O(N^2) matrix
multiplication remain. General fast and high-order error-controlled
methods remain issue #40. Primary background [104,106,109].

## Spin-weighted spherical harmonics — [`spin_spherical_harmonic_transform.m`](../spin_spherical_harmonic_transform.m)

The spin-weighted sphere basis is normalized with the **repo's**
integer-spin Wigner `wigner_rotation` convention:

```math
{}_sY_{\ell m}(\theta,\phi)=(-1)^s
\sqrt{\frac{2\ell+1}{4\pi}}
d^\ell_{m,-s}(\theta)e^{im\phi}.
```

Here `ell>=|s|`, `|m|<=ell`, and Wigner d uses the real generator
`K=(diag(c,1)-diag(c,-1))/2` with `d=expm(theta*K)`.
For fixed `(s,m)`, `ell0=max(|s|,|m|)` and
`a=|m+s|`, `b=|m-s|`, the normalized Wigner latitude basis
has the exact Jacobi representation

```math
G_\ell(x)=\sqrt{\frac{2\ell+1}{2}}d^\ell_{m,-s}(\arccos x)
=\frac{(-1)^{\max(0,m+s)}}{\sqrt{\mu_{a,b}}}
(1-x)^{a/2}(1+x)^{b/2}p_{\ell-\ell_0}^{(a,b)}(x),
```

where `mu_(a,b)=integral[-1,1] (1-x)^a(1+x)^b dx` and
`p_k` are positive-leading Jacobi polynomials orthonormal under
the **normalized** Jacobi weight. This gives a legitimate
**square** real orthogonal `(L-ell0)`-by-`(L-ell0)` Gauss–Jacobi
fixed-order block, with the spin/Condon–Shortley phase
`(-1)^(s+max(0,m+s))`. Its finite eigenphase spectrum is
**unknown**: the real Jacobi zeros are quadrature nodes, not
the eigenvalues of this orthogonal evaluation matrix.

**Full sphere sampling is fundamentally different.** At bandlimit
`L`, the spin-s coefficient space has `L^2-s^2` dimensions.
The `L`-point Gauss–Legendre colatitude nodes `x_j`,
with `2L-1` equispaced azimuths, give

```math
A_{(j,k),(\ell,m)}=
\sqrt{\frac{2\pi w_j}{2L-1}}\;{}_sY_{\ell m}
(\arccos x_j,2\pi k/(2L-1)).
```

Consequently `A` has size `L(2L-1)` by `L^2-s^2` and
`A^*A=I` exactly in mathematical arithmetic. Each azimuthal
mode is separated by discrete Fourier orthogonality and the
remaining Wigner d product is an ordinary polynomial in x
of degree at most `2L-2`, integrated exactly by Gauss–Legendre.
The grid *oversamples*: it is NOT a square unitary except
at `L=1,s=0`. The `sampling` interface explicitly returns
this rectangular matrix, while `sampling_operator` returns
forward/adjoint FFT-separated actions without forming it.

For multiple RHS, the angular blocks use `O(L^3)` stored
coefficients, `O(L^3+L^2 log L)` per-vector arithmetic; the
registered *square* Jacobi block remains dense with generic
`O(M^3)` Gaussian-eigensolver setup. This is a derived
Gauss–Legendre scheme, **not** the equiangular MW sampling
theorem or Huffenberger–Wandelt's original fast Wigner method.
The Condon–Shortley identity
`conj({}_sY_(ell,m))=(-1)^(m+s) {}_-sY_(ell,-m)`
is tested for spins `s=0,+/-1,+/-2`, and the local polar
spin-frame convention is not asserted globally smooth.
Background references [104,107,165].

## Peter–Weyl Fourier transform on SO(3) — [`so3_peter_weyl_transform.m`](../so3_peter_weyl_transform.m)

Unlike `wigner_rotation(ell,...)`, which is a square Wigner
representation at **one prescribed rotation**, this construction
analyzes functions on the entire 3D rotation group. In the
z-y-z Euler convention, use **normalized Haar measure**

```math
d\mu(R)=\frac{d\alpha\,\sin\beta\,d\beta\,d\gamma}{8\pi^2},
\quad\alpha,\gamma\in[0,2\pi),\quad\beta\in[0,\pi].
```

The L2(dmu)-normalized Peter–Weyl basis is

```math
F_{\ell mn}(\alpha,\beta,\gamma)
=\sqrt{2\ell+1}\,e^{-im\alpha}
d^\ell_{mn}(\beta)e^{-in\gamma}.
```

For bandlimit `L`, the coefficient dimension is
`sum_(ell=0..L-1)(2ell+1)^2=L(4L^2-1)/3`.
A Gauss–Legendre rule of `L` colatitude nodes and
`M=2L-1` equally spaced values of each outer Euler angle
gives `L*M^2` distinct sample locations. The **weighted**
sampling matrix is

```math
A_{(j,k,h),(\ell,m,n)}=
\sqrt{\frac{w_j}{2M^2}}\;
\sqrt{2\ell+1}\,D^\ell_{mn}(\alpha_k,\beta_j,\gamma_h).
```

Here rows are beta-fastest, then alpha, then gamma;
columns are ordered ell, m=-ell:ell, n=-ell:ell.
Both finite Fourier rules resolve differences of azimuthal
frequencies up to 2L-2. At fixed m,n, products of
Wigner-d latitude factors are polynomials in x=cos(beta)
of degree at most 2L-2, so Gauss–Legendre integrates
them exactly. Hence `A^*A=I` mathematically: exact
finite Plancherel/Parseval for the bandlimited space.
Except when L=1, **A is rectangular**, NOT a unitary.

The registered **strict square orthogonal generator** is a
fixed-Wigner-index latitude block. Set
`ell0=max(|m|,|n|)`, `a=|m-n|`, `b=|m+n|`.
The corresponding `L-ell0` Gauss–Jacobi modes yield

```math
Q_{mn}=(-1)^{\max(0,m-n)}Q_{\mathrm{Gauss\text{-}Jacobi}}
(L-\ell_0,a,b),\qquad Q_{mn}^*Q_{mn}=I.
```

Its *finite matrix eigenphases are unknown*; the
Jacobi zeros are **real quadrature nodes**, not those phases.
The full `sampling` rectangle is explicitly materialized only
for L<=7, while `sampling_operator` performs an exact two-angle
FFT separation with fixed-(m,n) Wigner–Jacobi latitude blocks
in `O(L^4+L^3 log L)` work and `O(L^4)` cached storage per
coefficient map, compared to `O(L^6)` naive dense evaluation.
The adjoint is the weighted group Fourier analysis.

The implementation independently reconstructs the
Gauss–Legendre scheme related to [167]; it does NOT
reproduce the distinct McEwen et al. 2015 three-torus
extension [166] or claim optimized high-bandlimit S2FFT
Wigner recurrences. The normalization, left/right z-axis
rotation covariance, and finite Haar Plancherel are tested.
Background [104,108,166,167].

## Bibliography

1. P. J. Davis, *Circulant Matrices*, Wiley, 1979. ISBN 9780471057710. [Library record](https://openlibrary.org/books/OL4407364M/Circulant_matrices).
2. N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, Holt, Rinehart and Winston, 1976. ISBN 9780030839931. [Book record](https://books.google.com/books?id=1C9HAQAAIAAJ).
3. S. Mataigne and K. A. Gallivan, “The eigenvalue decomposition of normal matrices by the skew-symmetric part,” *Electron. J. Linear Algebra* **42** (2026), 349–386. [DOI](https://doi.org/10.13001/ela.2026.9957).
4. J. L. Aurentz, T. Mach, L. Robol, R. Vandebril and D. S. Watkins, *Core-Chasing Algorithms for the Eigenvalue Problem*, SIAM, 2018. [DOI](https://doi.org/10.1137/1.9781611975345).
5. B. Simon, *Orthogonal Polynomials on the Unit Circle*, AMS Colloquium Publications **54**, 2005: Part 1, *Classical Theory* [DOI](https://doi.org/10.1090/coll054.1); Part 2, *Spectral Theory* [DOI](https://doi.org/10.1090/coll/054.2).
6. M. J. Cantero, L. Moral and L. Velázquez, “Five-diagonal matrices and zeros of orthogonal polynomials on the unit circle,” *Linear Algebra Appl.* **362** (2003), 29–56. [DOI](https://doi.org/10.1016/S0024-3795(02)00457-3).
7. B. Simon, “CMV matrices: Five years after,” *J. Comput. Appl. Math.* **208** (2007), 120–154. [DOI](https://doi.org/10.1016/j.cam.2006.10.033).
8. J. L. Aurentz, T. Mach, R. Vandebril and D. S. Watkins, “Fast and stable unitary QR algorithm,” *Electron. Trans. Numer. Anal.* **44** (2015), 327–341. [ETNA](https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341).
9. A. Bultheel, M. J. Cantero and R. Cruz-Barroso, “Matrix methods for quadrature formulas on the unit circle: A survey,” *J. Comput. Appl. Math.* **284** (2015), 78–100. [DOI](https://doi.org/10.1016/j.cam.2014.11.002).
10. W. B. Gragg, “The QR algorithm for unitary Hessenberg matrices,” *J. Comput. Appl. Math.* **16** (1986), 1–8. [DOI](https://doi.org/10.1016/0377-0427(86)90169-X).
11. B. Simon, “OPUC on one foot,” *Bull. Amer. Math. Soc.* **42** (2005), 431–460. [DOI](https://doi.org/10.1090/S0273-0979-05-01075-X); [full text](https://web.ma.utexas.edu/mp_arc/c/05/05-78.pdf).
12. NIST Digital Library of Mathematical Functions, [§18.33, Orthogonal Polynomials on the Unit Circle](https://dlmf.nist.gov/18.33).
13. R. Killip and I. Nenciu, “Matrix models for circular ensembles,” *Int. Math. Res. Not.* **2004**(50), 2665–2701. [DOI](https://doi.org/10.1155/S1073792804141597).
14. F. Mezzadri, “How to generate random matrices from the classical compact groups,” *Notices Amer. Math. Soc.* **54**(5) (2007), 592–604. [arXiv preprint DOI](https://doi.org/10.48550/arXiv.math-ph/0609050) (DataCite, not a DOI for the journal article); [preprint](https://arxiv.org/abs/math-ph/0609050).
15. I. Daubechies, *Ten Lectures on Wavelets*, SIAM, 1992. [DOI](https://doi.org/10.1137/1.9781611970104).
16. T. Kitagawa, M. S. Rudner, E. Berg and E. Demler, “Exploring topological phases with quantum walks,” *Phys. Rev. A* **82** (2010), 033429. [DOI](https://doi.org/10.1103/PhysRevA.82.033429).
17. A. Bisio, G. M. D’Ariano and A. Tosini, “Quantum field as a quantum cellular automaton: The Dirac free evolution in one dimension,” *Ann. Phys.* **354** (2015), 244–264. [DOI](https://doi.org/10.1016/j.aop.2014.12.016).
18. A. Pérez, “Asymptotic properties of the Dirac quantum cellular automaton,” *Phys. Rev. A* **93** (2016), 012328. [DOI](https://doi.org/10.1103/PhysRevA.93.012328).
19. L. K. Grover, “A fast quantum mechanical algorithm for database search,” *Proc. STOC* (1996), 212–219. [DOI](https://doi.org/10.1145/237814.237866).
20. N. L. Balazs and A. Voros, “The quantized Baker's transformation,” *Ann. Phys.* **190** (1989), 1–31. [DOI](https://doi.org/10.1016/0003-4916(89)90259-5).
21. M. Saraceno, “Classical structures in the quantized baker transformation,” *Ann. Phys.* **199** (1990), 37–60. [DOI](https://doi.org/10.1016/0003-4916(90)90367-W).
22. R. Artuso, G. Casati and D. L. Shepelyansky, “Fractal spectrum and anomalous diffusion in the kicked Harper model,” *Phys. Rev. Lett.* **68** (1992), 3826–3829. [DOI](https://doi.org/10.1103/PhysRevLett.68.3826).
23. J. L. Aurentz, T. Mach, R. Vandebril and D. S. Watkins, “Computing the eigenvalues of symmetric tridiagonal matrices via a Cayley transformation,” *Electron. Trans. Numer. Anal.* **46** (2017), 447–459. [ETNA](https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46).
24. N. J. Higham, “Computing the polar decomposition—with applications,” *SIAM J. Sci. Statist. Comput.* **7** (1986), 1160–1174. [DOI](https://doi.org/10.1137/0907079).
25. B. Sz.-Nagy, *Unitary Dilations of Hilbert Space Operators and Related Topics*, CBMS 19, AMS, 1974. [DOI](https://doi.org/10.1090/CBMS/019).
26. M. Szegedy, “Quantum speed-up of Markov chain based algorithms,” *Proc. FOCS* (2004), 32–41. [DOI](https://doi.org/10.1109/FOCS.2004.53).
27. P. Bourgade, A. Nikeghbali and A. Rouault, “Circular Jacobi ensembles and deformed Verblunsky coefficients,” *Int. Math. Res. Not.* **2009**(23) (2009), 4357–4394. [DOI](https://doi.org/10.1093/imrn/rnp092).
28. F. Haake, *Quantum Signatures of Chaos*, 3rd ed., Springer, 2010. [DOI](https://doi.org/10.1007/978-3-642-05428-0).
29. V. M. Bastidas, P. Pérez-Fernández, M. Vogl and T. Brandes, “Quantum criticality and dynamical instability in the kicked-top model,” *Phys. Rev. Lett.* **112** (2014), 140408. [DOI](https://doi.org/10.1103/PhysRevLett.112.140408).
30. E. Lieb, T. Schultz and D. Mattis, “Two soluble models of an antiferromagnetic chain,” *Ann. Phys.* **16** (1961), 407–466. [DOI](https://doi.org/10.1016/0003-4916(61)90115-4).
31. P. Pfeuty, “The one-dimensional Ising model with a transverse field,” *Ann. Phys.* **57** (1970), 79–90. [DOI](https://doi.org/10.1016/0003-4916(70)90270-8).
32. B. Bertini, P. Kos and T. Prosen, “Exact spectral form factor in a minimal model of many-body quantum chaos,” *Phys. Rev. Lett.* **121** (2018), 264101. [DOI](https://doi.org/10.1103/PhysRevLett.121.264101).
33. D. J. Yates and A. Mitra, “Entanglement properties of the time-periodic Kitaev chain,” *Phys. Rev. B* **96** (2017), 115108. [DOI](https://doi.org/10.1103/PhysRevB.96.115108).
34. M. van Caspel, S. E. Tapias Arze and I. Pérez Castillo, “Dynamical signatures of topological order in the driven-dissipative Kitaev chain,” *SciPost Phys.* **6** (2019), 026 (closed-system Floquet discussion in §4). [DOI](https://doi.org/10.21468/SciPostPhys.6.2.026).
35. F. M. Izrailev, “Simple models of quantum chaos: Spectrum and eigenfunctions,” *Phys. Rep.* **196** (1990), 299–392. [DOI](https://doi.org/10.1016/0370-1573(90)90067-C).
36. B. Bertini, P. Kos and T. Prosen, “Exact correlation functions for dual-unitary lattice models in 1+1 dimensions,” *Phys. Rev. Lett.* **123** (2019), 210601. [DOI](https://doi.org/10.1103/PhysRevLett.123.210601).
37. B. Bertini, P. Kos and T. Prosen, “Random matrix spectral form factor of dual-unitary quantum circuits,” *Commun. Math. Phys.* **387** (2021), 597–620. [DOI](https://doi.org/10.1007/s00220-021-04139-2).
38. E. Vernier, H.-C. Yeh, L. Piroli and A. Mitra, “Strong zero modes in integrable quantum circuits,” *Phys. Rev. Lett.* **133** (2024), 050606. [DOI](https://doi.org/10.1103/PhysRevLett.133.050606).
39. J.-S. Caux, “Another example: the trigonometric R-matrix and the anisotropic S=1/2 antiferromagnet (XXZ model),” in *The Bethe Ansatz*, online notes, section a.xxz (2024). [Permanent section URL](https://integrability.org/a_xxz.html).
40. M. P. Strzys, E.-M. Graefe and H. J. Korsch, “Kicked Bose–Hubbard systems and kicked tops—destruction and stimulation of tunneling,” *New J. Phys.* **10** (2008), 013024. [DOI](https://doi.org/10.1088/1367-2630/10/1/013024).
41. C. Khripkov, D. Cohen and A. Vardi, “Coherence dynamics of kicked Bose-Hubbard dimers: Interferometric signatures of chaos,” *Phys. Rev. E* **87** (2013), 012910. [DOI](https://doi.org/10.1103/PhysRevE.87.012910).
42. C. Liang, Y. Zhang and S. Chen, “Statistical and dynamical aspects of quantum chaos in a kicked Bose-Hubbard dimer,” *Phys. Rev. A* **109** (2024), 033316. [DOI](https://doi.org/10.1103/PhysRevA.109.033316).
43. E. T. Jaynes and F. W. Cummings, “Comparison of quantum and semiclassical radiation theories with application to the beam maser,” *Proc. IEEE* **51** (1963), 89–109. [DOI](https://doi.org/10.1109/PROC.1963.1664).
44. J. T. Chalker and P. D. Coddington, “Percolation, quantum tunnelling and the integer Hall effect,” *J. Phys. C: Solid State Phys.* **21** (1988), 2665–2679. [DOI](https://doi.org/10.1088/0022-3719/21/14/008).
45. C.-M. Ho and J. T. Chalker, “Models for the integer quantum Hall effect: The network model, the Dirac equation, and a tight-binding Hamiltonian,” *Phys. Rev. B* **54** (1996), 8708–8713. [DOI](https://doi.org/10.1103/PhysRevB.54.8708).
46. A. Joye, “Unitary and open scattering quantum walks on graphs,” *Rev. Math. Phys.* (2026), 2650005. [DOI](https://doi.org/10.1142/S0129055X26500054); [arXiv](https://arxiv.org/abs/2409.08428) (see §3.1, Eq. (3.3)).
47. R. N. Bracewell, "Discrete Hartley transform," *J. Opt. Soc. Amer.* **73** (1983), 1832–1835. [DOI](https://doi.org/10.1364/JOSA.73.001832).
48. D. C. Chu, "Polyphase codes with good periodic correlation properties," *IEEE Trans. Inform. Theory* **18** (1972), 531–532. [DOI](https://doi.org/10.1109/TIT.1972.1054840).
49. X. Ouyang and J. Zhao, “Orthogonal chirp division multiplexing,” *IEEE Trans. Commun.* **64**(9) (2016), 3946–3957. [DOI](https://doi.org/10.1109/TCOMM.2016.2594792); [preprint](https://arxiv.org/abs/1601.06576).
50. A. Bemani, N. Ksairi and M. Kountouris, “AFDM: A full diversity next generation waveform for high mobility communications,” *IEEE International Conference on Communications Workshops (ICC Workshops)* (2021), 1–6. [DOI](https://doi.org/10.1109/ICCWorkshops50388.2021.9473655); [preprint](https://arxiv.org/abs/2104.11331).
51. H. Bölcskei and F. Hlawatsch, "Discrete Zak transforms, polyphase transforms, and applications," *IEEE Trans. Signal Process.* **45** (1997), 851–866. [DOI](https://doi.org/10.1109/78.564174).
52. M. Püschel and J. M. F. Moura, “Algebraic signal processing theory: Foundation and 1-D time,” *IEEE Trans. Signal Process.* **56**(8) (2008), 3572–3585. [DOI](https://doi.org/10.1109/TSP.2008.925261); [author PDF](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf).
53. K. Gröchenig, *Foundations of Time-Frequency Analysis*, Birkhäuser, 2001. [DOI](https://doi.org/10.1007/978-1-4612-0003-1).
54. V. Britanak, P. C. Yip and K. R. Rao, *Discrete Cosine and Sine Transforms: General Properties, Fast Algorithms and Integer Approximations*, Academic Press, 2007. [DOI](https://doi.org/10.1016/b978-0-12-373624-6.x5000-0); [publisher](https://shop.elsevier.com/books/discrete-cosine-and-sine-transforms/britanak/978-0-12-373624-6).
55. K. G. Beauchamp, *Applications of Walsh and Related Functions: With an Introduction to Sequency Theory*, Academic Press, 1984. ISBN 9780120841806. [Library record](https://openlibrary.org/works/OL2903378W/Applications_of_Walsh_and_related_functions_with_an_introduction_to_sequency_theory); [book record](https://books.google.com/books?id=BvZSAAAAMAAJ).
56. D. Slepian, "Prolate spheroidal wave functions, Fourier analysis, and uncertainty—V: The discrete case," *Bell Syst. Tech. J.* **57** (1978), 1371–1430. [DOI](https://doi.org/10.1002/j.1538-7305.1978.tb02104.x).
57. D. I. Shuman, S. K. Narang, P. Frossard, A. Ortega and P. Vandergheynst, “The emerging field of signal processing on graphs: Extending high-dimensional data analysis to networks and other irregular domains,” *IEEE Signal Process. Mag.* **30**(3) (2013), 83–98. [DOI](https://doi.org/10.1109/MSP.2012.2235192).

58. T. Kottos and U. Smilansky, “Periodic orbit theory and spectral statistics for quantum graphs,” *Ann. Phys.* **274** (1999), 76–124. [DOI](https://doi.org/10.1006/aphy.1999.5904).
59. D. Damanik, A. Pushnitski and B. Simon, “The analytic theory of matrix orthogonal polynomials,” *Surveys Approx. Theory* **4** (2008), 1–85. [arXiv](https://arxiv.org/abs/0711.2703).
60. D. Rocchesso and J. O. Smith, “Circulant and elliptic feedback delay networks for artificial reverberation,” *IEEE Trans. Speech Audio Process.* **5** (1997), 51–63. [DOI](https://doi.org/10.1109/89.554269).
61. D. Alpay, P. Jorgensen and I. Lewkowicz, “Characterizations of rectangular (para)-unitary rational functions,” *Opuscula Math.* **36** (2016), 695–716. [DOI](https://doi.org/10.7494/OpMath.2016.36.6.695).
62. A. Bultheel and M.-J. Cantero, “A matricial computation of rational quadrature formulas on the unit circle,” *Numer. Algorithms* **52** (2009), 47–68. [DOI](https://doi.org/10.1007/s11075-008-9257-9).
63. T. Betcke, N. J. Higham, V. Mehrmann, C. Schröder and F. Tisseur, “NLEVP: A collection of nonlinear eigenvalue problems,” *ACM Trans. Math. Software* **39** (2013), article 7. [DOI](https://doi.org/10.1145/2427023.2427024); [NLEVP MATLAB source](https://github.com/ftisseur/nlevp).
64. L. Lu, T. Wang, Y.-C. Kuo, R.-C. Li and W.-W. Lin, “A fast algorithm for fast train palindromic quadratic eigenvalue problems,” *SIAM J. Sci. Comput.* **38** (2016), A3410–A3429. [DOI](https://doi.org/10.1137/16M1063563).
65. T.-M. Huang, T. Li, W.-W. Lin and C.-T. Wu, “Numerical studies on structure-preserving algorithms for surface acoustic wave simulations,” *J. Comput. Appl. Math.* **244** (2013), 140–154. [DOI](https://doi.org/10.1016/j.cam.2012.11.020).
66. P. Dörfler and G. Schmeisser, “Construction of unitary and normal companion matrices,” *Linear Algebra Appl.* **202** (1994), 193–220. [DOI](https://doi.org/10.1016/0024-3795(94)90190-2).
67. E. Fuentes and L. E. Garza, “CMV block matrices for symmetric matrix measures on the unit circle,” *Linear Algebra Appl.* **587** (2020), 314–333. [DOI](https://doi.org/10.1016/j.laa.2019.11.015).
68. B. Gutkin and U. Smilansky, “Can one hear the shape of a graph?” *J. Phys. A* **34** (2001), 6061–6068. [DOI](https://doi.org/10.1088/0305-4470/34/31/301).
69. A. Bunse-Gerstner and L. Elsner, “Schur parameter pencils for the solution of the unitary eigenproblem,” *Linear Algebra Appl.* **154–156** (1991), 741–778. [DOI](https://doi.org/10.1016/0024-3795(91)90402-I).
70. D. Damanik, M. Embree and J. Fillman, “Gap labels for zeros of the partition function of the 1D Ising model via the Schwartzman homomorphism,” *Indag. Math.* **35** (2024), 813–836. [DOI](https://doi.org/10.1016/j.indag.2023.05.004); [preprint](https://arxiv.org/abs/2211.01558).
71. L. Marin and H. Schulz-Baldes, “Scattering zippers and their spectral theory,” *J. Spectr. Theory* **3** (2013), 47–82. [DOI](https://doi.org/10.4171/JST/37); [preprint](https://arxiv.org/abs/1112.4959).
72. C. Cedzich, J. Fillman and D. C. Ong, “Almost Everything About the Unitary Almost Mathieu Operator,” *Commun. Math. Phys.* **403** (2023), 745–794. [DOI](https://doi.org/10.1007/s00220-023-04808-4); [preprint](https://arxiv.org/abs/2112.03216).
73. D. Emms, S. Severini, R. C. Wilson and E. R. Hancock, “Coined quantum walks lift the cospectrality of graphs and trees,” *Pattern Recognition* **42** (2009), 1988–2002. [DOI](https://doi.org/10.1016/j.patcog.2008.10.025).
74. S. Baid, Á. Sáiz, L. Lamata, P. Pérez-Fernández, A. M. Romero, A. Rios, J. M. Arias and J. E. García-Ramos, “Extended Lipkin model: Proposal for implementation in a quantum platform and machine learning analysis of its phase diagram,” *Phys. Rev. C* **110** (2024), 044318. [DOI](https://doi.org/10.1103/PhysRevC.110.044318); [preprint](https://arxiv.org/abs/2404.15558).
75. P. Pérez-Fernández, J.-M. Arias, J.-E. García-Ramos and L. Lamata, “A digital quantum simulation of the Agassi model,” *Phys. Lett. B* **829** (2022), 137133. [DOI](https://doi.org/10.1016/j.physletb.2022.137133); [preprint](https://arxiv.org/abs/2105.02834).
76. J. Khalouf-Rivera, Q. Wang, L. F. Santos, J.-E. García-Ramos, M. Carvajal and F. Pérez-Bernal, “Degeneracy in excited-state quantum phase transitions of two-level bosonic models and its influence on system dynamics,” *Phys. Rev. A* **109** (2024), 062219. [DOI](https://doi.org/10.1103/PhysRevA.109.062219); [preprint](https://arxiv.org/abs/2303.16551).


77. R. Kukulski, M. Stępniak, K. Hendzel, Ł. Pawela, B. Gardas and Z. Puchała, “Local certification of unitary operations,” *Scientific Reports* **14** (2024), 26588. [DOI](https://doi.org/10.1038/s41598-024-75148-z).
78. M. Markiewicz, Ł. Pawela and Z. Puchała, “Single-shot comparison of random quantum channels and measurements,” *Scientific Reports* (2026). [DOI](https://doi.org/10.1038/s41598-026-70695-z); [preprint](https://arxiv.org/abs/2501.17317).
79. M. Markiewicz, Ł. Pawela and Z. Puchała, “Choi-level twirling of quantum channels: finite constructions and non-compact transformations,” *Journal of Physics A* (2026). [DOI](https://doi.org/10.1088/1751-8121/ae965b); [preprint](https://arxiv.org/abs/2512.23586).
80. Ł. Pawela and Z. Puchała, “IntegrateUnitary.jl: A Julia package for symbolic integration over Haar measures,” arXiv:2605.23830 (2026). [Preprint](https://arxiv.org/abs/2605.23830).
81. L. Marcinkowska, Ł. Pawela, M. Markiewicz and Z. Puchała, “Asymptotic distinguishability of Haar-averaged measurement models,” arXiv:2605.31168 (2026). [Preprint](https://arxiv.org/abs/2605.31168).
82. A. Bílek, R. Kukulski, P. Lewandowska, Ł. Pawela and Z. Puchała, “Noise Resilience of Quantum Key Distribution Protocols Secured Against Independent Attacks With One-Way Communication,” arXiv:2607.11857 (2026). [Preprint](https://arxiv.org/abs/2607.11857). Context only: no new generator.
83. T. Śmierzchalski and collaborators, “Hybrid quantum-classical computation for automatic guided vehicles scheduling,” *Scientific Reports* **14** (2024), 21809. [DOI](https://doi.org/10.1038/s41598-024-72101-y). Not a unitary-matrix source.
84. S. Nishigaki, “Eigenphase distributions of unimodular circular ensembles,” *Progress of Theoretical and Experimental Physics* **2024** (2024), 021A02. [DOI](https://doi.org/10.1093/ptep/ptae018); [preprint](https://arxiv.org/abs/2401.09045).
85. M. Fasi and L. Robol, “Sampling the eigenvalues of random orthogonal and unitary matrices,” *Linear Algebra and its Applications* **620** (2021), 297–321. [DOI](https://doi.org/10.1016/j.laa.2021.02.031); [preprint](https://arxiv.org/abs/2009.11515).
86. M. A. Lohe, “Non-Abelian Kuramoto models and synchronization,” *J. Phys. A: Math. Theor.* **42** (2009), 395101. [DOI](https://doi.org/10.1088/1751-8113/42/39/395101).
87. H. Choi, S.-Y. Ha and H. Park, “Emergent behaviors of discrete Lohe aggregation flows,” *Discrete Contin. Dyn. Syst. B* **27** (2022), 6083–6123. [DOI](https://doi.org/10.3934/dcdsb.2021308); [preprint](https://arxiv.org/abs/2110.12366).

88. L. Velázquez, “Spectral methods for orthogonal rational functions,” *J. Funct. Anal.* **254** (2008), 954–986. [DOI](https://doi.org/10.1016/j.jfa.2007.11.004); [arXiv](https://arxiv.org/abs/0704.3456).
89. T. Trogdon, “On spectral and numerical properties of random butterfly matrices,” *Appl. Math. Lett.* **95** (2019), 48–58. [DOI](https://doi.org/10.1016/j.aml.2019.03.024); [arXiv](https://arxiv.org/abs/1710.00087).

90. P. J. Forrester and E. M. Rains, “Jacobians and rank 1 perturbations relating to unitary Hessenberg matrices,” *Int. Math. Res. Not.* **2006** (2006), 48306. [DOI](https://doi.org/10.1155/IMRN/2006/48306); [preprint](https://arxiv.org/abs/math/0505552).
91. P. Diţă, “Some results on the parametrization of complex Hadamard matrices,” *J. Phys. A* **37** (2004), 5355–5374. [DOI](https://doi.org/10.1088/0305-4470/37/20/008).
92. W. Tadej and K. Życzkowski, “A concise guide to complex Hadamard matrices,” *Open Syst. Inf. Dyn.* **13** (2006), 133–177. [DOI](https://doi.org/10.1007/s11080-006-8220-2); [preprint](https://arxiv.org/abs/quant-ph/0512154).

93. J. H. Hannay and M. V. Berry, “Quantization of linear maps on a torus—Fresnel diffraction by a periodic grating,” *Physica D* **1** (1980), 267–290. [DOI](https://doi.org/10.1016/0167-2789(80)90026-3).
94. P. Feinsilver and R. Schott, “Krawtchouk transforms and convolutions,” *Advances in Applied Clifford Algebras* (2018). [DOI](https://doi.org/10.1007/s13373-018-0132-2).

95. R. Koekoek, P. A. Lesky and R. F. Swarttouw, *Hypergeometric Orthogonal Polynomials and Their q-Analogues*, Springer Monographs in Mathematics, Springer, 2010. [DOI](https://doi.org/10.1007/978-3-642-05014-5); [DLMF duality](https://dlmf.nist.gov/18.21).
96. P.-A. Absil, R. Mahony and R. Sepulchre, *Optimization Algorithms on Matrix Manifolds*, Princeton University Press, 2008. [DOI](https://doi.org/10.1515/9781400830244); [author's book website](https://sites.uclouvain.be/absil/amsbook/).
97. A. Edelman, T. A. Arias and S. T. Smith, “The geometry of algorithms with orthogonality constraints,” *SIAM J. Matrix Anal. Appl.* **20** (1998), 303–353. [DOI](https://doi.org/10.1137/S0895479895290954).
98. K. Crane, M. Desbrun and P. Schröder, “Trivial connections on discrete surfaces,” *Comput. Graph. Forum* **29** (2010), 1525–1533. [DOI](https://doi.org/10.1111/j.1467-8659.2010.01761.x).
99. Y. Ge, P. van Goor and R. Mahony, “A geometric perspective on fusing Gaussian distributions on Lie groups,” *IEEE Control Systems Letters* **8** (2024), 844–849. [DOI](https://doi.org/10.1109/LCSYS.2024.3405485); [preprint](https://arxiv.org/abs/2403.16411).
100. S. Fiori, “Manifold Calculus in System Theory and Control—Fundamentals and First-Order Systems,” *Symmetry* **13** (2021), 2092. [DOI](https://doi.org/10.3390/sym13112092).
101. E. Hairer, C. Lubich and G. Wanner, *Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations*, 2nd ed., Springer, 2006. [DOI](https://doi.org/10.1007/3-540-30666-8).
102. E. Hairer, C. Lubich and G. Wanner, “Geometric numerical integration illustrated by the Störmer–Verlet method,” *Acta Numerica* **12** (2003), 399–450. [DOI](https://doi.org/10.1017/S0962492902000144).
103. B. Cockburn, S. Du and M. A. Sánchez, “Combining finite element space-discretizations with symplectic time-marching schemes for linear Hamiltonian systems,” *Frontiers in Applied Mathematics and Statistics* **9** (2023). [DOI](https://doi.org/10.3389/fams.2023.1165371).
104. G. H. Golub and J. H. Welsch, “Calculation of Gauss quadrature rules,” *Math. Comp.* **23** (1969), 221–230. [DOI](https://doi.org/10.1090/S0025-5718-69-99647-1).
105. N. Hale and A. Townsend, “A fast FFT-based discrete Legendre transform,” *IMA J. Numer. Anal.* **36** (2016), 1670–1684. [DOI](https://doi.org/10.1093/imanum/drv060).
106. S. Olver, R. M. Slevinsky and A. Townsend, “Fast algorithms using orthogonal polynomials,” *Acta Numerica* **29** (2020), 573–699. [Publisher](https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D).
107. R. M. Slevinsky, “Fast and backward stable transforms between spherical harmonic expansions and bivariate Fourier series,” *Appl. Comput. Harmon. Anal.* **47** (2019). [Preprint](https://arxiv.org/abs/1705.05448).
108. Ž. Gimbutas and L. Greengard, “A fast and stable method for rotating spherical harmonic expansions,” *J. Comput. Phys.* **228** (2009), 5621–5627. [DOI](https://doi.org/10.1016/j.jcp.2009.05.014).
109. M. D. Webb, “Stable Hermite transforms via the Golub--Welsch algorithm,” arXiv:2604.02041 (2026). [Preprint](https://arxiv.org/abs/2604.02041).
110. N. Mukunda and R. Simon, “Quantum Kinematic Approach to the Geometric Phase. I. General Formalism,” *Ann. Phys.* **228** (1993), 205–268. [DOI](https://doi.org/10.1006/aphy.1993.1093).
111. T. Bendokat, R. Zimmermann and P.-A. Absil, “A Grassmann manifold handbook: Basic geometry and computational aspects,” *Adv. Comput. Math.* **50** (2024), article 6. [DOI](https://doi.org/10.1007/s10444-023-10090-8). See published correction: [DOI](https://doi.org/10.1007/s10444-024-10115-w).
112. P. J. Kostelec and D. N. Rockmore, “FFTs on the Rotation Group,” *J. Fourier Anal. Appl.* **14** (2008), 145–179. [DOI](https://doi.org/10.1007/s00041-008-9013-5).
113. T. Needham, *Visual Differential Geometry and Forms: A Mathematical Drama in Five Acts*, Princeton University Press, 2021, especially its parallel-transport and holonomy treatment. [DOI](https://doi.org/10.1515/9780691219899); [publisher](https://press.princeton.edu/books/paperback/9780691203706/visual-differential-geometry-and-forms).
114. J. Voight, “Hyperbolic plane,” Chapter 33 of *Quaternion Algebras*, Springer, 2021; Section 33.6 on hyperbolic area and Gauss–Bonnet. [DOI](https://doi.org/10.1007/978-3-030-56694-4_33).
115. A. Ipp and D. Müller, “Implicit schemes for real-time lattice gauge theory,” *European Physical Journal C* **78** (2018), 884. [DOI](https://doi.org/10.1140/epjc/s10052-018-6323-x); [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC6315017/).

116. S. Bravyi, “Lagrangian representation for fermionic linear optics,” *Quantum Information and Computation* **5** (2005), 216–238. [DOI](https://doi.org/10.26421/QIC5.3-3); [preprint](https://arxiv.org/abs/quant-ph/0404180).

117. S. Scheel, “Permanents in linear optical networks,” arXiv:quant-ph/0406127 (2004). [Preprint](https://arxiv.org/abs/quant-ph/0406127).

118. M. Reck, A. Zeilinger, H. J. Bernstein and P. Bertani, “Experimental realization of any discrete unitary operator,” *Phys. Rev. Lett.* **73** (1994), 58–61. [DOI](https://doi.org/10.1103/PhysRevLett.73.58).
119. W. R. Clements, P. C. Humphreys, B. J. Metcalf, W. S. Kolthammer and I. A. Walmsley, “Optimal design for universal multiport interferometers,” *Optica* **3** (2016), 1460–1465. [DOI](https://doi.org/10.1364/OPTICA.3.001460).

120. P. P. Vaidyanathan, *Multirate Systems and Filter Banks*, Prentice-Hall, 1993. [Publisher](https://www.informit.com/store/multirate-systems-and-filter-banks-9780136057185).

121. R. E. A. C. Paley, “On Orthogonal Matrices,” *Journal of Mathematics and Physics* **12** (1933), 311–320. [DOI](https://doi.org/10.1002/sapm1933121311).

122. E. Celledoni, A. Marthinsen and B. Owren, “Commutator-free Lie group methods,” *Future Generation Computer Systems* **19** (2003), 341–352. [DOI](https://doi.org/10.1016/S0167-739X(02)00161-9).

123. A. Singer and H.-T. Wu, “Vector diffusion maps and the connection Laplacian,” *Communications on Pure and Applied Mathematics* **65** (2012), 1067–1144. [DOI](https://doi.org/10.1002/cpa.21395).

124. F. Wilczek and A. Zee, “Appearance of Gauge Structure in Simple Dynamical Systems,” *Physical Review Letters* **52** (1984), 2111–2114. [DOI](https://doi.org/10.1103/PhysRevLett.52.2111).

125. J. Qian and W.-W. Lin, “A numerical method for quadratic eigenvalue problems of gyroscopic systems,” *Journal of Sound and Vibration* **306** (2007), 284–296. [DOI](https://doi.org/10.1016/j.jsv.2007.05.009).
126. T. Ohsawa, “Preservation of Quadratic Invariants by Semiexplicit Symplectic Integrators for Nonseparable Hamiltonian Systems,” *SIAM Journal on Numerical Analysis* **61** (2023), 1293–1315. [DOI](https://doi.org/10.1137/22M1517718).
127. I. Kuzmanović Ivičić and S. Miodragović, “Perturbation bounds for stable gyroscopic systems in motion about an unstable equilibrium position,” *Journal of Computational and Applied Mathematics* **451** (2024), 116061. [DOI](https://doi.org/10.1016/j.cam.2024.116061).


128. J. Dehaene and B. De Moor, “Clifford group, stabilizer states, and linear and quadratic operations over GF(2),” *Physical Review A* **68** (2003), 042318. [DOI](https://doi.org/10.1103/PhysRevA.68.042318).
129. D. Gross, S. Nezami and M. Walter, “Schur–Weyl duality for the Clifford group with applications: Property testing, a robust Hudson theorem, and de Finetti representations,” *Communications in Mathematical Physics* **385** (2021), 1325–1393. [DOI](https://doi.org/10.1007/s00220-021-04118-7).


130. S. Trebst, M. Troyer, Z. Wang and A. W. W. Ludwig, “A short introduction to Fibonacci anyon models,” *Progress of Theoretical Physics Supplement* **176** (2008), 384–407. [DOI](https://doi.org/10.1143/PTPS.176.384); [preprint](https://arxiv.org/abs/0902.3275).
131. C. Carnahan, D. Zeuch and N. E. Bonesteel, “Systematically generated two-qubit anyon braids,” *Physical Review A* **93** (2016), 052328. [DOI](https://doi.org/10.1103/PhysRevA.93.052328); [preprint](https://arxiv.org/abs/1511.00719).


132. C. W. J. Beenakker, “Random-matrix theory of quantum transport,” *Reviews of Modern Physics* **69** (1997), 731–808. [DOI](https://doi.org/10.1103/RevModPhys.69.731); [preprint](https://arxiv.org/abs/cond-mat/9612179).
133. P. W. Brouwer, K. M. Frahm and C. W. J. Beenakker, “Quantum Mechanical Time-Delay Matrix in Chaotic Scattering,” *Physical Review Letters* **78** (1997), 4737–4740. [DOI](https://doi.org/10.1103/PhysRevLett.78.4737).
134. L. Erdős, T. Krüger and Y. Nemish, “Scattering in Quantum Dots via Noncommutative Rational Functions,” *Annales Henri Poincaré* **22** (2021), 4205–4269. [DOI](https://doi.org/10.1007/s00023-021-01085-6); [preprint](https://arxiv.org/abs/1911.05112).


135. A. C. Potter, T. Morimoto and A. Vishwanath, “Classification of Interacting Topological Floquet Phases in One Dimension,” *Physical Review X* **6** (2016), 041001. [DOI](https://doi.org/10.1103/PhysRevX.6.041001); [preprint](https://arxiv.org/abs/1602.05194).
136. D. V. Else and C. Nayak, “Classification of topological phases in periodically driven interacting systems,” *Physical Review B* **93** (2016), 201103(R). [DOI](https://doi.org/10.1103/PhysRevB.93.201103); [preprint](https://arxiv.org/abs/1602.04804).
137. X. Zhang, W. Jiang, J. Deng et al., “Digital quantum simulation of Floquet symmetry-protected topological phases,” *Nature* **607** (2022), 468–473. [DOI](https://doi.org/10.1038/s41586-022-04854-3); [open access](https://pmc.ncbi.nlm.nih.gov/articles/PMC9300455/).


138. J. Moser and A. P. Veselov, “Discrete versions of some classical integrable systems and factorization of matrix polynomials,” *Communications in Mathematical Physics* **139** (1991), 217–243. [DOI](https://doi.org/10.1007/BF02352494); [full text](https://projecteuclid.org/journals/communications-in-mathematical-physics/volume-139/issue-2/Discrete-versions-of-some-classical-integrable-systems-and-factorization-of/cmp/1104203302.full).
139. J. R. Cardoso and P. Miraldo, “Solving the discrete Euler–Arnold equations for the generalized rigid body motion,” *Journal of Computational and Applied Mathematics* **402** (2022), 113814. [DOI](https://doi.org/10.1016/j.cam.2021.113814); [preprint](https://arxiv.org/abs/2109.00505).
140. X. Shen, K. Tran and M. Leok, “High-order symplectic Lie group methods on SO(n) using the polar decomposition,” *Journal of Computational Dynamics* **9** (2022), 529–551. [DOI](https://doi.org/10.3934/jcd.2022003); [preprint](https://arxiv.org/abs/2201.10768).


141. M. Loveland, E. Valseth, M. Lukac and C. Dawson, “Extending FEniCS to work in higher dimensions using tensor product finite elements,” *Journal of Computational Science* **64** (2022), 101831. [DOI](https://doi.org/10.1016/j.jocs.2022.101831); [preprint](https://arxiv.org/abs/2202.00762).


142. Y. Nagano, S. Yamaguchi, K. Fujita and M. Koyama, “A Wrapped Normal Distribution on Hyperbolic Space for Gradient-Based Learning,” *Proceedings of Machine Learning Research* **97** (ICML 2019), 4693–4702. [Open proceedings](https://proceedings.mlr.press/v97/nagano19a.html); [preprint](https://arxiv.org/abs/1902.02992).


143. R. G. Campos, J. Rico-Melgoza and E. Chávez, “A New Formulation of the Fast Fractional Fourier Transform,” *SIAM Journal on Scientific Computing* **34** (2012), A1110–A1125. [DOI](https://doi.org/10.1137/100812677).
144. C. Candan, M. A. Kutay and H. M. Ozaktas, “The Discrete Fractional Fourier Transform,” *IEEE Transactions on Signal Processing* **48** (2000), 1329–1337. [DOI](https://doi.org/10.1109/78.839980).
145. M. J. Colbrook, A. Horning and A. Townsend, “Computing Spectral Measures of Self-Adjoint Operators,” *SIAM Review* **63** (2021), 489–524. [DOI](https://doi.org/10.1137/20M1330944); [preprint](https://arxiv.org/abs/2006.01766).
146. M. J. Colbrook and A. Townsend, “Rigorous Data-Driven Computation of Spectral Properties of Koopman Operators for Dynamical Systems,” *Communications on Pure and Applied Mathematics* **77** (2024), 221–283. [DOI](https://doi.org/10.1002/cpa.22125); [preprint](https://arxiv.org/abs/2111.14889).
147. M. J. Colbrook, “Computing Spectral Measures and Spectral Types,” *Communications in Mathematical Physics* **384** (2021), 433–501. [DOI](https://doi.org/10.1007/s00220-021-04072-4); [author PDF](https://www.damtp.cam.ac.uk/user/mjc249/pdfs/SpecMeasuresMJC.pdf).
148. M. J. Colbrook, “The mpEDMD Algorithm for Data-Driven Computations of Measure-Preserving Dynamical Systems,” *SIAM Journal on Numerical Analysis* **61** (2023), 1585–1608. [DOI](https://doi.org/10.1137/22M1521407); [preprint](https://arxiv.org/abs/2209.02244).


149. C. Cedzich and A. H. Werner, “Anderson Localization for Electric Quantum Walks and Skew-Shift CMV Matrices,” *Communications in Mathematical Physics* **387** (2021), 1257–1279. [DOI](https://doi.org/10.1007/s00220-021-04204-w); [preprint](https://arxiv.org/abs/1906.11931).
150. F. Yang, “Interval Spectrum for Electric Quantum Walk and Related Skew-Shift CMV Matrices,” *Letters in Mathematical Physics* **115** (2025), article 49. [DOI](https://doi.org/10.1007/s11005-025-01934-z); [preprint](https://arxiv.org/abs/2408.12724).

151. G. S. Ammar and W. B. Gragg, *Schur Flows for Orthogonal Hessenberg Matrices* (1993), Naval Postgraduate School report. [Full text](https://calhoun.nps.edu/bitstream/handle/10945/30447/Ammar_Gragg_Schur.pdf?sequence=3).
152. L. B. Golinskii, “Schur flows and orthogonal polynomials on the unit circle,” *Sbornik: Mathematics* **197** (2006), 1145–1165. [DOI](https://doi.org/10.1070/SM2006v197n08ABEH003792).
153. W. W. Symes, “The QR algorithm and scattering for the finite nonperiodic Toda lattice,” *Physica D* **4** (1982), 275–280. [DOI](https://doi.org/10.1016/0167-2789(82)90069-0).
154. A. Mukaihira and Y. Nakamura, “Schur flow for orthogonal polynomials on the unit circle and its integrable discretization,” *J. Comput. Appl. Math.* **139** (2002), 75–94. [DOI](https://doi.org/10.1016/S0377-0427(01)00388-0).
155. M. J. Cantero et al., “Darboux transformations for CMV matrices,” *Adv. Math.* (2016). [DOI](https://doi.org/10.1016/j.aim.2016.03.042).

156. F. J. Simons, F. A. Dahlen and M. A. Wieczorek, “Spatiospectral Concentration on a Sphere,” *SIAM Review* **48** (2006), 504–536. [DOI](https://doi.org/10.1137/S0036144504445765); [preprint](https://arxiv.org/abs/math/0408424).

157. B. Leistedt and J. D. McEwen, “Exact Wavelets on the Ball,” *IEEE Transactions on Signal Processing* **60** (2012). [DOI](https://doi.org/10.1109/TSP.2012.2215030); [preprint](https://arxiv.org/abs/1205.0792).
158. J. D. McEwen and B. Leistedt, “Fourier-Laguerre Transform, Convolution and Wavelets on the Ball,” *Sampling Theory and Applications (SampTA)* (2013), 329–332. [Preprint](https://arxiv.org/abs/1307.1307).

159. I. Nenciu, “Lax pairs for the Ablowitz-Ladik system via orthogonal polynomials on the unit circle,” *International Mathematics Research Notices* **2005** (2005), 647–686. [DOI](https://doi.org/10.1155/IMRN.2005.647); [preprint](https://arxiv.org/abs/math-ph/0412047).
160. M. Gekhtman and I. Nenciu, “Multi-Hamiltonian structure for the finite defocusing Ablowitz-Ladik equation,” *Proceedings of the American Mathematical Society* **136** (2008), 2785–2792. [Preprint](https://arxiv.org/abs/0706.2428).

161. Ç. Candan, M. A. Kutay and H. M. Ozaktas, “The discrete fractional Fourier transform,” *IEEE Transactions on Signal Processing* **48** (2000), 1329–1337. [DOI](https://doi.org/10.1109/78.839980); [full paper](https://www.ee.bilkent.edu.tr/~haldun/publications/ozaktas166.pdf).

162. N. Van Buggenhout, M. Van Barel and R. Vandebril, “Generation of orthogonal rational functions by procedures for structured matrices,” *Numerical Algorithms* **89** (2022). [DOI](https://doi.org/10.1007/s11075-021-01125-6); [preprint](https://arxiv.org/abs/2103.04788).

163. D. Camps, T. Mach, R. Vandebril and D. S. Watkins, *Pole-Swapping Algorithms for the Eigenvalue Problem*, SIAM Spotlights (2025), Chapters 4–5. [DOI](https://doi.org/10.1137/1.9781611978377).
164. D. Camps, T. Mach, R. Vandebril and D. S. Watkins, “On pole-swapping algorithms for the eigenvalue problem,” *Electronic Transactions on Numerical Analysis* **52** (2020), 480–508. [DOI](https://doi.org/10.1553/etna_vol52s480); [preprint](https://arxiv.org/abs/1906.08672).

165. B. P. Huffenberger and B. D. Wandelt, “Fast and Exact Spin-s Spherical Harmonic Transforms,” *Astrophysical Journal Supplement Series* **189** (2010), 255–260. [DOI](https://doi.org/10.1088/0067-0049/189/2/255); [preprint](https://arxiv.org/abs/1007.3514).

166. J. D. McEwen, M. Büttner, B. Leistedt, H. V. Peiris and Y. Wiaux, “A Novel Sampling Theorem on the Rotation Group,” *IEEE Signal Processing Letters* **22** (2015), 2425–2429. [DOI](https://doi.org/10.1109/LSP.2015.2490676); [preprint](https://arxiv.org/abs/1508.03101).
167. Z. Khalid, S. Durrani, R. A. Kennedy, Y. Wiaux and J. D. McEwen, “Gauss–Legendre Sampling on the Rotation Group,” *IEEE Signal Processing Letters* **23** (2016), 207–211. [DOI](https://doi.org/10.1109/LSP.2015.2503295); [preprint](https://arxiv.org/abs/1508.03353).

168. D. Nguyen and S. Sommer, “Parallel Transport on Matrix Manifolds and Exponential Action,” *SIAM J. Matrix Anal. Appl.* **46** (2025). [DOI](https://doi.org/10.1137/24M1684189); [preprint](https://arxiv.org/abs/2408.06054).
