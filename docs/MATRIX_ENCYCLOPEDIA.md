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

For a unitary $U\in\mathbb C^{n\times n}$, eigenvalues are written
$z_j=e^{i\theta_j}$. The **empirical spectral measure** and its moments are

$$
\mu_U^{\mathrm{trace}}=\frac1n\sum_{j=1}^n\delta_{z_j},\qquad
m_k(U)=\frac1n\mathop{\mathrm{tr}}(U^k).
$$

For a unit vector $v$, the different **vector spectral measure** is

$$
\mu_{U,v}=\sum_{j=1}^n |q_j^*v|^2\delta_{z_j},
$$

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

**Order:** $n$. **Structure:** real, orthogonal, sparse permutation,
circulant. **Spectrum:** precisely the $n$th roots of unity.

$$
P_ne_j=e_{j+1\pmod n},\qquad P_n^n=I,\qquad
\frac1n\mathop{\mathrm{tr}}(P_n^k)=\mathbf1_{n\mid k}.
$$

Useful as an exact uniform-phase reference and to detect phase aliasing in
Fourier moment methods. It is a single-cycle permutation rather than a
generic permutation with multiple cycle lengths. **Reference:** [1].

### 2. Phase-weighted cyclic shift — [`weighted_shift.m`](../weighted_shift.m)

**Order:** $n$. **Structure:** sparse monomial unitary; becomes a real signed
permutation only for phases $0$ or $\pi$. For real weights $\phi_j$,

$$
Ue_j=e^{i\phi_j}e_{j+1\pmod n},\qquad
U^n=e^{i\sum_j\phi_j}I.
$$

Its spectrum comprises the roots of $z^n=e^{i\sum_j\phi_j}$. A useful
counterexample to the idea that changing the entries of a sparse unitary
necessarily changes its relative eigenvalue spacings: changing the phases
only rotates the perfectly spaced cycle spectrum. **Reference:** [1].

### 3. Prescribed diagonal phases — [`diagonal_phase.m`](../diagonal_phase.m)

**Order:** $n$. **Structure:** sparse diagonal unitary. For real phase inputs
$\theta_j$,

$$
D=\mathop{\mathrm{diag}}(e^{i\theta_1},\ldots,e^{i\theta_n}).
$$

The spectrum, algebraic multiplicities, and trace moments are known exactly.
Use clustered phases, open gaps, isolated atoms, or an irregular deterministic
phase mesh to test spectral estimators independently of the eigensolver.
The default gives equally spaced roots of unity. **Reference:** [1].

### 4. Unitary circulant with prescribed phases — [`phase_circulant.m`](../phase_circulant.m)

**Order:** $n$. **Structure:** dense circulant with FFT-diagonalization.
For the normalized DFT $F$ and real phases $\theta_j$,

$$
U=F^*\mathop{\mathrm{diag}}(e^{i\theta_j})F.
$$

It has **the same eigenvalues as the diagonal control**, but generally
nontrivial eigenvectors and dense entries. This pairs well with entrywise or
matvec algorithms, exposing accidental dependence on eigenvector localization.
Our implementation uses an inverse FFT of the Fourier eigenvalue samples.
**Reference:** [1].

### 5. Periodic tight-binding propagator — [`cycle_propagator.m`](../cycle_propagator.m)

**Order:** $n$. **Structure:** dense circulant. Write $H=P_n+P_n^*$; then

$$
U(t)=e^{-itH},\quad
z_k=e^{-2it\cos(2\pi k/n)},\quad k=0,\ldots,n-1.
$$

The infinite-chain energy law is arcsine on $[-2,2]$, and its time-scaled
pushforward onto the circle produces a known nonuniform phase law. Large $|t|$
wraps the energy bands and can create overlapping singularities. This is a
strong test for circular spectral densities with analytic reference phases.
The $n=1,2$ cases retain the literal $H=P+P^*$ definition. **Reference:** [2].

### 6. Real orthogonal rotation blocks — [`orthogonal_blocks.m`](../orthogonal_blocks.m)

**Order:** $n$. **Structure:** sparse real block diagonal, orthogonal and
normal. Each $2\times2$ block is

$$
R(\theta)=\begin{pmatrix}\cos\theta&-\sin\theta\\
\sin\theta&\cos\theta\end{pmatrix},\qquad
\sigma(R)=\{e^{i\theta},e^{-i\theta}\}.
$$

An odd-dimensional matrix also has a chosen $+1$ or $-1$ terminal block.
Prescribe clusters, repeated conjugate pairs, determinant parity and eigenvalue
multiplicities. This is a controlled test for real-normal eigenvalue
algorithms and conjugation-symmetric phase density. **References:** [3, 4].

## B. OPUC, Szegő measures, and core-factorized representations

### 7. Finite CMV — [`finite_cmv.m`](../finite_cmv.m)

**Order:** $n$. **Structure:** sparse five-diagonal unitary. From
$\alpha_0,\ldots,\alpha_{n-2}$ in the open unit disk and a terminal
$\beta\in\mathbb T$, use the local core

$$
\Theta(\alpha)=
\begin{pmatrix}\overline\alpha&\sqrt{1-|\alpha|^2}\\
\sqrt{1-|\alpha|^2}&-\alpha\end{pmatrix}
$$

and alternating even/odd block-diagonal factors $\mathcal C=\mathcal L\mathcal M$.
The finite matrix closes with the terminal boundary phase. Its eigenvalues
are the zeros of a paraorthogonal polynomial on the circle. Real interior
coefficients and a real closure produce an orthogonal CMV matrix.
**References:** [5–7].

### 8. General Verblunsky multiplication matrix — [`opuc_matrix.m`](../opuc_matrix.m)

**Order:** $n$. **Representations:** `hessenberg`, `inverse_hessenberg`,
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

**Order:** $n$. **Source measure:** the Jacobi probability measure
$(1-t)^a(1+t)^b\,dt$ on $[-1,1]$, $a,b>-1$, lifted by $t=\cos\theta$
into a conjugation-invariant circle measure. The explicit coefficients are

$$
\alpha_{2k}=\frac{b-a}{2k+a+b+2},\qquad
\alpha_{2k+1}=-\frac{a+b+1}{2k+a+b+3}.
$$

Take the first $n-1$ and append a unit-modulus closing parameter. The
same measure can produce Hessenberg, inverse Hessenberg, CMV or twisted
matrices. Classical choices include Chebyshev I ($a=b=-1/2$), Legendre
($a=b=0$), Chebyshev II ($a=b=1/2$), and Gegenbauer
($a=b=\lambda-1/2$). The finite Gauss–Szegő nodes and weights are not
the original continuous measure. **References:** [5, 7, 11, 12].

### 10. Bernstein–Szegő / Poisson measure — [`bernstein_szego.m`](../bernstein_szego.m)

**Order:** $n$. **Source:** an explicit absolutely continuous probability
measure on the unit circle, with parameter $|a|<1$ and density (with respect
to $d\theta$)

$$
\frac{d\mu}{d\theta}=
\frac{1-|a|^2}{2\pi|e^{i\theta}-\overline a|^2}.
$$

The coefficients in the stated OPUC convention are $\alpha_0=a$ and
$\alpha_k=0$ for $k\geq1$, giving the exact moment
$\int z^k\,d\mu(z)=\overline a^{\,k}$ for $k\geq0$.
This provides an analytically transparent **nonuniform** weight with
rapidly terminating Verblunsky data, a useful dual check of the measure
and multiplication matrix. All four existing representation choices work.
**References:** [5, 12].

### 11. Circular-$\beta$ random CMV — [`circular_beta.m`](../circular_beta.m)

**Order:** $n$. **Structure:** sparse five-diagonal random unitary.
The Killip–Nenciu construction chooses independent interior coefficients
with uniform angle and

$$
|\alpha_{j}|^2\sim\mathop{\mathrm{Beta}}\!\left(
1,\frac{\beta(n-j-1)}{2}\right),\qquad 0\leq j\leq n-2,
$$

and an independent uniform terminal phase. The eigenangles have joint
density proportional to $\prod_{j<k}|e^{i\theta_j}-e^{i\theta_k}|^\beta$.
At $\beta=2$ the **eigenvalue law** is CUE, although the matrix is not
Haar-distributed on $U(n)$. This tests level repulsion while preserving
$O(n)$ sparse structure. The RNG-seeded option restores the caller's state.
No Statistics Toolbox is required. **Reference:** [13].

## C. Random and orthogonal signal-processing matrices

### 12. Haar unitary — [`haar_unitary.m`](../haar_unitary.m)

**Order:** $n$. **Structure:** dense, random. Phase-adjusted QR of a
complex Gaussian matrix gives Haar distribution on $U(n)$.
For the circular unitary ensemble (CUE), the **expected** one-point phase
density is uniform, but finite samples exhibit correlated phases with
repulsion; individual draws are not equally spaced. This is an ensemble
benchmark rather than a deterministic known-spectrum matrix. **Reference:** [14].

### 13. Haar orthogonal — [`haar_orthogonal.m`](../haar_orthogonal.m)

**Order:** $n$. **Structure:** dense, real orthogonal. Sign-adjusted QR of
an $n\times n$ real Gaussian matrix produces Haar distribution on $O(n)$.
Its eigenvalues obey complex-conjugate pairing and determinant/parity
constraints. **Haar on $O(n)$ is not Dyson's COE.** This distinction matters
when testing random matrix universality claims. **Reference:** [14].

### 14. Haar wavelet transform — [`haar_wavelet.m`](../haar_wavelet.m)

**Order:** $n=2^p$. **Structure:** sparse, real orthogonal, multiresolution.
The first analysis step forms adjacent sums/differences divided by
$\sqrt2$; subsequent levels act only on the coarse half. The `levels`
parameter ranges from zero (identity) to $p$ (full decomposition).
It is a signal-processing test of real orthogonality, highly nonuniform
sparsity and multiscale eigenvectors. **The eigenphases of a finite analysis
matrix are not the wavelet frequencies.** **Reference:** [15].

## D. Quantum maps, walks, and algorithms

### 15. Split-step quantum walk — [`split_step_walk.m`](../split_step_walk.m)

**Order:** $2n$ for $n$ sites. **Structure:** sparse, real orthogonal in the
selected rotation-coin convention. The periodic one-step operator is

$$
U=S_-C(\theta_2)S_+C(\theta_1).
$$

It has a two-state internal coin and periodic conditional shifts; the
basis alternates coin states site by site. Parameters change the two
quasienergy branches and potential gaps. The current periodic uniform-coin
version has no physical edge, so **edge states require a separate open or
interface construction**. **Reference:** [16].

### 16. Periodic Dirac quantum cellular automaton — [`dirac_automaton.m`](../dirac_automaton.m)

**Order:** $2n$ for $n$ sites. **Structure:** sparse complex unitary.
With mass $m\in[0,1]$, $\nu=\sqrt{1-m^2}$ and periodic cyclic shift $P_n$,

$$
U=\begin{pmatrix}\nu P_n&-imI\\-imI&\nu P_n^*\end{pmatrix},\qquad
\cos\omega_k=\nu\cos(2\pi k/n).
$$

The exact eigenvalues are $e^{\pm i\omega_k}$. For $0<m<1$, the large-$n$
phase density has gapped support and integrable band-edge singularities:

$$
\rho(\theta)=\begin{cases}
\dfrac{|\sin\theta|}{2\pi\sqrt{\nu^2-\cos^2\theta}},
  &|\cos\theta|<\nu,\\
0,&\text{otherwise}.
\end{cases}
$$

Thus this is a particularly useful nonuniform, exact-density benchmark.
At $m=0$ it reduces to two oppositely directed shifts; at $m=1$ its
spectrum is concentrated at $\pm i$. **References:** [17, 18].

### 17. Grover search iterate — [`grover_operator.m`](../grover_operator.m)

**Order:** $n\geq2$. **Structure:** dense real orthogonal product of two
reflections. Mark the first $M$ basis states, $1\leq M<n$, and let
$|s\rangle=n^{-1/2}(1,\ldots,1)^T$. Then

$$
G=(2|s\rangle\langle s|-I)(I-2\Pi_M),\qquad
\sin^2\alpha=M/n.
$$

The exact spectrum is $e^{\pm2i\alpha}$ (one each), $+1$ with multiplicity
$M-1$, and $-1$ with multiplicity $n-M-1$. Two nontrivial eigenphases
may have only $2/n$ empirical spectral mass, so a density estimator can
miss them. **Reference:** [19].

### 18. Quantum baker map — [`quantum_baker.m`](../quantum_baker.m)

**Order:** even $n$. **Structure:** dense, product of Fourier transforms.
The Balazs–Voros quantization is

$$
B_n=F_n^*\mathop{\mathrm{diag}}(F_{n/2},F_{n/2}).
$$

The optional Saraceno variant uses half-integer/antiperiodic Fourier phases
and a different parity-symmetry convention. Eigenphase statistics are
parameter/size/convention dependent; this matrix **does not inherit the
four-atom eigenvalue spectrum of $F_n$**. **References:** [20, 21].

### 19. Kicked Harper torus map — [`kicked_harper.m`](../kicked_harper.m)

**Order:** $n$. **Structure:** dense explicit matrix, FFT-factorized
mathematical action. With $D_x(K)$ and $D_p(L)$ defined by diagonal
phases $e^{-iK\cos(2\pi j/n)}$ and $e^{-iL\cos(2\pi j/n)}$,

$$
U=D_p(L)F_n^*D_x(K)F_n.
$$

This fixed finite-torus convention is exactly unitary. Quasienergy
statistics, clustering and band structure depend on $K,L,n$ and the torus
quantization; statements about multifractality in irrational limits must
not be transferred indiscriminately to every finite test matrix.
**Reference:** [22].

## E. Matrix functions, polar factors, and dilations

### 20. Cayley image of a Jacobi matrix — [`cayley_jacobi.m`](../cayley_jacobi.m)

**Order:** $n$. **Structure:** dense complex unitary, derived from a real
symmetric tridiagonal Jacobi matrix $J_n$ with Jacobi weight parameters
$a,b>-1$. For scale $t>0$,

$$
U=(J_n-itI)(J_n+itI)^{-1},\qquad
\lambda_U=\frac{x-it}{x+it}\quad(x\in\sigma(J_n)).
$$

This is an explicit Möbius pushforward of a classical real-line spectrum,
providing controlled phase clustering near $+1$ and $-1$ as the scale
changes. **It is not the same finite matrix as the Jacobi–Szegő lift.**
**Reference:** [23].

### 21. Nonsingular polar factor — [`polar_unitary.m`](../polar_unitary.m)

**Order:** $n$. **Structure:** generally dense. If $A=L\Sigma R^*$ is a
full SVD of a nonsingular matrix, the uniquely determined unitary polar
factor is $U=LR^*$. For real $A$, $U$ is orthogonal. It is also a
nearest unitary matrix to $A$ in every unitarily invariant norm.

We restrict the generator to numerically nonsingular inputs: for singular
$A$, an SVD-based square unitary extension is not unique and depends on
choices in nullspaces. Useful for polar algorithms, structure restoration,
Procrustes problems and perturbation analysis. **Reference:** [24].

### 22. Julia/Halmos unitary extension — [`julia_unitary.m`](../julia_unitary.m)

**Order:** $2n$ for an $n\times n$ contraction $T$, $\lVert T\rVert_2\leq1$. 
**Structure:** generally dense, block unitary. With defect operators
$D_T=(I-T^*T)^{1/2}$ and $D_{T^*}=(I-TT^*)^{1/2}$,

```math
\mathcal{U}(T)=
\begin{pmatrix}
T & D_{T^*} \\
D_T & -T^*
\end{pmatrix}.
```

This gives a finite unitary containing the nonnormal contraction $T$ in
its leading block. **It is not automatically a unitary power dilation**:
compressing arbitrary powers of $\mathcal U(T)$ need not yield $T^k$.
Useful for conservative system matrices, passive control and dilation theory.
**Reference:** [25].


## F. Further real and complex OPUC matrix generators

### 23. Generalized circular Jacobi — [`circular_jacobi.m`](../circular_jacobi.m)

**Order:** $n$. **Structure:** unitary, sparse five-diagonal by default.
For real $a>-1/2$ and $\eta\in\mathbb R$, the positive circle weight is
proportional to $e^{(\pi-\theta)\eta}|1-e^{i\theta}|^{2a}$ for
$0<\theta<2\pi$. With $b=a+i\eta$, the Verblunsky sequence is

$
\alpha_k=-\frac{(b)_{k+1}}{(1+\overline b)_{k+1}}.
$

At $\eta=0$ all coefficients are real; otherwise they are generally
complex. Supports an explicitly prescribed nonuniform *vector* measure,
but does not promise the same limiting unweighted phase density.
The finite boundary $\beta$ is independent of the measure.
**References:** [5, 12, 27].

### 24. Rogers--Szego wrapped Gaussian — [`rogers_szego.m`](../rogers_szego.m)

**Order:** $n$. **Structure:** real orthogonal if the terminal boundary is
real. For $0<q<1$, the circle measure is a wrapped Gaussian with variance
$-\log q$, and its Verblunsky coefficients satisfy

$
\alpha_k=(-1)^k q^{(k+1)/2}.
$

Its positive circular density is smooth and the $k$th Fourier moment is
$q^{k^2/2}$. This is a nontrivial parameter-dependent real-coefficient
OPUC benchmark with a known cyclic-vector measure. **References:** [5, 12].

### 25. Constant Geronimus CMV — [`geronimus_cmv.m`](../geronimus_cmv.m)

**Order:** $n$. **Structure:** unitary five-diagonal by default;
orthogonal for real $a$ and real $\beta$. Use $\alpha_k=a$ for
all $k\ge0$, $|a|<1$. In the infinite model the measure can have a
proper absolutely continuous arc and a point mass. This generator
tests finite approximations to arc gaps with complex coefficients.
**References:** [5, 7].

### 26. Periodic Verblunsky CMV — [`periodic_cmv.m`](../periodic_cmv.m)

**Order:** $n$. **Structure:** sparse unitary. Repeat a nonempty
complex or real coefficient pattern $a_0,\ldots,a_{p-1}$ with all
$|a_j|<1$; close the first $n-1$ entries using a terminal unit
phase. The infinite periodic theory has bands and gaps governed by
a transfer-matrix discriminant. Boundary effects and finite eigenvalues
must be kept distinct from the infinite spectral measure.
**References:** [5, 7].

### 27. Finite generalized Bernstein--Szego — [`bernstein_szego_finite.m`](../bernstein_szego_finite.m)

**Order:** $n$. **Structure:** unitary, sparse five-diagonal by default.
Specify a finite complex coefficient list $a_0,\ldots,a_{m-1}$
strictly inside the disk and set $\alpha_k=0$ for $k\ge m$.
With the associated reversed polynomial $\Phi_m^*$, the positive
weight is

$
\frac{d\mu}{d\theta}=
\frac{\prod_{j=0}^{m-1}(1-|a_j|^2)}
{2\pi|\Phi_m^*(e^{i\theta})|^2}.
$

The previous one-coefficient Poisson family is a special case.
These rational weights give arbitrarily chosen real or complex
Schur parameters with an exact and evaluable cyclic-vector density.
**References:** [5, 12].

## G. Families reserved for future versions (not implemented)

- **ETNA core-factor stress tests:** the five experimental constructions in
  Aurentz–Mach–Vandebril–Watkins (2015), including the nearly reducible Type
  III case with heavily clustered eigenvalues. Will use shared Verblunsky/core
  machinery but separate benchmark presets. [8]
- **Rational and isometric Arnoldi:** pole sequences, orthonormal rational
  basis conventions, unitary closures and appropriate generalized matrices/
  pencils. Not every finite compression is unitary. [9, 10]
- **Other OPUC measures:** positive Pastro and generalized $q$-weight
  specializations, inverse spectral finite atomic measures, Aleksandrov
  families, skew-shift parameters, and matrix-valued/block CMV. [5, 12]
- **Additional signal-processing matrices:** normalized DCT/DST variants,
  Daubechies wavelets, paraunitary polyphase banks and lossless scattering.
  Boundary and frequency conventions must be specified before admitting a
  finite matrix as unitary. [15]
- **Larger quantum models:** Szegedy graph walks, quantum-graph scattering,
  Floquet Ising models, kicked tops, cat maps and interferometer meshes.
  Different symmetry sectors and truncations cannot be conflated. [26]

## Bibliography

1. P. J. Davis, *Circulant Matrices*, Wiley, 1979.
2. N. W. Ashcroft and N. D. Mermin, *Solid State Physics*, Holt, Rinehart and Winston, 1976.
3. S. Mataigne and K. A. Gallivan, “The eigenvalue decomposition of normal matrices by the skew-symmetric part,” *Electron. J. Linear Algebra* **42** (2026), 349–386. [DOI](https://doi.org/10.13001/ela.2026.9957).
4. J. L. Aurentz, T. Mach, L. Robol, R. Vandebril and D. S. Watkins, *Core-Chasing Algorithms for the Eigenvalue Problem*, SIAM, 2018. [DOI](https://doi.org/10.1137/1.9781611975345).
5. B. Simon, *Orthogonal Polynomials on the Unit Circle*, Vols. 1–2, AMS, 2005. [AMS](https://www.ams.org/books/coll/054.1/).
6. M. J. Cantero, L. Moral and L. Velázquez, “Five-diagonal matrices and zeros of orthogonal polynomials on the unit circle,” *Linear Algebra Appl.* **362** (2003), 29–56. [DOI](https://doi.org/10.1016/S0024-3795(02)00457-3).
7. B. Simon, “CMV matrices: Five years after,” *J. Comput. Appl. Math.* **208** (2007), 120–154. [DOI](https://doi.org/10.1016/j.cam.2006.10.033).
8. J. L. Aurentz, T. Mach, R. Vandebril and D. S. Watkins, “Fast and stable unitary QR algorithm,” *Electron. Trans. Numer. Anal.* **44** (2015), 327–341. [ETNA](https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341).
9. A. Bultheel, M. J. Cantero and R. Cruz-Barroso, “Matrix methods for quadrature formulas on the unit circle: A survey,” *J. Comput. Appl. Math.* **284** (2015), 78–100. [DOI](https://doi.org/10.1016/j.cam.2014.11.002).
10. W. B. Gragg, “The QR algorithm for unitary Hessenberg matrices,” *J. Comput. Appl. Math.* **16** (1986), 1–8. [DOI](https://doi.org/10.1016/0377-0427(86)90169-X).
11. B. Simon, “OPUC on one foot,” *Bull. Amer. Math. Soc.* **42** (2005), 431–460. [Full text](https://web.ma.utexas.edu/mp_arc/c/05/05-78.pdf).
12. NIST Digital Library of Mathematical Functions, [§18.33, Orthogonal Polynomials on the Unit Circle](https://dlmf.nist.gov/18.33).
13. R. Killip and I. Nenciu, “Matrix models for circular ensembles,” *Int. Math. Res. Not.* **2004**(50), 2665–2701. [DOI](https://doi.org/10.1155/S1073792804141597).
14. F. Mezzadri, “How to generate random matrices from the classical compact groups,” *Notices Amer. Math. Soc.* **54** (2007), 592–604. [Preprint](https://arxiv.org/abs/math-ph/0609050).
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
27. P. Bourgade, A. Nikeghbali and A. Rouault, “Circular Jacobi ensembles and deformed Verblunsky coefficients,” *Int. Math. Res. Not.* (2009). [DOI](https://doi.org/10.1093/imrn/rnp092).
