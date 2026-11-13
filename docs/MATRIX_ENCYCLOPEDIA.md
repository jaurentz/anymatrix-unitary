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

```math
\mu_U^{\mathrm{trace}}=\frac1n\sum_{j=1}^n\delta_{z_j},\qquad
m_k(U)=\frac1n\mathop{\mathrm{tr}}(U^k).
```

For a unit vector $v$, the different **vector spectral measure** is

```math
\mu_{U,v}=\sum_{j=1}^n |q_j^*v|^2\delta_{z_j},
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

**Order:** $n$. **Structure:** real, orthogonal, sparse permutation,
circulant. **Spectrum:** precisely the $n$th roots of unity.

```math
P_ne_j=e_{j+1}\quad(e_{n+1}=e_1),\qquad P_n^n=I,\qquad
\frac1n\mathop{\mathrm{tr}}(P_n^k)=\mathbf{1}_{n\mid k}.
```

Useful as an exact uniform-phase reference and to detect phase aliasing in
Fourier moment methods. It is a single-cycle permutation rather than a
generic permutation with multiple cycle lengths. **Reference:** [1].

### 2. Phase-weighted cyclic shift — [`weighted_shift.m`](../weighted_shift.m)

**Order:** $n$. **Structure:** sparse monomial unitary; becomes a real signed
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

**Order:** $n$. **Structure:** sparse diagonal unitary. For real phase inputs
$\theta_j$,

```math
D=\mathop{\mathrm{diag}}(e^{i\theta_1},\ldots,e^{i\theta_n}).
```

The spectrum, algebraic multiplicities, and trace moments are known exactly.
Use clustered phases, open gaps, isolated atoms, or an irregular deterministic
phase mesh to test spectral estimators independently of the eigensolver.
The default gives equally spaced roots of unity. **Reference:** [1].

### 4. Unitary circulant with prescribed phases — [`phase_circulant.m`](../phase_circulant.m)

**Order:** $n$. **Structure:** dense circulant with FFT-diagonalization.
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

**Order:** $n$. **Structure:** dense circulant. Write $H=P_n+P_n^*$; then

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

**Order:** $n$. **Structure:** sparse real block diagonal, orthogonal and
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

### 7. Finite CMV — [`finite_cmv.m`](../finite_cmv.m)

**Order:** $n$. **Structure:** sparse five-diagonal unitary. From
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

**Order:** $n$. **Source:** an explicit absolutely continuous probability
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

### 11. Circular-$\beta$ random CMV — [`circular_beta.m`](../circular_beta.m)

**Order:** $n$. **Structure:** sparse five-diagonal random unitary.
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

```math
U=S_-C(\theta_2)S_+C(\theta_1).
```

It has a two-state internal coin and periodic conditional shifts; the
basis alternates coin states site by site. Parameters change the two
quasienergy branches and potential gaps. The current periodic uniform-coin
version has no physical edge, so **edge states require a separate open or
interface construction**. **Reference:** [16].

### 16. Periodic Dirac quantum cellular automaton — [`dirac_automaton.m`](../dirac_automaton.m)

**Order:** $2n$ for $n$ sites. **Structure:** sparse complex unitary.
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

**Order:** $n\geq2$. **Structure:** dense real orthogonal product of two
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

**Order:** even $n$. **Structure:** dense, product of Fourier transforms.
The Balazs–Voros quantization is

```math
B_n=F_n^*\mathop{\mathrm{diag}}(F_{n/2},F_{n/2}).
```

The optional Saraceno variant uses half-integer/antiperiodic Fourier phases
and a different parity-symmetry convention. Eigenphase statistics are
parameter/size/convention dependent; this matrix **does not inherit the
four-atom eigenvalue spectrum of $F_n$**. **References:** [20, 21].

### 19. Kicked Harper torus map — [`kicked_harper.m`](../kicked_harper.m)

**Order:** $n$. **Structure:** dense explicit matrix, FFT-factorized
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

**Order:** $n$. **Structure:** dense complex unitary, derived from a real
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
compressing arbitrary powers of $\mathcal{U}(T)$ need not yield $T^k$.
Useful for conservative system matrices, passive control and dilation theory.
**Reference:** [25].


## F. Further real and complex OPUC matrix generators

### 23. Generalized circular Jacobi — [`circular_jacobi.m`](../circular_jacobi.m)

**Order:** $n$. **Structure:** unitary, sparse five-diagonal by default.
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

**Order:** $n$. **Structure:** real orthogonal if the terminal boundary is
real. For $0<q<1$, the circle measure is a wrapped Gaussian with variance
$-\log q$, and its Verblunsky coefficients satisfy

```math
\alpha_k=(-1)^k q^{(k+1)/2}.
```

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

```math
\frac{d\mu}{d\theta}=
\frac{\prod_{j=0}^{m-1}(1-|a_j|^2)}
{2\pi|\Phi_m^*(e^{i\theta})|^2}.
```

The previous one-coefficient Poisson family is a special case.
These rational weights give arbitrarily chosen real or complex
Schur parameters with an exact and evaluable cyclic-vector density.
**References:** [5, 12].

## G. Floquet spin chains and driven fermionic systems

### 28. Quantum kicked top — [`kicked_top.m`](../kicked_top.m)

**Order:** $n=2j+1$ (integer or half-integer spin $j>0$). **Structure:**
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
exponential is generally dense. No uniformly $O(n)$-time action is claimed.
**References:** [28, 29].

### 29. Integrable kicked Ising chain — [`floquet_ising.m`](../floquet_ising.m)

**Order:** $2^n$ for $n$ spins. **Structure:** a diagonal Ising layer
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

**Order:** $2^n$ for $n\ge3$ periodic spins, capped at $n=11$
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

**Order:** $2n$ for $n$ sites in Nambu basis
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

**Order:** $n$, with normalized negative-exponent Fourier matrix $F_n$.
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

**Order:** $2^n$ for $n$ qubits; dense matrices are capped at $n=10$.
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

**Order:** $2^n$ for $n$ spin-$1/2$ sites, capped at $n=10$ for
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

## I. Bosonic, cavity-QED, and quantum Hall network operators

### 35. Kicked Bose–Hubbard dimer — [`kicked_bose_hubbard.m`](../kicked_bose_hubbard.m)

**Order:** $n+1$ for $n$ identical bosons on **two sites**, in the
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

**Order:** $2(n+1)$ for a photon cutoff at occupation $n$.
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

**Order:** $4n_xn_y$, with directed links $(x,y)$ on the periodic
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

## J. Non-quantum signal-processing transforms

This section gives independent finite matrix constructions. They are distinct
from the FFT, wavelet and paraunitary filter-bank families discussed earlier;
where a new generator is a specialized Fourier conjugation, the spectral
redundancy is explicit. Frequencies of a **signal** are not to be mistaken for
eigenphases of the **transform matrix**.

### Orthonormal cosine transforms — [`orthogonal_dct.m`](../orthogonal_dct.m)

**Order:** $n$. The parameter `kind` selects DCT-I through DCT-IV.
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

**Order:** $n$. The orthonormal kernel is

```math
H_{jk}=\frac{\cos(2\pi jk/n)+\sin(2\pi jk/n)}{\sqrt n}.
```

This is real, symmetric, orthogonal, and $H^2=I$.
The entire spectrum is in $\{-1,+1\}$; exact multiplicities may be
determined from $n$ and the trace. Useful as a two-atom spectrum in
real signal processing, independent of the Hadamard construction. [47]

### Chu/CAZAC unitary circulants — [`cazac_circulant.m`](../cazac_circulant.m)

**Order:** $n$. Choose an integer root $r$ coprime to $n$, and set

```math
c_k=\begin{cases}
\exp[-\pi i r k^2/n],&n\text{ even},\\
\exp[-\pi i r k(k+1)/n],&n\text{ odd}.
\end{cases}
\qquad U_{jk}=\frac{c_{(j-k)\bmod n}}{\sqrt n}.
```

Chu's perfect periodic autocorrelation yields $U^*U=I$,
with **constant-modulus entries** $|U_{jk}|=n^{-1/2}$.
The exact eigenvalues, including multiplicities, are
$z_\ell=n^{-1/2}\mathop{\mathrm{fft}}(c)_\ell$.
This is a restricted, application-derived subfamily of
unitary circulants (not a new general matrix class).
Applications: channel synchronization, radar/sonar pulse coding. [48]

### Discrete Fresnel transform — [`fresnel_transform.m`](../fresnel_transform.m)

**Order:** $n$. For $d=(j-k)\bmod n$ and
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

**Order:** $n$. For real parameters $c_1,c_2$, define

```math
A=D(c_2)F_nD(c_1),\qquad
D(c)=\mathop{\mathrm{diag}}_{0\leq k<n}
(e^{-2\pi i c k^2}),
```

where $F_n$ is the normalized negative-exponent DFT.
This is an exactly unitary chirp-DFT-chirp multiplication.
At $c_1=c_2=0$ it reduces to the DFT, whose spectrum has
at most four atoms; away from that point the phase law can differ.
Application: affine frequency-division multiplexing in
delay-Doppler channels. Do not conflate modulation $A$ with the
generally **nonunitary** wireless channel. [50]

### Finite Zak transforms — [`finite_zak.m`](../finite_zak.m)

**Order:** $n=pq$, $p\mid n$. Use the indexing $j+pk$
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

**Order:** $n$. Define $T e_j=e_{j+1\bmod n}$ and
$M e_j=e^{2\pi i j/n}e_j$, with zero-based indices.
The generator returns $U=M^pT^q$ for integer $p,q$.
It is a sparse monomial unitary, and its exact eigenvalues
are determined by the cycles of translation $q$ and the
product of modulation weights around each cycle. Uses:
Gabor time-frequency shifts and Weyl–Heisenberg systems. [51,53]

### Periodic Hilbert transform completion — [`periodic_hilbert.m`](../periodic_hilbert.m)

**Order:** $n$. Start from the DFT multiplier $-i$ at
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

**Order:** $n$. Define the negative-exponent DFT by

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

**Order:** $n$. Type II uses the kernel

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

**Order:** $n=2^m$. Natural Sylvester order satisfies

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
butterfly and row permutations, with $O(n)$ storage.
Applications: spread-spectrum signals, coding and image
processing. [55]

### Slepian/DPSS orthogonal basis — [`dpss_transform.m`](../dpss_transform.m)

**Order:** $n$, with half-bandwidth $0<W<1/2$.
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
fast $O(n)$-storage transform merely by wrapping its dense
eigenvectors in a function handle. Applications: multitaper
estimation, finite-band concentration and leakage suppression. [56]

### Graph Fourier transform — [`graph_fourier.m`](../graph_fourier.m)

**Order:** $n$. For a symmetric, nonnegative zero-diagonal
adjacency $A$, let $L=\mathop{\mathrm{diag}}(A\mathbf{1})-A$.
If $L=V\Lambda V^T$, sorted by nondecreasing eigenvalues,
the transform is $U=V^T$. It is real orthogonal. The generator
uses a path graph by default and fixes repeated eigenspaces
by vertex-order projection/orthonormalization. **Graph
Laplacian frequencies are not the eigenvalues of $U$**.
This is **explicit only**: a generic weighted graph has no
automatic fast graph-Fourier transform. Applications: graph
filtering and irregular-domain data analysis. [57]

## K. Families reserved for future versions (not implemented)

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
