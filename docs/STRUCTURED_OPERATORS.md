# Sparse and Efficient Structures for Unitary and Orthogonal Matrices

**Reference and implementation taxonomy — 21 November 2026**

This note classifies finite matrices and operator representations that use **less storage than a conventional dense $O(N^2)$ array** and/or permit **matrix–vector multiplication faster than the conventional dense $O(N^2)$ algorithm**. The primary setting is a square $U\in\mathbb C^{N\times N}$ with $U^*U=I$, or its real orthogonal specialization. Some structural formats also apply to nonunitary matrices; those cases are explicitly distinguished.

It complements [MATRIX_FREE.md](MATRIX_FREE.md) (the actual MATLAB function-handle API), [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md) (existing implementations and limitations), and [MATRIX_ENCYCLOPEDIA.md](MATRIX_ENCYCLOPEDIA.md) (application-derived families). **Being a mathematically structured matrix does not imply that the current repository implements a fast action.**


**Application-reference policy:** Every structural classification below includes linked primary application sources. Citations distinguish (a) a proven/implemented **fast operator structure**, (b) an **application of an exactly unitary** member, and (c) an application of the **broader nonunitary** structural class. No paper cited here is asserted to implement an Anymatrix generator unless explicitly stated.

## 1. What counts as efficient?

Use **$N$ for the matrix order** in every $N\times N$ unitary/orthogonal representation and its complexity bounds; physical counts and generator input parameters retain separate symbols (such as $n$, $q$, or $L$). Complexity figures below refer to one arbitrary dense vector after any stated setup; multiple right-hand sides usually scale with their count. Let $s=\mathop{\mathrm{nnz}}(U)$, $b$ be a fixed bandwidth/block size, $r$ a low-rank parameter, $g$ a number of local gates/cores, and $L$ a number of tensor factors.

Four separate attributes must be recorded:

1. **Explicit sparsity:** the matrix itself has substantially fewer than the $N^2$ entries of a general dense matrix; for example, a matrix with bandwidth $b$ has $s=O(Nb)$ nonzeros, which is linear in $N$ when $b$ is fixed.
2. **Compact exact representation:** a possibly dense matrix has a short exact description, such as $O(N)$ Fourier phases or $O(N)$ adjacent unitary cores, rather than requiring a full dense array.
3. **Fast action:** a specified algorithm computes $Ux$ and $U^*x$ with a documented cost such as $O(N)$ or $O(N\log N)$ instead of conventional dense $O(N^2)$ multiplication, taking into account factor lengths, rank growth, preprocessing, and tolerance.
4. **Approximate data sparsity:** a compressed approximation $\widetilde U$ meets a stated error tolerance. Unless constructed from unitary factors, it may **not** satisfy $\widetilde U^*\widetilde U=I$ exactly.

These concepts are independent. An unreduced unitary upper-Hessenberg matrix is generally **dense above the first subdiagonal** but has a short core factorization. By contrast, the manifold $U(N)$ has real dimension $N^2$, so a regular exact parametrization of unrestricted $N\times N$ unitaries cannot use only $O(N)$ or $O(N\log N)$ independent real parameters. There is no universal subquadratic-parameter exact description for arbitrary unitary matrices.

A function handle wrapping a stored dense array is **not** a compressed representation. Likewise, sparsity of a generator $H$ does not imply sparsity of $\exp(iH)$.

## 2. Summary of important structures

Entries marked **conditional** require hypotheses such as bounded ranks, particular transforms, or specialized algorithms; they are not blanket complexity guarantees. The storage column describes **operator data**, excluding the $N$-entry input/output vectors and any separately reported setup workspace.

| Class / representation | Operator storage | Action per vector | Exact unitary construction / qualification |
| --- | --- | --- | --- |
| Diagonal unit phases | $O(N)$ | $O(N)$ | Diagonal entries on the unit circle |
| Signed/phase-weighted permutation (monomial) | $O(N)$ | $O(N)$ | One unimodular nonzero per row and column |
| Banded unitary, fixed bandwidth $b$ | $O(Nb)$ | $O(Nb)$ | Must satisfy global orthogonality constraints; arbitrary banded matrices do not |
| Finite CMV / fixed-block CMV | $O(Nb)$, scalar $b=O(1)$ | $O(Nb)$ | Alternating unitary local blocks and valid terminal closure |
| Block diagonal / bounded-size direct sums | $O(Nb)$ | $O(Nb)$ | Unitary blocks; $b$ denotes maximum dense block size |
| Sparse quantum walks / scattering networks | $O(s)$ | $O(s)$ | Closed lossless network, phases, permutations and local unitary coins |
| Twisted/snake core orderings | $O(N)$ for one core per adjacent index | $O(N)$ | Same finite OPUC core data with alternative insertion order; explicit sparsity differs |
| Literal X-pattern paired unitary blocks | $O(N)$ | $O(N)$ | Diagonal + antidiagonal supports; independent mirrored $2\times2$ unitaries |
| X-shaped two-branch **product representation** | Depends on rank-structured core count | Depends on core count | Van Barel pull-through form; **not** the literal sparse X nonzero pattern |
| Chain of $g$ two-coordinate core transforms | $O(g)$ | $O(g)$ | Every core is unitary; $g=O(N)$ for a Hessenberg chain |
| $k$ Householder reflectors / rank-one phase reflections | $O(kN)$ | $O(kN)$ | Exact factor products; a full set of $k$ proportional to $N$ typically loses the advantage |
| Unitary rank-$r$ update of fast base $U_0$ | $O(Nr+r^2)+S(U_0)$ | $O(Nr+r^2)+T(U_0)$ | Constrained update, not an arbitrary additive low-rank perturbation |
| FFT/circulant unitary | $O(N)$ phase/eigenvalue data | $O(N\log N)$ | Fourier-diagonal phases unimodular |
| DFT / DCT / DST / Walsh-Hadamard | $O(1)$–$O(N)$ metadata | $O(N\log N)$ (Walsh: $O(N\log N)$) | Correct normalization, ordering and finite boundary model |
| Compactly supported orthonormal wavelets | $O(1)$ taps plus schedule for fixed filters | $O(N)$ for standard multilevel FWT | Exactly matched analysis, synthesis and boundary rules |
| Kronecker $U_1\otimes\cdots\otimes U_L$ | $O(\sum_j d_j^2)$ | $O(N\sum_j d_j)$ via mode products | $N=\prod_jd_j$, each $U_j$ unitary; can be faster when factors have structure |
| Product of sparse stages / circuit | Sum of stage descriptions | Sum of stage-application costs | Exact if every stage is unitary; gate count matters |
| Unitarily factored butterfly / multilevel network | Often $O(N\log N)$ | Often $O(N\log N)$ | Exact when each sparse stage is unitary; not every butterfly approximation is |
| Bounded-order quasiseparable / semiseparable | $O(Nr^2)$ conservative generator bound | $O(Nr^2)$, representation-dependent | Includes core-represented unitary Hessenberg; low ranks alone do not imply unitarity |
| HODLR, bounded off-diagonal rank $r$ | $O(Nr\log N)$ | $O(Nr\log N)$ | Usually an **approximation**; ordinary compression need not preserve unitarity |
| HSS / nested-basis hierarchical | $O(Nr+(N/b)r^2)$ typical | Similar order in a standard representation | Rank/leaf-size dependent; usually approximate unless specially constrained |
| Low displacement rank (Toeplitz-/Cauchy-like) | Often $O(Nr)$ generators | **Conditional** fast algorithms | Only particular members are unitary; fast multiplication depends on displacement operators |
| Complementary-low-rank / interpolative butterfly | Often $O(Nr\log N)$ | Often $O(Nr\log N)$ | Usually approximate, not automatically norm-preserving |
| Tensor train / matrix-product operator | $O(Ld^2\chi^2)$ for uniform local size $d$, bond $\chi$ | Depends on vector format and bond dimensions | Compact operator description does not reduce the $N$ entries of an arbitrary dense input vector |
| Generic unstructured dense unitary/orthogonal | $O(N^2)$ conventional storage | $O(N^2)$ standard action | Haar, generic polar/Procrustes; no general subquadratic exact representation (§5A) |
| Dense complete transform without verified fast action | Typically $O(N^2)$ if explicitly stored; defining parameters may be much fewer | $O(N^2)$ direct; fast complete action unverified | Full DPSS and arbitrary graph Fourier; **not** an impossibility claim (§5A) |
| Rational action of a sparse Hermitian matrix (Cayley) | Sparse factors / factorizations | Sparse solve cost, sometimes $O(N)$ | Exact mathematical unitary; cost depends on fill and factorization |
| Exponential of sparse skew-Hermitian generator | $O(\mathop{\mathrm{nnz}}(H))$ generator | **Accuracy- and time-dependent** action | $\exp(H)$ is exactly unitary mathematically but generally dense; applied approximately |

For a sequence of families, **subquadratic** requires the ranks, bandwidths, gate counts, tensor dimensions or other parameters to grow slowly enough with $N$. Big-O claims made for fixed parameters must not be applied when those parameters scale with $N$.

## 3. Literal sparse unitary matrices

### 3.1 Diagonal, monomial and block-diagonal

If $D=\mathop{\mathrm{diag}}(e^{i\theta_j})$, then $D$ is unitary and both its storage and multiplication are linear. A **monomial unitary** $U=DP$, with $P$ a permutation, is equally economical. Weighted cycles have explicitly computable eigenphases, making them particularly useful spectral benchmarks.

A direct sum $U=\bigoplus_j U_j$ is unitary if and only if each block is unitary. Bounded block sizes give linear storage/action, while a single growing dense block may dominate the cost. Relevant repository families: [diagonal_phase](../diagonal_phase.m), [weighted_shift](../weighted_shift.m), [cycle_shift](../cycle_shift.m), [orthogonal_blocks](../orthogonal_blocks.m), and [jaynes_cummings](../jaynes_cummings.m).


**Applications and primary references.** Phase-weighted permutations model finite measure-preserving/Koopman dynamics, transport along discrete cycles, and the shift stages in discrete quantum walks; [Davis, *Circulant Matrices* (1979)](https://books.google.com/books?q=Philip+J.+Davis+Circulant+Matrices+1979) explains cyclic/circulant models, while the database's [finite Koopman-torus generator](../finite_koopman_torus.m) gives an explicitly finite permutation application. Complete two-state excitation blocks occur in the Jaynes–Cummings cavity-QED model: [Jaynes and Cummings, *Proc. IEEE* (1963)](https://doi.org/10.1109/PROC.1963.1664). This block-diagonal example is an exactly unitary *finite direct sum of complete excitation-sector propagators*, not an arbitrary cutoff of a larger Hamiltonian.

### 3.2 Banded, CMV and block-CMV

For scalar finite CMV matrices,
```math
\Theta(\alpha)=
\begin{pmatrix}\overline\alpha&\rho\\ \rho&-\alpha\end{pmatrix},
\quad \rho=(1-|\alpha|^2)^{1/2},\quad |\alpha|<1,\qquad
\mathcal C=\mathcal L\mathcal M,
```
where $\mathcal L,\mathcal M$ are alternating direct sums of these $2\times2$ cores, with an appropriate **unit-modulus terminal boundary parameter**. This constructs a five-diagonal finite unitary. For fixed-size matrix-valued Verblunsky coefficients, block-CMV analogues have bounded block bandwidth. They remain $O(Nb)$ in stored entries for fixed block size $b$, but the boundary and contractivity constraints are essential [1–3].

Do **not** call a general unitary upper-Hessenberg matrix sparse merely because it is represented by the same core data: it usually has many nonzeros above its first subdiagonal. See [finite_cmv](../finite_cmv.m), [opuc_matrix](../opuc_matrix.m), [block_cmv_pencil](../block_cmv_pencil.m), and [MATRIX_FREE.md](MATRIX_FREE.md).


**Applications and primary references.** (i) **OPUC and spectral quadrature:** [Cantero–Moral–Velázquez (2003)](https://doi.org/10.1016/S0024-3795(02)00457-3) construct the five-diagonal representation, and [Simon, *CMV matrices: Five years after* (2007)](https://doi.org/10.1016/j.cam.2006.10.033) surveys its operator-theoretic uses. (ii) **Discrete-time quantum walks:** [Grünbaum, Velázquez, Werner and Werner, *Matrix-valued Szegő polynomials and quantum random walks* (2010)](https://doi.org/10.1002/cpa.20312) connects CMV/matrix-valued OPUC with quantum-walk dynamics; [Anderson localization for electric quantum walks and skew-shift CMV matrices](https://doi.org/10.1007/s00220-021-04204-w) uses the unitary band-matrix correspondence. (iii) **Random matrix theory:** [Killip–Nenciu (2004)](https://doi.org/10.1155/S1073792804141597) realize circular β eigenvalue laws using sparse CMV models. The classical asymptotic or bi-infinite theorems must not be asserted as exact finite truncation spectra.

### 3.3 Graphs, local walks and sparse scattering

Coin-plus-shift quantum walks, closed directed scattering networks, and block phase permutations are sparse when each vertex/link has bounded degree and each local coin/scattering matrix is small. Unitarity comes from the **product of compatible unitary stages**, not merely from the graph's sparsity. Open absorbing boundaries can invalidate it. Implemented examples include [split_step_walk](../split_step_walk.m), [dirac_automaton](../dirac_automaton.m), [chalker_coddington](../chalker_coddington.m), [metric_quantum_graph](../metric_quantum_graph.m), and [connection_cycle](../connection_cycle.m).


**Applications and primary references.** **Markov-chain search** motivates the reflection-product unitary walks of [Szegedy, FOCS (2004)](https://doi.org/10.1109/FOCS.2004.53); **graph search** is treated in [Magniez, Nayak, Roland and Santha, *Search via Quantum Walk* (2011)](https://doi.org/10.1137/090745854). **Metric quantum-graph scattering and wave propagation** are studied by [Kottos–Smilansky, *Ann. Phys.* (1999)](https://doi.org/10.1006/aphy.1999.5904). **Discrete relativistic propagation** motivates the one-dimensional Dirac quantum cellular automaton of [Bisio–D'Ariano–Tosini (2015)](https://doi.org/10.1016/j.aop.2014.12.016). In graph scattering the *fixed-wavenumber bond matrix* can be unitary while its secular equation defines a separate physical spectrum.

## 4. Dense matrices with exact fast representations

### 4.1 Adjacent cores, unitary Hessenberg, CMV snakes

An upper-Hessenberg unitary can be represented as
```math
U=G_1G_2\cdots G_{N-1}D,
```
where $G_j$ acts nontrivially on adjacent coordinates and $D$ is a suitable terminal diagonal phase. This representation has $O(N)$ parameters and an $O(N)$ forward/adjoint action **without constructing the usually dense Hessenberg matrix**. Reordering the same local cores gives CMV, inverse-Hessenberg and snake/twisted representations, with potentially different explicit nonzero patterns [1–4].

A chain of $O(N)$ cores is genuinely linear-time, whereas an arbitrary full Reck/Clements-style optical mesh of $O(N^2)$ cores is not. See [etna_hessenberg](../etna_hessenberg.m), [js_core_product](../js_core_product.m), and [um_core_apply](../um_core_apply.m). Distinguish **storage and matvec cost** from the complexity of a particular structured QR/QZ eigenvalue algorithm.


**Applications and primary references.** **Unitary eigenvalue algorithms:** [Gragg, *QR algorithm for unitary Hessenberg matrices* (1986)](https://doi.org/10.1016/0377-0427(86)90169-X) and [Aurentz–Mach–Vandebril–Watkins, *Fast and stable unitary QR* (2015)](https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341). **Harmonic analysis and signal frequency estimation:** [Gragg, *positive definite Toeplitz, isometric Arnoldi and quadrature* (1993)](https://doi.org/10.1016/0377-0427(93)90294-L) relates recurrence/core structure to spectral measures; [Ammar–Gragg–Reichel, *Determination of Pisarenko Frequency Estimates as Eigenvalues of an Orthogonal Matrix*, SPIE **826**, 143–145 (1987)](https://stage.eurekamag.com/research/109/661/109661106.php) gives the application directly. **Gauss–Szegő quadrature:** [Bultheel–Cantero–Cruz-Barroso (2015)](https://doi.org/10.1016/j.cam.2014.11.002) surveys unit-circle quadrature realizations. The quadrature vector measure must not be mistaken for the normalized trace measure.

### 4.1a Twisted / snake-shaped adjacent-core factorizations

The important generalization of a *single* forward core chain is an **arbitrary admissible insertion order** of successive adjacent unitary cores. Let $G_k$ act only on coordinates $k,k+1$. An ordering that inserts the next $G_k$ at one end or the other of the existing partial product yields a *snake-shaped factorization*. The conventional Hessenberg, inverse-Hessenberg and even/odd CMV products are special orders; other choices produce **twisted** representations of the same multiplication operator in different orthogonal Laurent polynomial bases.

For a finite order-$N$ unitary closure the factors consist of $N-1$ embedded $2\times2$ unitary cores and the correct terminal unimodular factor. All such factorizations require $O(N)$ scalar parameters and permit **$O(N)$ matvecs**, provided one keeps the product/factor schedule rather than expanding the generally less-sparse matrix. Their **explicit entry patterns vary**: the fact that an equivalent CMV representation is five-diagonal does *not* make all twisted forms explicitly five-diagonal.

The literature proves important spectral-equivalence results: for compatible coefficient data and the same boundary closure, the finite unitary snake realizations have the same eigenvalues and the same spectral measure relative to the distinguished starting vector; do not confuse this measure with the equal-weight trace measure. Core insertion sequences also occur naturally in unitary Arnoldi, Laurent-polynomial orthogonalization and structure-preserving eigenvalue methods.

**Applications / primary sources.**
- **Isometric Lanczos and circle spectral quadrature:** [D. S. Watkins, *Some Perspectives on the Eigenvalue Problem*, *SIAM Review* **35** (1993), 430–471](https://doi.org/10.1137/1035090) develops the unit-circle Lanczos/Laurent-polynomial perspective, and [Bultheel–Cantero–Cruz-Barroso, *Matrix methods for quadrature formulas on the unit circle: a survey* (2015)](https://doi.org/10.1016/j.cam.2014.11.002) treats the quadrature framework.
- **Snake-shaped factorizations and OPUC:** [Cruz-Barroso–Delvaux, *Orthogonal Laurent polynomials on the unit circle and snake-shaped matrix factorizations*, *J. Approx. Theory* (2009)](https://doi.org/10.1016/j.jat.2008.08.004) gives the basis-order/factor-order construction and the finite spectral relationships; see also [Cantero–Moral–Velázquez, *Minimal representations of unitary operators and orthogonal polynomials on the unit circle* (2005)](https://doi.org/10.1016/j.laa.2005.04.025).
- **Numerical eigensolvers:** [Aurentz–Mach–Robol–Vandebril–Watkins, *Core-Chasing Algorithms for the Eigenvalue Problem* (2018)](https://doi.org/10.1137/1.9781611975345) provides factor-chasing algorithms based on adjacent cores.

**Repository status.** The [MATHEMATICS.md](MATHEMATICS.md) discussion and the [\`um_core_order\` kernel](../um_core_order.m) encode admissible snake orderings, while [\`um_core_apply\`](../um_core_apply.m) supplies forward/adjoint products. More general *block* twists require separate core scheduling and block-dimension validation; arbitrary layouts should not be labeled implemented without tests.

### 4.1b X-shaped sparse unitary matrices versus X-shaped **factorizations**

There are **two distinct meanings** of “X-shaped” that must not be conflated:

**(A) Literal X-shaped nonzero pattern.** Define $J$ as the reversal permutation, $Je_j=e_{N+1-j}$. An X-pattern matrix has
```math
U_{ij}=0 \quad\text{unless}\quad j=i\ \text{or}\ j=N+1-i.
```
In a basis ordered by mirrored pairs $(e_j,e_{N+1-j})$, it is a direct sum of $2\times2$ blocks (plus a possible central $1\times1$ block for odd $N$). It is **unitary if and only if** every mirrored-pair block is unitary and the central entry (if present) has modulus one. It has at most $2N$ nonzeros, and therefore **$O(N)$ storage and multiply**. The entire matrix class is exactly analyzable as pairwise $2\times2$ scattering/rotation blocks. Applications: elementary beam-splitter pairings, paired reflections, mirror couplers and exactly soluble sparse orthogonal spectral benchmarks. This is a different sparse pattern from finite five-diagonal CMV.

**(B) X-shaped *Givens-product factorization*.** In [Van Barel, *Unitary rank structured matrices*, *J. Comput. Appl. Math.* (2008)](https://doi.org/10.1016/j.cam.2007.03.020), the **X-shaped factorization** is a diagram of two crossing Givens-product branches produced during the *block pull-through* process (in particular Figure 5.3), not a claim that the **entries of the resulting matrix** occur only on the diagonal and antidiagonal. It describes efficient representations and reordering of **unitary rank-structured matrices**, and the factor count/rank profile determines the action complexity. This should be cross-indexed with **quasiseparable unitary** representations in §5.1.

**Attribution note.** [Watkins (1993)](https://doi.org/10.1137/1035090) is the foundational cited review for unit-circle Lanczos and sparse Laurent representations. The **specific phrase “X-shaped factorization” and the crossing-branches figure** have been verified in *Van Barel (2008)*, not in the accessible abstract of Watkins (1993). Until the original Watkins figure is checked, do not state that Watkins introduced this precise diagram or that its product factor is a literally X-sparse matrix. Preserve both constructions as distinct database structure tags.

### 4.1c Strict core-factored unitary Hessenberg pencils versus general RQR

`unitary_hessenberg_pencil` is a genuinely O(N)-parameter
pair A,B of independently unitary upper-Hessenberg matrices.
A=diag(d)H_a and B=H_b each use N-1 embedded 2-by-2 unitary
cores and a terminal phase. The **quotient U=B'*A** therefore
is mathematically exactly unitary and has O(N) matvec/adjoint
complexity via `um_hessenberg_pencil_apply`. A and B generally
have dense upper triangles; the core representation, not explicit
entry sparsity, supplies the efficiency.

The proper Hessenberg pencil pole at index j is the projective
subdiagonal ratio `xi_j=A(j+1,j)/B(j+1,j)`. Independent
core-sine choices realize any finite nonzero pole magnitude
(with finite precision range restrictions), while a left row
phase supplies its argument. Exactly zero/infinite poles use
one vanishing core sine and force a REDUCIBLE factor; one
cannot assert both subdiagonals nonzero and a zero or infinite
pole simultaneously.

The auxiliary `um_hessenberg_unitary_generalized` deforms
the numerator on the diagonal, preserving poles and upper
Hessenberg form while losing numerator/quotient unitarity.
It is filed as a **nonunitary generalized pencil** in
`metadata/eigenvalue_problems.json`, not counted among the
strict unitary generator classes. This is distinct from the
full rational pole-swapping algorithms in
[Camps et al. (2020)](https://doi.org/10.1553/etna_vol52s480)
and their [2025 monograph](https://doi.org/10.1137/1.9781611978377).

### 4.2 Reflections and constrained low-rank updates

A unitary rank-one phase update is
```math
R(v,\phi)=I+(e^{i\phi}-1)vv^*,\qquad v^*v=1.
```
The Householder reflector is $\phi=\pi$. Its action is $x+(e^{i\phi}-1)v(v^*x)$, requiring $O(N)$ work. More generally, with $V^*V=I_r$, $K\in U(r)$, and a fast unitary $U_0$,
```math
U=U_0\bigl[I+V(K-I_r)V^*\bigr]
```
is **exactly unitary**, stores $V,K$ plus $U_0$, and costs $T(U_0)+O(Nr+r^2)$ per vector. The Hermitian projector $VV^*$ makes this an invariant-subspace unitary update. An arbitrary $U_0+XY^*$ is **not** automatically unitary.

Relevant examples: [grover_operator](../grover_operator.m), [unitary_rank_one_update](../unitary_rank_one_update.m), [entangled_reflection](../entangled_reflection.m), and geometrically motivated low-rank transport maps. [References 4, 5.]


**Applications and primary references.** **Quantum search and amplitude amplification:** [Grover, STOC (1996)](https://doi.org/10.1145/237814.237866) uses products of low-rank reflections about marked and uniform states. **Orthogonal transformations in numerical linear algebra:** [Householder, *J. ACM* (1958)](https://doi.org/10.1145/320941.320947) gives the classical reflector application to triangularization. **Geometric optimization:** [Edelman–Arias–Smith, *Geometry of algorithms with orthogonality constraints* (1998)](https://doi.org/10.1137/S0895479895290954) develops methods on orthogonality-constrained manifolds, where low-dimensional rotations/retractions are useful; this is an *application of orthogonality*, not a claim that every manifold update has the particular rank-$r$ formula above. **Rootfinding distinction:** [Aurentz et al., *Fast and backward stable computation of roots of polynomials* (2015)](https://doi.org/10.1137/140983434) exploits a **unitary-plus-rank-one companion matrix**, whose *full matrix is generally not unitary*.

### 4.3 FFT-diagonal, circulant, Fourier and chirp factorizations

Let $F$ be the normalized DFT and $D=\mathop{\mathrm{diag}}(e^{i\phi_j})$. Then
```math
U=F^*DF
```
is an exactly unitary **circulant** matrix, generally dense, with $O(N)$ stored phases and $O(N\log N)$ matvecs. Multiplying several diagonal and Fourier stages (e.g., $D_1F^*D_2F$) retains unitary structure; complexity is proportional to the number of FFT stages. Not every Toeplitz matrix is unitary, and the fact that a Toeplitz matrix is FFT-multipliable is **not** a proof of unitarity.

Related exactly normalized Fourier, DCT/DST, Hartley, Zak, Weyl–Heisenberg and Walsh constructions exploit FFT, permutations, diagonal chirps, or fast butterfly recursions. Conventional DCT/DST type and endpoints matter. Examples: [phase_circulant](../phase_circulant.m), [kicked_harper](../kicked_harper.m), [discrete_fourier](../discrete_fourier.m), [orthogonal_dct](../orthogonal_dct.m), [discrete_sine](../discrete_sine.m), [walsh_hadamard](../walsh_hadamard.m), and [affine_fourier](../affine_fourier.m). [References 6, 7.]


**Applications and primary references.** **Periodic convolution and Fourier spectral discretization:** [Davis, *Circulant Matrices* (1979)](https://books.google.com/books?q=Philip+J.+Davis+Circulant+Matrices+1979) and [Cooley–Tukey, FFT (1965)](https://doi.org/10.1090/S0025-5718-1965-0178586-1) underpin fast Fourier-diagonal actions; only unimodular Fourier multipliers define a unitary circulant. **Quantum-chaotic Floquet models:** [Izrailev, *Simple models of quantum chaos* (1990)](https://doi.org/10.1016/0370-1573(90)90067-C) discusses kicked-rotor models, and [Artuso–Casati–Shepelyansky, kicked Harper (1992)](https://doi.org/10.1103/PhysRevLett.68.3826) motivates diagonal-phase/FFT compositions. **Fast trigonometric and algebraic signal processing:** [Püschel–Moura (2008)](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf) supplies transform-algebra context.

### 4.4 Wavelets and sparse-stage factorizations

Orthonormal compactly supported wavelet transforms have a fast multilevel filter-bank algorithm. For a fixed filter length, the ordinary decimated discrete wavelet transform typically requires $O(N)$ operations even though its complete matrix is not usually banded. A finite transform is orthogonal only for **compatible normalization and boundary extension**.

More generally, if
```math
U=S_mS_{m-1}\cdots S_1,\qquad S_j^*S_j=I,
```
then $U$ is exactly unitary. Storage and action cost are respectively the **sum** of stored stage data and work across stages. This covers local gate circuits, Fourier butterflies, optical networks, lossless filter banks and many Floquet models. It does not imply fast application without a stage-count bound [6–8]. See [daubechies_wavelet](../daubechies_wavelet.m), [optical_mesh](../optical_mesh.m), [dual_unitary_circuit](../dual_unitary_circuit.m), and [orthogonal_butterfly](../orthogonal_butterfly.m).


**Applications and primary references.** **Multiresolution analysis, compression and denoising:** [Mallat (1989)](https://doi.org/10.1109/34.192463) and [Daubechies, *Ten Lectures on Wavelets* (1992)](https://doi.org/10.1137/1.9781611970104). **Lossless multirate filter banks:** [Soman–Vaidyanathan, *On orthonormal wavelets and paraunitary filter banks* (1993)](https://doi.org/10.1109/78.205722) establishes the link between orthonormal wavelet bases and paraunitary polyphase networks. **Programmable linear optics:** [Reck et al. (1994)](https://doi.org/10.1103/PhysRevLett.73.58) and [Clements et al. (2016)](https://doi.org/10.1364/OPTICA.3.001460) use beam-splitter/phase meshes to realize arbitrary lossless finite unitaries; the *universal meshes generally use $O(N^2)$ cores*, so their universal parametrization is not a subquadratic method. **Local many-body circuits:** [Bertini–Kos–Prosen (2019)](https://doi.org/10.1103/PhysRevLett.123.210601) gives dual-unitary dynamics, a circuit structure distinct from an arbitrary butterfly approximation.

### 4.5 Tensor and Kronecker products

For $U=\bigotimes_{j=1}^LU_j$, each $U_j\in U(d_j)$, $N=\prod_jd_j$. Only $\sum_jd_j^2$ scalar entries need be stored. Reshape a length-$N$ vector into a tensor and apply the small matrices to each mode. With dense $d_j\times d_j$ factors the action takes $O(N\sum_jd_j)$, not $O(N^2)$ when the local dimensions remain small. An **arbitrary dense input still consists of $N$ entries that must be stored and read**. The Kronecker matrix need never be materialized.

See [haar_tensor_power](../haar_tensor_power.m), [haar_block_tensor_power](../haar_block_tensor_power.m), [tensor_permutation](../tensor_permutation.m), [tensor_fem_wave](../tensor_fem_wave.m), and [um_tensor_apply](../um_tensor_apply.m). Extra fast local factors can reduce the modewise cost.


**Applications and primary references.** **Quantum information:** [Nielsen–Chuang, *Quantum Computation and Quantum Information* (2010)](https://doi.org/10.1017/CBO9780511976667) treats local gates acting on tensor-product Hilbert spaces. **Quantum channels:** [Watrous, *The Theory of Quantum Information* (2018)](https://cs.uwaterloo.ca/~watrous/TQI/) describes superoperators, including the exact unitary-conjugation representation $\overline V\otimes V$, whose order is $d^2$. **Separable multidimensional PDE discretizations:** tensor-product sine transforms diagonalize particular constant-coefficient or separable wave/elliptic operators, and the repository's [tensor_fem_wave](../tensor_fem_wave.m) applies this structure to a mass-weighted finite-element wave propagator. A *general nonseparable PDE coefficient field* cannot be assumed to retain tensor diagonalization.

### 4.6 Cayley transforms and rational implicit operators

For Hermitian $H$, the Cayley transform
```math
U=(H-iI)(H+iI)^{-1}
```
is exactly unitary and is generally dense even when $H$ is sparse. Applying it through an already factored shifted sparse system can be inexpensive. For Hermitian tridiagonal $H$, banded solves are $O(N)$ after $O(N)$ setup; higher-dimensional sparse problems may have substantial fill. Never claim $O(N)$ for arbitrary sparse matrices or ill-conditioned solves. [References 4, 9.] See [cayley_jacobi](../cayley_jacobi.m), [nonuniform_fem_wave](../nonuniform_fem_wave.m), and [um_cayley_apply](../um_cayley_apply.m).

Likewise, $\exp(iH)$ is unitary for Hermitian $H$, but applying it with Krylov, Taylor, Chebyshev or other exponential-action methods is usually **approximate** and iteration/degree dependent. Store the sparse generator; do not infer that the exponential is sparse or that its action has a uniform linear-time bound [10]. The existing [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md) documents these distinctions.


**Applications and primary references.** **Symmetric-tridiagonal eigenvalue algorithms:** [Aurentz–Mach–Vandebril–Watkins, *Computing eigenvalues ... via a Cayley transformation* (2017)](https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46) exploits the circle-valued image of a Hermitian problem. **Geometric and energy-preserving time integration:** [Hairer–Lubich–Wanner, *Geometric Numerical Integration*](https://doi.org/10.1007/3-540-30666-8) describes structure-preserving numerical integration, for which Cayley/implicit midpoint steps are important instances. **Large-scale exponential propagation:** [Al-Mohy–Higham, *Computing the Action of the Matrix Exponential* (2011)](https://doi.org/10.1137/100788860) derives accuracy-aware actions without assembling the full exponential. A *time-stepping approximation* to the exponential should be identified separately from the exact unitary exponential.

## 5. Low-rank and hierarchical matrix structure

These formats are important to numerical linear algebra even though a generic member is **not** unitary. They should enter this database as a fast unitary representation only after an exact structural construction or a validated unitary-preserving approximation is supplied.

### 5.1 Quasiseparable and semiseparable matrices

A matrix is $r$-**quasiseparable** if off-diagonal submatrices crossing any index split have rank at most $r$ (upper and lower orders may differ). Generators encode low-rank factors and their transitions. With fixed generator order, matrix–vector products are linear or near-linear in $N$, with constants/powers depending on the chosen generator format. The class includes unreduced unitary Hessenberg matrices and other low-rank off-diagonal patterns. A quasiseparable matrix need not be unitary: **rank geometry and orthogonality are separate constraints** [4, 11, 12].

Semiseparable terminology varies across sources; record the precise rank definition (strict upper/lower off-diagonal blocks, inclusion of diagonal, generator convention), not only a broad label.


**Applications and primary references.** **Fast unitary eigenvalue and eigenphase counting methods:** [Aurentz–Mach–Robol–Vandebril–Watkins, *Core-Chasing Algorithms* (2018)](https://doi.org/10.1137/1.9781611975345) analyzes the structured algorithms, while [Eidelman–Haimovici, *Bisection for unitary Hessenberg matrices via quasiseparable structure* (2023)](https://cris.tau.ac.il/en/publications/the-bisection-eigenvalue-method-for-unitary-hessenberg-matrices-v/) is directly concerned with eigenvalue counting. **Structured linear systems and recurrence-based computations:** [Vandebril–Van Barel–Mastronardi, *Matrix Computations and Semiseparable Matrices*, vols. 1–2](https://www.press.jhu.edu/books/title/9714/matrix-computations-and-semiseparable-matrices) provides the broader generator framework. The second application category also includes nonunitary matrices; the unitary subclass needs explicit constraints.

### 5.2 HODLR, HSS and related nested-basis formats

**HODLR** recursively approximates sibling off-diagonal blocks by low-rank products. Uniform rank $r$ yields typical $O(Nr\log N)$ storage/matvec cost. **HSS** adds nested bases, potentially improving asymptotic constants and reaching linear complexity for fixed ranks and leaf sizes. An HSS representation with dense leaf blocks of size $b$ has an additional dependence on $b$ and on transfer matrices, roughly $O(Nr+(N/b)r^2+Nb)$ in a common accounting. Rank, tree geometry and accuracy must be reported [13, 14].

**$\mathcal H$, $\mathcal H^2$**, fast multipole and nested interpolation formats address related operator/kernel compression problems, with application-dependent admissibility and ranks. These representations **do not automatically enforce $U^*U=I$**. Postcompressing a dense unitary may produce a matrix that is only approximately unitary; applying an exact polar correction need not preserve low ranks or fast action.


**Applications and primary references.** **Discretized elliptic PDEs and integral equations:** [Xia–Chandrasekaran–Gu–Li, *Fast algorithms for hierarchically semiseparable matrices* (2010)](https://doi.org/10.1002/nla.691) develops HSS methods for large rank-structured problems. **Compressed direct solvers and kernel matrices:** [Xia–Chandrasekaran–Gu–Li, *Superfast Multifrontal Method for Large Structured Linear Systems of Equations* (2010)](https://doi.org/10.1137/09074543X) applies structured low-rank compression within a sparse direct solver, and [the HSS reference software](https://github.com/fastsolvers/HSS) provides numerical examples. These works concern hierarchical matrix operations and solvers, **not** a general construction of exactly unitary HSS matrices. For the database, unitary scattering/propagation applications would require a separately verified unitary parametrization or explicit defect/error bounds.

### 5.3 Toeplitz, Hankel, Vandermonde, Cauchy and displacement rank

Low displacement rank means that an operator such as
```math
\nabla_{A,B}(M)=AM-MB
```
(or the Stein analogue) has low rank, so that the structured matrix can be described by generator matrices plus the displacement pair. Toeplitz-like, Cauchy-like and Vandermonde-like classes are key examples. Suitable displacement operators support fast or superfast matvec algorithms, but their complexity and conditioning depend on the operators, nodes, representation and preprocessing [15, 16].

**An arbitrary Toeplitz, Hankel, Cauchy, or Vandermonde matrix is not unitary.** In a unitary database this heading is a *candidate structural representation*, not an independent guarantee of unitary admissibility. Fourier, orthogonal-polynomial and specially normalized transform matrices can exhibit overlapping displacement and fast-transform structure.


**Applications and primary references.** **Time-series covariance, prediction and harmonic estimation:** [Gragg, *Toeplitz matrices, isometric Arnoldi and Gaussian quadrature* (1993)](https://doi.org/10.1016/0377-0427(93)90294-L) connects Toeplitz moments to a **derived unitary Hessenberg representation**; the covariance matrix itself is not unitary. **Fast structured linear systems:** [Kailath–Kung–Morf, *Displacement ranks of matrices and linear equations* (1979)](https://www.sciencedirect.com/science/article/pii/0022247X79901240) supplies the classical framework, and [Xia et al., *Superfast structured solver for Toeplitz linear systems* (2012)](https://doi.org/10.1137/110831982) applies displacement/HSS ideas. **Polynomial coefficient conversion:** [Townsend–Webb–Olver, *Fast polynomial transforms based on Toeplitz and Hankel matrices* (2018)](https://arxiv.org/abs/1604.07486) exploits Toeplitz–Hankel structure. These conversions are **not automatically unitary matrices**, even when the source/target polynomial bases are orthogonal under related but distinct inner products.

### 5.4 Complementary low rank and butterfly compression

Fourier integral operators and some high-frequency/special-function transforms exhibit the **complementary low-rank property**: blocks indexed by fine groups in one variable and coarse groups in the other have low numerical rank. Interpolative butterfly factorizations exploit this to achieve approximately $O(Nr\log N)$ storage and action for bounded numerical ranks, with setup costs that can be much larger [17, 18].

Distinguish:
- **Exact unitary butterfly circuits:** products of sparse *unitary* stages; unitarity is guaranteed by construction.
- **Numerical butterfly factorizations:** low-rank approximations to a matrix; usually neither exact nor automatically unitary.
- **Fast polynomial-transform algorithms:** may use asymptotic expansions, recurrences, Fourier transforms, or butterfly-like compression; individual algorithms require their own stability and complexity assumptions.

The distinction matters for tests that need eigenvalues **exactly** on the unit circle, or for repeated unitary-power applications.


**Applications and primary references.** **Wave-equation propagation, high-frequency acoustics, reflection seismology and tomography:** [Candès–Demanet–Ying, *Fast Computation of Fourier Integral Operators* (2007)](https://doi.org/10.1137/060671139) and [Li–Yang–Ying, *Multiscale Butterfly Algorithm for Multidimensional Fourier Integral Operators* (2015)](https://doi.org/10.1137/140997658). **Reusable black-box operator compression:** [Li–Yang–Martin–Ho–Ying, *Butterfly Factorization* (2015)](https://doi.org/10.1137/15M1007173) constructs compact approximations for complementary-low-rank matrices. FIO discretizations may have amplitudes, damping, quadrature weights or boundary treatments that make the matrix **nonunitary**; being amenable to a butterfly algorithm does not establish a unitary spectral benchmark.

### 5.5 Orthogonal-polynomial and spherical-harmonic transforms

Properly weighted *square* orthogonal-polynomial evaluation matrices at a complete Gaussian quadrature rule are orthogonal in exact arithmetic, but **do not automatically have an $O(N\log N)$ transform**. Chebyshev transforms have FFT reductions; faster Jacobi and spherical-harmonic conversion algorithms require more specialized representations and often nontrivial precomputation. Slevinsky's fast Chebyshev–Jacobi transforms and spherical-harmonic conversion work are specific references to consult [19, 20]. Townsend and Slevinsky's [FastTransforms project](https://mikaelslevinsky.github.io/FastTransforms/index.html) is a useful implementation/literature entry point.

In the current MATLAB implementation, [gauss_jacobi_transform](../gauss_jacobi_transform.m) has an $O(N\log N)$ operator only for the exact Gauss--Chebyshev I and II cases ($\alpha=\beta=-1/2$ and $\alpha=\beta=1/2$). Other Jacobi parameters reject the operator mode. [gauss_hermite_transform](../gauss_hermite_transform.m), [gauss_laguerre_transform](../gauss_laguerre_transform.m), and [spherical_harmonic_rotation](../spherical_harmonic_rotation.m) remain explicitly dense. Specialized literature algorithms must not be attributed to these generators without matching normalization, parameters, and actual implemented forward/adjoint kernels.


**Applications and primary references.** **Spectral PDE solvers and rapid conversion between polynomial bases:** [Hale–Townsend, *Fast, Simple, and Stable Chebyshev–Legendre Transform* (2014)](https://doi.org/10.1137/130932223), [Slevinsky, *Chebyshev–Jacobi transform* (2018)](https://doi.org/10.1093/imanum/drw070), and [Townsend–Webb–Olver, *Fast polynomial transforms based on Toeplitz and Hankel matrices* (2018)](https://arxiv.org/abs/1604.07486). **Spherical harmonic analysis and geometric computation:** [Slevinsky, *Fast and backward stable transforms between spherical harmonic expansions and bivariate Fourier series* (2019)](https://doi.org/10.1016/j.acha.2017.11.001). **Software/literature collection:** [FastTransforms](https://github.com/MikaelSlevinsky/FastTransforms). Distinguish the *unitary/orthogonal, properly quadrature-weighted finite evaluation transform* from a **general change-of-basis matrix**, which may be nonunitary; moreover, a fast transform from coefficients to values does not automatically give a fast action for a differently normalized database generator.

### 5.5a Fast transforms in other orthogonal/unitary families

Beyond Jacobi/Chebyshev–Legendre, there are many application-specific fast transforms. A **mathematically orthogonal matrix** and a **fast implementation of that exact finite matrix** are separate claims. The table below uses $N$ for the **order of the complete square transform** where it exists; in spherical settings, formulas quoted in terms of a bandlimit $L$ must be converted using the relevant $N\asymp L^2$ (sphere) or $N\asymp L^3$ (rotation group), as applicable.

| Family and physical use | Literature-supported fast mechanism | Is the finite transform unitary/orthogonal? | Important source / status |
| --- | --- | --- | --- |
| **DFT, DCT-I–IV, DST-I–IV, discrete Hartley**; spectral PDEs, image/audio processing | FFT/chirp/even–odd extension, typically $O(N\log N)$ | Yes after correct **square** normalization, types and boundary conventions | [Püschel–Moura (2008)](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf); repository fast kernels implemented |
| **Walsh–Hadamard and finite dyadic group transforms**; coding, randomized projections | Recursive butterflies, $O(N\log N)$ | Yes for normalized Sylvester Walsh; not all Hadamard orderings imply the same eigenphases | [Püschel–Moura (2008)](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf); implemented |
| **Haar/Daubechies orthonormal wavelets and wavelet packets**; multiresolution and denoising | Finite-impulse-response filter bank, $O(N)$ for fixed tap length and standard decimation | Yes with matched finite boundary extension and appropriate complete basis | [Mallat (1989)](https://doi.org/10.1109/34.192463); [Daubechies (1992)](https://doi.org/10.1137/1.9781611970104); implemented variants |
| **Orthogonal Legendre/Jacobi/Gegenbauer/Chebyshev transforms**; spectral PDE and quadrature | FFT-based asymptotics, hierarchical rank structure, Toeplitz–Hankel products, often quasi-linear **under assumptions** | A square *quadrature-weighted evaluation matrix* is orthogonal; **change-of-polynomial-basis matrices generally are not** | [Hale–Townsend (2014)](https://doi.org/10.1137/130932223); [Slevinsky (2018)](https://doi.org/10.1093/imanum/drw070); [Townsend–Webb–Olver (2018)](https://arxiv.org/abs/1604.07486) |
| **Gaussian Hermite and Laguerre transforms**; quantum oscillators, tomography and radial PDEs | Specialized fast recurrence/phase-function/hierarchical methods for selected grids and parameter regimes | Square weighted Gauss-node evaluation matrix can be orthogonal; fast **continuous integral** approximation need not equal this matrix | [*A Fast Hermite Transform* (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2630232/); [Tygert, *An algorithm for the rapid evaluation of special function transforms* (2009)](https://www.sciencedirect.com/science/article/pii/S1063520309000888); do **not** tag every Gauss–Hermite/Laguerre generator fast |
| **Spherical harmonic and associated-Legendre transforms**; geoscience, climate, electromagnetics | FFT in longitude plus fast Legendre/low-rank conversion; often subquadratic in total coefficient order | Properly normalized weighted **square** sampling may be unitary; common equiangular sampling is **rectangular**, with exact analysis/synthesis only under sampling theorem | [Mohlenkamp (1999)](https://doi.org/10.1007/BF01261607); [Slevinsky (2019)](https://doi.org/10.1016/j.acha.2017.11.001); [Price–McEwen (2024)](https://www.sciencedirect.com/science/article/pii/S0021999124003589) |
| **Wigner-$D$, SO(3), and spherical-harmonic rotations**; molecular scattering, rigid rotations, 3D imaging | Block representation and fast harmonic/Wigner recurrences; specific fast algorithms with bandlimit-dependent cost | Wigner rotation **within each irreducible angular momentum block** is unitary; sampled SO(3) transforms require weights and dimension checks | [Kostelec–Rockmore, fast SO(3) transforms](https://citeseerx.ist.psu.edu/document?doi=1916f38d0e0ca8d670a6516f5342057b982cd07b&repid=rep1&type=pdf); [fast spherical-harmonic rotation (2012)](https://www.sciencedirect.com/science/article/abs/pii/S0021999111005468); code status depends on requested operation |
| **Fourier–Laguerre / wavelets on the ball**; 3D astronomy, tomography | Separable radial Laguerre and angular spherical-harmonic fast transforms | Usually a **sampling theorem and analysis/synthesis pair**; do not assume each implementation gives a square Euclidean-unitary matrix | [Leistedt–McEwen, *Exact Wavelets on the Ball* (2012)](https://arxiv.org/abs/1205.0792); [FLAG implementation](https://astro-informatics.github.io/flag-archive/index.html) |
| **DPSS / Slepian transforms**; time–frequency localization and compressed sensing | FFT plus low-rank correction for **leading DPSS subspace** and selected transforms | The **full** DPSS eigenvector matrix is orthogonal; a fast leading-subspace projection is *not automatically* a fast action on that complete matrix | [Karnik et al., *The Fast Slepian Transform* (2019)](https://doi.org/10.1016/j.acha.2017.07.005); [Zhu et al., *ROAST: Rapid Orthogonal Approximate Slepian Transform*](https://arxiv.org/search/?query=ROAST%3A+Rapid+Orthogonal+Approximate+Slepian+Transform&searchtype=all) distinguishes a fast **orthonormal approximation to a subspace** from a general full DPSS transform |
| **Krawtchouk**; imaging, discrete harmonic analysis and finite spin models | Some exact finite-length sparse-factor formulas and short-length algorithms; generic near-linear complete transform **not established here** | Normalized complete discrete-polynomial evaluation matrix is orthogonal | [Fast discrete Krawtchouk algorithms for short lengths (2025)](https://doi.org/10.3390/electronics14193958); [finite families in repository](../krawtchouk_transform.m) |
| **Hahn / dual Hahn / discrete Chebyshev / Racah**; orthogonal moments and image features | Fast evaluation of basis coefficients/moments for selected orders; **not necessarily** a fast arbitrary-vector matvec with the full dense orthogonal transform | With positive weights and a complete finite support, properly normalized polynomial evaluation matrix is orthogonal | [Fast Hahn polynomial computation (2021)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8356550/); [discrete orthogonal-polynomial toolbox (2026)](https://doi.org/10.1016/j.softx.2026.102545); repository [hahn](../hahn_transform.m) and [dual Hahn](../dual_hahn_transform.m) |
| **Graph Fourier**; graph signal processing and spectral clustering | FFT when graph Laplacian is explicitly diagonalized by Fourier/sine/tensor bases (cycles, paths, Cartesian products) | Yes for a full Hermitian graph-Laplacian eigenbasis; no generic fast eigencoordinate transform on arbitrary graphs | [Shuman–Narang–Frossard–Ortega–Vandergheynst, graph signal processing (2013)](https://doi.org/10.1109/MSP.2012.2235192); [graph_fourier](../graph_fourier.m) intentionally dense |
| **General compact exact orthogonal-product transforms**; adaptive PCA, rotations, compressed projections | Products of a controlled number of Givens, reflectors, structured butterfly factors | Yes if each retained stage is orthogonal even when the target matrix is only approximated | [*Fast approximation of orthogonal matrices and application to PCA* (2022)](https://www.sciencedirect.com/science/article/abs/pii/S0165168421004886); factor count and approximation error need recording |

**What the fast algorithm actually computes matters.** A fast basis-conversion algorithm $c\mapsto \widetilde c$, a quadrature projection $f\mapsto c$, an approximate integral transform, and the action of one *particular square unitary matrix* are different maps. To register a family as fast in this database, specify the **exact finite matrix, scaling, ordering, domain/codomain inner products, asymptotic parameter and forward/adjoint implementation**. Complexity must include any offline factorization and any rank/tolerance dependence.

**Implementation classification recommended for this repository:** label each entry **verified fast complete operator**, **fast special case only**, **fast subspace/projector only**, **fast related transform but different matrix**, or **no demonstrated fast action**. In particular, [dpss_transform](../dpss_transform.m), [graph_fourier](../graph_fourier.m), [gauss_jacobi_transform](../gauss_jacobi_transform.m), [gauss_hermite_transform](../gauss_hermite_transform.m), [gauss_laguerre_transform](../gauss_laguerre_transform.m), [hahn_transform](../hahn_transform.m) and [dual_hahn_transform](../dual_hahn_transform.m) should not acquire an unconditional $O(N\log N)$ tag merely by citing a paper about a related transform.

For a broad numerical survey, see [Olver–Slevinsky–Townsend, *Fast algorithms using orthogonal polynomials*, *Acta Numerica* (2020)](https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D) and [FastTransforms](https://github.com/MikaelSlevinsky/FastTransforms), which covers Jacobi/Legendre/Laguerre, spherical and spin transforms, Proriol polynomials on triangles, and Zernike functions on disks. Many of those are orthogonal **only after specifying the correct discrete weights and square sampling operator**.

### 5.5b Verified FFT Gauss--Chebyshev subset (implemented)

The existing `gauss_jacobi_transform` now provides a **genuine fast
matrix-free action** for the exact positive-leading, ascending-node
Gauss-Chebyshev I and II matrices:
`Q_I=DCTIII*diag((-1).^k)` at (alpha,beta)=(-1/2,-1/2),
and `Q_II=DSTI*diag((-1).^k)` at (+1/2,+1/2).
Both share verified O(N log N) FFT kernels for forward and adjoint,
with $O(N)$ operator description/workspace; no $N\times N$ matrix is stored.
The full Gauss-Jacobi and fixed-order associated-Legendre cases
still use O(N^3) Golub-Welsch setup, O(N^2) multiplication, and
**do not claim** a fast operator. This specialization is a consequence
of exact Chebyshev trigonometric identities and standard DCT/DST,
*not* a claimed implementation of the general Hale--Townsend
or Slevinsky fast Jacobi-conversion algorithms.
Source background: [Hale--Townsend (2016)](https://doi.org/10.1093/imanum/drv060).

### 5.6 Tensor trains and matrix-product operators

An operator on $L$ local $d$-dimensional spaces ($N=d^L$) can sometimes be represented as a matrix-product operator with moderate bond dimension $\chi$, requiring $O(Ld^2\chi^2)$ stored parameters. Products and sums of local operators often admit such descriptions, but bond dimension may grow quickly under composition or time evolution [21].

For **a full arbitrary vector** of length $N$, applying the compact operator cannot avoid reading/writing $O(N)$ data, and the contraction cost depends on the bonds and algorithm. Applying an MPO to an **MPS-compressed input** is a different task, with extra compression errors and possible bond growth. Unitarity of an MPO requires global constraints or an exact product of unitary gates; unrestricted local tensors are not automatically a unitary operator.


**Applications and primary references.** **One-dimensional quantum many-body simulation and DMRG:** [Schollwöck, *The density-matrix renormalization group in the age of matrix product states* (2011)](https://arxiv.org/abs/1008.3477) describes MPS/MPO representations; [Paeckel et al., *Time-evolution methods for matrix-product states* (2019)](https://arxiv.org/abs/1901.05824) compares practical evolution algorithms. **Higher-dimensional numerical tensor approximation:** [Oseledets, *Tensor-Train Decomposition* (2011)](https://doi.org/10.1137/090752286) supplies a numerical linear algebra foundation. **Quantum-circuit evolutions:** short-depth local unitary gate sequences can have compact MPO descriptions, but repeated time evolution and truncation can increase bonds and spoil exact unitarity. The physically large state-space order is $N=d^L$, even if the operator's bond dimensions are small.

## 5A. Dense unitary and orthogonal families without a verified fast representation

A **genuinely dense or presently non-fast family** deserves a positive classification in the database. An explicit dense matrix alone is **not** enough: a unitary Hessenberg matrix may be dense but have an $O(N)$ core representation; a circulant unitary has an FFT action; a sparse-generator exponential can admit an implicit action. These belong in their respective *efficient* classes, not automatically here.

Distinguish two statuses, attached to a **specific family, parameter regime and implementation**:

1. **generic_dense_unstructured.** A generic member has no distinguished exploitable short exact representation or algorithm. The unitary manifold $U(N)$ has real dimension $N^2$ and the orthogonal group $O(N)$ has real dimension $N(N-1)/2$. Thus a regular exact parametrization for a generic member of these unrestricted classes needs a number of independent real parameters proportional to $N^2$, ruling out subquadratic parameter counts. Conventional dense storage and standard matrix–vector multiplication cost $O(N^2)$. This is *not* a proof that every conceivable algorithm or particular sample needs that much arithmetic.
2. **dense_no_verified_fast_action.** There may be a recurrence, compact mathematical definition, sparse auxiliary eigenproblem, or specialized partial transform, but **no validated subquadratic action for the complete finite square matrix and its adjoint** is currently available in this repository. This is a research/implementation status, not a mathematical impossibility result.

### 5A.1 Representative families, applications and references

| Family | Classification in the generic case | Application references and exact qualification |
| --- | --- | --- |
| **Haar unitary** and unrestricted circular-unitary entry ensembles | **generic_dense_unstructured**; $O(N^2)$ data/action in the standard model | Random-matrix theory, quantum-chaotic spectra: [Mezzadri (2007)](https://arxiv.org/abs/math-ph/0609050), [Dyson (1962)](https://doi.org/10.1063/1.1703863); [haar_unitary](../haar_unitary.m). A sparse CMV matrix can reproduce CUE **eigenvalues** but is *not* a Haar matrix in its entries: [Killip–Nenciu (2004)](https://doi.org/10.1155/S1073792804141597). |
| **Haar real orthogonal** and unstructured orthogonal rotations | **generic_dense_unstructured**; a generic full Householder/Givens factorization also takes quadratically many stored parameters or arithmetic | Classical compact groups and random-orthogonal testing: [Mezzadri (2007)](https://arxiv.org/abs/math-ph/0609050); [haar_orthogonal](../haar_orthogonal.m). A short product of reflectors belongs to the *efficient* class instead. |
| **Polar factor of a general dense nonsingular matrix** | **generic_dense_unstructured** in the absence of structure in its input | Nearest unitary, reorthogonalization, Procrustes: [Higham (1986)](https://doi.org/10.1137/0907079); [polar_unitary](../polar_unitary.m). A normal/circulant/fast-structured input is a separate subfamily and may have a fast polar representation. |
| **Unstructured Procrustes rotation** | **generic_dense_unstructured** when its cross-covariance is unrestricted | Registration/shape alignment: [Arun–Huang–Blostein (1987)](https://doi.org/10.1109/TPAMI.1987.4767965). Low-dimensional or low-rank special cases should be classified separately. |
| **Arbitrary Haar-conjugated prescribed-spectrum matrix** $QDQ^*$ | **generic_dense_unstructured** when $Q$ is generic dense | Synthetic inverse-spectrum and spectral-density benchmarks; the diagonal phases alone do **not** compress arbitrary eigenvectors. For Fourier/core/short-circuit $Q$, instead classify by that fast factorization. |
| **Universal lossless optical mesh** with $O(N^2)$ independent couplers | **generic_dense_unstructured** in its full universal parametrization despite being a product of unitary cores | Programmable interferometers and linear optics: [Reck et al. (1994)](https://doi.org/10.1103/PhysRevLett.73.58), [Clements et al. (2016)](https://doi.org/10.1364/OPTICA.3.001460); [optical_mesh](../optical_mesh.m). Sparse-depth, butterfly-wired or otherwise restricted meshes can be efficient. |
| **Complete DPSS/Slepian eigenvector matrix** | **dense_no_verified_fast_action** for the repo's general full-basis action | Time–frequency concentration and multitaper methods: [Slepian (1978)](https://doi.org/10.1002/j.1538-7305.1978.tb02104.x), [Karnik et al., *The Fast Slepian Transform* (2019)](https://doi.org/10.1016/j.acha.2017.07.005); [dpss_transform](../dpss_transform.m). Fast **leading-subspace projections** do not establish a fast complete orthogonal transform. |
| **Arbitrary graph Fourier eigenvector matrix** | **dense_no_verified_fast_action** for general graph topology | Graph signal processing: [Shuman et al. (2013)](https://doi.org/10.1109/MSP.2012.2235192); [graph_fourier](../graph_fourier.m). Sparse Laplacian storage does not make its full eigenvector matrix sparse or fast. Circulant, product-graph and special graph topologies may be fast. |
| **Square Gaussian-quadrature orthogonal-polynomial transforms** where no matching fast full action is implemented | **dense_no_verified_fast_action**, *not* intrinsically unstructured | Spectral methods/quadrature: [Gautschi (2004)](https://doi.org/10.1093/oso/9780198506720.001.0001); [Olver–Slevinsky–Townsend (2020)](https://doi.org/10.1017/S0962492920000045). Nodes/weights and a three-term recurrence can provide an $O(N)$-parameter mathematical specification, yet direct full evaluation costs $O(N^2)$. Specialized Hermite, Jacobi, Laguerre and other fast transforms are discussed in §§5.5–5.5a; do not blanket-label them dense. |
| **Generic dense multichannel scattering matrix** | **generic_dense_unstructured** when channels are unconstrained and no fast network decomposition is specified | Quantum transport: [Beenakker (1997)](https://doi.org/10.1103/RevModPhys.69.731). By contrast, *closed networks of bounded-size local scattering nodes* belong to the sparse-network class (§3.3). |

### 5A.2 Criteria and database metadata

For each **representation** of a family, distinguish whether its numeric matrix is genuinely sparse, whether it has a specified compact exact factorization, whether its actual MATLAB operator provides a fast forward/adjoint action, and whether it has a documented approximate compression. An optional operator mode is **not** by itself proof of fast multiplication (a universal optical mesh uses $O(N^2)$ couplers). If a representation is efficient, tag it with the appropriate class and its precise rank/bandwidth/stage/tolerance assumptions. Only when these checks fail use **generic_dense_unstructured** (quadratic intrinsic degrees of freedom for generic samples) or **dense_no_verified_fast_action** (compact definition or specialized literature, but no verified fast **full** multiply).

Do **not** claim that a missing implementation proves a lower bound. In particular, distinctions between **full transform versus leading eigenspace projection**, **exactly unitary versus approximately orthogonal**, and **setup time versus repeated matvec cost** are mandatory.

Proposed additional fields (not yet implemented in `metadata/catalog.json`): **efficiency_class**, **dense_reason**, **subquadratic_storage_known**, **subquadratic_full_action_known**, **partial_transform_only**, **fast_action_reference**, **representation_exactness**, **checked_date**, **setup_complexity**, and **unitarity_defect**. The current catalog instead records **matrix_free_supported**, **operator_kernel**, and **operator_note** per generator; those specify API coverage and its caveats, not a universal speedup. A family may have multiple representations with different classifications.

This separate class prevents the database from overpromising matrix-free efficiency while leaving room for future discoveries of fast structure.

## 6. Exactness, unitary preservation and numerical hazards

| Situation | Safe statement | Unsafe inference |
| --- | --- | --- |
| Unitary Hessenberg via cores | Dense matrix, $O(N)$ descriptor and action | “Hessenberg means sparse matrix” |
| Sparse Hermitian $H$ | $\exp(iH)$ exactly unitary as an operator; numerical action can be efficient | “The exponential is sparse” or “always $O(N)$” |
| Numerical HODLR/HSS/butterfly approximation | Compact to tolerance under rank assumptions | “The approximation is exactly unitary” |
| Ordinary Toeplitz/Cauchy or graph adjacency | Structured; may admit fast multiplication | “Structured means unitary” |
| Generic Haar/polar matrices | Usually dense; no general exact short descriptor | “Random butterfly has Haar distribution” |
| Quantum circuit on $q$ qubits | Gate schedule compact; full vector length $N=2^q$ | “Simulation costs polynomial time in $q$ for arbitrary dense states” |
| Approximate sparse solves | Applied inverse depends on factorization, fill and residual | “Every sparse Cayley transform is linear time” |
| Square orthogonal polynomial transform | Orthogonal for matched weights/nodes | “Every such transform has a published fast algorithm” |
| $U_0+XY^*$ | Unitary only under additional algebraic constraints | “Low rank perturbation preserves unitarity” |

**Validation protocol.** For each claimed structure, verify both forward and adjoint actions against explicit small instances; verify $U^*U\approx I$ in small explicit cases; test multiple right-hand sides and complex adjoints; compare storage/operation counts at growing $N$; document setup work and any numerical truncation; and test spectral invariants only in the physical/symmetry sectors actually represented. For approximate structures additionally measure $\|\widetilde U^*\widetilde U-I\|$ and $\|U-\widetilde U\|$ (when the latter is accessible).

## 7. Suggested structured metadata for the database

A family may have **multiple valid representations** (e.g., explicit sparse CMV, core product, implicit operator, dense diagnostic matrix). Record each representation separately rather than assigning a single unconditional label:

| Metadata field | Meaning / example |
| --- | --- |
| **matrix_order** | $N$, distinct from underlying physical size or number of factors |
| **exact_unitary_model** | true/false, with mathematical conditions and boundary convention |
| **basis_and_factor_order** | adjacent-core insertion/twists, Laurent monomial order, block boundaries, transform normalization |
| **fast_algorithm_scope** | full exact operator; special parameter range; subspace/projector only; related but distinct transform |
| **representation_kind** | explicit_sparse, explicit_dense, core_chain, fft, tensor, rational_solve, hierarchical, butterfly, MPO, other |
| **structural_parameters** | bandwidth, nnz, core count, block size, ranks, bond dimensions, stage depth |
| **representation_exactness** | exact, floating-point factorized, or tolerance-dependent approximation |
| **storage_complexity** | symbolic expression **with parameter dependence** |
| **matvec_complexity / adjoint_complexity** | symbolic expression, work per RHS, assumptions |
| **setup_complexity** | preprocessing, factorization, FFT plans, compressed construction |
| **workspace_complexity** | additional scratch beyond input/output, and RHS batching |
| **unitarity_defect** | exact by construction vs measured $\|\widetilde U^*\widetilde U-I\|$ |
| **operator_interface** | implemented flag and verified action modes/dimensions; does not imply a subquadratic algorithm |
| **available_spectral_reference** | exact finite spectrum, moments, limiting law, or unknown |
| **reference_ids / provenance** | primary literature and code implementation |
| **implementation_status** | implemented, planned, or theoretical reference only |
| **efficiency_class / dense_reason** | efficient representation or dense status (see §5A); include evidence, restrictions and whether a fast **full** action exists |

For the existing MATLAB contract and current kernel coverage, see [MATRIX_FREE.md](MATRIX_FREE.md) and [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md). The registered operator flags and per-family caveats live in [metadata/catalog.json](../metadata/catalog.json). In particular, the repository uses **'transp' for the conjugate transpose**, not the nonconjugating MATLAB **A.'** operation; this needs to remain consistent across new kernels.

## 8. Priorities for this repository

**Already represented:** CMV/Hessenberg core chains; sparse shifts, diagonals, walks and scattering; FFT/circulant/fast trigonometric transforms; Walsh and wavelets; local circuit gates; Kronecker mode products; constrained reflections; and structured Cayley/exponential actions. Check the current [implementation audit](IMPLEMENTATION_AUDIT.md) rather than treating this list as an exhaustive assertion of fast support for every generator.

**High-value follow-up work:**

1. Add **per-representation storage and matvec metadata** to the catalog; distinguish actual implemented costs from literature-possible costs.
2. Add rigorous **operation-count and adjoint/orthogonality tests** for the existing core, FFT, wavelet, tensor, butterfly and Cayley operators.
3. Survey **twisted/snake core orderings, X-pattern sparse blocks and X-shaped pull-through factors**, then **bounded-rank quasiseparable unitary and block-CMV variants** with explicit unitary generators; connect to fast QR/QZ rather than merely using general HSS approximations.
4. Evaluate **structure-preserving hierarchical or butterfly approximations of unitary operators**, including norm-preserving factorizations, approximation error and repeated-power drift.
5. Extend the implemented **Gauss--Chebyshev I/II** FFT operator subset to additional Jacobi/associated-Legendre and spherical cases only after verifying the *same finite transform matrix*, including normalization, order, and adjoint.
6. Keep the **generic Haar, polar, DPSS and arbitrary graph eigenvector matrices** labeled dense unless an exact, family-specific fast representation is demonstrated.

## References

Numbers below are local to this note. Existing family-specific references remain in the [encyclopedia](MATRIX_ENCYCLOPEDIA.md) and [mathematical notes](MATHEMATICS.md).

1. M. J. Cantero, L. Moral, and L. Velázquez, “Five-diagonal matrices and zeros of orthogonal polynomials on the unit circle,” *Linear Algebra Appl.* **362** (2003), 29–56. [DOI](https://doi.org/10.1016/S0024-3795(02)00457-3).
2. B. Simon, *Orthogonal Polynomials on the Unit Circle*, Parts 1–2, AMS, 2005. [AMS](https://www.ams.org/books/coll/054.1/).
3. R. Killip and I. Nenciu, “Matrix models for circular ensembles,” *Int. Math. Res. Not.* (2004), 2665–2701. [DOI](https://doi.org/10.1155/S1073792804141597).
4. J. L. Aurentz, T. Mach, R. Vandebril, and D. S. Watkins, “Fast and stable unitary QR algorithm,” *Electron. Trans. Numer. Anal.* **44** (2015), 327–341. [ETNA](https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341); J. L. Aurentz, T. Mach, L. Robol, R. Vandebril, and D. S. Watkins, *Core-Chasing Algorithms for the Eigenvalue Problem*, SIAM, 2018. [DOI](https://doi.org/10.1137/1.9781611975345).
5. A. S. Householder, “Unitary triangularization of a nonsymmetric matrix,” *J. ACM* **5** (1958), 339–342. [DOI](https://doi.org/10.1145/320941.320947).
6. J. W. Cooley and J. W. Tukey, “An algorithm for the machine calculation of complex Fourier series,” *Math. Comp.* **19** (1965), 297–301. [DOI](https://doi.org/10.1090/S0025-5718-1965-0178586-1).
7. M. Püschel and J. M. F. Moura, “Algebraic signal processing theory: Foundation and 1-D time,” *IEEE Trans. Signal Process.* **56** (2008), 3572–3585. [Author PDF](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf).
8. I. Daubechies, *Ten Lectures on Wavelets*, SIAM, 1992. [DOI](https://doi.org/10.1137/1.9781611970104); W. R. Clements et al., “Optimal design for universal multiport interferometers,” *Optica* **3** (2016), 1460–1465. [DOI](https://doi.org/10.1364/OPTICA.3.001460).
9. J. L. Aurentz, T. Mach, R. Vandebril, and D. S. Watkins, “Computing the eigenvalues of symmetric tridiagonal matrices via a Cayley transformation,” *Electron. Trans. Numer. Anal.* **46** (2017), 447–459. [ETNA](https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46).
10. A. H. Al-Mohy and N. J. Higham, “Computing the action of the matrix exponential, with an application to exponential integrators,” *SIAM J. Sci. Comput.* **33** (2011), 488–511. [DOI](https://doi.org/10.1137/100788860).
11. R. Vandebril, M. Van Barel, and N. Mastronardi, *Matrix Computations and Semiseparable Matrices*, Johns Hopkins University Press, Vols. 1–2, 2007–2008. [Publisher](https://www.press.jhu.edu/books/title/9714/matrix-computations-and-semiseparable-matrices).
12. “Nested Product Decomposition of Quasiseparable Matrices,” *SIAM J. Matrix Anal. Appl.* (2013). [DOI](https://doi.org/10.1137/12086354X).
13. J. Xia, S. Chandrasekaran, M. Gu, and X. S. Li, “Fast algorithms for hierarchically semiseparable matrices,” *Numer. Linear Algebra Appl.* (2010). [DOI](https://doi.org/10.1002/nla.691).
14. “Mixed Precision HODLR Matrices,” *SIAM J. Sci. Comput.* (2025). [DOI](https://doi.org/10.1137/24M1683925).
15. T. Kailath, S.-Y. Kung, and M. Morf, “Displacement ranks of matrices and linear equations,” *J. Math. Anal. Appl.* **68** (1979), 395–407. [Publisher](https://www.sciencedirect.com/science/article/pii/0022247X79901240).
16. I. Gohberg and V. Olshevsky, “Complexity of multiplication with vectors for structured matrices,” [Author abstract](https://www2.math.uconn.edu/~olshevsky/abstracts/compl.html).
17. Y. Li, H. Yang, E. R. Martin, K. L. Ho, and L. Ying, “Butterfly factorization,” *Multiscale Model. Simul.* **13** (2015). [DOI](https://doi.org/10.1137/15M1007173).
18. Y. Li and H. Yang, “Interpolative butterfly factorization,” *SIAM J. Sci. Comput.* (2017). [DOI](https://doi.org/10.1137/16M1074941).
19. R. M. Slevinsky, “On the use of Hahn’s asymptotic formula and stabilized recurrence for a fast, simple and stable Chebyshev–Jacobi transform,” *IMA J. Numer. Anal.* **38** (2018), 102–124. [DOI](https://doi.org/10.1093/imanum/drw070).
20. R. M. Slevinsky, “Fast and backward stable transforms between spherical harmonic expansions and bivariate Fourier series,” *Appl. Comput. Harmon. Anal.* (2019). [DOI](https://doi.org/10.1016/j.acha.2017.11.001).
21. I. V. Oseledets, “Tensor-train decomposition,” *SIAM J. Sci. Comput.* **33** (2011), 2295–2317. [DOI](https://doi.org/10.1137/090752286).

*Scope note:* This is a taxonomy of efficient representations and literature **not** a claim that each general matrix structure defines a unitary matrix family, nor that each technique is implemented in MATLAB here. Exactness and complexity are representation- and parameter-dependent.

### Verified blockwise SO(3) real/complex spherical-harmonic actions

`spherical_harmonic_rotation` and the new
`real_spherical_harmonic_rotation` both expose a block direct-sum
`'operator'` on degrees 0,...,L-1. The latter uses an exact
Condon–Shortley complex-to-real harmonic basis transform,
so every degree block is real orthogonal rather than the
nonorthogonal entrywise real part of complex Wigner D.
For N=L^2, per-matvec work/storage is O(L^3)=O(N^(3/2)).
The important caveat is O(L^4) *dense* per-degree setup and
O(ell^2) per-degree application. This is NOT the asymptotically
fast Wigner or spherical rotation algorithm of
[Gimbutas–Greengard (2009)](https://doi.org/10.1016/j.jcp.2009.05.014)
or the later S2FFT family. The new matrix-free label refers
to the full block direct sum only, not to a sparse or fast
single-irrep Wigner factor. Issue #41 retains that research task.

### Hermite/Laguerre stability improvement without a fake fast method

`um_gauss_basis` now orients eigenvectors of irreducible positive
Jacobi recurrences using the more reliable FIRST or LAST endpoint
component and the Sturm parity rule. This prevents underflow in
the tiny first Gaussian eigenvector component from determining
an arbitrary sign. The optional `quadrature` mode in the
Hermite/Laguerre public generators independently evaluates
normalized log Christoffel weights via rescaled polynomial
recurrences (O(N^2)), with physical Hermite/Laguerre function
row scaling kept distinct from the orthogonal weighted polynomial
evaluation. The quadrature nodes and weights do NOT provide
an FFT-fast action of the full coefficient transform.
The explicit Q remains O(N^3) setup, O(N^2) storage and apply,
in contrast with genuinely fast specialized transforms
([Webb 2026](https://arxiv.org/abs/2604.02041), and [Olver et al.](https://www.cambridge.org/core/journals/acta-numerica/article/abs/fast-algorithms-using-orthogonal-polynomials/4FAD8C7C28EC20EE7583465C1A89AA3D)).


### Spin-weighted harmonics: Jacobi square blocks vs FFT sampling

For fixed spin s and azimuth m, Wigner d^(ell)_(m,-s)
reduces to a signed positive-leading Jacobi basis with
a=abs(m+s), b=abs(m-s). This yields a strict
(L-max(abs(m),abs(s)))-square Gauss–Jacobi **orthogonal**
weighted evaluation, the `spin_spherical_harmonic_transform`
default. Its Gauss–Jacobi nodes are **not** matrix eigenphases.

The full s-spin spherical transform has L²-s² coefficients
and L*(2L-1) Gauss–Legendre/equispaced-azimuth samples,
so the weighted evaluation has A'*A=I but is rectangular.
`um_spin_sampling_apply` exploits Fourier separation in
azimuth to apply the full map and adjoint in O(L³+L² log L)
rather than O(L^4) dense sampling arithmetic. Exact
discrete quadrature follows from degree at most 2L-2
Wigner-d products and the L-node Legendre rule. It is
a derived exact GL implementation, not a copy of
[Huffenberger–Wandelt (2010)](https://doi.org/10.1088/0067-0049/189/2/255)
or the MW equiangular grid. Pole local frames and
spin +/- conjugation are documented explicitly.


### Peter–Weyl SO(3) FFT isometry versus fixed Wigner rotations

The unitary `wigner_rotation(ell,...)` acts at a single
rotation on a 2ell+1-dimensional representation.
The SO(3) group Fourier **sampling** map acts on the
full coefficient space of L*(4L^2-1)/3 dimensions.
Gauss–Legendre beta nodes and two 2L-1-point
equiangular outer Euler grids provide L*(2L-1)^2
Haar-weighted rows. Two exact Fourier orthogonalities
and Legendre exactness for products of Wigner d functions
of degree at most 2L-2 prove A'*A=I. The full map
is rectangular except L=1. The separate square
fixed-(m,n) Wigner–Jacobi latitude block is an exact
orthogonal matrix with UNKNOWN eigenphases.
FFT-separated full synthesis/adjoint uses O(L^4)
latitude cache and O(L^4+L^3 log L) apply; it does
not reproduce the original McEwen torus-extension
algorithm ([DOI](https://doi.org/10.1109/LSP.2015.2490676)),
but realizes a Gauss–Legendre scheme of the type
studied by Khalid et al. ([DOI](https://doi.org/10.1109/LSP.2015.2503295)).
