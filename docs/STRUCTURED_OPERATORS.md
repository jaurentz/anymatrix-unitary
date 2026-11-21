# Sparse and Efficient Structures for Unitary and Orthogonal Matrices

**Reference and implementation taxonomy — 21 November 2026**

This note classifies finite matrices and operator representations that reduce **storage below \(\Theta(N^2)\)** and/or provide **matrix–vector multiplication faster than \(\Theta(N^2)\)**. The primary setting is a square \(U\in\mathbb C^{N\times N}\) with \(U^*U=I\), or its real orthogonal specialization. Some structural formats also apply to nonunitary matrices; those cases are explicitly distinguished.

It complements [MATRIX_FREE.md](MATRIX_FREE.md) (the actual MATLAB function-handle API), [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md) (existing implementations and limitations), and [MATRIX_ENCYCLOPEDIA.md](MATRIX_ENCYCLOPEDIA.md) (application-derived families). **Being a mathematically structured matrix does not imply that the current repository implements a fast action.**


**Application-reference policy:** Every structural classification below includes linked primary application sources. Citations distinguish (a) a proven/implemented **fast operator structure**, (b) an **application of an exactly unitary** member, and (c) an application of the **broader nonunitary** structural class. No paper cited here is asserted to implement an Anymatrix generator unless explicitly stated.

## 1. What counts as efficient?

Let \(N\) denote the **matrix order** (not the number of qubits, grid points before doubling, tensor modes, or physical sites). Complexity figures below refer to one arbitrary dense vector after any stated setup; multiple right-hand sides usually scale with their count. Let \(s=\operatorname{nnz}(U)\), \(b\) be a fixed bandwidth/block size, \(r\) a low-rank parameter, \(g\) a number of local gates/cores, and \(L\) a number of tensor factors.

Four separate attributes must be recorded:

1. **Explicit sparsity:** the matrix itself has \(s=o(N^2)\) nonzeros.
2. **Compact exact representation:** a possibly dense matrix is specified with \(o(N^2)\) parameters, e.g., Fourier phases or a chain of cores.
3. **Fast action:** a specified algorithm computes \(Ux\) and \(U^*x\) in \(o(N^2)\) operations, taking into account factor lengths, rank growth, preprocessing and tolerance.
4. **Approximate data sparsity:** a compressed approximation \(\widetilde U\) meets a stated error tolerance. Unless constructed from unitary factors, it may **not** satisfy \(\widetilde U^*\widetilde U=I\) exactly.

These concepts are independent. An unreduced unitary upper-Hessenberg matrix is generally **dense above the first subdiagonal** but has a short core factorization; a general \(N\times N\) Haar unitary, on the other hand, requires \(\Theta(N^2)\) real parameters to specify exactly. In fact the manifold \(U(N)\) has real dimension \(N^2\). There is no universal subquadratic exact description for arbitrary unitary matrices.

A function handle wrapping a stored dense array is **not** a compressed representation. Likewise, sparsity of a generator \(H\) does not imply sparsity of \(\exp(iH)\).

## 2. Summary of important structures

Entries marked **conditional** require hypotheses such as bounded ranks, particular transforms, or specialized algorithms; they are not blanket complexity guarantees. The storage column describes **operator data**, excluding the \(N\)-entry input/output vectors and any separately reported setup workspace.

| Class / representation | Operator storage | Action per vector | Exact unitary construction / qualification |
| --- | --- | --- | --- |
| Diagonal unit phases | \(O(N)\) | \(O(N)\) | Diagonal entries on the unit circle |
| Signed/phase-weighted permutation (monomial) | \(O(N)\) | \(O(N)\) | One unimodular nonzero per row and column |
| Banded unitary, fixed bandwidth \(b\) | \(O(Nb)\) | \(O(Nb)\) | Must satisfy global orthogonality constraints; arbitrary banded matrices do not |
| Finite CMV / fixed-block CMV | \(O(Nb)\), scalar \(b=O(1)\) | \(O(Nb)\) | Alternating unitary local blocks and valid terminal closure |
| Block diagonal / bounded-size direct sums | \(O(Nb)\) | \(O(Nb)\) | Unitary blocks; \(b\) denotes maximum dense block size |
| Sparse quantum walks / scattering networks | \(O(s)\) | \(O(s)\) | Closed lossless network, phases, permutations and local unitary coins |
| Chain of \(g\) two-coordinate core transforms | \(O(g)\) | \(O(g)\) | Every core is unitary; \(g=O(N)\) for a Hessenberg chain |
| \(k\) Householder reflectors / rank-one phase reflections | \(O(kN)\) | \(O(kN)\) | Exact factor products; full \(k=\Theta(N)\) typically loses the advantage |
| Unitary rank-\(r\) update of fast base \(U_0\) | \(O(Nr+r^2)+S(U_0)\) | \(O(Nr+r^2)+T(U_0)\) | Constrained update, not an arbitrary additive low-rank perturbation |
| FFT/circulant unitary | \(O(N)\) phase/eigenvalue data | \(O(N\log N)\) | Fourier-diagonal phases unimodular |
| DFT / DCT / DST / Walsh-Hadamard | \(O(1)\)–\(O(N)\) metadata | \(O(N\log N)\) (Walsh: \(O(N\log N)\)) | Correct normalization, ordering and finite boundary model |
| Compactly supported orthonormal wavelets | \(O(1)\) taps plus schedule for fixed filters | \(O(N)\) for standard multilevel FWT | Exactly matched analysis, synthesis and boundary rules |
| Kronecker \(U_1\otimes\cdots\otimes U_L\) | \(O(\sum_j d_j^2)\) | \(O(N\sum_j d_j)\) via mode products | \(N=\prod_jd_j\), each \(U_j\) unitary; can be faster when factors have structure |
| Product of sparse stages / circuit | Sum of stage descriptions | Sum of stage-application costs | Exact if every stage is unitary; gate count matters |
| Unitarly factored butterfly / multilevel network | Often \(O(N\log N)\) | Often \(O(N\log N)\) | Exact when each sparse stage is unitary; not every butterfly approximation is |
| Bounded-order quasiseparable / semiseparable | \(O(Nr^2)\) conservative generator bound | \(O(Nr^2)\), representation-dependent | Includes core-represented unitary Hessenberg; low ranks alone do not imply unitarity |
| HODLR, bounded off-diagonal rank \(r\) | \(O(Nr\log N)\) | \(O(Nr\log N)\) | Usually an **approximation**; ordinary compression need not preserve unitarity |
| HSS / nested-basis hierarchical | \(O(Nr+(N/b)r^2)\) typical | Similar order in a standard representation | Rank/leaf-size dependent; usually approximate unless specially constrained |
| Low displacement rank (Toeplitz-/Cauchy-like) | Often \(O(Nr)\) generators | **Conditional** fast algorithms | Only particular members are unitary; fast multiplication depends on displacement operators |
| Complementary-low-rank / interpolative butterfly | Often \(O(Nr\log N)\) | Often \(O(Nr\log N)\) | Usually approximate, not automatically norm-preserving |
| Tensor train / matrix-product operator | \(O(Ld^2\chi^2)\) for uniform local size \(d\), bond \(\chi\) | Depends on vector format and bond dimensions | Compact operator description does not mean an arbitrary full vector can be read in \(o(N)\) |
| Rational action of a sparse Hermitian matrix (Cayley) | Sparse factors / factorizations | Sparse solve cost, sometimes \(O(N)\) | Exact mathematical unitary; cost depends on fill and factorization |
| Exponential of sparse skew-Hermitian generator | \(O(\operatorname{nnz}(H))\) generator | **Accuracy- and time-dependent** action | \(\exp(H)\) is exactly unitary mathematically but generally dense; applied approximately |

For a sequence of families, **subquadratic** requires the ranks, bandwidths, gate counts, tensor dimensions or other parameters to grow slowly enough with \(N\). Big-O claims made for fixed parameters must not be applied when those parameters scale with \(N\).

## 3. Literal sparse unitary matrices

### 3.1 Diagonal, monomial and block-diagonal

If \(D=\operatorname{diag}(e^{i\theta_j})\), then \(D\) is unitary and both its storage and multiplication are linear. A **monomial unitary** \(U=DP\), with \(P\) a permutation, is equally economical. Weighted cycles have explicitly computable eigenphases, making them particularly useful spectral benchmarks.

A direct sum \(U=\bigoplus_j U_j\) is unitary if and only if each block is unitary. Bounded block sizes give linear storage/action, while a single growing dense block may dominate the cost. Relevant repository families: [diagonal_phase](../diagonal_phase.m), [weighted_shift](../weighted_shift.m), [cycle_shift](../cycle_shift.m), [orthogonal_blocks](../orthogonal_blocks.m), and [jaynes_cummings](../jaynes_cummings.m).


**Applications and primary references.** Phase-weighted permutations model finite measure-preserving/Koopman dynamics, transport along discrete cycles, and the shift stages in discrete quantum walks; [Davis, *Circulant Matrices* (1979)](https://books.google.com/books?q=Philip+J.+Davis+Circulant+Matrices+1979) explains cyclic/circulant models, while the database's [finite Koopman-torus generator](../finite_koopman_torus.m) gives an explicitly finite permutation application. Complete two-state excitation blocks occur in the Jaynes–Cummings cavity-QED model: [Jaynes and Cummings, *Proc. IEEE* (1963)](https://doi.org/10.1109/PROC.1963.1664). This block-diagonal example is an exactly unitary *finite direct sum of complete excitation-sector propagators*, not an arbitrary cutoff of a larger Hamiltonian.

### 3.2 Banded, CMV and block-CMV

For scalar finite CMV matrices,
\[
\Theta(\alpha)=
\begin{pmatrix}\overline\alpha&\rho\\ \rho&-\alpha\end{pmatrix},
\quad \rho=(1-|\alpha|^2)^{1/2},\quad |\alpha|<1,\qquad
\mathcal C=\mathcal L\mathcal M,
\]
where \(\mathcal L,\mathcal M\) are alternating direct sums of these \(2\times2\) cores, with an appropriate **unit-modulus terminal boundary parameter**. This constructs a five-diagonal finite unitary. For fixed-size matrix-valued Verblunsky coefficients, block-CMV analogues have bounded block bandwidth. They remain \(O(Nb)\) in stored entries for fixed block size \(b\), but the boundary and contractivity constraints are essential [1–3].

Do **not** call a general unitary upper-Hessenberg matrix sparse merely because it is represented by the same core data: it usually has many nonzeros above its first subdiagonal. See [finite_cmv](../finite_cmv.m), [opuc_matrix](../opuc_matrix.m), [block_cmv_pencil](../block_cmv_pencil.m), and [MATRIX_FREE.md](MATRIX_FREE.md).


**Applications and primary references.** (i) **OPUC and spectral quadrature:** [Cantero–Moral–Velázquez (2003)](https://doi.org/10.1016/S0024-3795(02)00457-3) construct the five-diagonal representation, and [Simon, *CMV matrices: Five years after* (2007)](https://doi.org/10.1016/j.cam.2006.10.033) surveys its operator-theoretic uses. (ii) **Discrete-time quantum walks:** [Grünbaum, Velázquez, Werner and Werner, *Matrix-valued Szegő polynomials and quantum random walks* (2010)](https://doi.org/10.1002/cpa.20312) connects CMV/matrix-valued OPUC with quantum-walk dynamics; [Anderson localization for electric quantum walks and skew-shift CMV matrices](https://doi.org/10.1007/s00220-021-04204-w) uses the unitary band-matrix correspondence. (iii) **Random matrix theory:** [Killip–Nenciu (2004)](https://doi.org/10.1155/S1073792804141597) realize circular-\(\beta\) eigenvalue laws using sparse CMV models. The classical asymptotic or bi-infinite theorems must not be asserted as exact finite truncation spectra.

### 3.3 Graphs, local walks and sparse scattering

Coin-plus-shift quantum walks, closed directed scattering networks, and block phase permutations are sparse when each vertex/link has bounded degree and each local coin/scattering matrix is small. Unitarity comes from the **product of compatible unitary stages**, not merely from the graph's sparsity. Open absorbing boundaries can invalidate it. Implemented examples include [split_step_walk](../split_step_walk.m), [dirac_automaton](../dirac_automaton.m), [chalker_coddington](../chalker_coddington.m), [metric_quantum_graph](../metric_quantum_graph.m), and [connection_cycle](../connection_cycle.m).


**Applications and primary references.** **Markov-chain search** motivates the reflection-product unitary walks of [Szegedy, FOCS (2004)](https://doi.org/10.1109/FOCS.2004.53); **graph search** is treated in [Magniez, Nayak, Roland and Santha, *Search via Quantum Walk* (2011)](https://doi.org/10.1137/090745854). **Metric quantum-graph scattering and wave propagation** are studied by [Kottos–Smilansky, *Ann. Phys.* (1999)](https://doi.org/10.1006/aphy.1999.5904). **Discrete relativistic propagation** motivates the one-dimensional Dirac quantum cellular automaton of [Bisio–D'Ariano–Tosini (2015)](https://doi.org/10.1016/j.aop.2014.12.016). In graph scattering the *fixed-wavenumber bond matrix* can be unitary while its secular equation defines a separate physical spectrum.

## 4. Dense matrices with exact fast representations

### 4.1 Adjacent cores, unitary Hessenberg, CMV snakes

An upper-Hessenberg unitary can be represented as
\[
U=G_1G_2\cdots G_{N-1}D,
\]
where \(G_j\) acts nontrivially on adjacent coordinates and \(D\) is a suitable terminal diagonal phase. This representation has \(O(N)\) parameters and an \(O(N)\) forward/adjoint action **without constructing the usually dense Hessenberg matrix**. Reordering the same local cores gives CMV, inverse-Hessenberg and snake/twisted representations, with potentially different explicit nonzero patterns [1–4].

A chain of \(O(N)\) cores is genuinely linear-time, whereas an arbitrary full Reck/Clements-style optical mesh of \(O(N^2)\) cores is not. See [etna_hessenberg](../etna_hessenberg.m), [js_core_product](../js_core_product.m), and [um_core_apply](../um_core_apply.m). Distinguish **storage and matvec cost** from the complexity of a particular structured QR/QZ eigenvalue algorithm.


**Applications and primary references.** **Unitary eigenvalue algorithms:** [Gragg, *QR algorithm for unitary Hessenberg matrices* (1986)](https://doi.org/10.1016/0377-0427(86)90169-X) and [Aurentz–Mach–Vandebril–Watkins, *Fast and stable unitary QR* (2015)](https://etna.ricam.oeaw.ac.at/volumes/2011-2020/vol44/abstract.php?pages=327-341). **Harmonic analysis and signal frequency estimation:** [Gragg, *positive definite Toeplitz, isometric Arnoldi and quadrature* (1993)](https://doi.org/10.1016/0377-0427(93)90294-L) relates recurrence/core structure to spectral measures; [Ammar–Gragg–Reichel, *Determination of Pisarenko Frequency Estimates as Eigenvalues of an Orthogonal Matrix*, SPIE **826**, 143–145 (1987)](https://stage.eurekamag.com/research/109/661/109661106.php) gives the application directly. **Gauss–Szegő quadrature:** [Bultheel–Cantero–Cruz-Barroso (2015)](https://doi.org/10.1016/j.cam.2014.11.002) surveys unit-circle quadrature realizations. The quadrature vector measure must not be mistaken for the normalized trace measure.

### 4.2 Reflections and constrained low-rank updates

A unitary rank-one phase update is
\[
R(v,\phi)=I+(e^{i\phi}-1)vv^*,\qquad v^*v=1.
\]
The Householder reflector is \(\phi=\pi\). Its action is \(x+(e^{i\phi}-1)v(v^*x)\), requiring \(O(N)\) work. More generally, with \(V^*V=I_r\), \(K\in U(r)\), and a fast unitary \(U_0\),
\[
U=U_0\bigl[I+V(K-I_r)V^*\bigr]
\]
is **exactly unitary**, stores \(V,K\) plus \(U_0\), and costs \(T(U_0)+O(Nr+r^2)\) per vector. The Hermitian projector \(VV^*\) makes this an invariant-subspace unitary update. An arbitrary \(U_0+XY^*\) is **not** automatically unitary.

Relevant examples: [grover_operator](../grover_operator.m), [unitary_rank_one_update](../unitary_rank_one_update.m), [entangled_reflection](../entangled_reflection.m), and geometrically motivated low-rank transport maps. [References 4, 5.]


**Applications and primary references.** **Quantum search and amplitude amplification:** [Grover, STOC (1996)](https://doi.org/10.1145/237814.237866) uses products of low-rank reflections about marked and uniform states. **Orthogonal transformations in numerical linear algebra:** [Householder, *J. ACM* (1958)](https://doi.org/10.1145/320941.320947) gives the classical reflector application to triangularization. **Geometric optimization:** [Edelman–Arias–Smith, *Geometry of algorithms with orthogonality constraints* (1998)](https://doi.org/10.1137/S0895479895290954) develops methods on orthogonality-constrained manifolds, where low-dimensional rotations/retractions are useful; this is an *application of orthogonality*, not a claim that every manifold update has the particular rank-\(r\) formula above. **Rootfinding distinction:** [Aurentz et al., *Fast and backward stable computation of roots of polynomials* (2015)](https://doi.org/10.1137/140983434) exploits a **unitary-plus-rank-one companion matrix**, whose *full matrix is generally not unitary*.

### 4.3 FFT-diagonal, circulant, Fourier and chirp factorizations

Let \(F\) be the normalized DFT and \(D=\operatorname{diag}(e^{i\phi_j})\). Then
\[
U=F^*DF
\]
is an exactly unitary **circulant** matrix, generally dense, with \(O(N)\) stored phases and \(O(N\log N)\) matvecs. Multiplying several diagonal and Fourier stages (e.g., \(D_1F^*D_2F\)) retains unitary structure; complexity is proportional to the number of FFT stages. Not every Toeplitz matrix is unitary, and the fact that a Toeplitz matrix is FFT-multipliable is **not** a proof of unitarity.

Related exactly normalized Fourier, DCT/DST, Hartley, Zak, Weyl–Heisenberg and Walsh constructions exploit FFT, permutations, diagonal chirps, or fast butterfly recursions. Conventional DCT/DST type and endpoints matter. Examples: [phase_circulant](../phase_circulant.m), [kicked_harper](../kicked_harper.m), [discrete_fourier](../discrete_fourier.m), [orthogonal_dct](../orthogonal_dct.m), [discrete_sine](../discrete_sine.m), [walsh_hadamard](../walsh_hadamard.m), and [affine_fourier](../affine_fourier.m). [References 6, 7.]


**Applications and primary references.** **Periodic convolution and Fourier spectral discretization:** [Davis, *Circulant Matrices* (1979)](https://books.google.com/books?q=Philip+J.+Davis+Circulant+Matrices+1979) and [Cooley–Tukey, FFT (1965)](https://doi.org/10.1090/S0025-5718-1965-0178586-1) underpin fast Fourier-diagonal actions; only unimodular Fourier multipliers define a unitary circulant. **Quantum-chaotic Floquet models:** [Izrailev, *Simple models of quantum chaos* (1990)](https://doi.org/10.1016/0370-1573(90)90067-C) discusses kicked-rotor models, and [Artuso–Casati–Shepelyansky, kicked Harper (1992)](https://doi.org/10.1103/PhysRevLett.68.3826) motivates diagonal-phase/FFT compositions. **Fast trigonometric and algebraic signal processing:** [Püschel–Moura (2008)](https://users.ece.cmu.edu/~moura/papers/t-sp-aug08-part1--pueschelmoura-ieeexplore.pdf) supplies transform-algebra context.

### 4.4 Wavelets and sparse-stage factorizations

Orthonormal compactly supported wavelet transforms have a fast multilevel filter-bank algorithm. For a fixed filter length, the ordinary decimated discrete wavelet transform typically requires \(O(N)\) operations even though its complete matrix is not usually banded. A finite transform is orthogonal only for **compatible normalization and boundary extension**.

More generally, if
\[
U=S_mS_{m-1}\cdots S_1,\qquad S_j^*S_j=I,
\]
then \(U\) is exactly unitary. Storage and action cost are respectively the **sum** of stored stage data and work across stages. This covers local gate circuits, Fourier butterflies, optical networks, lossless filter banks and many Floquet models. It does not imply fast application without a stage-count bound [6–8]. See [daubechies_wavelet](../daubechies_wavelet.m), [optical_mesh](../optical_mesh.m), [dual_unitary_circuit](../dual_unitary_circuit.m), and [orthogonal_butterfly](../orthogonal_butterfly.m).


**Applications and primary references.** **Multiresolution analysis, compression and denoising:** [Mallat (1989)](https://doi.org/10.1109/34.192463) and [Daubechies, *Ten Lectures on Wavelets* (1992)](https://doi.org/10.1137/1.9781611970104). **Lossless multirate filter banks:** [Soman–Vaidyanathan, *On orthonormal wavelets and paraunitary filter banks* (1993)](https://doi.org/10.1109/78.205722) establishes the link between orthonormal wavelet bases and paraunitary polyphase networks. **Programmable linear optics:** [Reck et al. (1994)](https://doi.org/10.1103/PhysRevLett.73.58) and [Clements et al. (2016)](https://doi.org/10.1364/OPTICA.3.001460) use beam-splitter/phase meshes to realize arbitrary lossless finite unitaries; the *universal meshes generally contain \(\Theta(N^2)\) cores*, so their universal parametrization is not a subquadratic method. **Local many-body circuits:** [Bertini–Kos–Prosen (2019)](https://doi.org/10.1103/PhysRevLett.123.210601) gives dual-unitary dynamics, a circuit structure distinct from an arbitrary butterfly approximation.

### 4.5 Tensor and Kronecker products

For \(U=\bigotimes_{j=1}^LU_j\), each \(U_j\in U(d_j)\), \(N=\prod_jd_j\). Only \(\sum_jd_j^2\) scalar entries need be stored. Reshape a length-\(N\) vector into a tensor and apply the small matrices to each mode. With dense \(d_j\times d_j\) factors the action takes \(O(N\sum_jd_j)\), not \(O(N^2)\) when the local dimensions remain small. An **arbitrary dense input vector still uses \(\Omega(N)\) memory and arithmetic**. The Kronecker matrix need never be materialized.

See [haar_tensor_power](../haar_tensor_power.m), [haar_block_tensor_power](../haar_block_tensor_power.m), [tensor_permutation](../tensor_permutation.m), [tensor_fem_wave](../tensor_fem_wave.m), and [um_tensor_apply](../um_tensor_apply.m). Extra fast local factors can reduce the modewise cost.


**Applications and primary references.** **Quantum information:** [Nielsen–Chuang, *Quantum Computation and Quantum Information* (2010)](https://doi.org/10.1017/CBO9780511976667) treats local gates acting on tensor-product Hilbert spaces. **Quantum channels:** [Watrous, *The Theory of Quantum Information* (2018)](https://cs.uwaterloo.ca/~watrous/TQI/) describes superoperators, including the exact unitary-conjugation representation \(\overline V\otimes V\), whose order is \(d^2\). **Separable multidimensional PDE discretizations:** tensor-product sine transforms diagonalize particular constant-coefficient or separable wave/elliptic operators, and the repository's [tensor_fem_wave](../tensor_fem_wave.m) applies this structure to a mass-weighted finite-element wave propagator. A *general nonseparable PDE coefficient field* cannot be assumed to retain tensor diagonalization.

### 4.6 Cayley transforms and rational implicit operators

For Hermitian \(H\), the Cayley transform
\[
U=(H-iI)(H+iI)^{-1}
\]
is exactly unitary and is generally dense even when \(H\) is sparse. Applying it through an already factored shifted sparse system can be inexpensive. For Hermitian tridiagonal \(H\), banded solves are \(O(N)\) after \(O(N)\) setup; higher-dimensional sparse problems may have substantial fill. Never claim \(O(N)\) for arbitrary sparse matrices or ill-conditioned solves. [References 4, 9.] See [cayley_jacobi](../cayley_jacobi.m), [nonuniform_fem_wave](../nonuniform_fem_wave.m), and [um_cayley_apply](../um_cayley_apply.m).

Likewise, \(\exp(iH)\) is unitary for Hermitian \(H\), but applying it with Krylov, Taylor, Chebyshev or other exponential-action methods is usually **approximate** and iteration/degree dependent. Store the sparse generator; do not infer that the exponential is sparse or that its action has a uniform linear-time bound [10]. The existing [IMPLEMENTATION_AUDIT.md](IMPLEMENTATION_AUDIT.md) documents these distinctions.


**Applications and primary references.** **Symmetric-tridiagonal eigenvalue algorithms:** [Aurentz–Mach–Vandebril–Watkins, *Computing eigenvalues ... via a Cayley transformation* (2017)](https://www.maths.tcd.ie/EMIS/journals/ETNA/volumes/2011-2020/vol46/abstract_vol46_pp447-459.html?pages=447-459&vol=46) exploits the circle-valued image of a Hermitian problem. **Geometric and energy-preserving time integration:** [Hairer–Lubich–Wanner, *Geometric Numerical Integration*](https://doi.org/10.1007/3-540-30666-8) describes structure-preserving numerical integration, for which Cayley/implicit midpoint steps are important instances. **Large-scale exponential propagation:** [Al-Mohy–Higham, *Computing the Action of the Matrix Exponential* (2011)](https://doi.org/10.1137/100788860) derives accuracy-aware actions without assembling the full exponential. A *time-stepping approximation* to the exponential should be identified separately from the exact unitary exponential.

## 5. Low-rank and hierarchical matrix structure

These formats are important to numerical linear algebra even though a generic member is **not** unitary. They should enter this database as a fast unitary representation only after an exact structural construction or a validated unitary-preserving approximation is supplied.

### 5.1 Quasiseparable and semiseparable matrices

A matrix is \(r\)-**quasiseparable** if off-diagonal submatrices crossing any index split have rank at most \(r\) (upper and lower orders may differ). Generators encode low-rank factors and their transitions. With fixed generator order, matrix–vector products are linear or near-linear in \(N\), with constants/powers depending on the chosen generator format. The class includes unreduced unitary Hessenberg matrices and other low-rank off-diagonal patterns. A quasiseparable matrix need not be unitary: **rank geometry and orthogonality are separate constraints** [4, 11, 12].

Semiseparable terminology varies across sources; record the precise rank definition (strict upper/lower off-diagonal blocks, inclusion of diagonal, generator convention), not only a broad label.


**Applications and primary references.** **Fast unitary eigenvalue and eigenphase counting methods:** [Aurentz–Mach–Robol–Vandebril–Watkins, *Core-Chasing Algorithms* (2018)](https://doi.org/10.1137/1.9781611975345) analyzes the structured algorithms, while [Eidelman–Haimovici, *Bisection for unitary Hessenberg matrices via quasiseparable structure* (2023)](https://cris.tau.ac.il/en/publications/the-bisection-eigenvalue-method-for-unitary-hessenberg-matrices-v/) is directly concerned with eigenvalue counting. **Structured linear systems and recurrence-based computations:** [Vandebril–Van Barel–Mastronardi, *Matrix Computations and Semiseparable Matrices*, vols. 1–2](https://www.press.jhu.edu/books/title/9714/matrix-computations-and-semiseparable-matrices) provides the broader generator framework. The second application category also includes nonunitary matrices; the unitary subclass needs explicit constraints.

### 5.2 HODLR, HSS and related nested-basis formats

**HODLR** recursively approximates sibling off-diagonal blocks by low-rank products. Uniform rank \(r\) yields typical \(O(Nr\log N)\) storage/matvec cost. **HSS** adds nested bases, potentially improving asymptotic constants and reaching linear complexity for fixed ranks and leaf sizes. An HSS representation with dense leaf blocks of size \(b\) has an additional dependence on \(b\) and on transfer matrices, roughly \(O(Nr+(N/b)r^2+Nb)\) in a common accounting. Rank, tree geometry and accuracy must be reported [13, 14].

**\(\mathcal H\), \(\mathcal H^2\)**, fast multipole and nested interpolation formats address related operator/kernel compression problems, with application-dependent admissibility and ranks. These representations **do not automatically enforce \(U^*U=I\)**. Postcompressing a dense unitary may produce a matrix that is only approximately unitary; applying an exact polar correction need not preserve low ranks or fast action.


**Applications and primary references.** **Discretized elliptic PDEs and integral equations:** [Xia–Chandrasekaran–Gu–Li, *Fast algorithms for hierarchically semiseparable matrices* (2010)](https://doi.org/10.1002/nla.691) develops HSS methods for large rank-structured problems. **Compressed direct solvers and kernel matrices:** [Xia–Chandrasekaran–Gu–Li, *Superfast Multifrontal Method for Large Structured Linear Systems of Equations* (2010)](https://doi.org/10.1137/09074543X) applies structured low-rank compression within a sparse direct solver, and [the HSS reference software](https://github.com/fastsolvers/HSS) provides numerical examples. These works concern hierarchical matrix operations and solvers, **not** a general construction of exactly unitary HSS matrices. For the database, unitary scattering/propagation applications would require a separately verified unitary parametrization or explicit defect/error bounds.

### 5.3 Toeplitz, Hankel, Vandermonde, Cauchy and displacement rank

Low displacement rank means that an operator such as
\[
\nabla_{A,B}(M)=AM-MB
\]
(or the Stein analogue) has low rank, so that the structured matrix can be described by generator matrices plus the displacement pair. Toeplitz-like, Cauchy-like and Vandermonde-like classes are key examples. Suitable displacement operators support fast or superfast matvec algorithms, but their complexity and conditioning depend on the operators, nodes, representation and preprocessing [15, 16].

**An arbitrary Toeplitz, Hankel, Cauchy, or Vandermonde matrix is not unitary.** In a unitary database this heading is a *candidate structural representation*, not an independent guarantee of unitary admissibility. Fourier, orthogonal-polynomial and specially normalized transform matrices can exhibit overlapping displacement and fast-transform structure.


**Applications and primary references.** **Time-series covariance, prediction and harmonic estimation:** [Gragg, *Toeplitz matrices, isometric Arnoldi and Gaussian quadrature* (1993)](https://doi.org/10.1016/0377-0427(93)90294-L) connects Toeplitz moments to a **derived unitary Hessenberg representation**; the covariance matrix itself is not unitary. **Fast structured linear systems:** [Kailath–Kung–Morf, *Displacement ranks of matrices and linear equations* (1979)](https://www.sciencedirect.com/science/article/pii/0022247X79901240) supplies the classical framework, and [Xia et al., *Superfast structured solver for Toeplitz linear systems* (2012)](https://doi.org/10.1137/110831982) applies displacement/HSS ideas. **Polynomial coefficient conversion:** [Townsend–Webb–Olver, *Fast polynomial transforms based on Toeplitz and Hankel matrices* (2018)](https://arxiv.org/abs/1604.07486) exploits Toeplitz–Hankel structure. These conversions are **not automatically unitary matrices**, even when the source/target polynomial bases are orthogonal under related but distinct inner products.

### 5.4 Complementary low rank and butterfly compression

Fourier integral operators and some high-frequency/special-function transforms exhibit the **complementary low-rank property**: blocks indexed by fine groups in one variable and coarse groups in the other have low numerical rank. Interpolative butterfly factorizations exploit this to achieve approximately \(O(Nr\log N)\) storage and action for bounded numerical ranks, with setup costs that can be much larger [17, 18].

Distinguish:
- **Exact unitary butterfly circuits:** products of sparse *unitary* stages; unitarity is guaranteed by construction.
- **Numerical butterfly factorizations:** low-rank approximations to a matrix; usually neither exact nor automatically unitary.
- **Fast polynomial-transform algorithms:** may use asymptotic expansions, recurrences, Fourier transforms, or butterfly-like compression; individual algorithms require their own stability and complexity assumptions.

The distinction matters for tests that need eigenvalues **exactly** on the unit circle, or for repeated unitary-power applications.


**Applications and primary references.** **Wave-equation propagation, high-frequency acoustics, reflection seismology and tomography:** [Candès–Demanet–Ying, *Fast Computation of Fourier Integral Operators* (2007)](https://doi.org/10.1137/060671139) and [Li–Yang–Ying, *Multiscale Butterfly Algorithm for Multidimensional Fourier Integral Operators* (2015)](https://doi.org/10.1137/140997658). **Reusable black-box operator compression:** [Li–Yang–Martin–Ho–Ying, *Butterfly Factorization* (2015)](https://doi.org/10.1137/15M1007173) constructs compact approximations for complementary-low-rank matrices. FIO discretizations may have amplitudes, damping, quadrature weights or boundary treatments that make the matrix **nonunitary**; being amenable to a butterfly algorithm does not establish a unitary spectral benchmark.

### 5.5 Orthogonal-polynomial and spherical-harmonic transforms

Properly weighted *square* orthogonal-polynomial evaluation matrices at a complete Gaussian quadrature rule are orthogonal in exact arithmetic, but **do not automatically have an \(O(N\log N)\) transform**. Chebyshev transforms have FFT reductions; faster Jacobi and spherical-harmonic conversion algorithms require more specialized representations and often nontrivial precomputation. Slevinsky's fast Chebyshev–Jacobi transforms and spherical-harmonic conversion work are specific references to consult [19, 20]. Townsend and Slevinsky's [FastTransforms project](https://mikaelslevinsky.github.io/FastTransforms/index.html) is a useful implementation/literature entry point.

The repository's [gauss_jacobi_transform](../gauss_jacobi_transform.m), [gauss_hermite_transform](../gauss_hermite_transform.m), [gauss_laguerre_transform](../gauss_laguerre_transform.m), and [spherical_harmonic_rotation](../spherical_harmonic_rotation.m) should **not** be assigned fast-matvec tags simply because specialized fast algorithms exist somewhere in the literature. The algorithm actually used by each generator must be checked.


**Applications and primary references.** **Spectral PDE solvers and rapid conversion between polynomial bases:** [Hale–Townsend, *Fast, Simple, and Stable Chebyshev–Legendre Transform* (2014)](https://doi.org/10.1137/130932223), [Slevinsky, *Chebyshev–Jacobi transform* (2018)](https://doi.org/10.1093/imanum/drw070), and [Townsend–Webb–Olver, *Fast polynomial transforms based on Toeplitz and Hankel matrices* (2018)](https://arxiv.org/abs/1604.07486). **Spherical harmonic analysis and geometric computation:** [Slevinsky, *Fast and backward stable transforms between spherical harmonic expansions and bivariate Fourier series* (2019)](https://doi.org/10.1016/j.acha.2017.11.001). **Software/literature collection:** [FastTransforms](https://github.com/MikaelSlevinsky/FastTransforms). Distinguish the *unitary/orthogonal, properly quadrature-weighted finite evaluation transform* from a **general change-of-basis matrix**, which may be nonunitary; moreover, a fast transform from coefficients to values does not automatically give a fast action for a differently normalized database generator.

### 5.6 Tensor trains and matrix-product operators

An operator on \(L\) local \(d\)-dimensional spaces (\(N=d^L\)) can sometimes be represented as a matrix-product operator with moderate bond dimension \(\chi\), requiring \(O(Ld^2\chi^2)\) stored parameters. Products and sums of local operators often admit such descriptions, but bond dimension may grow quickly under composition or time evolution [21].

For **a full arbitrary vector** of length \(N\), applying the compact operator cannot avoid reading/writing \(O(N)\) data, and the contraction cost depends on the bonds and algorithm. Applying an MPO to an **MPS-compressed input** is a different task, with extra compression errors and possible bond growth. Unitarity of an MPO requires global constraints or an exact product of unitary gates; unrestricted local tensors are not automatically a unitary operator.


**Applications and primary references.** **One-dimensional quantum many-body simulation and DMRG:** [Schollwöck, *The density-matrix renormalization group in the age of matrix product states* (2011)](https://arxiv.org/abs/1008.3477) describes MPS/MPO representations; [Paeckel et al., *Time-evolution methods for matrix-product states* (2019)](https://arxiv.org/abs/1901.05824) compares practical evolution algorithms. **Higher-dimensional numerical tensor approximation:** [Oseledets, *Tensor-Train Decomposition* (2011)](https://doi.org/10.1137/090752286) supplies a numerical linear algebra foundation. **Quantum-circuit evolutions:** short-depth local unitary gate sequences can have compact MPO descriptions, but repeated time evolution and truncation can increase bonds and spoil exact unitarity. The physically large state-space order is \(N=d^L\), even if the operator's bond dimensions are small.

## 6. Exactness, unitary preservation and numerical hazards

| Situation | Safe statement | Unsafe inference |
| --- | --- | --- |
| Unitary Hessenberg via cores | Dense matrix, \(O(N)\) descriptor and action | “Hessenberg means sparse matrix” |
| Sparse Hermitian \(H\) | \(\exp(iH)\) exactly unitary as an operator; numerical action can be efficient | “The exponential is sparse” or “always \(O(N)\)” |
| Numerical HODLR/HSS/butterfly approximation | Compact to tolerance under rank assumptions | “The approximation is exactly unitary” |
| Ordinary Toeplitz/Cauchy or graph adjacency | Structured; may admit fast multiplication | “Structured means unitary” |
| Generic Haar/polar matrices | Usually dense; no general exact short descriptor | “Random butterfly has Haar distribution” |
| Quantum circuit on \(q\) qubits | Gate schedule compact; full vector length \(N=2^q\) | “Simulation costs polynomial time in \(q\) for arbitrary dense states” |
| Approximate sparse solves | Applied inverse depends on factorization, fill and residual | “Every sparse Cayley transform is linear time” |
| Square orthogonal polynomial transform | Orthogonal for matched weights/nodes | “Every such transform has a published fast algorithm” |
| \(U_0+XY^*\) | Unitary only under additional algebraic constraints | “Low rank perturbation preserves unitarity” |

**Validation protocol.** For each claimed structure, verify both forward and adjoint actions against explicit small instances; verify \(U^*U\approx I\) in small explicit cases; test multiple right-hand sides and complex adjoints; compare storage/operation counts at growing \(N\); document setup work and any numerical truncation; and test spectral invariants only in the physical/symmetry sectors actually represented. For approximate structures additionally measure \(\|\widetilde U^*\widetilde U-I\|\) and \(\|U-\widetilde U\|\) (when the latter is accessible).

## 7. Suggested structured metadata for the database

A family may have **multiple valid representations** (e.g., explicit sparse CMV, core product, implicit operator, dense diagnostic matrix). Record each representation separately rather than assigning a single unconditional label:

| Metadata field | Meaning / example |
| --- | --- |
| **matrix_order** | \(N\), distinct from underlying physical size or number of factors |
| **exact_unitary_model** | true/false, with mathematical conditions and boundary convention |
| **representation_kind** | explicit_sparse, explicit_dense, core_chain, fft, tensor, rational_solve, hierarchical, butterfly, MPO, other |
| **structural_parameters** | bandwidth, nnz, core count, block size, ranks, bond dimensions, stage depth |
| **representation_exactness** | exact, floating-point factorized, or tolerance-dependent approximation |
| **storage_complexity** | symbolic expression **with parameter dependence** |
| **matvec_complexity / adjoint_complexity** | symbolic expression, work per RHS, assumptions |
| **setup_complexity** | preprocessing, factorization, FFT plans, compressed construction |
| **workspace_complexity** | additional scratch beyond input/output, and RHS batching |
| **unitarity_defect** | exact by construction vs measured \(\|\widetilde U^*\widetilde U-I\|\) |
| **operator_interface** | yes/no; verified action modes and dimensions |
| **available_spectral_reference** | exact finite spectrum, moments, limiting law, or unknown |
| **reference_ids / provenance** | primary literature and code implementation |
| **implementation_status** | implemented, planned, or theoretical reference only |

For the existing MATLAB contract, see [MATRIX_FREE.md](MATRIX_FREE.md). In particular, the repository uses **'transp' for the conjugate transpose**, not the nonconjugating MATLAB **A.'** operation; this needs to remain consistent across new kernels.

## 8. Priorities for this repository

**Already represented:** CMV/Hessenberg core chains; sparse shifts, diagonals, walks and scattering; FFT/circulant/fast trigonometric transforms; Walsh and wavelets; local circuit gates; Kronecker mode products; constrained reflections; and structured Cayley/exponential actions. Check the current [implementation audit](IMPLEMENTATION_AUDIT.md) rather than treating this list as an exhaustive assertion of fast support for every generator.

**High-value follow-up work:**

1. Add **per-representation storage and matvec metadata** to the catalog; distinguish actual implemented costs from literature-possible costs.
2. Add rigorous **operation-count and adjoint/orthogonality tests** for the existing core, FFT, wavelet, tensor, butterfly and Cayley operators.
3. Survey **bounded-rank quasiseparable unitary and block-CMV variants** with explicit unitary generators; connect to fast QR/QZ rather than merely using general HSS approximations.
4. Evaluate **structure-preserving hierarchical or butterfly approximations of unitary operators**, including norm-preserving factorizations, approximation error and repeated-power drift.
5. Examine **fast polynomial transforms** (especially Jacobi/associated Legendre and spherical harmonics) for cases where the *specific orthogonal transform matrix* in the collection can actually use a verified fast algorithm.
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
